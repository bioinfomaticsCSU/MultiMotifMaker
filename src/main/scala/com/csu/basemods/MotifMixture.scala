// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Patrick Marks (pmarks@pacificbiosciences.com)

/**
  2017, Central South University.
  Modified by Tao Li(litao_csu@csu.edu.cn)
 */
package com.csu.basemods

import scala.io.Source
import scala._
import scala.collection.mutable
import java.util.zip.GZIPInputStream
import java.io.{FileInputStream, FileReader, InputStreamReader}
import java.util.concurrent._

import scala.math.Numeric.Implicits._
import java.util.Date
import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer
object Constants {
  // Context size is hard-coded.  It's 41bp -- -20bp to +20 around the cognate base
  val CONTEXT_SIZE: Int = 41
  val CONTEXT_OFFSET: Int = 20
  var CONTEXT_POOLTHREADS:Int = 17
}


// A single position in a motif
case class Cut(
    pos : Int,        // The position of the constrained base, relative to the modification
    selBase : Char)      // The identity of the constrained base.

case class ScorePoint(
  score : Float,
  pos : Int)

  // Define a sort order on GenomePointers -- sort by reference, then by site, then by strand
object ScorePointOrder extends Ordering[ScorePoint]
{
  def compare(a : ScorePoint, that : ScorePoint) : Int =
  {
    var c1 = a.score compare that.score
    if(c1 != 0)
    {
      c1
    }
    else
    {
      a.pos compare that.pos
    }
  }
}


// State of the motif search
case class SearchNode(
    max : Float,      // Maximum possible score of any child node
    nodeScore : Float,      // Actual score of the current node
    penalty : Float,
    nMods : Int,            // Number of modification hits matching the current node
    nGenome : Int,          // Number of genome positions matching the current node
    fractionHit: Float,     // Fraction of genome positions that are modified (nMods / nGenome)
    cuts : Set[Cut])    // The set of cuts that make up the motif

// Modification record -- one modification detected on the genome. Generally this is read in from the modifications.gff of
// the PacBio Base Modification Detection workflow.
case class Modification(
    hitScore : Float,       // Score of hit
    context : Array[Byte],     // Context snippet from gff
    coverage : Int,         // Coverage from gff
    ipdRatio : Float,      // Kinetic IPD ratio from gff
    modType : String,
    refId : String,          // Reference Id string from gff
    genomePointer : GenomePointer)


// Holds an array of tuples containing the fwd and reverse-complement reference strings
case class GenomeInfo(size : Int, contigs : Array[Tuple2[String, String]], nameToId : Map[String, Int])
{
  def get(gp : GenomePointer) =
  {
    var ctg = contigs(gp.refId)

    if(gp.strand == 0)
      ctg._1(gp.site)
    else
    {
      var str = ctg._2
      str(str.length - 1 - gp.site)
    }
  }
}

// Pointer to a position in the genome. RefId is an index into GenomeInfo.contigs, site is the position in that contig, strand is 0 for forward, 1 for reverse
case class GenomePointer(refId : Int, site : Int, strand : Byte)

// Define a sort order on GenomePointers -- sort by reference, then by site, then by strand
object PointerCompare extends Ordering[GenomePointer]
{
  def compare(a : GenomePointer, that : GenomePointer) : Int =
  {
    var c1 = a.refId compare that.refId
    if(c1 != 0)
    {
      c1
    }
    else
    {
      var c2 = a.site compare that.site

      if(c2 != 0)
      {
        c2
      }
      else
      {
          a.strand compare that.strand
      }
    }
  }
}

// Motif helper functions
object Motif {

  def cutString(cuts : Set[Cut]) =
  {
    var startPos = (cuts.minBy { _.pos }).pos
    var endPos = (cuts.maxBy { _.pos }).pos

    def pickChar(idx : Int) =
    {
       cuts find (c => c.pos == idx) match  {
        case Some(cut) => cut.selBase
        case None => 'N'
     }
    }

    var motifSeq = ((startPos to endPos) map pickChar).mkString
    motifSeq
  }

}

// Helper functions for reading raw data
object Reader {

  // Read an (optionally gzipped) text file as an iterator of strings
  def getReader(fn : String) : InputStreamReader = {
    var is =
      if(fn.endsWith("gz"))
        new GZIPInputStream(new FileInputStream(fn), 1e6.toInt)
      else
        new FileInputStream(fn)

    new InputStreamReader(is)
  }


  // Read an (optionally gzipped) text file as an iterator of strings
  def read(fn : String) : Iterator[String] = {
    var is =
      if(fn.endsWith("gz"))
        new GZIPInputStream(new FileInputStream(fn), 1e6.toInt)
      else
        new FileInputStream(fn)

  var lines = Source.fromInputStream(is).getLines()
    lines
  }


  // Reverse complement a string
  def reverseComplement(s : String) : String =
  {
    var r = scala.StringBuilder.newBuilder

    def complement(c : Char) =
      c match
      {
        case 'A' => 'T'
        case 'a' => 'T'
        case 'C' => 'G'
        case 'G' => 'C'
        case 'T' => 'A'
        case 'K' => 'M'
        case 'M' => 'K'
        case 'R' => 'Y'
        case 'Y' => 'R'
        case 'S' => 'S'
        case 'W' => 'W'
        case 'B' => 'V'
        case 'D' => 'H'
        case 'H' => 'D'
        case 'V' => 'B'
        case a  =>  a
      }


    for (i <- 0 until s.length)
    {
      r += complement(s(s.length - 1 - i))
    }

    r.toString()
  }

  def reverseComplement2(s : String) : String =
  {
    def complement(c : Char) =
    {
      c match {
        case 'A' => 'T'
        case 'C' => 'G'
        case 'G' => 'C'
        case 'T' => 'A'
      }
    }

    s.reverseIterator.map(complement).toString
  }

  // Read in the modifications.gff file emitted by the modification detection phase
  def loadModificationsFromGff(f : String, genome : GenomeInfo) =
  {
    // Parse a line of GFF - pull out the fields needed to populate a Modification record
    def parseLine(l : String) =
    {
      if(l(0) == '#')
      {
        None
      }
      else
      {
        var cols = l.split("\t")

        var kvp = cols(8)
        var pairs = kvp.split(";") map (i => { var r = i.split("="); (r(0).trim(), r(1))}) toMap


        var refId = cols(0).split("\\s+")(0)
        //var refNum = refId.substring(3).toInt - 1
        var refNum = genome.nameToId(refId)

        var startPos = cols(3).toInt - 1

        var strand = (if(cols(6)(0) == '+') 0 else 1).toByte

        var modificationType = cols(2)

        var genomePointer = new GenomePointer(refNum, startPos, strand)

        // The 'context' snippet is hard-coded to be a 41 base string covering the cognate base and -20bp to +20bp inclusive
        var context = (pairs("context") map MotifMixture.charToCode).toArray

        // Per-strand coverage
        var coverage = pairs("coverage").toDouble.toInt
        var ipdRatio = pairs("IPDRatio").toFloat

        // The score is the max of the standard 'detection' qv and the identificationQv (if it exists)
        var score =
          pairs.get("identificationQv") match {
              case Some(score) => scala.math.max(cols(5).toFloat, pairs("identificationQv").toFloat)
              case None => cols(5).toFloat
          }

        Some(new Modification(score, context, coverage, ipdRatio, modificationType, refId, genomePointer))
      }
    }

    var lines = read(f)
    (lines map parseLine).flatten
  }


  // Make the special GenomeInfo class from fasta data
  def makeGenomeInfo(fastaFile : String) =
  {
    var fastaData = readFastaContigs(fastaFile)
    var contigs = fastaData.map { case(n,v)  => (v, reverseComplement(v)) }.toArray
    var contigIds = fastaData.map { case(n,_) => n }.zipWithIndex.toMap
    var size = contigs.map(c => c._1.length).sum * 2

    new GenomeInfo(size, contigs, contigIds)
  }


  // Read in raw fasta data
  def readFastaContigs(f : String) : List[(String, String)] =
  {
    var contigs = List[(String,String)]()

    var lines = read(f)

    var currentHeader : Option[String] = None

    var buf : StringBuffer = null

    for(l <- lines)
    {

      if(l.length > 0 && l(0) == '>')
      {
        currentHeader match
        {
          case Some(h) => contigs = contigs ++ List((h, buf.toString()))
          case _ => ()
        }

        currentHeader = Some(l.substring(1).split("\\s+")(0))
        buf = new StringBuffer()
      }
      else
      {
        buf.append(l.toUpperCase())
      }
    }

    // Put in final contig
    currentHeader match
      {
        case Some(h) => contigs = contigs ++ List((h, buf.toString()))
        case _ => ()
      }

    contigs
  }
}

object MotifMixture {

  import Constants._
  val n = CONTEXT_SIZE
  val c = CONTEXT_OFFSET
  var searchBreadthStarting:Int = 16
  var searchBreathEnding:Int = 6
  var depthToSwitching:Int = 3
  var layer:Int = 1
  val poolThreads:Int = CONTEXT_POOLTHREADS


  val threadPool = Executors.newFixedThreadPool(poolThreads*2)
  val branchThreadPool = Executors.newFixedThreadPool(poolThreads)

  def submitTask[T](f : () => T) : Future[T] =
  {
    threadPool.submit(f)
  }

  def submitBranch[T](f : () => T) : Future[T] =
  {
    branchThreadPool.submit(f)
  }

  var parallelize = false

  implicit def callable[T](f: () => T): Callable[T] =  new Callable[T]() { def call() = f() }

  // Intepret IUPAC wobble bases
  def compareBases(b : Char, selection : Char) =
  {
    if(b == selection)
    {
      true
    }
    else
    {
      selection match {
        case 'K' if b == 'G' || b == 'T' => true
        case 'M' if b == 'A' || b == 'C' => true
        case 'R' if b == 'A' || b == 'G' => true
        case 'Y' if b == 'C' || b == 'T' => true
        case 'S' if b == 'C' || b == 'G' => true
        case 'W' if b == 'A' || b == 'T' => true
        case 'B' if b != 'A' => true
        case 'D' if b != 'C' => true
        case 'H' if b != 'G' => true
        case 'V' if b != 'T' => true
        case _ => false
      }
    }
  }


  // Char of DNA base to unsigned byte
  def charToCode(c : Char) : Byte =
  {
    c match {
      case 'N' => 4
      case 'A' => 0
      case 'C' => 1
      case 'G' => 2
      case 'T' => 3
      case  _  => 4
    }
  }

  // Unsigned byte to Char of DNA base
  def codeToChar(b : Byte) : Char =
  {
    b match {
      case 0 => 'A'
      case 1 => 'C'
      case 2 => 'G'
      case 3 => 'T'
      case 4 => 'N'
    }
  }

  // What can you get by adding bases to a wobble
    def promotions = List(
                        ('A', List('M', 'R', 'W')),
                        ('C', List('M', 'Y', 'S')),
                        ('G', List('K', 'R', 'S')),
                        ('T', List('K', 'Y', 'W')),
                        ('K', List('B', 'D')),
                        ('M', List('V', 'H')),
                        ('R', List('V', 'D')),
                        ('Y', List('H', 'B')),
                        ('S', List('V', 'B')),
                        ('W', List('H', 'D')),
                        ('B', List()),
                        ('D', List()),
                        ('H', List()),
                        ('V', List())) toMap


  // Does context match the motif 'cut'
  def passCut(cut : Cut, context : Array[Byte]) : Boolean =
  {
    compareBases(codeToChar(context(cut.pos + c)), cut.selBase)
  }
  // Does context match the motif 'cuts'
  def passCuts(base:Char,cuts : Set[Cut]) : Boolean =
  {
    cuts forall (cut =>compareBases(base, cut.selBase))
  }

  // Test if a whole motif is present in the modification context
  def hasCuts(cuts : Set[Cut]) (example : Modification) : Boolean =
  {
    cuts forall (c => passCut(c, example.context))
  }



  // Take a list of modifications and summarize the number and total score of the remaining modifications after all possible 1 base cuts are applied
    // returns (int[,] * float[,]) containing the count and sum-of-scores respectively.  The arrays are laid indexed as [baseCode, positionOfCut], so the
    // number of remaining mods for a cut 'C' at position -2  would be in counts.[(charToCode 'C'), -2 + !c]
  def countByColumn(examples : Array[Modification]) : Tuple2[Array[Array[Int]], Array[Array[Float]]] =
  {
    var counts = Array.fill[Int](5,n)(0)
    var scores = Array.fill[Float](5,n)((0.0).toFloat)

    var i = 0
    var numExamples = examples.length

    while(i < numExamples)
    //for (ex <- examples)
    {
      val ex = examples(i)
      var ctx = ex.context
      var score = ex.hitScore
      var p = 0

      while(p < n)
      //for (p <- 0 to n - 1)
      {
        var b = ctx(p)
        counts(b)(p) += 1
        scores(b)(p) += score
        p += 1
      }

      i += 1
    }

    (counts, scores)
  }


    // tests for whether the 'cuts' are all present at a given position in the string genome
    def testGenomeForMotif(genome : String, cuts : Set[Cut], position : Int) =
    {
      cuts forall (c => compareBases(genome(position + c.pos), c.selBase))
    }

    // Find all the positions in 'genomeInfo' that match 'cuts'
    def findMotifInstances (genomeInfo : GenomeInfo, cuts : Set[Cut]) =
    {
      var cutPositions = cuts map (c => c.pos)
      var hits = new scala.collection.mutable.ArrayBuffer[GenomePointer]()

      var cutArr = cuts.toArray

      if(cutPositions.size > 0)
        {
          var pre = - (cutPositions min)
          var post = cutPositions max

          for(i <- 0 to genomeInfo.contigs.length - 1)
          {
            var (fwdSeq, revSeq) = genomeInfo.contigs(i)
            var len = fwdSeq.length()

            var first = pre
            var last = fwdSeq.length - post - 1

            var pos = first
            while(pos <= last)
            {
              //if( testGenomeForMotif(fwdSeq, cuts, pos))
              if(cutArr forall (c => compareBases(fwdSeq(pos + c.pos), c.selBase)))
                hits += new GenomePointer(i, pos, 0)

              //if( testGenomeForMotif(revSeq, cuts, pos))
              if(cutArr forall (c => compareBases(revSeq(pos + c.pos), c.selBase)))
                hits += new GenomePointer(i, len - 1 - pos, 1)

              pos += 1
            }
          }
        }

      hits
    }

        // Find all the positions in 'genomeInfo' that match 'cuts'
    def findMotifInstances2 (genomeInfo : GenomeInfo, cuts : Set[Cut]) =
    {
      var cutPositions = cuts map (c => c.pos)
      var cutArr = cuts.toArray


      def cutsForContig(seqs : ((String, String), Int)) =
        {
          var ((fwdSeq, revSeq), i) = seqs
          var len = fwdSeq.length()

          var hitsFwd = new scala.collection.mutable.ArrayBuffer[GenomePointer]()
          var hitsRev = new scala.collection.mutable.ArrayBuffer[GenomePointer]()

          var pre = - (cutPositions min)
          var post = cutPositions max
          var first = pre
          var last = fwdSeq.length - post - 1

          var pos = first
          while(pos <= last)
          {
              if(cutArr forall (c => compareBases(fwdSeq(pos + c.pos), c.selBase)))
                  hitsFwd += new GenomePointer(i, pos, 0)

              if(cutArr forall (c => compareBases(revSeq(pos + c.pos), c.selBase)))
                hitsRev += new GenomePointer(i, pos, 1)

              pos += 1
          }
          (hitsFwd.toArray, hitsRev.toArray)
        }

      genomeInfo.contigs.zipWithIndex.map(cutsForContig)
    }

    // Pick a strand from the contig
    def selectStrand(strand : Byte, contig : Tuple2[String, String]) = if(strand == 0) contig._1 else contig._2


    // Count the number of given genome instances that have any possible next cut
      // Return the cut count as an int[,]
    def countSubCuts (genomeInfo : GenomeInfo, instances : Array[Tuple2[Array[GenomePointer], Array[GenomePointer]]]) =
    {
      var counts = Array.ofDim[Int](4, n)

      var w = n - c - 1

      def innerCount (sequence : String, hits : Array[GenomePointer]) =
      {
          var j = 0
          var l = hits.length
          while(j < l)
        {
            var hit = hits(j)
            var start = Math.max(0, hit.site - w)
            var last = Math.min(sequence.length - 1, hit.site + w)

            var offset = -hit.site + c
            var p = start

            //for(p <- start to last)
            while(p <= last)
            {
              var idx = charToCode(sequence(p))
              if(idx < 4)
                  counts(idx)(p + offset) += 1

                p += 1
            }

            j += 1
        }
      }

      for (i <- 0 to instances.length - 1)
      {
          var (fwdSeq, revSeq) = genomeInfo.contigs(i)
          var (fwdHits, revHits) = instances(i)

          innerCount(fwdSeq, fwdHits)
          innerCount(revSeq, revHits)
      }

      counts
    }


    // Count the number of given genome instances that have any possible next cut
      // Return the cut count as an int[,]
    def countSubCutsPar (genomeInfo : GenomeInfo, instances : Array[Tuple2[Array[GenomePointer], Array[GenomePointer]]]) =
    {

      def innerCount (sequence : String, hits : Array[GenomePointer]) : Array[Array[Int]] =
      {
          var w = n - c - 1
          var counts = Array.ofDim[Int](4, n)
          var j = 0
          var l = hits.length
          while(j < l)
        {
            var hit = hits(j)
            var start = Math.max(0, hit.site - w)
            var last = Math.min(sequence.length - 1, hit.site + w)

            var offset = -hit.site + c
            var p = start

            //for(p <- start to last)
            while(p <= last)
            {
              var idx = charToCode(sequence(p))
              if(idx < 4)
                  counts(idx)(p + offset) += 1

                p += 1
            }

            j += 1
        }

          counts
      }


      var finalCounts = Array.ofDim[Int](4, n)

      // Accumulate the counts we've got
      def addCounts(addArr : Array[Array[Int]]) =
      {
        var r = 0
        while(r < 4)
        {
          var i = 0
          while(i < n )
          {
            finalCounts(r)(i) += addArr(r)(i)
            i += 1
          }

          r += 1
        }
      }

      // Spawn a job for counting each strand of each contig
      var futures =
        (0 to instances.length - 1).flatMap(i =>
        {
          val (fwdSeq, revSeq) = genomeInfo.contigs(i)
          val (fwdHits, revHits) = instances(i)

          val ff = submitTask(() => innerCount(fwdSeq, fwdHits) )
          val fr = submitTask(() => innerCount(revSeq, revHits) )
          List(ff,fr)
      })

      // Sum up the counts from the different contigs
      for (f <- futures)
      {
        addCounts(f.get())
      }

      finalCounts
    }

    // GenomePointer in a layout the mirrors GenomeInfo.contig
    // We use this special representation for performance
    type PointerSet = Array[Tuple2[Array[GenomePointer], Array[GenomePointer]]]

    // Filter down a list of genomePointer based on '
    def filterGenomePointers(genomeInfo : GenomeInfo, instances : PointerSet, newCut : Cut) =
    {
      var w = n - c
      def filterHits ( sequence : String, hits : Array[GenomePointer] ) =
      {
        var hitBuf = hits filter (h => {
          var testPos = h.site + newCut.pos
          testPos >= 0 && testPos < sequence.length && compareBases(sequence(testPos), newCut.selBase)
        })

        hitBuf
      }

      var r = instances.zipWithIndex map (arg => {
        var ((fwdHits, revHits), idx) = arg
        var (fwdSeq, revSeq) = genomeInfo.contigs(idx)
        (filterHits(fwdSeq, fwdHits), filterHits(revSeq, revHits))
        })

      r
    }

    def filterGenomePointers2(genomeInfo : GenomeInfo, instances : PointerSet, cuts : Set[Cut]) =
    {
      var w = n - c
      def filterHits ( sequence : String, hits : Array[GenomePointer] ) =
      {
        var hitBuf = hits filter (h => {
          cuts forall(newCut => {
            var testPos = h.site + newCut.pos
            testPos >= 0 && testPos < sequence.length && compareBases(sequence(testPos), newCut.selBase)
          })
        })

        hitBuf
      }

      var r = instances.zipWithIndex map (arg => {
        var ((fwdHits, revHits), idx) = arg
        var (fwdSeq, revSeq) = genomeInfo.contigs(idx)
        (filterHits(fwdSeq, fwdHits), filterHits(revSeq, revHits))
      })

      r
    }

    val canonicalBases = List('A', 'C', 'G', 'T')
    val iupacBases = List('A','C','G','T','K','M','R','Y','S','W','B','D','H','V')

    val wobbleMaps = List(
                        ('A', List('A')),
                        ('C', List('C')),
                        ('G', List('G')),
                        ('T', List('T')),
                        ('K', List('G', 'T')),
                        ('M', List('A', 'C')),
                        ('R', List('A', 'G')),
                        ('Y', List('C', 'T')),
                        ('S', List('G', 'C')),
                        ('W', List('A', 'T')),
                        ('B', List('C' ,'G', 'T')),
                        ('D', List('A', 'G', 'T')),
                        ('H', List('A', 'C', 'T')),
                        ('V', List('A', 'C', 'G')))

    val baseIdx = canonicalBases.zipWithIndex toMap
    //val baseRow = (bases.zipWithIndex map (c => (c._2, c._1))) toMap
      val wobbleIdx = (wobbleMaps map { case(b, l) => (b, l map baseIdx toArray) }).toMap

    // Generate the list of all possible cuts -- we can cut on any canonical base anywhere in the 41bp context
    def enumerateCuts(currentCuts : Set[Cut], baseSet : List[Char]) =
    {
      var l = new mutable.ArrayBuffer[Cut]()

      for(b <- baseSet)
        l.append(Cut(0, b))

      for(i <- 1 to c; b <- baseSet)
        {
        l.append(Cut(-i, b))
        l.append(Cut(i, b))
        }

      var currentPos = currentCuts map (c => c.pos)
      l  filter (c => ! (currentPos.contains(c.pos)))
    }

    // This converts the 'methylated fraction' to a score
    // The branch-and-bound search approach only works if this score is guaranteed to less than 1
    def fractionScore(f : Double) = Math.pow(f, 0.9).toFloat

      // Limits on how small / rare the motif can be
    def maxContextSize = 9

    // Only call motifs with at least this number of hits
    def minExamples = 5

    // Motifs must show up at least this often in the genome
    def minMotifRateInGenome = 2e-5

    // Motifs with less that this number of genome instances get a more stringent handling -- they aren't allowed to have too many Ns in them
    def rareMotifThreshold = 1000

    // Maximum number of Ns that a rare motif is allowed to have
    def maxNsRareMotifs = 9


    // Rare motifs (below minMotifRateInGenome) will be called if they are strongly modified
    def rateMotifMinModifiedFraction = 0.75


    // Build a function that scores all possible sub-cuts of a motif
    // A sub-cut of a motif is any new motif that has one added cut positions
    // This function tabulates the number of genome instance of each sub-cut, and the total number and score of matching modifications.
    def makeScorer(genomeInfo : GenomeInfo, baseSet : List[Char], leftToRight : Boolean) =
    {
        // all possible cuts
      var allCuts = enumerateCuts(Set(), baseSet)

      // Compute the score of all possible motifs that have one additional cut added to 'currentCut'
          // Only returns new cuts whose maximum possible score exceeds 'minScore.
          // mods are the modifications matching currentCut
          // genomeInstances are the genomePointers matching currentCut
      def scoreSubCuts(minScore : Double, currentCut : Set[Cut], mods : Array[Modification], genomeInstances : PointerSet) =
      {
          // How many genome instances are present for each sub cut
        var genomeCounts =
            if(parallelize)
              countSubCutsPar(genomeInfo, genomeInstances)
            else
              countSubCuts(genomeInfo, genomeInstances)

        // How many modifications are present for each sub cut (and their total scores)
        var (counts, scores) = countByColumn(mods)

        val currentCutPositions = currentCut map { _.pos }
        val maxCutPosition = currentCutPositions.fold(Int.MinValue)(Math.max)

        def allowedCutPosition(p : Int) =
          if (! leftToRight)
            ! currentCutPositions.contains(p)
          else
              p > maxCutPosition

        // find the smallest two score entries not on existing cut positions
        // to achieve an improvement over the current node score, a further cut will have to discard
        // at least that much score.  This lets us reduce our max estimate for a sub cut.

        var pq = new scala.collection.mutable.PriorityQueue[ScorePoint]()(ScorePointOrder)


        for(b <- List('A', 'C', 'G', 'T'))
        {
         for(cIdx <- new Range(-c, c, 1))
         {
           if(allowedCutPosition(cIdx))
           {
             pq.enqueue(ScorePoint(-scores(baseIdx(b))(cIdx + c), cIdx))
             pq.drop(10)
           }
         }
        }

        // Find the smallest possible delta in total score
        def minNextCut (cutPos : Int) =
        {
          pq.find { pt => pt.pos != cutPos } match
          {
            case Some(p) => -p.score
            case None => (0).toFloat
          }
        }

        // Implement wobble cuts by summing up the appropriate rows of the counts matrix
        def sumRows[T] (arr : Array[Array[T]],  rows : Array[Int],  col : Int)(implicit numeric: scala.math.Numeric[T]) =
        {
          var i = 0
          var sum = numeric.zero
          val n = rows.length

          while(i < n)
          {
            sum = sum + arr(rows(i))(col)
            i+=1
          }
          sum
        }

        // Cut extents
        var firstCut = currentCut.map(c => c.pos).fold(Int.MaxValue)((a,b) => Math.min(a,b))
        var lastCut = currentCut.map(c => c.pos).fold(Int.MinValue)((a,b) => Math.max(a,b))
        var nCuts = currentCut.size


        // Possibly construct SearchNode for a sub-cut.  This SearchNode summarizes the cuts, the current score and maximum possible of any descendant of this cut
        // This is where the 'bound' of branch and bound happens.  Sub-cuts whose maxScore is less than minScore are eliminated - they can never beat the current best motif.
        def scoreCut(newCut : Cut) =
        {
          var tryCut = allowedCutPosition(newCut.pos)

          if(tryCut)
          {
              var pos = newCut.pos + c
              var rows = wobbleIdx(newCut.selBase)

              // Compute the raw material for the score by summing over the rows that make up the wobble
              var nPass = sumRows(counts, rows, pos)

              var genomeCount = sumRows(genomeCounts, rows, pos)

              var fracHit = Math.min(1, nPass.toFloat/genomeCount.toFloat)

              // Best possible sub-cut score
              var totalScore = sumRows(scores, rows, pos)

                // Score of the current node
              var nodeScore = totalScore * fractionScore(fracHit)

              // minimum possible penalty
              var minCutPenalty = minNextCut(newCut.pos)

              // Best possible current score
              var maxScore = Math.max(totalScore - minCutPenalty, nodeScore)

              var firstBase = Math.min(firstCut, newCut.pos)
              var lastBase = Math.max(lastCut, newCut.pos)
              var numNs = lastBase - firstBase + 1 - (nCuts + 1)


              if( // bound condition
                  maxScore > minScore &&

                  // Minimum abundance of motif in genome, or a good hit
                  ((genomeCount.toDouble / genomeInfo.size) >= minMotifRateInGenome || (fracHit > rateMotifMinModifiedFraction && numNs <= maxNsRareMotifs)) &&

                  // If the motif is 'rare', then it can't have too many Ns
                  (genomeCount > rareMotifThreshold || numNs <= maxNsRareMotifs) &&

                  // Minimum number of hits
                  nPass >= minExamples
              )
              {
                  Some(
                    new SearchNode(
                      maxScore,
                      nodeScore,
                      minCutPenalty,
                      nPass,
                      genomeCount,
                      fracHit,
                      currentCut + newCut), newCut)
              }
              else
                None
          }
          else
          {
             None
          }
        }

        // Return nodes sorted in descending current score order
        var scoredCuts = (allCuts map scoreCut).flatten
        scoredCuts.sortBy { -_._1.nodeScore }
      }

      scoreSubCuts _
    }


    // Score a cut the hard way -- go through the entire genome looking for instances
    def scoreCut(genomeInfo : GenomeInfo, modifications : Array[Modification], cuts : Set[Cut]) =
    {
        var genomeHits =
          if (cuts.size == 0)
          {
            genomeInfo.size
          }
          else
          {
            findMotifInstances(genomeInfo, cuts).length
          }

        var matchingExamples = modifications filter (hasCuts(cuts))
        var nMatch = matchingExamples.length

        var maxScore = matchingExamples map { _.hitScore } sum

        var fracHit = Math.min(1.0, nMatch.toFloat / genomeHits.toFloat).toFloat

        var nodeScore = maxScore * fractionScore(fracHit)

        var res = new SearchNode(maxScore, nodeScore,0, nMatch, genomeHits, fracHit, cuts)
        res
    }


    def emptyNode = new SearchNode(0, 0,0, 0, 0, 0, Set.empty)

    def status(s : String) = println(s)

    // The overall search function finds and refines the best motif that exists, given the genomeInfo, and
    // 'allExamples' the list of modification hits.
    def searchAndRefine(genomeInfo : GenomeInfo, allExamples : Array[Modification], allPointers : PointerSet, leftToRight : Boolean) =
    {
        // This function will be used to score sub-cuts of each SearchNode we decide to explore
        var scorer = makeScorer(genomeInfo, canonicalBases, leftToRight)

        // seacher parameterizes the main search function with a searchBreadth
        // searchBreadth is the maximum number of child motifs to test.  If we want to be sure to find the best possible node, we would have to
        // explore all children that are not excluded by the bound condition - that leads to a prohibitively large search space, so we do a best-first
        // exploration of the children, up to a maximum of 'searchBreadth'.  For strong motifs the bounding condition will usually exclude almost all the child nodes.
        def searcher(searchBreadthStart : Int, searchBreadthEnd : Int, depthToSwitch : Int) =
        {
            // For multi-base motifs there are many paths through the search DAG that lead to the same place - this is because we can add positions to the
            // motif in different orders.  To protect against re-exploring some of the same paths, everytime we finish the search down from a SearchNode, we
            // cache the result for reuse if we ever reach this motif again by a different path.
          var currentNodeCache = mutable.Map[Set[Cut], SearchNode]()
          var currentBestNode = new SearchNode(0, 0,0, 0, 0, 0, Set.empty)
          //The first iteration data is divided into 16 groups
          def parallelSearch(initBestNode:SearchNode,initCurrentNode : SearchNode, initExamples : Array[Modification], initGenomeInstances : PointerSet):SearchNode =
          {
            // bail if we are too deep - don't explore absurdly long motifs
            if(initCurrentNode.cuts.size >= maxContextSize || initExamples.length < minExamples)
              return initCurrentNode

            // What is the current best nodeScore in the overall search?
            // To stay in the search, a child node must have a maxScore > branchMinScore
            var initBranchMinScore = Math.max(initBestNode.nodeScore, initCurrentNode.nodeScore)

            // Get the scores of all feasible child nodes, sorted by their nodeScore. Take the top 'searchBreadth' child nodes
            var initAllBranches = scorer(initBranchMinScore, initCurrentNode.cuts, initExamples, initGenomeInstances)

            val initSearchBreadth = if (initCurrentNode.cuts.size >= depthToSwitch) searchBreadthEnd else searchBreadthStart
            var initBranches = initAllBranches.slice(0, Math.min(initSearchBreadth, initAllBranches.length))

            if(initBranches.size == 0)
              return initBestNode

            currentBestNode = initBestNode
            //calculate subBranches from branches of first layer
            if(layer > 1) {
              val temp = initBranches.sortBy(-_._1.nodeScore)
              currentBestNode = temp(0)._1
              var count = 2
              val branchFutures = initBranches.flatMap(branch =>{
                val bf = submitBranch(() =>branchSplit(count,currentBestNode,branch,initExamples,initGenomeInstances))
                List(bf)
              })

              var resultBranches = new ArrayBuffer[(ArrayBuffer[(SearchNode,Cut)],SearchNode)]()

              for(b <- branchFutures){
                try {
                  if(b.get()!= null){
                    resultBranches.append(b.get())
                  }else{
                    System.out.println("FutureTask.get() didn't get the result.")
                  }
                } catch{
                  case e:Exception =>e.printStackTrace()
                }
              }
              var input = resultBranches.flatMap(r => r._1)                        //branches
              var tempNodes = resultBranches.map(r => r._2).sortBy(-_.nodeScore)   //subBest

              if(tempNodes(0).nodeScore > currentBestNode.nodeScore)
                currentBestNode = tempNodes(0)
              initBranches = input
            }

            if(initBranches.size == 0)
              return currentBestNode

            if(leftToRight)
              initBranches = initBranches.sortBy { case(n,c) => c.pos }

            var resultNodesSet = new ArrayBuffer[SearchNode]()

//            println("layer:" + layer + "branchSize:" + initBranches.size)
            //achieve in multiThreads
            val nodeFutures = initBranches.flatMap(branch =>{
              val nf = submitBranch(() =>filterSearch(currentBestNode,branch,initExamples,initGenomeInstances,currentNodeCache))
              List(nf)
            })

            for (b <- nodeFutures)
            {
              try {
                if(b.get()!=null){
                  resultNodesSet.append(b.get())
                }else{
                  System.out.println("FutureTask.get() didn't get the result.")
                }
              } catch{
                case e:Exception =>e.printStackTrace()
              }

            }

            //compare bestSubNodes of different batches to get the candidate best motifSearchNode
            val results = resultNodesSet.sortBy(-_.nodeScore)
            if(currentBestNode.nodeScore < results(0).nodeScore)
              currentBestNode = results(0)

            currentBestNode
          }

          //calculate the subBranches
          def branchSplit(count:Int,BestNode:SearchNode,branch : Tuple2[SearchNode, Cut],examples : Array[Modification],
                          GenomeInstances : PointerSet):(ArrayBuffer[(SearchNode,Cut)],SearchNode)={
            var filterNode = branch._1
            var filterCut = branch._2

            var b = new mutable.ArrayBuffer[(SearchNode,Cut)]()
            var currentBest =
              if(BestNode.nodeScore > filterNode.nodeScore) BestNode else filterNode

            if(filterNode.max < BestNode.nodeScore){
              b.append(branch)
              (b,currentBest)
            }else{
              // Filter down the list of matching modification hits
              var filterExamples = examples filter (ex => passCut(filterCut, ex.context))

              // Filter down the list of matching genome positions
              var filterGenomeInstances = filterGenomePointers(genomeInfo, GenomeInstances, filterCut)


              if(filterNode.cuts.size >= maxContextSize || filterExamples.length < minExamples){
                b.append(branch)
                if(currentBest.nodeScore > currentBestNode.nodeScore){
                  currentBestNode = currentBest
                }
                return (b,currentBest)
              }
              // What is the current best nodeScore in the overall search?
              // To stay in the search, a child node must have a maxScore > branchMinScore
              var BranchMinScore = Math.max(BestNode.nodeScore, filterNode.nodeScore)

              // Get the scores of all feasible child nodes, sorted by their nodeScore. Take the top 'searchBreadth' child nodes
              var allBranches = scorer(BranchMinScore, filterNode.cuts, filterExamples, filterGenomeInstances)

              val searchBreadth = if (filterNode.cuts.size >= depthToSwitch) searchBreadthEnd else searchBreadthStart
              val num = Math.min(searchBreadth, allBranches.length)
              var initBranches = allBranches.slice(0,num)

              if(currentBest.nodeScore > currentBestNode.nodeScore){
                currentBestNode = currentBest
              }

              if(initBranches.size==0){
                b.append(branch)
                return (b,currentBest)
              }

              val nextCount = count + 1
              if(nextCount <= layer) {
                val results = initBranches.map(branch => branchSplit(nextCount,currentBest,branch,filterExamples,filterGenomeInstances))
                var temp = new mutable.ArrayBuffer[(SearchNode,Cut)]()
                var temp1 = new mutable.ArrayBuffer[SearchNode]()
                temp = results.map(r => r._1).flatten
                temp1 = results.map(r => r._2).sortBy(-_.nodeScore)
                if(temp1(0).nodeScore > currentBest.nodeScore)
                  currentBest = temp1(0)
                initBranches = temp
              }

              (initBranches,currentBest)
            }
          }

          //filter work before the iterate function of branchAndBoundMotifSearch
          def filterSearch(initNode : SearchNode, branch : Tuple2[SearchNode, Cut],
                           currentExamples : Array[Modification], currentGenomeInstances : PointerSet,
                           initNodeCache:mutable.Map[Set[Cut], SearchNode]):SearchNode =
          {
            //              System.out.println(initNode + ";" + branch + ";" + currentExamples.size + ";" + currentGenomeInstances.length + ";")
            var nodeCache = mutable.Map[Set[Cut], SearchNode]()
            nodeCache = initNodeCache.clone()

            var filterNode = branch._1
            var filterCuts = branch._1.cuts

            var currentBest =
              if(initNode.nodeScore > filterNode.nodeScore) initNode else filterNode

            if(filterNode.max < initNode.nodeScore)
            {
              initNode
            }else {
              // Filter down the list of matching modification hits
              var filterExamples = currentExamples filter (ex => hasCuts(filterCuts)(ex))

              // Filter down the list of matching genome positions
              var filterGenomeInstances = filterGenomePointers2(genomeInfo, currentGenomeInstances, filterCuts)

              if(filterNode.cuts.size >= maxContextSize || filterExamples.length < minExamples){
                if(currentBest.nodeScore > currentBestNode.nodeScore){
                  currentBestNode = currentBest
                }
                return currentBest
              }

              //run branchAndBoundMotifSearch to find subBest motif
              var motifNode = branchAndBoundMotifSearch(currentBest, filterNode, filterExamples, filterGenomeInstances, nodeCache)

              if(poolThreads != searchBreadthStarting && layer!=1){
                currentNodeCache.synchronized{
                  for (entry <- nodeCache.entrySet) {
                    currentNodeCache.get(entry.getKey) match {
                      case Some(v) => ()
                      case None => currentNodeCache.put(entry.getKey(), entry.getValue())
                    }
                  }
                }
              }

              if (motifNode.nodeScore > currentBest.nodeScore)
              {
                if(motifNode.nodeScore > currentBestNode.nodeScore){
                  currentBestNode = motifNode
                }

                motifNode
              } else {
                if(currentBest.nodeScore > currentBestNode.nodeScore){
                  currentBestNode = currentBest
                }
                currentBest
              }
            }
          }

          //run branchAndBoundMotifSearch to find subBest motif after data filter
          def branchAndBoundMotifSearch(bestOverallNode : SearchNode, currentNode : SearchNode, examples : Array[Modification], genomeInstances : PointerSet,nodeCache:mutable.Map[Set[Cut], SearchNode]) : SearchNode =
          {
            // Bail if node is cached -- we already know the answer
            nodeCache.get(currentNode.cuts) match {
              case Some(n) => return n
              case None => ()
            }

            // bail if we are too deep - don't explore absurdly long motifs
            if(currentNode.cuts.size >= maxContextSize || examples.length < minExamples)
              return currentNode

            // What is the current best nodeScore in the overall search?
            // To stay in the search, a child node must have a maxScore > branchMinScore
            var branchMinScore = Math.max(bestOverallNode.nodeScore, currentNode.nodeScore)

            // Get the scores of all feasible child nodes, sorted by their nodeScore. Take the top 'searchBreadth' child nodes
            var allBranches = scorer(branchMinScore, currentNode.cuts, examples, genomeInstances)

            val searchBreadth = if (currentNode.cuts.size >= depthToSwitch) searchBreadthEnd else searchBreadthStart
            var trialBranches = allBranches.slice(0, Math.min(searchBreadth, allBranches.length))
            //              System.out.println("branch"+trialBranches.size)
            if(leftToRight)
              trialBranches = trialBranches.sortBy { case(n,c) => c.pos }

            // Explore a childNode, while updating the best node if we find a new winner
            def nodeFold (bestLocalBranch : SearchNode, trial : Tuple2[SearchNode, Cut]) =
            {
              // Node to explore
              var trialNode = trial._1

              // newly added cut that turns the current node in the child node we are about to test
              var newCut = trial._2

              // update the best score -- it could come from bestOverallNode, or hae been discovered during the fold
              var currentBest =
                if(bestOverallNode.nodeScore > bestLocalBranch.nodeScore) bestOverallNode else bestLocalBranch

              // Don't bother exploring this one if it can't win
              if(trialNode.max < currentBest.nodeScore)
              {
                bestLocalBranch
              }
              else
              {
                // Explore the child node!

                // Filter down the list of matching modification hits
                var cutExamples = examples filter (ex => passCut(newCut, ex.context))

                // Filter down the list of matching genome positions
                var cutGenomeInstances = filterGenomePointers(genomeInfo, genomeInstances, newCut)

                // Recurse, then see if we have a new winner
                var resultNode = branchAndBoundMotifSearch(currentBest, trialNode, cutExamples, cutGenomeInstances,nodeCache:mutable.Map[Set[Cut], SearchNode])
                if (resultNode.nodeScore > bestLocalBranch.nodeScore) resultNode else bestLocalBranch
              }
            }

            var bestSubNode = trialBranches.foldLeft(currentNode)(nodeFold)

            // Cache the best child of this node.
            nodeCache(currentNode.cuts) = bestSubNode

            bestSubNode
          }

          parallelSearch _
        }


        def findWobbleSubCut(node : SearchNode) : SearchNode =
        {
          // Don't find more than 9
          if (node.cuts.size >= 9)
          {
            node
          }
          else
          {
            var motifInstances = findMotifInstances2(genomeInfo, node.cuts)
            var currentHits = allExamples.filter(hasCuts(node.cuts))

            var scorer = makeScorer(genomeInfo, iupacBases, false)
            var opts = scorer(node.nodeScore, node.cuts, currentHits, motifInstances)

            var newOpts = opts.map(v => v._1).filter(newNode => newNode.nodeScore > node.nodeScore)
            newOpts = newOpts.sortBy(v => -v.nodeScore)

            newOpts.headOption  match {
              case Some(n) => findWobbleSubCut(n)
              case None => node
            }
            }
        }

        // Try a series of edits to the motif, greedily accepting any that improve the score, and iterating to a fix-point
          // Promote -- change an existing cut base into a more general 'wobble' base
          // Drop -- remove an existing cut
        def refineNode(startNode: SearchNode) : SearchNode =
        {
          var node = findWobbleSubCut(startNode)
          var cuts = node.cuts

          var alternateMotifs = new mutable.ArrayBuffer[Set[Cut]]()
          for (c <- cuts)
          {
            var cDrop = cuts - c

            // Don't drop the last cut base!
            if(cDrop.size > 0)
              alternateMotifs += cDrop

            for(p <- promotions(c.selBase))
            {
              alternateMotifs += (cDrop + new Cut(c.pos, p))
            }
          }

          var alternateNodes = alternateMotifs map (m => scoreCut(genomeInfo, allExamples, m))
          var alternate = alternateNodes find (n => n.nodeScore > node.nodeScore)

          alternate match {
            case Some(alt) => refineNode(alt)
            case None => node
          }
        }

        status("Motif Searching...")

        // Start with all GenomePositions and no cuts.  Currently 16 seems to be a good search breadth
        var greedy = searcher(searchBreadthStarting, searchBreathEnding, depthToSwitching)(emptyNode, emptyNode, allExamples, allPointers)

//        System.out.println("greedy:" + greedy)

        status("Motif Refinement...")
        if(greedy.cuts.size > 0)
        {
          refineNode(greedy)
        }
        else
        {
            greedy
        }
    }

    // Find the best motif, remove the matching examples from the list of hits, then search again
    def loopSearch(genomeInfo : GenomeInfo, examples : Array[Modification]) : List[SearchNode] =
    {

      // Construct the complete set of genomePointers
      def makeAllPointers(genomeInfo : GenomeInfo) : PointerSet =
      {
        genomeInfo.contigs.zipWithIndex map (c => {

          var (contigs, idx) = c
          var (fwd, rev) = contigs

          var fHits = new mutable.ArrayBuffer[GenomePointer]()
          var rHits = new mutable.ArrayBuffer[GenomePointer]()

          var first = 0
          var last = fwd.length

          for (pos <- first to last)
          {
            fHits.append( new GenomePointer(idx, pos, 0) )
            rHits.append( new GenomePointer(idx, pos, 1) )
          }

          (fHits.toArray, rHits.toArray)
        })
      }


      val allPointers = makeAllPointers(genomeInfo)

      def _go(examples : Array[Modification], motifAcc : List[SearchNode]) : List[SearchNode] =
      {
//        println("examples.size: " + examples.size)
        println("currentTime；" + new Date())
        if (motifAcc.length > 12 && motifAcc.head.fractionHit < 0.25)
        {
          status("Reached maximum number of motifs")
          motifAcc.reverse
        }
        else if (examples.length < 4)
        {
          status("Not enough modifications left to cluster")
          motifAcc.reverse
        }
        else
        {
          var r = searchAndRefine(genomeInfo, examples, allPointers, false)
          r match {
            // bail if we get an 'empty' cut
            case m if m.cuts.size == 0 => motifAcc.reverse

            // stop searching last motif has less than 200 hits and less than 12.5% hit
            case m if m.nMods < 500 && m.fractionHit < 0.15 => motifAcc.reverse

            // stop searching when we have a very low fraction motif
            case m if m.nMods < 1000 && m.fractionHit < 0.015 => motifAcc.reverse

            // found a good motif - continue
            case m =>
                status("Found motif: " + Motif.cutString(m.cuts))
                println("motifFoundTime；" + new Date())
                var remainingExamples = examples filter (ex => ! hasCuts(m.cuts)(ex))
                _go(remainingExamples, m :: motifAcc)
          }
        }
      }

      var hits = _go(examples, List())
      hits
    }



}
