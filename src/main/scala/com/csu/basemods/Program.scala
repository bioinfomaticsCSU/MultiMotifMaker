// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

import java.io._
import java.util.Date

import au.com.bytecode.opencsv.{CSVReader, CSVWriter => JCSVWriter}
import com.beust.jcommander.{JCommander, Parameter, Parameters}
import com.csu.basemods.MotifMixture._

import scala.collection.JavaConversions._
import scala.collection.mutable
import scala.io.Source
import scala.util.{Failure, Success, Try}
import scala.xml.XML
import scopt.OptionParser._

abstract class MergeCase[A,B]
case class Left[A,B](l: A) extends MergeCase[A,B]
case class Right[A,B](r: B) extends MergeCase[A,B]
case class Both[A,B](l: A, r: List[B]) extends MergeCase[A,B]

class HetMerge[A,B,C](s1: Seq[A], f1: A=>C, s2: Seq[B], f2: B=>C)(implicit ord: Ordering[C]) extends Iterator[MergeCase[A,B]]
{

  val e1 = s1.iterator
  var e2 = s2.iterator

  var vLeft: Option[A] = None
  var vRight: Option[B] = None

  def incrLeft() = vLeft = if (e1.hasNext) Some(e1.next()) else None
  def incrRight() = vRight = if(e2.hasNext) Some(e2.next()) else None

  var started = false

  var nextElement: Option[MergeCase[A,B]] = None

  def hasNext: Boolean =
  {
    if(!started)
    {
      nextElement = getNext()
    }

    nextElement match { case Some(_) => true; case None => false }
  }

  def next(): MergeCase[A,B] =
  {
    if(!started)
    {
      nextElement = getNext()
    }

    nextElement match {
      case Some(v) =>
        nextElement = getNext()
        v
      case None => throw new Exception("iterator is finished")
      }
  }


  def getNext(): Option[MergeCase[A,B]] =
  {
    if(!started)
    {
      incrLeft()
      incrRight()
      started = true
    }

    (vLeft, vRight) match
    {
      case (None, None) => None
      case (Some(i), None) => incrLeft(); Some(Left(i))
      case (None, Some(i)) => incrRight(); Some(Right(i))
      case (Some(i), Some(j)) =>

          val c = ord.compare(f1(i), f2(j))
          if(c < 0)
          {
            incrLeft()
            Some(Left(i))
          }
          else if (c > 0)
          {
            incrRight()
            Some(Right(j))
          }
          else
          {
            var lv = f1(i)

            def loop(j: List[B]): List[B] =
            {
              incrRight()
              vRight match {
                case Some(v) if ord.compare(lv, f2(v)) == 0 => loop(v :: j)
                case _ => j
              }
            }

            incrLeft()
            val matching = loop(List(j))
            Some(Both(i, matching))
          }
    }
  }
}


class CSVWriter(writer: Writer) {

  def this(file: File, enc: String = "UTF-8") =
    this(new OutputStreamWriter(new FileOutputStream(file), enc))

  val csvWriter = new JCSVWriter(writer)

  def close() = csvWriter.close()

  def writeAll(allLines: List[List[String]]) =
    csvWriter.writeAll(allLines.map(_.toArray))

  // XXX this is a little gross - in opencsv 2.4 (the only version in the maven
  // repo), the call signature is writeNext(String... nextLine)
//  def writeRow(fields: List[String]) = csvWriter.writeNext(fields:_*)
  def writeRow(fields: List[String]) = csvWriter.writeNext(fields:_*)

}


object MotifUtils
{
  import Constants._

  // Reconstruct the context snippet string of a genome pointer
  def getContextSnippets(genomeInfo: GenomeInfo, genomePointer: GenomePointer): String =
  {
    def getStrand(s: Byte, ctg: Tuple2[String, String]) = if (s == 0) ctg._1 else ctg._2

    val ctg = getStrand(genomePointer.strand, genomeInfo.contigs(genomePointer.refId))

    val center = if (genomePointer.strand == 0) genomePointer.site else ctg.length - 1 - genomePointer.site

    (for (idx <- 0 to CONTEXT_SIZE) yield {
        val gIdx = center - CONTEXT_OFFSET + idx
        if ( gIdx < 0 || gIdx >= ctg.length) 'N' else ctg(gIdx)
      }).toList.mkString
  }


  def loadMotifCsv(fn: String, minFraction: Double) =
  {
    var input = new FileInputStream(fn)
    var lines = Source.fromInputStream(input).getLines().toArray

    var headers = lines(0).split(",").map(s => s.replaceAll("\"", ""))
    var headerMap = headers.zipWithIndex toMap

    def readLine(s : String) =
    {
      var cols = s.split(",")
      def gc(name: String) = cols(headerMap(name)).replaceAll("\"", "")

      var groupTag =
        headerMap.get("groupTag").map(i => cols(i)) match
        {
          case Some(s) if s.length > 0 => Some(s.replaceAll("\"", ""))
          case _ => None
        }

      var partnerMotifString =
        headerMap.get("partnerMotifString").map(i => cols(i)) match
        {
          case Some(s) if s.length > 0 => Some(s.replaceAll("\"", ""))
          case _ => None
        }

      var modificationType =
        headerMap.get("modificationType").map(i => cols(i)) match
        {
          case Some(s) if s.length > 0 => s.replaceAll("\"", "")
          case _ => "modified_base"
        }

      MotifSummary(
          gc("objectiveScore").toFloat,
          gc("nDetected").toInt,
          gc("nGenome").toInt,
          gc("fraction").toFloat,
          gc("meanScore").toFloat,
          gc("meanIpdRatio").toFloat,
          gc("meanCoverage").toFloat,
          gc("motifString"),
          gc("centerPos").toInt,
          modificationType,
          groupTag,
          partnerMotifString
      )
    }

    lines.drop(1).map(readLine).filter(v => v.fraction >= minFraction)
  }


  def loadRawModsData(genomeInfo: GenomeInfo, fn: String) = {
    var reader = new CSVReader(Reader.getReader(fn), ',', '"', '\0');
    // Read first line to get csv headers
    var nextLine = reader.readNext()
    var hm = nextLine.zipWithIndex toMap
    var dictSet = Array.fill(genomeInfo.contigs.length, 2)(new mutable.HashMap[Int, KineticSummary]())
    def grabRow(items: Array[String]) = {
      var refName = items(hm("refName"))
      var refId = genomeInfo.nameToId(refName)
      //var refId = items(hm("refId")).substring(3).toInt - 1
      var strand = items(hm("strand")).toInt
      var tpl = items(hm("tpl")).toInt - 1
      var summary = KineticSummary(
        items(hm("ipdRatio")).toFloat,
        Math.min(32000, items(hm("coverage")).toInt).toShort,
        items(hm("score")).toShort)
      dictSet(refId)(strand)(tpl) = summary
    }
    println("Loading raw modification data")
    nextLine = reader.readNext()
    // Read remaining rows and store raw kinetics
    while(nextLine != null){
      grabRow(nextLine)
      nextLine = reader.readNext()
    }
    println("Done loading raw mods")
    def find(gp: GenomePointer) = {
      var d = dictSet(gp.refId)(gp.strand)
      d.get(gp.site)
    }
    find _
  }

  // alternative to CSV parsing - just extract the necessary data from the GFF
  def loadRawModsDataGff(genomeInfo: GenomeInfo, gff: GffData) = {
    var dictSet = Array.fill(genomeInfo.contigs.length, 2)(new mutable.HashMap[Int, KineticSummary]())
    for (row <- gff.rows) {
      var refId = genomeInfo.nameToId(row.seqId)
      var strand = if (row.strand == '+') 0 else 1
      var tpl = row.startPos
      var summary = KineticSummary(
        row.attributes("ipdRatio").toFloat,
        Math.min(32000, row.attributes("coverage").toInt).toShort,
        row.attributes("score").toShort)
      dictSet(refId)(strand)(tpl) = summary
    }
    println("Done extracting raw kinetics data")
    def find(gp: GenomePointer) = {
      var d = dictSet(gp.refId)(gp.strand)
      d.get(gp.site)
    }
    find _
  }


  // Make the 'master file' of modification hits and any associated motif data
  // Pull in the raw kinetics from the csv for sites that match a motif but do not appear in the gff file
  def makeCombinedMotifsAndModsGff(genomeInfo: GenomeInfo, motifs: Array[MotifSummary], gffData: GffData, kineticLookup: GenomePointer => Option[KineticSummary], outGff: String) =
  {
    def getCutSet(motif: MotifSummary ) =
    {
      motif.motifString.zipWithIndex map { case (c, idx) => Cut(pos = idx - motif.centerPos, selBase = c) } filter { _.selBase != 'N'} toSet
    }

    var motifsAndCuts = motifs map (m => (m, getCutSet(m)))

    // Find all the hits of all the motifs, sort by genomePointer
    println("Finding motif hits...")
    var motifInstances = motifsAndCuts map { case (m, c) => (m, findMotifInstances(genomeInfo, c)) }
    var motifHits = motifInstances.flatMap { case (m, hits) => (hits map (h => (h, m))) }
    motifHits = motifHits.sortBy(c => c._1)(PointerCompare)

    println("Getting modifications...")
    // Load all the gffRows (and a converted version), sort by genomePointer
    var modRows = (gffData.rows map (r => (r, r.toModification(genomeInfo)))).toArray
    modRows = modRows.sortBy({ case (r, mod) => mod.genomePointer })(PointerCompare)

    println("Merging...")
      // Merge the sorted mods and motif hits according to their genomePointer
      var merged = new HetMerge[(GffRow, Modification), (GenomePointer, MotifSummary), GenomePointer](
          modRows, { case (_,m) => m.genomePointer },
          motifHits, { _._1 }) (PointerCompare)


      var contigNames = genomeInfo.nameToId.map { case(k,v) => (v,k )}

        // Convert the merged modification / motif hit data into the final gff rows
      def createFinalRow(r: MergeCase[(GffRow, Modification), (GenomePointer, MotifSummary)]): Option[GffRow] =
      {
      r match {

        // Matches a motif but was not detected - pull in the raw single site kinetics
        case Right((gp, motif)) =>
          kineticLookup(gp) match
          {
            case Some(k) =>

              var ctx = MotifUtils.getContextSnippets(genomeInfo, gp)

              Some(GffRow(
                  seqId = contigNames(gp.refId),
                  source = "kinModCall",
                  entryType = ".",
                  startPos = gp.site + 1,
                  endPos = gp.site + 1,
                  strand = if(gp.strand == 0) '+' else '-',
                  score = k.score.toFloat,
                  phase = ".",
                  attributes =
                    Map("IPDRatio" -> k.ipdRatio.toString,
                        "context" -> ctx,
                        "coverage" -> k.coverage.toString,
                        "motif" -> motif.motifString,
                        "id" -> motif.groupTag.getOrElse(""))))
            case None => None
          }

          // Detected & not a motif - goes through unchanged
        case Left((gffRow, modification)) => Some(gffRow)

        // Detected & matches motif
        case Both((gffRow, modification), motifList) =>

          var motifStringArray = motifList map { case (_, motif) => motif.motifString }
          var motifString = motifStringArray.mkString(",")

          var idStringArray = motifList map { case (_, motif) => motif.groupTag } flatten
          var idString = idStringArray.mkString(",")

          var attrs = gffRow.attributes + ("motif" -> motifString) + ("id" -> idString)

          var editRow = gffRow.copy(attributes = attrs)
          Some(editRow)
      }
      }

    var rows = (merged map createFinalRow).flatten.toIndexedSeq
    var finalGff = GffData(gffData.headers, rows)

    println("Writing output file: " + outGff)
    // Write the final result
    finalGff.writeToFile(outGff)
  }



  def makeMotifsOnlyGff(genomeInfo: GenomeInfo, motifs: Array[MotifSummary], gffData: GffData, kineticLookup: GenomePointer => Option[KineticSummary], outGff: String) =
  {
    def getCutSet( motif: MotifSummary ) =
    {
      motif.motifString.zipWithIndex map { case (c, idx) => Cut(pos = idx - motif.centerPos, selBase = c) } filter { _.selBase != 'N'} toSet
    }

    var motifsAndCuts = motifs map (m => (m, getCutSet(m)))

    // Find all the hits of all the motifs, sort by genomePointer
    println("Finding motif hits...")
    var motifInstances = motifsAndCuts map { case (m, c) => (m, findMotifInstances(genomeInfo, c)) }
    var motifHits = motifInstances.flatMap { case (m, hits) => (hits map (h => (h, m))) }

    // Don't sort -- leave in motif order
    motifHits = motifHits.sortBy(c => c._1)(PointerCompare)

    // We don't even look at the original gff now
    //var modRows = (gffData.rows map (r => (r, r.toModification(genomeInfo)))).toArray

      // Convert the merged modification / motif hit data into the final gff rows
      def createFinalRow(r: (GenomePointer, MotifSummary)): Option[GffRow] =
      {
      var (gp, motif) = r;

      kineticLookup(gp) match
      {
        case Some(k) =>
          var modType = modMap(genomeInfo.get(gp))
          var attrs =  Map(
                  "IPDRatio" -> k.ipdRatio.toString,
                  "coverage" -> k.coverage.toString,
                  "motif" -> motif.motifString,
                  "modification" -> modType
              )

          attrs = motif.groupTag match {
            case Some(s) => attrs + ("id" -> s)
            case _ => attrs
          }

          Some(GffRow(
              seqId = "ref%06d" format (gp.refId + 1),
              source = "MotifMaker",
              entryType = "DNA_motif",
              startPos = gp.site + 1,
              endPos = gp.site + 1,
              strand = if(gp.strand == 0) '+' else '-',
              score = k.score.toFloat,
              phase = ".",
              attributes = attrs))

        case None => None
      }
      }

    var rows = (motifHits map createFinalRow).flatten.toIndexedSeq
    var finalGff = GffData(gffData.headers, rows)

    println("Writing output file: " + outGff)
    // Write the final result
    finalGff.writeToFile(outGff)
  }


    // sample function for converting from generic modified base to a specific type
  def modMap(base: Char) =
  {
    base match
    {
      case 'A' => "m6A"
      case 'C' => "m4C"
      case _ => "modified_base"
    }
  }

   // sample function for converting from generic modified base to a specific type
  def newModType(gffRow: GffRow, mod: Modification) =
  {
    codeToChar(mod.context(CONTEXT_OFFSET)) match
    {
      case 'A' => "methylated_A"
      case 'C' => "methylated_C"
      case _ => gffRow.entryType
    }
  }

}

case class KineticSummary(
    ipdRatio: Float,
    coverage: Short,
    score: Short
)

object FloatExt {
  implicit def rFloat(f: Float) = new {
    def fmt(d: Int) = ("%." + d + "f").format(f)
  }
}

import com.csu.basemods.FloatExt._

case class MotifSummary(
    objectiveScore: Float,
    nDetected: Int,
    nGenome: Int,
    fraction: Float,
    meanScore: Float,
    meanIpdRatio: Float,
    meanCoverage: Float,
    motifString: String,
    centerPos: Int,
    modificationType: String,
    groupTag: Option[String] = None,
    partnerMotifString: Option[String] = None
)
 {

  def isReverseComplement(other: MotifSummary): Boolean =
  {
    return Reader.reverseComplement(this.motifString) == other.motifString
  }

  def toXml = {

    var motifEntry =
      if (motifString == "Not Clustered")
      {
        <td><i>{motifString}</i></td>
      }
      else
      {
        if (centerPos < 0 || centerPos >= motifString.length())
        {
          <td>{motifString}</td>
        }
        else
        {
          var pre = motifString.substring(0, centerPos)
          var cb = motifString.substring(centerPos, centerPos+1)
          var post = motifString.substring(centerPos + 1, motifString.length())
          <td>{pre}<b>{cb}</b>{post}</td>
        }
      }

    var modString = if(modificationType == "modified_base") "unknown" else modificationType

    <tr>
      {motifEntry}
      <td>{centerPos+1}</td>
      <td>{modString}</td>
      <td>{(fraction*100).fmt(2)}</td>
      <td>{nDetected}</td>
      <td>{nGenome}</td>
      <td>{meanScore.fmt(1)}</td>
      <td>{meanCoverage.fmt(1)}</td>
      <td>{partnerMotifString.getOrElse("")}</td>
  </tr>
  }


}

object MotifSummary {
    def xmlHeader = {
        <tr>
      <th>Motif</th>
      <th>Modified Position</th>
        <th>Modification Type</th>
      <th>% Motifs Detected</th>
      <th># of Motifs Detected</th>
      <th># of Motifs in Genome</th>
      <th>Mean Modification QV</th>
      <th>Mean Motif Coverage</th>
      <th>Partner Motif</th>
    </tr>
  }

  def makeXmlTable(motifs: Seq[MotifSummary], unclustered: MotifSummary, filename: String) =
  {
    var tbl =
        <report>
      <attributes>
      </attributes>
          <title>Motifs</title>

          <table>
            <thead>{xmlHeader}</thead>
            <tbody>
              {motifs.map(m => m.toXml)}
              {unclustered.toXml}
            </tbody>
          </table>
      <graphGroup>
          <title>Modification QV Histogram By Motif</title>
          <graph>
            <title>Modification QV Histogram</title>
            <image>motifHistogram.png</image>
            <caption>Motif and base modification type identification is dependent on sequence coverage, set threshold levels, degree of DNA methylation, and other factors. Results should be interpreted in consideration of these factors.</caption>
          </graph>
      </graphGroup>
        </report>

    XML.save(filename, tbl, "UTF-8", true, null)
  }
}


object Program extends App {


  def loadModificationsGff(gffFile:String, minScore: Float, genome: GenomeInfo) =
  {
    var mods = (Reader.loadModificationsFromGff(gffFile, genome) filter { _.hitScore > minScore }).toArray
    mods
  }


  // Go back over the genome and count motif instances.
  // Summarize the kinetic evidence in the matching modification hits.
  def summarizeEvidence(genomeInfo: GenomeInfo, allMods: Array[Modification], node: SearchNode) =
  {
    var examples = allMods filter (hasCuts(node.cuts) _)

    var exampleCount = examples.groupBy(m => m.modType).mapValues(c => c.length)
    var (bestModificationType,_) = exampleCount.maxBy { case (x,count) => count }

    var startPos = (node.cuts.minBy { _.pos }).pos

    var motifSeq = Motif.cutString(node.cuts)

    var fraction = examples.length.toFloat / node.nGenome.toFloat

    var meanScore = examples.map(e => e.hitScore).sum / examples.length
    var meanCoverage = examples.map(e => e.coverage.toFloat).sum / examples.length
    var meanIpdRatio = examples.map(e => e.ipdRatio).sum / examples.length

    new MotifSummary(
        node.nodeScore,
        examples.length,
        node.nGenome,
        fraction,
        meanScore,
        meanIpdRatio,
        meanCoverage,
        motifSeq,
        -startPos, bestModificationType, None, None)
  }



abstract class MotifGroup[A]
case class Lone[A](v1: A) extends MotifGroup[A]
case class Self[A](v1: A) extends MotifGroup[A]
case class Pair[A](v1: A, v2: A) extends MotifGroup[A]

  /*
   * Group motifs that are mutually or self reverse-complementary
   * Name each grouping as "<motifString1>/<motifString2>"
   * and add this name to the MotifSummary
   */
  def groupMotifs(motifs: List[MotifSummary]) =
  {
    var groups = IndexedSeq[MotifGroup[MotifSummary]]()
    var leftMotifs = motifs.toSet

    while(leftMotifs.size > 0)
    {
      var test = leftMotifs.maxBy(motif => motif.fraction)
      leftMotifs = leftMotifs - test

      // Check if we are self-complementary
      if (test.isReverseComplement(test))
      {
        groups = groups :+ Self(test)
      }
      else
      {
      // Scan the remaining motifs for rc motifs
      leftMotifs.find(other => test.isReverseComplement(other)) match
      {
        case Some(other) =>
          leftMotifs = leftMotifs - other
          groups = groups :+ Pair(test, other)

        case None =>
          groups = groups :+ Lone(test)
      }
      }
    }

    // Emit the motifs in grouped order
    // Enrich each motif hit with a groupId and possibly a paired motif string
    groups.flatMap(g =>
      g match {
      case Lone(m) => List(m.copy(groupTag=Some(m.motifString), partnerMotifString=None))
      case Self(m) => List(m.copy(groupTag=Some(m.motifString), partnerMotifString=Some(m.motifString)))
      case Pair(m1, m2) =>
          var groupName = Some(m1.motifString + "/" + m2.motifString)
          List( m1.copy(groupTag=groupName, partnerMotifString=Some(m2.motifString)),
                m2.copy(groupTag=groupName, partnerMotifString=Some(m1.motifString)))
    })
  }

  def computeUnclustered(allMods: Array[Modification], genomeInfo: GenomeInfo, motifsFound: List[SearchNode]) =
  {

    // Count mods that didn't make it into a motif:
    var cuts = motifsFound map { _.cuts }

    var examples = allMods.filter(m => cuts.forall(c => !hasCuts(c)(m)))
    var nonClustered = examples.length

    var nGenome = genomeInfo.size - motifsFound.map(m => m.nGenome).sum
    var fraction = nonClustered.toFloat / nGenome

    var meanScore = examples.map(e => e.hitScore).sum / examples.length
    var meanCoverage = examples.map(e => e.coverage.toFloat).sum / examples.length
    var meanIpdRatio = examples.map(e => e.ipdRatio).sum / examples.length

    MotifSummary(-1, nonClustered, nGenome, fraction,
        meanScore,
    meanIpdRatio,
    meanCoverage,
    "Not Clustered",
    -1,
    "",
    None,
    None)
  }


  def findMotifs(allMods: Array[Modification], genomeInfo: GenomeInfo) =
  {
    var motifsFound = loopSearch(genomeInfo, allMods)

    var goodMotifs = motifsFound.filter(m => m.cuts.exists(c => c.pos == 0 && MotifMixture.canonicalBases.contains(c.selBase)))
    //summaries = summaries.filter(m => m.centerPos >= 0 && m.centerPos < m.motifString.length() && MotifMixture.canonicalBases.contains(m.motifString(m.centerPos)))

    // Remove motifs that don't have their center base on a canonical base
    // These are false positives by definition.
    var summaries = goodMotifs map (m => summarizeEvidence(genomeInfo, allMods, m))

    var unclustered = computeUnclustered(allMods, genomeInfo, goodMotifs)
    (summaries.sortBy(c => - c.objectiveScore), unclustered)
  }

  def writeMotifCsv(fn: File, motifs: List[MotifSummary]) =
    {
      var w = new CSVWriter(fn, "UTF-8")

      def ts(s : Option[String]) = s.getOrElse("")

      w.writeRow(
        List(
          "motifString",
          "centerPos",
          "modificationType",

          "fraction",
          "nDetected",
          "nGenome",

          "groupTag",
          "partnerMotifString",

          "meanScore",
          "meanIpdRatio",
          "meanCoverage",

          "objectiveScore"
          ))

      for (s <- motifs) {
        val ln = List[String](

          s.motifString.toString,
          s.centerPos.toString,
          s.modificationType,

          s.fraction.toString,
          s.nDetected.toString,
          s.nGenome.toString,

          ts(s.groupTag),
          ts(s.partnerMotifString),

          s.meanScore.toString,
          s.meanIpdRatio.toString,
          s.meanCoverage.toString,

          s.objectiveScore.toString
          )

        w.writeRow(ln)
      }

      w.close()
    }


  def runMotifReport(gff: String, fasta: String, output: String, minScore: Float, xmlOut: String) =
  {
    var genome = Reader.makeGenomeInfo(fasta)
    var mods = loadModificationsGff(gff, minScore, genome)
    var (motifs, unclustered) = findMotifs(mods, genome)
    // Group motifs into partners
    motifs = groupMotifs(motifs).toList
    println("Got results:")
    for(m <- motifs)
      println(m)
    writeMotifCsv(new File(output), motifs.sortBy(m => - m.objectiveScore))
    // Optionally write output xml
    xmlOut match {
      case null => ()
      case "" => ()
      case fn =>  MotifSummary.makeXmlTable(motifs, unclustered, fn)
    }
  }


  def runGffReprocess() = {
    var genomeInfo = Reader.makeGenomeInfo(ReprocessCommand.fasta)
    var gff = Gff.loadGff(ReprocessCommand.gff)
    def noRawKinetics(gp: GenomePointer): Option[KineticSummary] = None
    var rawKinData = ReprocessCommand.csv match {
      case null => Try {
        MotifUtils.loadRawModsDataGff(genomeInfo, gff)
      } match {
        case Success(d) => d
        case Failure(_) => noRawKinetics _
      }
      case s => MotifUtils.loadRawModsData(genomeInfo, s)
    }
    println("Loading motifs: " + ReprocessCommand.motifs)
    var motifs = MotifUtils.loadMotifCsv(ReprocessCommand.motifs, ReprocessCommand.minFraction)
    MotifUtils.makeCombinedMotifsAndModsGff(genomeInfo, motifs, gff, rawKinData, ReprocessCommand.output)
  }



  // Members declared as var because JCommander assigns a new collection declared
  // as java.util.List because that's what JCommander will replace it with.
  // It'd be nice if JCommander would just use the provided List so this
  // could be a val and a Scala LinkedList.

  object Args {
    @Parameter(names = Array("-h", "--help"), help=true)
    var help: Boolean = false
  }

  // Parameters for motif finding
  @Parameters(commandDescription="Run motif finding")
  object FindCommand
  {

    @Parameter(names = Array("-m", "--minScore"), description="Minimum Qmod score to use in motif finding")
    var minScore: Double = 30

    @Parameter(names = Array("-g", "--gff"), description="modifications.gff or .gff.gz file", required=true)
    var gff: String = null

    @Parameter(names = Array("-f", "--fasta"), description="Reference fasta file", required=true)
    var fasta: String = null

    @Parameter(names = Array("-o", "--output"), description="Output motifs csv file", required=true)
    var csv: String = null

    @Parameter(names = Array("-x", "--xml"), description="Output motifs xml file", required=false)
    var xml: String = null

    @Parameter(names = Array("-p", "--parallelize"), description="Parallelize motif finder", required=false)
    var parallelize: Boolean = true

    @Parameter(names = Array("-l", "--layer"), description="Search depth used to Parallelize motif finder", required=false)
    var layer: Int = 1

    @Parameter(names = Array("-s", "--searchBreadthStart"), description="searchBreadthStart for motif finder", required=false)
    var searchBreadthStart: Int = 16

    @Parameter(names = Array("-t", "--thread"), description="The concurrency to parallelize motif finder", required=false)
    var poolThreads: Int = searchBreadthStart

    @Parameter(names = Array("-e", "--searchBreadthEnd"), description="searchBreadthEnd for motif finder", required=false)
    var searchBreadthEnd: Int = 6

    @Parameter(names = Array("-d", "--depthToSwitch"), description="depthToSwitch for motif finder", required=false)
    var depthToSwitch: Int = 3
  }

  // Parameters for code the annotates the gff with motif information.
  @Parameters(commandDescription="Reprocess gff file with motif information")
  object ReprocessCommand
  {

    @Parameter(names = Array("--minFraction"), description="Only use motifs above this methylated fraction")
    var minFraction: Double = 0.75

    @Parameter(names = Array("-m", "--motifs"), description="motifs csv", required=true)
    var motifs: String = null

    @Parameter(names = Array("-g", "--gff"), description="original modifications.gff or .gff.gz file", required=true)
    var gff: String = null

    @Parameter(names = Array("-f", "--fasta"), description="Reference fasta file", required=true)
    var fasta: String = null

    @Parameter(names = Array("-c", "--csv"), description="Raw modifications.csv file", required=false)
    var csv: String = null

    @Parameter(names = Array("-o", "--output"), description="Reprocessed modifications.gff file", required=true)
    var output: String = null

  }

  // Main function -- use JCommander (http://www.jcommander.org) for sub-command style argument parsing
  override def main(args: Array[String]) {
    println("startTime；" + new Date())

    // Some magic to call a particular constructor of JCommander, they appear ambiguous to Scala.
    // For details see http://stackoverflow.com/questions/3313929/how-do-i-disambiguate-in-scala-between-methods-with-vararg-and-without
    val mkAmbig2 = classOf[JCommander].getConstructors.filter(_.getParameterTypes.length==1)
    var jc = mkAmbig2.head.newInstance(Args).asInstanceOf[JCommander]

    // Register the sub-commands
    jc.addCommand("find", FindCommand)
    jc.addCommand("reprocess", ReprocessCommand)
    jc.setProgramName("MultiMotifMaker")

    // Turn an array into a java varags
    jc.parse(args.toArray : _*)

    if(Args.help || jc.getParsedCommand() == null)
    {
      jc.usage()
    }
    else
    {
      try {
        jc.getParsedCommand() match
        {
          // Run motif finding
          case "find" =>
            if(FindCommand.layer == 1 && FindCommand.poolThreads > FindCommand.searchBreadthStart)
              Constants.CONTEXT_POOLTHREADS = FindCommand.searchBreadthStart
            else
              Constants.CONTEXT_POOLTHREADS = FindCommand.poolThreads
            MotifMixture.parallelize = FindCommand.parallelize
            MotifMixture.layer = FindCommand.layer
            MotifMixture.searchBreadthStarting = FindCommand.searchBreadthStart
            MotifMixture.searchBreathEnding = FindCommand.searchBreadthEnd
            MotifMixture.depthToSwitching = FindCommand.depthToSwitch
//            println("layer:"+MotifMixture.layer)
//            println("poolThreads:" + Constants.CONTEXT_POOLTHREADS)
            runMotifReport(FindCommand.gff, FindCommand.fasta, FindCommand.csv, FindCommand.minScore.toFloat, FindCommand.xml)

          // Annotate the modifications.gff with motif tags
          case "reprocess" => runGffReprocess()
        }
      }
      finally
      {
        // Shut down our thread pool
        MotifMixture.threadPool.shutdown()
        MotifMixture.branchThreadPool.shutdown()
      }
    }
    println("endTime；" + new Date())
  }
}
