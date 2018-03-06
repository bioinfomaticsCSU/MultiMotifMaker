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

package com.csu.basemods

import java.io._


case class GffRow(
	seqId : String,
	source : String,
	
	entryType : String,
	
	startPos : Int,
	endPos : Int,
	
	score : Float,
	strand : Char,
	phase : String,
	attributes : Map[String, String]
)
{
	override def toString() =
	{
	  var attrString = attributes.map({case (k,v) => k + "=" + v}).mkString(";")	  
	  var s = List(seqId, source, entryType, startPos.toString, endPos.toString, score.toInt.toString, strand.toString, phase, attrString).mkString("\t")
	  s
	}
	
	def toModification(genomeInfo : GenomeInfo) =
	{	  
	  var gp = genomePointer(genomeInfo)
	  
	  var ctxString = attributes("context")
	  var ctxArray = (ctxString map MotifMixture.charToCode).toArray
	  
	  Modification(score, ctxArray, attributes("coverage").toInt, attributes("IPDRatio").toFloat, entryType, seqId, gp)
	}
	
	def genomePointer(genomeInfo : GenomeInfo) =
	{
	  var seqNumber = genomeInfo.nameToId(seqId)
	  var str = if (strand == '+') 0.toByte else 1.toByte
	  GenomePointer(seqNumber, startPos - 1, str)
	}	
}


case class GffData(
    headers : IndexedSeq[String],
    rows : IndexedSeq[GffRow]
    )
{
  def writeToFile(f : String)
  {
    val writer = new PrintWriter(new File(f))
    
    for(h <- headers)
    	writer.println(h)
    	
    for(r <- rows)
      writer.println(r.toString())

    writer.close()
  }
}


object Gff {
  
  def parseGffRow(s : String) =
  {
    var cols = s.split("\t")
    
    var pairs = cols(8).split(";") map (i => { var r = i.split("="); (r(0).trim(), r(1))}) toMap
    
    new GffRow(
        seqId = cols(0).split(" ")(0),
        source = cols(1),
        entryType = cols(2),
        
        startPos = cols(3).toInt,
        endPos = cols(4).toInt,
        
        score = cols(5).toInt,
        strand = cols(6)(0),
        phase = cols(7),
        attributes = pairs        
    )
  }
  
  def loadGff(fn : String) =
  {
    var headers = Vector[String]()
    var rows = Vector[GffRow]()
    
    def parseLine(l : String) =
    {
      if(l(0) == '#')
        headers = headers :+ l
      else
        rows = rows :+ parseGffRow(l)
    }
    
    for (l <- Reader.read(fn))
      parseLine(l)
    
    GffData(headers = headers, rows = rows)
  }
  
}
