/* 
 * Copyright (c) 2015-2016 TU Delft, The Netherlands.
 * All rights reserved.
 * 
 * You can redistribute this file and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Authors: Hamid Mushtaq
 *
*/
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.log4j.Logger
import org.apache.log4j.Level

import sys.process._
import scala.sys.process.Process

import java.io._

import tudelft.utils.ChromosomeRange
import tudelft.utils.DictParser
import tudelft.utils.Configuration
import tudelft.utils.SAMRecordIterator

import htsjdk.samtools._

object DNASeqAnalyzer 
{
final val MemString = "-Xmx5120m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"
//////////////////////////////////////////////////////////////////////////////
def bwaRun (x: String, config: Configuration) : 
	Array[(Int, SAMRecord)] = 
{
	val refFolder = config.getRefFolder
	val tmpFolder = config.getTmpFolder
		
	// Create the command string (bwa mem...)and then execute it using the Scala's process package. More help about 
	//	Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package. 
	
	//bwa mem refFolder/RefFileName -p -t numOfThreads fastqChunk > outFileName
	val outFileName = config.getTmpFolder + x + ".output"
	val cmd = Seq(config.getToolsFolder+"bwa","mem",refFolder+RefFileName,"-p","-t",config.getNumInstances(),config.getInputFolder+x) #> new File(outFileName) !
//	val cmd = Seq(config.getToolsFolder+"bwa", "mem", "/data/spark/ref/ucsc.hg19.fasta","-p", "-t 4", "/home/gsouto/chunks/0.fq.gz"," home/gsouto/chunks/test")
//	cmd.lines

	val bwaKeyValues = new BWAKeyValues(outFileName)
	bwaKeyValues.parseSam()
	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()
	
	// Delete the temporary files
	
	return kvPairs // Replace this with return kvPairs
}
	 
def writeToBAM(fileName: String, samRecordsSorted: Array[SAMRecord], config: Configuration) : ChromosomeRange = 
{
	val header = new SAMFileHeader()
	header.setSequenceDictionary(config.getDict())
	val outHeader = header.clone()
	outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	val factory = new SAMFileWriterFactory();
	val writer = factory.makeBAMWriter(outHeader, true, new File(fileName));
	
	val r = new ChromosomeRange()
	val input = new SAMRecordIterator(samRecordsSorted, header, r)
	while(input.hasNext()) 
	{
		val sam = input.next()
		writer.addAlignment(sam);
	}
	writer.close();
	
	return r
}

def compareSam(s1: SAMRecord,s2: SAMRecord) : Boolean = 
{
	if(s1.getReferenceIndex == s2.getReferenceIndex)
		return s1.getAlignmentStart < s2.getAlignmentStart
	else 
		return s1.getReferenceIndex < s2.getReferenceIndex
}

def variantCall (chrRegion: Int, samRecordsSorted: Array[SAMRecord], config: Configuration) : 
	Array[(Integer, (Integer, String))] = 
{	
	val tmpFolder = config.getTmpFolder
	val toolsFolder = config.getToolsFolder
	val refFolder = config.getRefFolder
	val numOfThreads = config.getNumThreads
	
	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion). 
	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More 
	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
	//	Note that MemString here is -Xmx6144m, and already defined as a constant variable above, and so are reference files' names.
	

	var sortedArray = samRecordsSorted.sortWith((s1: SAMRecord, s2: SAMRecord) => compareSam(s1,s2))

	// SAM records should be sorted by this point
	val chrRange = writeToBAM(tmpFolder + "region" + chrRegion + "-p1.bam", sortedArray, config)
	
	// Picard preprocessing
	//	java MemString -jar toolsFolder/CleanSam.jar INPUT=tmpFolder/regionX-p1.bam OUTPUT=tmpFolder/regionX-p2.bam
	val picard1 = Seq("java",MemString,"-jar",toolsFolder + "CleanSam.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p1.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+"-p2.bam").!
	//	java MemString -jar toolsFolder/MarkDuplicates.jar INPUT=tmpFolder/regionX-p2.bam OUTPUT=tmpFolder/regionX-p3.bam	
	//		METRICS_FILE=tmpFolder/regionX-p3-metrics.txt 
	val picard2 = Seq("java",MemString,"-jar",toolsFolder + "MarkDuplicates.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p2.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+"-p3.bam", "METRICS_FILE="+tmpFolder+"region"+chrRegion+"-p3-metrics.txt").!
	//	java MemString -jar toolsFolder/AddOrReplaceReadGroups.jar INPUT=tmpFolder/regionX-p3.bam OUTPUT=tmpFolder/regionX.bam 
	//		RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1
	val picard3 = Seq("java",MemString,"-jar",toolsFolder + "AddOrReplaceReadGroups.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p3.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+".bam","RGID=GROUP1","RGLB=LIB1","RGPL=ILLUMINA","RGPU=UNIT1","RGSM=SAMPLE1").!
	// 	java MemString -jar toolsFolder/BuildBamIndex.jar INPUT=tmpFolder/regionX.bam
	val picard4 = Seq("java",MemString,"-jar",toolsFolder + "BuildBamIndex.jar","INPUT="+tmpFolder+"region"+chrRegion+".bam").!
	//	delete tmpFolder/regionX-p1.bam, tmpFolder/regionX-p2.bam, tmpFolder/regionX-p3.bam and tmpFolder/regionX-p3-metrics.txt
	val picard5 = Seq("rm",tmpFolder+"region"+chrRegion+"-p1.bam",tmpFolder+"region"+chrRegion+"-p2.bam",tmpFolder+"region"+chrRegion+"-p3.bam",tmpFolder+"region"+chrRegion+"-p3-metrics.txt").!


	// Make region file 
	val tmpBed = new File(tmpFolder+"tmp"+chrRegion+".bed")
	chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
	//	toolsFolder/bedtools intersect -a refFolder/ExomeFileName -b tmpFolder/tmpX.bed -header > tmpFolder/bedX.bed
	val reg1 = Seq(toolsFolder+"bedtools","intersect","-a",refFolder+ExomeFileName,"-b",tmpFolder+"tmp"+chrRegion+".bed","-header") #> new File(tmpFolder+"bed"+chrRegion+".bed") !
	//	delete tmpFolder/tmpX.bed
	val reg2 = Seq("rm",tmpFolder+"tmp"+chrRegion+".bed").!
	
	// Indel Realignment 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt numOfThreads -R refFolder/RefFileName 
	//		-I tmpFolder/regionX.bam -o tmpFolder/regionX.intervals -L tmpFolder/bedX.bed
	val indR1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","RealignerTargetCreator", "-nt", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+".bam","-o",tmpFolder+"region"+chrRegion+".intervals","-L",tmpFolder+"bed"+chrRegion+".bed").!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T IndelRealigner -R refFolder/RefFileName -I tmpFolder/regionX.bam 
	//		-targetIntervals tmpFolder/regionX.intervals -o tmpFolder/regionX-2.bam -L tmpFolder/bedX.bed
	val indR2 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","IndelRealigner", "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+".bam","-L",tmpFolder+"bed"+chrRegion+".bed", "-targetIntervals", tmpFolder+"region"+chrRegion+".intervals", "-o", tmpFolder+"region"+chrRegion+"-2.bam","-L",tmpFolder+"bed"+chrRegion+".bed").!
	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
	val indR3 = Seq("rm",tmpFolder+"region"+chrRegion+".bam", tmpFolder+"region"+chrRegion+".bai", tmpFolder+"region"+chrRegion+".intervals").!


	// Base quality recalibration 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T BaseRecalibrator -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX.table -L tmpFolder/bedX.bed --disable_auto_index_creation_and_locking_when_reading_rods
	//		-knownSites refFolder/SnpFileName
	val bqr1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","BaseRecalibrator", "-nct", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-2.bam","-o",tmpFolder+"region"+chrRegion+".table","-L",tmpFolder+"bed"+chrRegion+".bed", "--disable_auto_index_creation_and_locking_when_reading_rods", "-knownSites", refFolder+SnpFileName).!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T PrintReads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX-3.bam -BQSR tmpFolder/regionX.table -L tmpFolder/bedX.bed 
	val bqr2 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","PrintReads", "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-2.bam","-o",tmpFolder+"region"+chrRegion+"-3.bam", "-BQSR", tmpFolder+"region"+chrRegion+".table", "-L", tmpFolder+"bed"+chrRegion+".bed").!
	// delete tmpFolder/regionX-2.bam, tmpFolder/regionX-2.bai, tmpFolder/regionX.table
	val bqr3 = Seq("rm",tmpFolder+"region"+chrRegion+"-2.bam", tmpFolder+"region"+chrRegion+"-2.bai", tmpFolder+"region"+chrRegion+".table").!

	// Haplotype -> Uses the region bed file
	// java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T HaplotypeCaller -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-3.bam -o tmpFolder/regionX.vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L tmpFolder/bedX.bed 
	//		--no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods
	val hp1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","HaplotypeCaller", "-nct", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-3.bam","-o",tmpFolder+"region"+chrRegion+".vcf","-stand_call_conf","30.0","-stand_emit_conf","30.0","-L", tmpFolder+"bed"+chrRegion+".bed","--no_cmdline_in_header","--disable_auto_index_creation_and_locking_when_reading_rods").lines
	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
	val hp2 = Seq("rm",tmpFolder+"region"+chrRegion+"-3.bam", tmpFolder+"region"+chrRegion+"-3.bai", tmpFolder+"bed"+chrRegion+".bed").!



	// return the content of the vcf file produced by the haplotype caller.
	//	Return those in the form of <Chromsome number, <Chromosome Position, line>>
	return Array((0,(0,hp1.toString))) // Replace this with what is described in the above 2 lines
}

def main(args: Array[String]) 
{
	val config = new Configuration()
	config.initialize()
		 
	val conf = new SparkConf().setAppName("DNASeqAnalyzer")
	// For local mode, include the following two lines
	conf.setMaster("local[" + config.getNumInstances() + "]")
	conf.set("spark.cores.max", config.getNumInstances())
	//conf.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
	//conf.set("spark.kryo.classesToRegister", "org.apache.hadoop.io.LongWritable,org.apache.hadoop.io.Text,htsjdk.samtools.SAMRecord,import tudelft.utils.SAMRecordIterator")
	
	val sc = new SparkContext(conf)
	sc.broadcast(config)
	
	// Comment these two lines if you want to see more verbose messages from Spark
	Logger.getLogger("org").setLevel(Level.OFF);
	Logger.getLogger("akka").setLevel(Level.OFF);
		
	var t0 = System.currentTimeMillis
	
	// Rest of the code goes here

//	val initialRdd = sc.newAPIHadoopFile(config.getInputFolder(),classOf[TextInputFormat], classOf[LongWritable], classOf[Text])	

//	println(initialRdd.partitions.size)
  	val folder = new File(config.getInputFolder)
    val files = folder.listFiles.filter(_.getName.endsWith(".fq.gz")).toArray

    //val arraynerino: ListBuffer[Array[(Int, SAMRecord)]] = new ListBuffer[Array[(Int, SAMRecord)]]()

    //String
    val initRDD = sc.parallelize(files, config.getNumInstances.toInt)

    //(Chromossome, record)
    //(Int, Sam)
    val chromoRDD = initRDD.mapPartitions(files => {
    	 	for(file <- files)
    			yield bwaRun(file.getName, config)
    	}).flatMap(y => y)


    //(Chromossome, recordCount), sorted in descending order of SamRecords
    //(Int, int)
    val chromoCountRDD = chromoRDD.mapPartitions(
    	lines => 
    	for(line <- lines)
    		yield(line._1, 1)
    	).reduceByKey(_+_).sortBy(_._2, false)

    //after zipWithIndex
    //  ((Chromossome, recordCount), index)
    //  ((int,int), int)
    //after map
    // since the rdd was sorted, if we map and index and get the moduled of it by the number of tasks, we can map every chromossome to a region and it will be distributed
    // 
    // (Chromossome, Region)
    val regionChromoRDD = chromoCountRDD.zipWithIndex().mapPartitions(lines => 
																		for(line <- lines)
																			yield( line._1._1 ,line._2 % config.getNumInstances.toInt) )

    //Joins Chromossome region with its SAMRecords, mapping it to
    //(Region, SamRecord)
    //(Int, Array(SamRecord))
    val regionSamRDD = regionChromoRDD.join(chromoRDD).mapPartitions(lines => for(line <- lines) yield (line._2._1, Array[SAMRecord](line._2._2)))
    

    //Groups every SamRecord with the same region
    val variantCallRDD = regionSamRDD.reduceByKey(_ ++ _).mapPartitions(lines =>for(line <- lines) yield variantCall(line._1.toInt, line._2, config))

    println("begin_________")
  	variantCallRDD.take(10).foreach(println)
    //chromoCountRDD.collect.foreach(println)
   // println(chromoRDD.count)
	
	val et = (System.currentTimeMillis - t0) / 1000 
	println("|Execution time: %d mins %d secs|".format(et/60, et%60))
}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition
