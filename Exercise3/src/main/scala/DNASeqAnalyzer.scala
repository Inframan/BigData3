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

import org.apache.commons.io.FileUtils

import sys.process._
import scala.sys.process.Process

import java.io._
import java.text.DateFormat
import java.text.SimpleDateFormat
import java.util.Calendar

import scala.io.Source

import tudelft.utils.ChromosomeRange
import tudelft.utils.DictParser
import tudelft.utils.Configuration
import tudelft.utils.SAMRecordIterator

import htsjdk.samtools._

import hdsfFile.Operations

object DNASeqAnalyzer 
{
final val MemString = "-Xmx5120m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"
//////////////////////////////////////////////////////////////////////////////

def getTimeStamp() : String =
{
	return "[" + new SimpleDateFormat("HH:mm:ss").format(Calendar.getInstance().getTime()) + "] "
}


def bwaRun (x: String, config: Configuration) : 
	Array[(Int, SAMRecord)] = 
{
	val hdfs = new Operations()
	val refFolder = config.getRefFolder
	val tmpFolder = config.getTmpFolder 
		
	// Create the command string (bwa mem...)and then execute it using the Scala's process package. More help about 
	//	Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package. 

	new File(tmpFolder).mkdirs

	val toolsFolder = tmpFolder + "tools/"


	val outFileName = tmpFolder + x + ".output"

	val logFile = config.getOutputFolder + "log/bwa/log" + x.split("fq")(0) + "txt"


	//create log file
	hdfs.createFile(logFile)
	
	//download input file
	hdfs.downloadFile(config.getInputFolder+x,tmpFolder+x )



	//copy tools folder
	val copyTool = Seq("cp","-r",config.getToolsFolder, tmpFolder)
 	hdfs.writeData(logFile,getTimeStamp() + "\t"+ copyTool.mkString(" ")+"\n")
 	copyTool.!

	
	//bwa mem refFolder/RefFileName -p -t numOfThreads fastqChunk > outFileName
	hdfs.writeData(logFile,getTimeStamp() + "\tbwa mem started: "+ toolsFolder+"bwa mem " + refFolder+RefFileName + " -p -t " + config.getNumInstances() + " " +tmpFolder+x + " > " + outFileName + "\n")
	val cmd = Seq(toolsFolder+"bwa","mem",refFolder+RefFileName,"-p","-t",config.getNumInstances(),tmpFolder+x) #> new File(outFileName) !


	val bwaKeyValues = new BWAKeyValues(outFileName)
	bwaKeyValues.parseSam()
	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()

	// Delete the temporary files
	val delete = Seq("rm", outFileName, tmpFolder+x )
	hdfs.writeData(logFile,getTimeStamp() + delete.mkString(" ")+"\n")
	delete.!

	//Write to log/bwa
	hdfs.writeData(logFile, getTimeStamp() + "\tbwa mem completed for "+ x +". Number of key value pairs = " + kvPairs.length + "\n")
	


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
	//val toolsFolder = config.getToolsFolder
	val toolsFolder = tmpFolder + "tools/"
	val refFolder = config.getRefFolder
	val numOfThreads = config.getNumThreads
	
	val hdfs = new Operations()

	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion). 
	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More 
	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
	//	Note that MemString here is -Xmx6144m, and already defined as a constant variable above, and so are reference files' names.
	val outputFolder = config.getOutputFolder + "log/vc"
	new File(outputFolder).mkdirs

	new File(tmpFolder).mkdirs

	//val vcLog = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFolder+ "/log" +chrRegion+ ".txt"), "UTF-8"))
	hdfs.createFile(outputFolder+ "/log" +chrRegion+ ".txt")
	//var sortedArray = samRecordsSorted.sortWith((s1: SAMRecord, s2: SAMRecord) => compareSam(s1,s2))
	
	//copy tools folder
	val copyTool = Seq("cp","-r",config.getToolsFolder, tmpFolder)
 	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t"+ copyTool.mkString(" ")+"\n")
 	copyTool.!

	// SAM records should be sorted by this point
	//val chrRange = writeToBAM(tmpFolder + "region" + chrRegion + "-p1.bam", sortedArray, config)
	val chrRange = writeToBAM(tmpFolder + "region" + chrRegion + "-p1.bam", samRecordsSorted, config)
	
	// Picard preprocessing
	//	java MemString -jar toolsFolder/CleanSam.jar INPUT=tmpFolder/regionX-p1.bam OUTPUT=tmpFolder/regionX-p2.bam
	
	val picard1 = Seq("java",MemString,"-jar",toolsFolder + "CleanSam.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p1.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+"-p2.bam")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + picard1.mkString(" ") + "\n")
	picard1.!
	//	java MemString -jar toolsFolder/MarkDuplicates.jar INPUT=tmpFolder/regionX-p2.bam OUTPUT=tmpFolder/regionX-p3.bam	
	//		METRICS_FILE=tmpFolder/regionX-p3-metrics.txt 
	//vcLog.write("java",MemString,"-jar",toolsFolder + "MarkDuplicates.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p2.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+"-p3.bam", "METRICS_FILE="+tmpFolder+"region"+chrRegion+"-p3-metrics.txt")
	val picard2 = Seq("java",MemString,"-jar",toolsFolder + "MarkDuplicates.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p2.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+"-p3.bam", "METRICS_FILE="+tmpFolder+"region"+chrRegion+"-p3-metrics.txt")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + picard2.mkString(" ") + "\n")
	picard2.!
	//	java MemString -jar toolsFolder/AddOrReplaceReadGroups.jar INPUT=tmpFolder/regionX-p3.bam OUTPUT=tmpFolder/regionX.bam 
	//		RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1
	val picard3 = Seq("java",MemString,"-jar",toolsFolder + "AddOrReplaceReadGroups.jar","INPUT="+tmpFolder+"region"+chrRegion+"-p3.bam", "OUTPUT="+tmpFolder+"region"+chrRegion+".bam","RGID=GROUP1","RGLB=LIB1","RGPL=ILLUMINA","RGPU=UNIT1","RGSM=SAMPLE1")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + picard3.mkString(" ") + "\n")
	picard3.!
	// 	java MemString -jar toolsFolder/BuildBamIndex.jar INPUT=tmpFolder/regionX.bam
	val picard4 = Seq("java",MemString,"-jar",toolsFolder + "BuildBamIndex.jar","INPUT="+tmpFolder+"region"+chrRegion+".bam")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + picard4.mkString(" ") + "\n")
	picard4.!
	//	delete tmpFolder/regionX-p1.bam, tmpFolder/regionX-p2.bam, tmpFolder/regionX-p3.bam and tmpFolder/regionX-p3-metrics.txt
	val picard5 = Seq("rm",tmpFolder+"region"+chrRegion+"-p1.bam",tmpFolder+"region"+chrRegion+"-p2.bam",tmpFolder+"region"+chrRegion+"-p3.bam",tmpFolder+"region"+chrRegion+"-p3-metrics.txt").!

	// Make region file 
	val tmpBed = new File(tmpFolder+"tmp"+chrRegion+".bed")
	chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
	//	toolsFolder/bedtools intersect -a refFolder/ExomeFileName -b tmpFolder/tmpX.bed -header > tmpFolder/bedX.bed
	val reg1 = Seq(toolsFolder+"bedtools","intersect","-a",refFolder+ExomeFileName,"-b",tmpFolder+"tmp"+chrRegion+".bed","-header") 
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + reg1.mkString(" ") + " > " +tmpFolder+"bed"+chrRegion+".bed" + "\n")
	reg1 #> new File(tmpFolder+"bed"+chrRegion+".bed") !
	//	delete tmpFolder/tmpX.bed
	val reg2 = Seq("rm",tmpFolder+"tmp"+chrRegion+".bed").!
	
	// Indel Realignment 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt numOfThreads -R refFolder/RefFileName 
	//		-I tmpFolder/regionX.bam -o tmpFolder/regionX.intervals -L tmpFolder/bedX.bed
	val indR1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","RealignerTargetCreator", "-nt", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+".bam","-o",tmpFolder+"region"+chrRegion+".intervals","-L",tmpFolder+"bed"+chrRegion+".bed")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + indR1.mkString(" ") + "\n")
	indR1.!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T IndelRealigner -R refFolder/RefFileName -I tmpFolder/regionX.bam 
	//		-targetIntervals tmpFolder/regionX.intervals -o tmpFolder/regionX-2.bam -L tmpFolder/bedX.bed
	val indR2 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","IndelRealigner", "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+".bam","-L",tmpFolder+"bed"+chrRegion+".bed", "-targetIntervals", tmpFolder+"region"+chrRegion+".intervals", "-o", tmpFolder+"region"+chrRegion+"-2.bam","-L",tmpFolder+"bed"+chrRegion+".bed")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + indR2.mkString(" ") + "\n")
	indR2.!
	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
	val indR3 = Seq("rm",tmpFolder+"region"+chrRegion+".bam", tmpFolder+"region"+chrRegion+".bai", tmpFolder+"region"+chrRegion+".intervals").!


	// Base quality recalibration 
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T BaseRecalibrator -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX.table -L tmpFolder/bedX.bed --disable_auto_index_creation_and_locking_when_reading_rods
	//		-knownSites refFolder/SnpFileName
	val bqr1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","BaseRecalibrator", "-nct", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-2.bam","-o",tmpFolder+"region"+chrRegion+".table","-L",tmpFolder+"bed"+chrRegion+".bed", "--disable_auto_index_creation_and_locking_when_reading_rods", "-knownSites", refFolder+SnpFileName)
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + bqr1.mkString(" ") + "\n")
	bqr1.!
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T PrintReads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX-3.bam -BQSR tmpFolder/regionX.table -L tmpFolder/bedX.bed 
	val bqr2 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","PrintReads", "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-2.bam","-o",tmpFolder+"region"+chrRegion+"-3.bam", "-BQSR", tmpFolder+"region"+chrRegion+".table", "-L", tmpFolder+"bed"+chrRegion+".bed")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + bqr2.mkString(" ") + "\n")
	bqr2.!
	// delete tmpFolder/regionX-2.bam, tmpFolder/regionX-2.bai, tmpFolder/regionX.table
	val bqr3 = Seq("rm",tmpFolder+"region"+chrRegion+"-2.bam", tmpFolder+"region"+chrRegion+"-2.bai", tmpFolder+"region"+chrRegion+".table").!

	// Haplotype -> Uses the region bed file
	// java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T HaplotypeCaller -nct numOfThreads -R refFolder/RefFileName -I 
	//		tmpFolder/regionX-3.bam -o tmpFolder/regionX.vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L tmpFolder/bedX.bed 
	//		--no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods
	val hp1 = Seq("java",MemString,"-jar",toolsFolder + "GenomeAnalysisTK.jar","-T","HaplotypeCaller", "-nct", numOfThreads, "-R", refFolder+RefFileName, "-I", tmpFolder+"region"+chrRegion+"-3.bam","-o",tmpFolder+"region"+chrRegion+".vcf","-stand_call_conf","30.0","-stand_emit_conf","30.0","-L", tmpFolder+"bed"+chrRegion+".bed","--no_cmdline_in_header","--disable_auto_index_creation_and_locking_when_reading_rods")
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\t" + hp1.mkString(" ") + "\n")
	hp1.!
	hdfs.writeData(outputFolder+ "/log" +chrRegion+ ".txt",getTimeStamp() + "\tOutput written to vcf file")
	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
	val hp2 = Seq("rm",tmpFolder+"region"+chrRegion+"-3.bam", tmpFolder+"region"+chrRegion+"-3.bai", tmpFolder+"bed"+chrRegion+".bed").!


	var ret = new Array[(Integer, (Integer, String))](0)
	
		// return the content of the vcf file produced by the haplotype caller.
	//	Return those in the form of <Chromsome number, <Chromosome Position, line>>

	//Reads every line of the temporary vcf file into the return array
	for (line <- Source.fromFile(tmpFolder+"region"+chrRegion+".vcf").getLines; if(line.startsWith("chr"))) {
		var chrm = 0

		if(line.startsWith("chrX"))
			chrm = 23
		else  if(line.startsWith("chrY"))
			chrm = 24
		else
			chrm = line.split("\t")(0).split("chr")(1).toInt

  		ret = ret ++ Array[(Integer, (Integer, String))]((chrm, (line.split("\t")(1).toInt,line)))
	}
	
	hdfs.uploadFile(tmpFolder+"region"+chrRegion+".vcf", config.getOutputFolder+"region"+chrRegion+".vcf")

	//Removes 
	val remVcf = Seq("rm","-r",tmpFolder+"region"+chrRegion+".vcf").!

	return ret // Replace this with what is described in the above 2 lines
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
	val hdfs = new Operations()


	val sc = new SparkContext(conf)
	sc.broadcast(config)
	
	// Comment these two lines if you want to see more verbose messages from Spark
	Logger.getLogger("org").setLevel(Level.OFF);
	Logger.getLogger("akka").setLevel(Level.OFF);
		
	var t0 = System.currentTimeMillis
	
	// Rest of the code goes here

	hdfs.createDirectory(config.getOutputFolder)

	val logFolder = config.getOutputFolder + "/log"
	//new File(logFolder).mkdirs
	hdfs.createDirectory(logFolder)
	hdfs.createDirectory(logFolder + "/bwa")
	hdfs.createDirectory(logFolder + "/vc")

	hdfs.createFile(config.getOutputFolder+ "chromossome" + ".vcf")
	hdfs.createFile(logFolder + "/log.txt")

  //	val bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(config.getOutputFolder+ "chromossome" + ".vcf"), "UTF-8"))
  //	val bwLog = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(logFolder + "/log.txt"), "UTF-8"))

//	println(initialRdd.partitions.size)
  	//val folder = new File(config.getInputFolder)
    val files = hdfs.directoryList(config.getInputFolder);


    //val arraynerino: ListBuffer[Array[(Int, SAMRecord)]] = new ListBuffer[Array[(Int, SAMRecord)]]()

    //String    
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp()+ "Creating initalRDD\n")
    val initRDD = sc.parallelize(files, config.getNumInstances.toInt)

    //(Chromossome, record)
    //(Int, Sam)    
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Running BWA\n")
    val chromoRDD = initRDD.mapPartitions(files => {
    	for(file <- files )
    			yield bwaRun(file, config)
    	}).flatMap(y => y)

    
    //(Chromossome, recordCount), sorted in descending order of SamRecords
    //(Int, int)    
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Grouping SAMRecord by chromossome and sorting each chromossome by the number of SamRecords in descending order\n")
    val chromoCountRDD = chromoRDD.mapPartitions(
    	lines => 
    	for(line <- lines)
    		yield (line._1, 1)
    	).reduceByKey(_+_, config.getNumInstances().toInt).sortBy(_._2, false)

    //after zipWithIndex
    //  ((Chromossome, recordCount), index)
    //  ((int,int), int)
    //after map
    // since the rdd was sorted, if we map and index and get the moduled of it by the number of tasks, we can map every chromossome to a region and it will be distributed
    // 
    // (Chromossome, Region)
    //bwLog.write(getTimeStamp() + "Creating the regions and assigning the chromossomes to these regions\n")
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Creating the regions and assigning the chromossomes to these regions\n")
    val regionChromoRDD = chromoCountRDD.zipWithIndex().mapPartitions(lines => 
																		for(line <- lines)
																			yield( line._1._1 ,line._2 % config.getNumInstances.toInt) )

    //Joins Chromossome region with its SAMRecords, mapping it to
    //(Region, SamRecord)
    //(Int, Array(SamRecord))
    //bwLog.write(getTimeStamp() + "Assigning each region to its SamRecord value\n")
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Assigning each region to its SamRecord value\n")
    val regionSamRDD = regionChromoRDD.join(chromoRDD).mapPartitions(lines => for(line <- lines) yield (line._2._1, Array[SAMRecord](line._2._2)))
    
    
    //Groups every SamRecord with the same region
   // bwLog.write(getTimeStamp() + "Variant call\n")
    hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Variant call\n")
    val variantCallRDD = regionSamRDD.reduceByKey(_ ++ _, config.getNumInstances().toInt).mapPartitions(lines =>{
    	var line =lines.next
    	variantCall(line._1.toInt, line._2.sortWith((s1: SAMRecord, s2: SAMRecord) => compareSam(s1,s2)), config).toIterator
    	} ).distinct

	
	hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Sorting the results of variantCall\n")
  	val outputRDD = variantCallRDD.sortBy(line => (line._1, line._2._1))

  	
  	val bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(config.getTmpFolder+ "chromossome" + ".vcf"), "UTF-8"))


  	hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Writing the final vcf\n")
  	
  	outputRDD.collect.foreach(x => bw.write(x._2._2 + "\n"))
	bw.close()


	hdfs.writeData(logFolder + "/log.txt" , getTimeStamp() + "Uploading vcf\n")
	hdfs.uploadFile(config.getTmpFolder+ "chromossome" + ".vcf", config.getOutputFolder+"chromossome"+".vcf")
 

 
    val remTemp = Seq("rm", "-rf", config.getTmpFolder).!
  	new File(config.getTmpFolder).mkdirs
	
	//val et = (System.currentTimeMillis - t0) / 1000 
	//println("|Execution time: %d mins %d secs|".format(et/60, et%60))
}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition
