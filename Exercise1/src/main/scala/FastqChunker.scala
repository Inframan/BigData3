import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.log4j.Logger
import org.apache.log4j.Level

import org.apache.hadoop.io.compress.GzipCodec
import java.util.zip._

import java.io._
import java.nio.file.{Paths, Files}

object FastqChunker 
{
	def main(args: Array[String]) 
	{
		if (args.size < 3)
		{
			println("Not enough arguments!\nArg1 = number of parallel tasks = number of chunks\nArg2 = input folder\nArg3 = output folder")
			System.exit(1)
		}
		
		val prllTasks = args(0)
		val inputFolder = args(1)
		val outputFolder = args(2)
		
		if (!Files.exists(Paths.get(inputFolder)))
		{
			println("Input folder " + inputFolder + " doesn't exist!")
			System.exit(1)
		}
			 
		// Create output folder if it doesn't already exist
		new File(outputFolder).mkdirs
		
		println("Number of parallel tasks = number of chunks = " + prllTasks + "\nInput folder = " + inputFolder + "\nOutput folder = " + outputFolder)
		
		val conf = new SparkConf().setAppName("DNASeqAnalyzer")
		conf.setMaster("local[" + prllTasks + "]")
		conf.set("spark.cores.max", prllTasks)	
		

		val sc = new SparkContext(conf)
		
		// Comment these two lines if you want to see more verbose messages from Spark
		Logger.getLogger("org").setLevel(Level.OFF);
		Logger.getLogger("akka").setLevel(Level.OFF);	
			
		var t0 = System.currentTimeMillis
		
		// Rest of the code goes here 


		//Read the line
		//Associate every line with an id
		//Since every 4 lines should be associated with the same id, it's divided by 4 
		//And we map the number to the key
		//Then we reduce so that the 4 lines join each other becomes together
		//(Id, DNA sequence)
		//(Long, Text)
		val first_avroRdd =  sc.textFile(args(1) + "/fastq1.fq").zipWithIndex().mapPartitions(
			lines => 
			{
			for(line <- lines)
		 		yield (line._2 /4 , line._1)  
		 	}).reduceByKey(_ + "\n"+ _, prllTasks.toInt).setName("rdd_first_read")
		
		//Read the line
		//Associate every line with an id
		//Since every 4 lines should be associated with the same id, it's divided by 4 
		//And we map the number to the key
		//(Id, DNA sequence)
		//(Long, Text)
		val second_avroRdd = sc.textFile(args(1) + "/fastq2.fq").zipWithIndex().mapPartitions(lines => 
			{
			for(line <- lines)
		 		yield (line._2 /4 , line._1)  
		 	} ).reduceByKey(_+ "\n" +_, prllTasks.toInt).setName("rdd_second_read")

		//(Id, ())
		val joined = first_avroRdd.join(second_avroRdd).sortByKey(true).setName("rdd_joined")
		

		

		val outputRDD = joined.mapPartitionsWithIndex((index: Int, it: Iterator[(Long, (String, String))]) => {
			val bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFolder+ "/" + index + ".fq.gz")), "UTF-8"))
			
			for(line <- it)
			{
				bw.write(line._2._1 + "\n" + line._2._2)
				bw.newLine()
			}
		
			bw.close()

			it
			} ).setName("rdd_output")

		outputRDD.count
		val et = (System.currentTimeMillis - t0) / 1000 
		println("|Execution time: %d mins %d secs|".format(et/60, et%60))
	}
	//////////////////////////////////////////////////////////////////////////////




} // End of Class definition
