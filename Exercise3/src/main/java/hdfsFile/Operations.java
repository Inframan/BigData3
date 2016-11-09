package hdsfFile; 
import java.io.*;

import java.util.ArrayList;

import org.apache.hadoop.conf.Configuration;

import org.apache.hadoop.fs.*;

 

public class Operations {

	FileSystem hdfs;
	Path workingDir;

	public Operations() throws IOException
	{
		Configuration conf = new Configuration();
        conf.set("fs.defaultFS", "hdfs://127.0.0.1:9000/user/gsouto");
        conf.set("hadoop.job.ugi", "gsouto");

		hdfs = FileSystem.get(conf);
		workingDir=hdfs.getWorkingDirectory();
	}

	public boolean exists(String path)throws IOException
	{
		return hdfs.exists(new Path(path));
	}

	public void deleteDirectory(String path) throws IOException
	{
		Path folder= new Path(path);
		//folder = Path.mergePaths(workingDir, folder);
		if(hdfs.exists(folder))
		{

			//Delete existing Directory
			hdfs.delete(folder, true);
		}
	}

	public void createDirectory(String path) throws IOException
	{		

		Path newFolder= new Path(path);
		//newFolder = Path.mergePaths(workingDir, newFolder);
		deleteDirectory(path);

		hdfs.mkdirs(newFolder);     //Create new Directory
	}

	public void uploadFile(String localPath, String hdfsPath) throws IOException
	{
		//Copying File from local to HDFS

		Path localFile = new Path(localPath);

		Path hdfsFile= new Path(hdfsPath);
		//hdfsFile = Path.mergePaths(workingDir, hdfsFile);
		hdfs.copyFromLocalFile(localFile, hdfsFile);
	}

	public void downloadFile(String hdfsPath, String localPath) throws IOException
	{
		//Copying File from HDFS to local

		Path localFile =new Path(localPath);
		Path hdfsFile= new Path(hdfsPath);
		
		//hdfsFile = Path.mergePaths(workingDir, hdfsFile);
		hdfs.copyToLocalFile(hdfsFile, localFile);
	}

	public void createFile(String path) throws IOException
	{

		//Creating a file in HDFS

		Path newFile = new Path(path);

		//newFile = Path.mergePaths(workingDir, newFile);
		hdfs.createNewFile(newFile);

	}

	public void writeData(String path, String data) throws IOException
	{

		//Writing data to a HDFS file

		byte[] byt = data.getBytes();

		FSDataOutputStream writer = hdfs.append(new Path(path));

		writer.write(byt);

		writer.close();
	}

	public String[] directoryList(String path) throws IOException
	{
		FileStatus[] fileStatus = hdfs.listStatus(new Path(path));
		ArrayList<String> temp = new ArrayList<String>();

		for(int i = 0; i < fileStatus.length; i++)
		{
			if(fileStatus[i].getPath().getName().endsWith(".fq.gz"))
				temp.add(fileStatus[i].getPath().getName());
		}
		return temp.toArray(new String[temp.size()]);
	}


}

