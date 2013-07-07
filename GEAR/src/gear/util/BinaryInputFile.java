package gear.util;

import gear.ConstValues;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

public class BinaryInputFile
{
	public BinaryInputFile(String fileName, String fileType)
	{
		try
		{
			stream = new FileInputStream(fileName);
		}
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "Cannot find the file '" + fileName + "'.");
		}
		
		channel = stream.getChannel();
		buffer = ByteBuffer.allocateDirect(BUFFER_CAPACITY);
		
		this.fileName = fileName;
		this.fileType = fileType;
	}
	
	public void close()
	{
		try
		{
			stream.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when closing the file '" + fileName + "'.");
		}
	}
	
	public long available()
	{
		long size = 0;
		try
		{
			size = channel.size();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O excpetion occurred when getting the size of file '" + fileName + "'.");
			System.exit(1);
		}
		
		long position = 0;
		try
		{
			position = channel.position();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O excpetion occurred when reading the file '" + fileName + "'.");
		}
		
		return size - position;
	}
	
	public float readFloat()
	{
		buffer.clear();
		buffer.limit(ConstValues.FLOAT_SIZE);
		
		try
		{
			int nread = channel.read(buffer);
			
			if (nread != ConstValues.FLOAT_SIZE)
			{
				String msg = "";
				msg += "The " + fileType + " file '" + fileName + "' contains an error at byte " + channel.position() + ". ";
				if (nread == -1)
				{
					msg += "Unexpected end-of-file was reached.";
				}
				else
				{
					msg += "The program expected to read " + ConstValues.FLOAT_SIZE + " bytes, but actually read " + nread + ".";
				}
				Logger.printUserError(msg);
				System.exit(1);
			}
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O excpetion occurred when reading the file '" + fileName + "'.");
		}
		
		buffer.rewind();
		
		return buffer.getFloat();
	}
	
	public void setLittleEndian(boolean isLittleEndian)
	{
		buffer.order(isLittleEndian ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN);
	}
	
	private FileInputStream stream;
	private FileChannel channel;
	private ByteBuffer buffer;
	private static final int BUFFER_CAPACITY = 32;
	private String fileName;
	private String fileType;
}
