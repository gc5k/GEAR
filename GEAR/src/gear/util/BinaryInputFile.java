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
		// Little-Endian is default to false because Java itself uses big-endian
		this(fileName, fileType, /*littleEndian*/false);
	}
	
	public BinaryInputFile(String fileName, String fileType, boolean littleEndian)
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
		setLittleEndian(littleEndian);
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
	
	public int readInt()
	{
		read(ConstValues.INT_SIZE);
		return buffer.getInt();
	}
	
	public long readLong()
	{
		read(ConstValues.LONG_SIZE);
		return buffer.getLong();
	}
	
	public float readFloat()
	{
		read(ConstValues.FLOAT_SIZE);
		return buffer.getFloat();
	}
	
	public double readDouble()
	{
		read(ConstValues.DOUBLE_SIZE);
		return buffer.getDouble();
	}
	
	public void setLittleEndian(boolean isLittleEndian)
	{
		buffer.order(isLittleEndian ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN);
	}
	
	private void read(int size)
	{
		buffer.clear();
		buffer.limit(size);
		
		try
		{
			int nread = channel.read(buffer);
			
			if (nread != size)
			{
				String msg = "";
				msg += "The " + fileType + " file '" + fileName + "' contains an error at byte " + channel.position() + ". ";
				if (nread == -1)
				{
					msg += "Unexpected end-of-file was reached.";
				}
				else
				{
					msg += "The program expected to read " + size + " bytes, but actually read " + nread + ".";
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
	}
	
	private FileInputStream stream;
	private FileChannel channel;
	private ByteBuffer buffer;
	private static final int BUFFER_CAPACITY = 32;
	private String fileName;
	private String fileType;
}
