package gear.util;

import java.io.*;
import java.nio.*;

/**
 * @author Zhu Zhixiang, zzxiang21cn@hotmail.com
 */

public class LittleEndianDataInputStream
{
	public LittleEndianDataInputStream(DataInputStream inStream, int bufferSize)
	{
		inStream_ = inStream;
		buffer_ = ByteBuffer.allocate(bufferSize);
	}

	public float readFloat() throws IOException
	{
		buffer_.clear();
		buffer_.order(ByteOrder.BIG_ENDIAN);
		buffer_.putFloat(inStream_.readFloat());
		buffer_.flip();
		return buffer_.order(ByteOrder.LITTLE_ENDIAN).getFloat();
	}

	public int available() throws IOException
	{
		return inStream_.available();
	}

	private DataInputStream inStream_;
	private ByteBuffer buffer_;
}
