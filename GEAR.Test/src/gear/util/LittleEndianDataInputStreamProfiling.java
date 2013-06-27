package gear.util;

import java.io.*;

public class LittleEndianDataInputStreamProfiling
{

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
	{
		long startTime = System.currentTimeMillis();
		FileInputStream fileStream = new FileInputStream("data/LittleEndian.dat");
		DataInputStream bigEndianDataStream = new DataInputStream(fileStream);
		LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(bigEndianDataStream, Float.SIZE);
		int cnt = 0;
		while (littleEndianDataStream.available() > 0)
		{
			@SuppressWarnings("unused")
			float f = littleEndianDataStream.readFloat();
			++cnt;
		}
		long timeElapsed = (System.currentTimeMillis() - startTime) / 1000;
		System.out.println("It takes " + timeElapsed + " sec to read " + cnt + " float(s).");
	}

}
