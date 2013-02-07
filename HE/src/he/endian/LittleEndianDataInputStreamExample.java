package he.endian;

import java.io.*;

/**
 * @author Zhu Zhixiang, zzxiang21cn@hotmail.com
 */

public class LittleEndianDataInputStreamExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		FileInputStream fileStream = new FileInputStream("/Users/uqgchen5/Documents/workspace/Test/poly.grm.bin");
		DataInputStream bigEndianDataStream = new DataInputStream(fileStream);
		LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(bigEndianDataStream, Float.SIZE);
		
		for (int i = 0; i < 10; ++i) {
			if (littleEndianDataStream.available() > 0) {
				System.out.println(littleEndianDataStream.readFloat());
			}
		}
	}

}
