package gear.util;

import gear.ConstValues;

public class BinaryInputFileProfile
{
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
	{
		long startTime = System.currentTimeMillis();
		int cnt = 0;
		for (int i = 0; i < 5; ++i)
		{
			BinaryInputFile file = new BinaryInputFile("data/LittleEndianProfile.dat", "data");
			file.setLittleEndian(true);
			while (file.available() >= ConstValues.FLOAT_SIZE)
			{
				@SuppressWarnings("unused")
				float f = file.readFloat();
				++cnt;
			}
			file.close();
		}
		long timeElapsed = (System.currentTimeMillis() - startTime) / 1000;
		System.out.println("It takes " + timeElapsed + " sec to read " + cnt + " float(s).");
	}
}
