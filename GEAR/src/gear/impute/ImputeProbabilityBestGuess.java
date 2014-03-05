package gear.impute;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.BufferedReader;

public class ImputeProbabilityBestGuess
{
	private String batchFile = null;
	private ArrayList<String> batchList = NewIt.newArrayList();

	public ImputeProbabilityBestGuess()
	{
		batchFile = CmdArgs.INSTANCE.imputeBatchFile;
		BufferedReader br = BufferedReader.openTextFile(batchFile, "IMPUTE batch");

		String line = null;
		while( (line = br.readLine()) != null)
		{
			batchList.add(line);
		}
		Logger.printUserLog(batchList.size() + " files in '" + batchFile + "'." );
	}

	public void convert()
	{
		DataOutputStream bedout = null;
		PrintWriter bim = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(CmdArgs.INSTANCE.out + ".bed"));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out + ".bim")));
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, and .bim files.");
		}

		String line = null;
		int GenoCnt = 0;
		int SampleSize = 0;
		try {
			Logger.printUserLog("Converting IMPUTE genotype probabilities into 'best-guess' genotypes");
			for(int i = 0; i < batchList.size(); i++)
			{
				int gC = 0;
				BufferedReader br = BufferedReader.openTextFile(batchList.get(i), "batch list");
				while((line = br.readLine()) != null)
				{
					gC++;
					GenoCnt++;
					String [] dos= line.split("\\s+");
					bim.println(dos[0] + "\t" + dos[1] + "\t" + "0" + "\t" + dos[2] + "\t" + dos[3] + "\t" + dos[4]);
					if (GenoCnt == 1) 
					{
						SampleSize = (dos.length - 5)/3;
						Logger.printUserLog( SampleSize + " individuals in total.");
					}
					else
					{
						if ( (dos.length -5)/3 != SampleSize)
						{
							Logger.printUserLog("ERROR! " +  batchList.get(i) + " has " + ((dos.length - 5)/3) + " individuals.");
							System.exit(0);
						}
					}

					byte gbyte = 0;
					int idx = 0;
					for(int j = 0; j < SampleSize; j++)
					{
						int g = bestGuess(Float.parseFloat(dos[4 + j*3 + 1]), Float.parseFloat(dos[4 + j*3 + 2]), Float.parseFloat(dos[4 + j*3 + 3]));
						g <<= 2 * idx;
						gbyte |= g;
						idx++;

						if (j != (SampleSize - 1))
						{
							if (idx == 4)
							{
								bedout.writeByte(gbyte);
								gbyte = 0;
								idx = 0;
							}
						}
						else
						{
							bedout.writeByte(gbyte);
						}
					}
				}
				Logger.printUserLog("Read " + gC + " loci in '" + batchList.get(i) + "'.");
			}
			Logger.printUserLog( GenoCnt + " loci have been converted into best gusess genotypes, which have been saved in '" + CmdArgs.INSTANCE.out + ".bed' and " + "'" + CmdArgs.INSTANCE.out + ".bim'.");

			bedout.close();
			bim.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An exception occurred when writing the bed file.");					
		}
	}

	private int bestGuess(float g1, float g2, float g3)
	{

		int g = 0;
		if(g1 >= 0.5)
		{
			g = 0;
		}
		else if (g2 >= 0.5)
		{
			g = 1;
		}
		else if (g3 >= 0.5)
		{
			g = 2;
		} else {
			float t = g1;
			if (t < g2) {
				g = 1;
				t = g2;
			}
			if (t < g3) {
				g = 2;
				t = g3;
			}
		}
		switch (g)
		{
			case 0:
				g = 0;
				break;
			case 1:
				g = 2;
				break;
			case 2:
				g = 3;
				break;
			default:
				g = 1;
				break; // missing
		}

		return g;
	}
}
