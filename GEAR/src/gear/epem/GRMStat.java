package gear.epem;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.BinaryInputFile;
import gear.util.FileUtil;
import gear.util.Logger;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;

public class GRMStat
{

	private final String delim = "\\s+";

	private int M = -1;
	private double Nt = 0;
	private double N = 0;
	private double mean = 0;
	private double prod = 0;
	private double v = 0;
	private double Ne = 0;
	private double Me = 0;
	private double grmCutoff = 0;

	private StringBuffer sb = new StringBuffer();

	public GRMStat()
	{
		M = CmdArgs.INSTANCE.getHEArgs().getGrmM();
		if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
		{
			grmCutoff = CmdArgs.INSTANCE.getHEArgs().AbsGrmCutoff();
		} else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
		{
			grmCutoff = CmdArgs.INSTANCE.getHEArgs().GrmCutoff();
		}

	}

	public void GetGRMStats()
	{
		if (CmdArgs.INSTANCE.getHEArgs().isGrmBinary())
		{
			BinaryGRM();
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isGrmTxt())
		{
			txtGRM();
		}
		else
		{
			gzGRM();
		}

		double s = Nt - N;
		double r_sq = 0;
		sb.append("grm file: " + CmdArgs.INSTANCE.getHEArgs().getGrm()
				+ "\n");

		sb.append("Total lines in grm: " + Nt + "\n");
		sb.append("Read " + N + " lines." + "\n");
		sb.append("Individuals: " + s + "\n");
		sb.append("Mean is " + mean + "\n");
		sb.append("Variance is " + v + "\n");
		sb.append("The effective sample size is " + Ne + "\n");
		sb.append("Effective number of markers is " + Me + "\n");

		if (M > 0) 
		{
			r_sq = (M - Me) / (Me * (M-1));
			sb.append("Average r^2 LD between markers is " + r_sq + "\n");
		}

		Logger.printUserLog(sb.toString());

		StringBuilder fsb = new StringBuilder();
		fsb.append(CmdArgs.INSTANCE.out);
		fsb.append(".gs");
		PrintWriter pw = null;
		try
		{
			pw = new PrintWriter(fsb.toString());
		} 
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		pw.append(sb);
		pw.close();

	}

	private void BinaryGRM()
	{
		BufferedReader reader = FileUtil.FileOpen(CmdArgs.INSTANCE
				.getHEArgs().getGrmId());

		int size = 0;
		try
		{
			while (reader.readLine() != null)
			{
				size++;
			}
			reader.close();
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the ID file '"
							+ CmdArgs.INSTANCE.getHEArgs().getGrmId()
							+ "'.");
		}

		// *************************************read grm file	
		BinaryInputFile grmFile = new BinaryInputFile(CmdArgs.INSTANCE.getHEArgs().getGrm(), "GRM");
		grmFile.setLittleEndian(true);
		
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				Nt++;
				double g = 0;
				if (grmFile.available() >= ConstValues.FLOAT_SIZE)
				{
					g = grmFile.readFloat();
				}

				int id1 = i;
				int id2 = j;
				if (id1 == id2)
					continue;

				if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
				{
					if (Math.abs(g) > grmCutoff)
						continue;
				} 
				else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
				{
					if (g > grmCutoff)
						continue;
				}

				mean += g;
				prod += g * g;
				N++;
			}
		}
		
		grmFile.close();

		getEffectiveNumber();
	}

	private void gzGRM()
	{
		FileInputStream fin = null;
		try
		{
			fin = new FileInputStream(CmdArgs.INSTANCE.getHEArgs()
					.getGrm());
		} 
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "Cannot open the GRM file '"
					+ CmdArgs.INSTANCE.getHEArgs().getGrm() + "'.");
		}

		GZIPInputStream gzis = null;
		try
		{
			gzis = new GZIPInputStream(fin);
		} 
		catch (IOException e)
		{
			Logger.handleException(e, "Cannot open the GRM archive '"
					+ CmdArgs.INSTANCE.getHEArgs().getGrm() + "'.");
		}
		InputStreamReader xover = new InputStreamReader(gzis);

		BufferedReader grmFile = new BufferedReader(xover);

		String line;
		try
		{
			while ((line = grmFile.readLine()) != null)
			{
				Nt++;
				String[] s = line.split(delim);
				int id1 = Integer.parseInt(s[0]) - 1;
				int id2 = Integer.parseInt(s[1]) - 1;

				if (id1 == id2)// exclude diagonal if it is required
					continue;
				double g = Double.parseDouble(s[3]);

				if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
				{
					if (Math.abs(g) > grmCutoff)
						continue;
				} 
				else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
				{
					if (g > grmCutoff)
						continue;
				}

				mean += g;
				prod += g * g;
				N++;
			}
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the GRM archive '"
							+ CmdArgs.INSTANCE.getHEArgs().getGrm()
							+ "'.");
		}

		getEffectiveNumber();
	}

	private void txtGRM()
	{
		BufferedReader grmFile = FileUtil.FileOpen(CmdArgs.INSTANCE
				.getHEArgs().getGrm());

		String line;
		try
		{
			while ((line = grmFile.readLine()) != null)
			{
				Nt++;
				String[] s = line.split(delim);
				int id1 = Integer.parseInt(s[0]) - 1;
				int id2 = Integer.parseInt(s[1]) - 1;

				if (id1 == id2)// exclude diagonal if it is required
					continue;
				double g = Double.parseDouble(s[3]);

				if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
				{
					if (Math.abs(g) > grmCutoff)
						continue;
				} 
				else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
				{
					if (g > grmCutoff)
						continue;
				}

				mean += g;
				prod += g * g;
				N++;
			}
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the GRM file '"
							+ CmdArgs.INSTANCE.getHEArgs().getGrm()
							+ "'.");
		}

		getEffectiveNumber();
	}

	private void getEffectiveNumber()
	{
		mean /= N;
		prod /= N;
		v = prod - mean * mean;
		Ne = -1 / mean + 1;
		Me = 1 / v;
	}

}
