package gear.simulation;

import gear.CommandArguments;
import gear.CommandImpl;
import gear.ConstValues;
import gear.util.Logger;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;

public final class SimuFamilyCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		this.cmdArgs = (SimuFamilyCommandArguments)cmdArgs;
		
		init();

		for (int i = 0; i < this.cmdArgs.getNumberOfFamilies(); i++)
		{
			generateNuclearFamily(NKid[i], NAffKid[i], i);
		}
	}
	
	private void init()
	{
		rnd = new RandomDataImpl();
		
		if (cmdArgs.getSeed() == null)
		{
			rnd.reSeed();
		}
		else
		{
			rnd.reSeed(cmdArgs.getSeed());
		}

		NKid = new int[cmdArgs.getNumberOfFamilies()];
		Arrays.fill(NKid, 2);
		NAffKid = new int[cmdArgs.getNumberOfFamilies()];
		Arrays.fill(NAffKid, 1);

		maf = new double[cmdArgs.getNumberOfMarkers()];
		Arrays.fill(maf, 0.5);
		LD = new double[cmdArgs.getNumberOfMarkers()];
		Arrays.fill(LD, 0);
		rec = new double[cmdArgs.getNumberOfMarkers()];
		Arrays.fill(rec, 0.5);
		rec[0] = maf[0];

		gm = new int[cmdArgs.getNumberOfFamilies() * famSize][cmdArgs.getNumberOfMarkers()];
	}

	private void generateNuclearFamily(int nkid, int affKid, int famIdx)
	{
		int[][] p = sampleChromosome(famIdx, 0);
		int[][] m = sampleChromosome(famIdx, 1);
		
		for (int i = 0; i < nkid; i++)
		{
			generateBaby(p, m, famIdx, i + 2);
		}

		if (cmdArgs.getMakeBed())
		{
			writeBFile();
		}
		else
		{
			writeFile();
		}
	}

	private int[][] sampleChromosome(int famIdx, int shift)
	{
		int[][] v = new int[maf.length][2];
		for (int i = 0; i < maf.length; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				double r = rnd.nextUniform(0, 1);
				if (i == 0)
				{
					v[i][j] = r < maf[i] ? 0 : 1;
				}
				else
				{
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[i - 1][j];
					double f1 = a == 0 ? maf[i - 1] : (1 - maf[i - 1]);
					double f2 = a == 0 ? maf[i] : (1 - maf[i]);
					v[i][j] = d < (f1 * f2 + LD[i - 1]) / f1 ? v[i - 1][j]
							: (1 - v[i - 1][j]);
				}
			}
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}

		return v;
	}

	private void generateBaby(int[][] p, int[][] m, int famIdx, int shift)
	{
		int[][] v = new int[maf.length][2];

		for (int i = 0; i < 2; i++)
		{
			int[][] chr = i == 0 ? p : m;
			int idx = 1;
			for (int j = 0; j < maf.length; j++)
			{
				double r = rnd.nextUniform(0, 1);
				idx = r < rec[j] ? 1 - idx : idx;
				v[j][i] = chr[j][idx];
			}
		}

		for (int i = 0; i < maf.length; i++)
		{
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}
	}

	public void writeFile()
	{
		PrintWriter ped = null;
		PrintWriter map = null;

		try
		{
			ped = new PrintWriter(new BufferedWriter(new FileWriter(cmdArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(cmdArgs.getOutRoot() + ".map")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .ped and .map files.");
		}

		for (int h = 0; h < cmdArgs.getNumberOfFamilies(); h++)
		{
			int fid = (h + 1) * 10000;
			int pid = fid + 1;
			int mid = fid + 2;

			ped.print(fid + " ");
			ped.print(pid + " ");
			ped.print(0 + " ");
			ped.print(0 + " ");
			ped.print(1 + " ");
			ped.print(1 + " ");

			for (int i = 0; i < maf.length; i++)
			{
				StringBuilder sb = new StringBuilder();
				switch (gm[h * famSize][0])
				{
				case 0:
					sb.append(A[0] + " " + A[0]);
					break;
				case 1:
					sb.append(A[0] + " " + A[0]);
					break;
				case 2:
					sb.append(A[0] + " " + A[0]);
					break;
				default:
					break;
				}
				ped.print(sb.toString());
			}
			ped.print("\n");

			ped.print(fid + " ");
			ped.print(mid + " ");
			ped.print(0 + " ");
			ped.print(0 + " ");
			ped.print(2 + " ");
			ped.print(1 + " ");

			for (int i = 0; i < maf.length; i++)
			{
				StringBuilder sb = new StringBuilder();
				switch (gm[h * famSize + 1][0])
				{
				case 0:
					sb.append(A[0] + " " + A[0]);
					break;
				case 1:
					sb.append(A[0] + " " + A[0]);
					break;
				case 2:
					sb.append(A[0] + " " + A[0]);
					break;
				default:
					break;
				}
				ped.print(sb.toString());
			}
			ped.print("\n");

			for (int j = 0; j < 2; j++)
			{
				ped.print(fid + " ");
				ped.print((fid + 3 + j) + " ");
				ped.print(pid + " ");
				ped.print(mid + " ");

				try
				{
					ped.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
					ped.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				}
				catch (MathException e)
				{
					Logger.handleException(e, "Failed to generate the random values.");
				}

				for (int i = 0; i < maf.length; i++)
				{
					StringBuilder sb = new StringBuilder();
					switch (gm[h * famSize + 2 + j][0])
					{
					case 0:
						sb.append(A[0] + " " + A[0]);
						break;
					case 1:
						sb.append(A[0] + " " + A[0]);
						break;
					case 2:
						sb.append(A[0] + " " + A[0]);
						break;
					default:
						break;
					}
					ped.print(sb.toString());
				}
				ped.print("\n");
			}
		}
		
		for (int i = 0; i < maf.length; i++)
		{
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (maf.length * 1.0) + " ");
			map.println(i * 100);
		}

		ped.close();
		map.close();
	}

	public void writeBFile()
	{
		DataOutputStream bed = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		
		try
		{
			bed = new DataOutputStream(new FileOutputStream(cmdArgs.getOutRoot() + ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(cmdArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(cmdArgs.getOutRoot() + ".bim")));

		} 
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, .fam and .bim files.");
		}

		for (int i = 0; i < gm.length; i++)
		{
			fam.print("sample_" + i + " ");
			fam.print(1 + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			try
			{
				fam.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				fam.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
			}
			catch (MathException e)
			{
				Logger.handleException(e, "Failed to generate the random values.");
			}
		}

		try
		{
			bed.writeByte(ConstValues.PLINK_BED_BYTE1);
			bed.writeByte(ConstValues.PLINK_BED_BYTE2);
			bed.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < cmdArgs.getNumberOfMarkers(); i++)
			{
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < gm.length; j++)
				{
					int g = gm[j][i] + 1;
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

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (gm.length - 1))
					{
						if (idx == 4)
						{
							bed.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} 
					else
					{
						bed.writeByte(gbyte);
					}
				}
			}
			bed.close();
		} 
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing the .bed file.");
		}

		for (int i = 0; i < cmdArgs.getNumberOfMarkers(); i++)
		{
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (cmdArgs.getNumberOfMarkers() * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A[0] + " " + A[1]);
		}

		bim.close();
		fam.close();

	}

	private SimuFamilyCommandArguments cmdArgs;
	
	private RandomDataImpl rnd;

	private final String[] A = { "A", "C" };
	private int[] NKid = null;
	private int[] NAffKid = null;
	private double[] LD = null;
	private double[] rec = null;
	private double[] maf = null;

	private int[][] gm = null;
	private final int famSize = 4; 
}
