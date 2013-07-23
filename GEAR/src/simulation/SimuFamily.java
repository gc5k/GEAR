package simulation;

import gear.CmdArgs;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;

public class SimuFamily
{
	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;
	
	private RandomDataImpl rnd;
	private long seed = 2011;

	private final String[] A = { "A", "C" };
	private int NFam = 100;
	private int NMarker = 10;
	private int[] NKid = null;
	private int[] NAffKid = null;
	private double[] LD = null;
	private double[] rec = null;
	private double[] maf = null;

	private int[][] gm = null;
	private int famSize = 4;


	public SimuFamily(int NFam, int NMarker, long seed)
	{
		this.NFam = NFam;
		this.NMarker = NMarker;
		this.seed = seed;
		initial();
	}

	private void initial()
	{
		rnd = new RandomDataImpl();
		rnd.reSeed(seed);

		NKid = new int[NFam];
		Arrays.fill(NKid, 2);
		NAffKid = new int[NFam];
		Arrays.fill(NAffKid, 1);

		maf = new double[NMarker];
		Arrays.fill(maf, 0.5);
		LD = new double[NMarker];
		Arrays.fill(LD, 0);
		rec = new double[NMarker];
		Arrays.fill(rec, 0.5);
		rec[0] = maf[0];

		gm = new int[NFam * famSize][NMarker];
	}

	public void generateSample()
	{

		for (int i = 0; i < NFam; i++)
		{
			generateNuclearFamily(NKid[i], NAffKid[i], i);
		}

	}

	private void generateNuclearFamily(int nkid, int affKid, int famIdx)
	{

		int[][] p = sampleChromosome(famIdx, 0);
		int[][] m = sampleChromosome(famIdx, 1);
		for (int i = 0; i < nkid; i++)
		{
			generateBaby(p, m, famIdx, i + 2);
		}

		if (CmdArgs.INSTANCE.makebedFlag)
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
		PrintWriter pedout = null;
		PrintWriter map = null;

		try
		{
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(
					CmdArgs.INSTANCE.out + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(
					CmdArgs.INSTANCE.out + ".map")));
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (int h = 0; h < famSize; h++)
		{
			int fid = (h + 1) * 10000;
			int pid = fid + 1;
			int mid = fid + 2;

			pedout.print(fid + " ");
			pedout.print(pid + " ");
			pedout.print(0 + " ");
			pedout.print(0 + " ");
			pedout.print(1 + " ");
			pedout.print(1 + " ");

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
				pedout.print(sb.toString());
			}
			pedout.print("\n");

			pedout.print(fid + " ");
			pedout.print(mid + " ");
			pedout.print(0 + " ");
			pedout.print(0 + " ");
			pedout.print(2 + " ");
			pedout.print(1 + " ");

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
				pedout.print(sb.toString());
			}
			pedout.print("\n");

			for (int j = 0; j < 2; j++)
			{
				pedout.print(fid + " ");
				pedout.print((fid + 3 + j) + " ");
				pedout.print(pid + " ");
				pedout.print(mid + " ");

				try
				{
					pedout.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
					pedout.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				}
				catch (MathException e)
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
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
					pedout.print(sb.toString());
				}
				pedout.print("\n");
			}
		}
		
		for (int i = 0; i < maf.length; i++)
		{
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (maf.length * 1.0) + " ");
			map.println(i * 100);
		}

		pedout.close();
		map.close();
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		try
		{
			bedout = new DataOutputStream(new FileOutputStream(CmdArgs.INSTANCE.out + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".bim")));

		} 
		catch (IOException e)
		{
			e.printStackTrace();
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
				e.printStackTrace();
			}
		}

		try
		{
			bedout.writeByte(byte1);
			bedout.writeByte(byte2);
			bedout.writeByte(byte3);
			for (int i = 0; i < NMarker; i++)
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
			bedout.close();
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < NMarker; i++)
		{
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (NMarker * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A[0] + " " + A[1]);
		}

		bim.close();
		fam.close();

	}

}
