package gear.subcommands.simulationfammpheno;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.pop.PopStat;
import gear.util.BufferedReader;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.NonSquareMatrixException;
import org.apache.commons.math.linear.NotPositiveDefiniteMatrixException;
import org.apache.commons.math.linear.NotSymmetricMatrixException;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

public final class SimuFamilyMPCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		famMPArgs = (SimuFamilyMPCommandArguments) cmdArgs;

		init();
		try
		{
			ibdF = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".ibdo")));
			ibdSibF = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".ibd")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .ibdo file.");
			Logger.handleException(e, "An I/O exception occurred when creating the .ibd file.");
		}

		for (int i = 0; i < famMPArgs.getNumberOfFamilies(); i++)
		{
			generateNuclearFamily(NKid[i], NAffKid[i], i);
		}
		ibdF.close();
		ibdSibF.close();

		writePheno();

		if (famMPArgs.getMakeBed())
		{
			writeBFile();
		}
		else
		{
			writeFile();
		}
	}

	private void init()
	{
		rnd = new RandomDataImpl();
		rnd.reSeed(famMPArgs.getSeed());

		M = famMPArgs.getNumberOfMarkers();
		nullM = famMPArgs.getNullMarker();

		NKid = new int[famMPArgs.getNumberOfFamilies()];
		Arrays.fill(NKid, 2);
		NAffKid = new int[famMPArgs.getNumberOfFamilies()];
		Arrays.fill(NAffKid, 1);

		hsq = famMPArgs.getHsq();

		readCorMatrix();
		getEffect();
		getMAF();
		getDPrime();
		calLD();
		getRec();

		gm = new int[famMPArgs.getNumberOfFamilies() * famSize][M];
		phe = new double[famMPArgs.getNumberOfFamilies() * famSize][hsq.length];
	}

	private void calLD()
	{
		LD = PopStat.CalcLDfromDPrime(maf, DPrime);
	}

	private void getDPrime()
	{
		DPrime = new double[M - 1];
		if (famMPArgs.isPlainLD())
		{
			Arrays.fill(DPrime, famMPArgs.getLD());
		}
		else if (famMPArgs.isRandLD())
		{
			for (int i = 0; i < DPrime.length; i++)
			{
				DPrime[i] = rnd.nextUniform(-1, 1);
			}
		}
	}

	private void getRec()
	{
		recSex = new double[M][2];
		if (famMPArgs.isPlainRec())
		{
			for (int i = 0; i < recSex.length; i++)
			{
				recSex[i][0] = famMPArgs.getRec();
				recSex[i][1] = famMPArgs.getRec();
			}
		}
		else if (famMPArgs.isSexRec())
		{
			double[] rs = famMPArgs.getRecSex();
			for (int i = 0; i < recSex.length; i++)
			{
				recSex[i][0] = rs[0];
				recSex[i][1] = rs[1];
			}
		}
		else if (famMPArgs.isRandRec())
		{
			for (int i = 0; i < recSex.length; i++)
			{
				recSex[i][0] = recSex[i][1] = rnd.nextUniform(0.01, 0.5);
			}
		}
		recSex[0][0] = 0.5;
		recSex[0][1] = 0.5;	
	}

	private void getMAF()
	{
		maf = new double[M];
		if (famMPArgs.isPlainMAF())
		{
			Arrays.fill(maf, famMPArgs.getMAF());
		}
		else if (famMPArgs.isUnifMAF())
		{
			for (int i = 0; i < maf.length; i++)
			{
				maf[i] = rnd.nextUniform(0.01, 0.5);
			}
		}
	}
	
	private void getEffect()
	{
		int[] idx = Sample.SampleIndex(0, M-1, M-nullM);
		Arrays.sort(idx);

		effect = new double[M][hsq.length];
		if (famMPArgs.isPolyEffectFile())
		{
			java.io.BufferedReader reader = FileUtil.FileOpen(famMPArgs.getPolyEffectFile());

			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if (c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic effects.  Ignore the rest of the content in '" + famMPArgs.getPolyEffectFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < hsq.length)
					{
						Logger.printUserError("Incorrect columns for genetic effect file '" + famMPArgs.getPolyEffectFile() + "' in line " + (c+1) + ".");
						Logger.printUserLog("GEAR quited");
						System.exit(0);
					}
					for (int i = 0; i < hsq.length; i++)
					{
						effect[c][i] = Double.parseDouble(l[i]);						
					}
					c++;
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the effect file '"
								+ famMPArgs.getPolyEffectFile() + "'.");
			}
		}
		else
		{
			double[][] z = new double[M][hsq.length];
			for (int i = 0; i < hsq.length; i++)
			{
				for (int j = 0; j < idx.length; j++)
				{
					z[idx[j]][i] = rnd.nextGaussian(0, Math.sqrt(1/( idx.length * 1.0)));
				}
			}

			for (int i = 0; i < hsq.length; i++)
			{
				for (int j = 0; j < idx.length; j++)
				{
					for (int k = 0; k < hsq.length; k++)
					{
						effect[idx[j]][i] += z[idx[j]][k] * cholC[i][k];						
					}
				}
			}
		}

		for (int i = 0; i < hsq.length; i++)
		{
			for (int j = 0; j < effect.length; j++)
			{
				if (hsq[i] == 0)
				{
					effect[j][i] = 0;
				}
			}
		}
		
		
		PrintWriter addEff = null;
		try
		{
			addEff = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".rnd")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .rnd file.");
		}

		for (int i = 0; i < M; i++)
		{
			addEff.print("rs" + i + " " + A[0]);

			for (int j = 0; j < hsq.length; j++)
			{
				addEff.print(" " + effect[i][j]);
			}
			addEff.println();
		}

		addEff.close();
	}

	private void readCorMatrix()
	{
		corMat = new double[hsq.length][hsq.length];
		corEMat = new double[hsq.length][hsq.length];
		for(int i = 0; i < corMat.length; i++)
		{
			corMat[i][i] = 1;
			corEMat[i][i] = 1;
		}

		if (famMPArgs.getCMFile() != null)
		{
			BufferedReader bf = BufferedReader.openTextFile(famMPArgs.getCMFile(), "cm file.");
			Logger.printUserLog("Reading '" + famMPArgs.getCMFile() + "'.");

			String[] d = null;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != hsq.length)
				{
					Logger.printUserError("incorrect '" + famMPArgs.getCMFile() + "'.");
					System.exit(0);
				}
				int c = 0;
				for (int i = 0; i < d.length; i++)
				{
					corMat[cIdx][c++] = Double.parseDouble(d[i]);
				}
				cIdx++;
			}
			Logger.printUserLog(corMat.length + "X" + corMat.length + " correlation matrix has been read in.");	
		}

		if (famMPArgs.getCMEFile() != null)
		{
			BufferedReader bf = BufferedReader.openTextFile(famMPArgs.getCMEFile(), "cm file.");
			Logger.printUserLog("Reading '" + famMPArgs.getCMEFile() + "'.");

			String[] d = null;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != hsq.length)
				{
					Logger.printUserError("incorrect '" + famMPArgs.getCMEFile() + "'.");
					System.exit(0);
				}
				int c = 0;
				for (int i = 0; i < d.length; i++)
				{
					corEMat[cIdx][c++] = Double.parseDouble(d[i]);
				}
				cIdx++;
			}
			Logger.printUserLog(corMat.length + "X" + corMat.length + " environmental correlation matrix has been read in.");	
		}

		CholeskyDecompositionImpl cholCImpl = null;
		CholeskyDecompositionImpl cholCEImpl = null;

		try {
			cholCImpl = new CholeskyDecompositionImpl(new Array2DRowRealMatrix(corMat));
			cholCEImpl = new CholeskyDecompositionImpl(new Array2DRowRealMatrix(corEMat));
		} catch (NonSquareMatrixException e) {
			e.printStackTrace();
		} catch (NotSymmetricMatrixException e) {
			e.printStackTrace();
		} catch (NotPositiveDefiniteMatrixException e) {
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
		} catch (NullPointerException e) {
			e.printStackTrace();
		}
		RealMatrix L = cholCImpl.getL();
		cholC = L.getData();

		RealMatrix L1 = cholCEImpl.getL();
		cholCE = L1.getData();
	}

	private void generateNuclearFamily(int nkid, int affKid, int famIdx)
	{
		int[][] p = sampleChromosome(famIdx, 0);
		int[][] m = sampleChromosome(famIdx, 1);

			int[][] rc1 = generateBaby(p, m, famIdx, 2);
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < rc1.length; k++)
				{
					ibdF.print(rc1[k][j] + " ");
				}
				ibdF.println();
			}

			int[][] rc2 = generateBaby(p, m, famIdx, 3);
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < rc2.length; k++)
				{
					ibdF.print(rc2[k][j] + " ");
				}
				ibdF.println();
			}
		
		int[] ibd = new int[maf.length];
		for (int i = 0; i < maf.length; i++)
		{
			ibd[i] = 2 - Math.abs(rc1[i][0] - rc2[i][0]) - Math.abs(rc1[i][1] - rc2[i][1]);
			ibdSibF.print(ibd[i] + " ");
		}
		ibdSibF.println();
	}

	private int[][] sampleChromosome(int famIdx, int shift)
	{
		int[][] v = new int[maf.length][2];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < maf.length; j++)
			{
				double r = rnd.nextUniform(0, 1);
				if (j == 0)
				{
					v[j][i] = r < maf[j] ? 0 : 1;
				}
				else
				{
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[j - 1][i];
					double f1 = a == 0 ? maf[j - 1] : (1 - maf[j - 1]);
					double f2 = a == 0 ? maf[j] : (1 - maf[j]);
					v[j][i] = d < (f1 * f2 + LD[j - 1]) / f1 ? v[j - 1][i]
							: (1 - v[j - 1][i]);
				}
			}
		}

		for (int i = 0; i < maf.length; i++)
		{
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}

		for (int i = 0; i < maf.length; i++)
		{
			for (int j = 0; j < hsq.length; j++)
			{
				phe[famIdx * famSize + shift][j] += gm[famIdx * famSize + shift][i] * effect[i][j];				
			}
		}
		return v;
	}

	private int[][] generateBaby(int[][] p, int[][] m, int famIdx, int shift)
	{
		int[][] v = new int[maf.length][2];
		int[][] rc = new int[maf.length][2];

		for (int i = 0; i < 2; i++)
		{
			int[][] chr = i == 0 ? p : m;
			int idx = 1;
			try
			{
				idx = rnd.nextBinomial(1, recSex[0][i]);
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}

			for (int j = 0; j < maf.length; j++)
			{
				double r = rnd.nextUniform(0, 1);
				idx = r < recSex[j][i] ? 1 - idx : idx;
				rc[j][i] = idx;
				v[j][i] = chr[j][idx];
			}
		}

		for (int i = 0; i < maf.length; i++)
		{
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}

		for (int i = 0; i < maf.length; i++)
		{
			for (int j = 0; j < hsq.length; j++)
			{
				phe[famIdx * famSize + shift][j] += gm[famIdx * famSize + shift][i] * effect[i][j];				
			}
		}

		return rc;
	}

	public void writePheno()
	{
		PrintWriter pheno = null;

		for (int i = 0; i < hsq.length; i++)
		{
			try
			{
				pheno = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() +"." + (i+1) + ".phe")));
			}
			catch (IOException e)
			{
				Logger.handleException(e, "An I/O exception occurred when creating the .phe file.");
			}

			double[] p = new double[famMPArgs.getNumberOfFamilies() * 2];
			int cn=0;

			for (int h = 0; h < famMPArgs.getNumberOfFamilies(); h++)
			{
				for (int j = 0; j < famSize; j++)
				{
					if(j < 2)
					{
						p[cn++] = phe[h*famSize+j][i];
					}
				}
			}
			double vg = StatUtils.variance(p);
			Logger.printUserLog("Genetic variation is "+vg);
			double ve = hsq[i] == 0 ? 1:vg * (1 - hsq[i]) / hsq[i];

			for (int h = 0; h < famMPArgs.getNumberOfFamilies(); h++)
			{
				int fid = (h + 1) * 10000;
				
				for (int j = 0; j < famSize; j++)
				{
					double pv = phe[h*famSize+j][i] + rnd.nextGaussian(0, Math.sqrt(ve));
					int pid = fid + 1 + j;

					pheno.print(fid + " ");
					pheno.print(pid + " ");
					pheno.println(pv);
				}
			}

			pheno.close();
		}
	}

	public void writeFile()
	{
		PrintWriter ped = null;
		PrintWriter map = null;

		try
		{
			ped = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".map")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .ped and .map files.");
		}

		for (int h = 0; h < famMPArgs.getNumberOfFamilies(); h++)
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
				switch (gm[h * famSize][i])
				{
				case 0:
					sb.append(A[0] + " " + A[0]);
					break;
				case 1:
					sb.append(A[0] + " " + A[1]);
					break;
				case 2:
					sb.append(A[1] + " " + A[1]);
					break;
				default:
					break;
				}
				if (i == (maf.length - 1))
				{
					ped.print(sb.toString());			
				}
				else
				{
					ped.print(sb.toString() + " ");
				}
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
				switch (gm[h * famSize + 1][i])
				{
				case 0:
					sb.append(A[0] + " " + A[0]);
					break;
				case 1:
					sb.append(A[0] + " " + A[1]);
					break;
				case 2:
					sb.append(A[1] + " " + A[1]);
					break;
				default:
					break;
				}
				if (i == (maf.length - 1))
				{
					ped.print(sb.toString());						
				}
				else
				{
					ped.print(sb.toString() + " ");
				}
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
					switch (gm[h * famSize + 2 + j][i])
					{
					case 0:
						sb.append(A[0] + " " + A[0]);
						break;
					case 1:
						sb.append(A[0] + " " + A[1]);
						break;
					case 2:
						sb.append(A[1] + " " + A[1]);
						break;
					default:
						break;
					}
					if (i == (maf.length - 1))
					{
						ped.print(sb.toString());						
					}
					else
					{
						ped.print(sb.toString() + " ");
					}
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
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(famMPArgs.getOutRoot() + ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(famMPArgs.getOutRoot() + ".bim")));

		} 
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, .fam and .bim files.");
		}

		int fid = 0;
		int fc = 0;
		int cn = 0;
		int pid = 0;
		int mid = 0;
		for (int i = 0; i < gm.length; i++)
		{
			if (i % 4 == 0)
			{
				fc++;
				fid = fc * 10000;
				pid = fid + 1;
				mid = fid + 2;
				cn = 0;
			}
			cn++;
			fam.print(fid + " ");
			if (cn <= 2) 
			{
				fam.print((fid+cn) + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
			}
			else
			{
				fam.print((fid+cn) + " ");
				fam.print(pid + " ");
				fam.print(mid + " ");				
			}
			try
			{
				if (cn <= 2)
				{
					fam.print(cn + " ");
				}
				else 
				{
					fam.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				}
				fam.print((rnd.nextBinomial(1, 0.5) + 1) + "\n");

			}
			catch (MathException e)
			{
				Logger.handleException(e, "Failed to generate the random values.");
			}
		}

		try
		{
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < M; i++)
			{
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < gm.length; j++)
				{
					int g = (int) gm[j][i];
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
			Logger.handleException(e, "An I/O exception occurred when writing the .bed file.");
		}

		for (int i = 0; i < M; i++)
		{
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (M * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A[0] + " " + A[1]);
		}

		bim.close();
		fam.close();

	}

	private SimuFamilyMPCommandArguments famMPArgs;

	private double[][] corEMat;
	private double[][] cholCE;
	private double[][] corMat;
	private double[][] cholC;

	private RandomDataImpl rnd;
	private final String[] A = { "A", "C" };
	private int[] NKid = null;
	private int[] NAffKid = null;
	private double[] LD = null;
//	private double[] rec = null;
	private double[][] recSex = null;
	private double[] maf = null;
	private double[] DPrime = null;

	private int[][] gm = null;
	private final int famSize = 4;
	private double[][] phe = null;
//	private double[] h2;

	private int M;
	private int nullM;

	private double[] hsq;
	private double[][] effect;

	PrintWriter ibdF = null;
	PrintWriter ibdSibF = null;
	
}
