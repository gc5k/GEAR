package he;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import he.endian.LittleEndianDataInputStream;

public class HEMRead
{
	private final String delim = "\\s+";
	protected boolean[] flag;
	HashMap<String, Integer> ID;
	protected String grmListFile;
	protected String grmFile;
	protected String grmID;
	protected String keepFile;
	protected String phenoFile;
	protected String output;
	protected int[] mpheno;
	protected int perm;
	protected boolean permFlag = false;

	protected boolean reverse;
	protected boolean k_button;
	protected double k;
	protected gear.HEType heType;

	protected HashMap<String, Integer> ID2Idx;

	protected double yyProd;

	protected double[][] XtX;
	protected double[] XtY;
	protected double[][] y;
	protected int dim;
	protected double P;
	protected boolean isCC = false;

	protected ArrayList<BufferedReader> grmList;
	protected ArrayList<LittleEndianDataInputStream> binList;
	protected ArrayList<String> grmFileList;
	protected ArrayList<String> idFileList;
	protected HashMap<Double, Integer> cat = new HashMap<Double, Integer>();

	protected boolean isSingleGrm;
	protected boolean isBinGrm = false;

	protected long nRec = 0;
	protected Lambda lambda;
	StringBuffer sb = new StringBuffer();

	public HEMRead()
	{
		if (CmdArgs.INSTANCE.getHEArgs().isSingleGrm())
		{
			grmFile = CmdArgs.INSTANCE.getHEArgs().getGrm();
			grmID = CmdArgs.INSTANCE.getHEArgs().getGrmId();
			existGrm();
			isSingleGrm = true;
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isMultiGrm())
		{
			grmListFile = CmdArgs.INSTANCE.getHEArgs().getGrmList();
			existGrmList();
			isSingleGrm = false;
		}

		keepFile = CmdArgs.INSTANCE.keepFile;
		phenoFile = CmdArgs.INSTANCE.getHEArgs().getPheno();
		mpheno = CmdArgs.INSTANCE.getHEArgs().getMPheno();
		reverse = CmdArgs.INSTANCE.reverse;
		k_button = CmdArgs.INSTANCE.k_button;
		k = CmdArgs.INSTANCE.k;
		output = CmdArgs.INSTANCE.out;
		heType = CmdArgs.INSTANCE.getHEArgs().getType();
		permFlag = CmdArgs.INSTANCE.permFlag;
		perm = CmdArgs.INSTANCE.perm;

		XtX = new double[mpheno.length + 1][mpheno.length + 1];
		XtY = new double[mpheno.length + 1];

		String line;

		// *********************************** read grm id
		BufferedReader reader = FileProcessor.FileOpen(grmID);

		int i2 = 0;
		ID2Idx = new HashMap<String, Integer>();
		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				ID2Idx.put(sb.toString(), i2++);
			}
			reader.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}

		flag = new boolean[ID2Idx.size()];
		Arrays.fill(flag, false);

		nRec = flag.length * (flag.length + 1) / 2;

		// *********************************** read pheno file
		reader = FileProcessor.FileOpen(phenoFile);

		y = new double[flag.length][mpheno.length + 1];
		HashMap<Double, Integer> cat = new HashMap<Double, Integer>();

		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				if (ID2Idx.containsKey(sb.toString()))
				{
					int ii = ID2Idx.get(sb.toString());
					boolean f = true;
					y[ii][0] = 1;
					for (int j = 0; j < mpheno.length; j++)
					{
						if (ConstValues.isNA(s[1 + mpheno[j]]))
						{
							f = false;
							break;
						}
						else
						{
							y[ii][j + 1] = Double.parseDouble(s[1 + mpheno[j]]);
						}
					}
					flag[ii] = f;
					
					if (cat.containsKey(y[ii][1]))
					{
						Integer I = (Integer) cat.get(y[ii][1]);
						I++;
						cat.put(y[ii][1], I);
					}
					else
					{
						cat.put(y[ii][1], 1);
					}
				}
			}
			reader.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "Error in reading phenotype file '" + phenoFile + "'.");
		}

		// ************************keep
		if (keepFile != null)
		{
			reader = FileProcessor.FileOpen(keepFile);
			boolean[] ff = new boolean[flag.length];
			Arrays.fill(ff, false);
			try
			{
				while ((line = reader.readLine()) != null)
				{
					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder(s[0] + "." + s[1]);
					if (ID2Idx.containsKey(sb.toString()))
					{
						int ii = ID2Idx.get(sb.toString());
						ff[ii] = true;
					}
				}
				reader.close();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}

			for (int i = 0; i < ff.length; i++)
			{
				flag[i] &= ff[i];
			}
		}

		int Len = 0;
		for (int i = 0; i < flag.length; i++)
			if (flag[i])
				Len++;
		dim = Len * (Len - 1) / 2;

		// ************************************standardising
		double[] ss = new double[y[0].length - 1];
		double[] ssx = new double[y[0].length - 1];
		for (int i = 0; i < flag.length; i++)
		{
			if (!flag[i])
				continue;
			for (int j = 0; j < ss.length; j++)
			{
				ss[j] += y[i][j + 1];
				ssx[j] += y[i][j + 1] * y[i][j + 1];
			}
		}
		double[] sd = new double[ssx.length];
		for (int i = 0; i < sd.length; i++)
		{
			ss[i] /= Len;
			sd[i] = Math.sqrt((ssx[i] - Len * ss[i] * ss[i]) / (Len - 1));
		}

		if (CmdArgs.INSTANCE.scale)
		{
			Logger.printUserLog("Standardising phentoype.");
			for (int i = 0; i < flag.length; i++)
			{
				if (!flag[i])
					continue;
				for (int j = 1; j < y[i].length; j++)
				{
					y[i][j] = (y[i][j] - ss[j - 1]) / sd[j - 1];
				}
			}
		}
	}

	private void existGrm()
	{
		grmList = NewIt.newArrayList();
		grmFileList = NewIt.newArrayList();
		idFileList = NewIt.newArrayList();
		binList = NewIt.newArrayList();

		if (CmdArgs.INSTANCE.getHEArgs().isGrmTxt())
		{
			FileProcessor.exists(grmFile);
			FileProcessor.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			grmList.add(FileProcessor.FileOpen(grmFile
					.toString()));
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isGrm())
		{
			FileProcessor.exists(grmFile);
			FileProcessor.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			FileInputStream fin = null;
			try
			{
				fin = new FileInputStream(grmFile);
			}
			catch (FileNotFoundException e1)
			{
				Logger.handleException(e1,
										"Error in opening '" + grmFile + "'.");
			}
			GZIPInputStream gzis = null;
			try
			{
				gzis = new GZIPInputStream(fin);
			}
			catch (IOException e1)
			{
				Logger.handleException(e1,
										"Error in opening gz file '" + grmFile + "'.");
			}
			InputStreamReader xover = new InputStreamReader(gzis);

			BufferedReader grmFile = new BufferedReader(xover);
			grmList.add(grmFile);
		}
		else
		{
			FileProcessor.exists(grmFile);
			FileProcessor.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			FileInputStream fileStream = null;
			try
			{
				fileStream = new FileInputStream(grmFile);
			}
			catch (FileNotFoundException e)
			{
				Logger.handleException(	e,
										"Error in opening GRM bin file '" + grmFile + "'.");
			}
			DataInputStream bigEndianDataStream = new DataInputStream(
					fileStream);
			LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(
					bigEndianDataStream, Float.SIZE);
			binList.add(littleEndianDataStream);
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrm() || CmdArgs.INSTANCE
				.getHEArgs().isGrmTxt())
		{
			if (grmList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmBinary())
		{
			if (binList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
			isBinGrm = true;
		}
	}

	private void existGrmList()
	{
		grmList = NewIt.newArrayList();
		grmFileList = NewIt.newArrayList();
		idFileList = NewIt.newArrayList();
		binList = NewIt.newArrayList();

		BufferedReader listReader = FileProcessor.FileOpen(grmListFile);
		String line = null;
		try
		{
			while ((line = listReader.readLine()) != null)
			{
				String root = line.trim();
				StringBuilder tsb = new StringBuilder(root);

				if (CmdArgs.INSTANCE.getHEArgs().isGrmTxtList())
				{
					FileProcessor.exists(tsb.append(".grm.txt").toString());
					tsb = new StringBuilder(root);
					FileProcessor.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.txt").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmList.add(FileProcessor.FileOpen(tsb.append(".grm.txt")
							.toString()));
				}
				else if (CmdArgs.INSTANCE.getHEArgs().isGrmList())
				{
					FileProcessor.exists(tsb.append(".grm.gz").toString());
					tsb = new StringBuilder(root);
					FileProcessor.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.gz").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					FileInputStream fin = null;
					try
					{
						tsb = new StringBuilder(root);
						fin = new FileInputStream(tsb.append(".grm.gz")
								.toString());
					}
					catch (FileNotFoundException e1)
					{
						Logger.handleException(	e1,
												"Error in opening '" + line
														.trim() + "' in " + CmdArgs.INSTANCE
														.getHEArgs()
														.getGrmList() + "'.");
					}
					GZIPInputStream gzis = null;
					try
					{
						gzis = new GZIPInputStream(fin);
					}
					catch (IOException e1)
					{
						Logger.handleException(	e1,
												"Error in opening gz file '" + line
														.trim() + "' in " + CmdArgs.INSTANCE
														.getHEArgs()
														.getGrmList() + "'.");
					}
					InputStreamReader xover = new InputStreamReader(gzis);

					BufferedReader grmFile = new BufferedReader(xover);
					grmList.add(grmFile);
				}
				else
				{
					tsb = new StringBuilder(root);
					FileProcessor.exists(tsb.append(".grm.bin").toString());
					tsb = new StringBuilder(root);
					FileProcessor.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.bin").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					FileInputStream fileStream = null;
					try
					{
						tsb = new StringBuilder(root);
						fileStream = new FileInputStream(tsb.append(".grm.bin")
								.toString());
					}
					catch (FileNotFoundException e)
					{
						Logger.handleException(	e,
												"Error in opening GRM bin file '" + line
														.trim() + "' in " + CmdArgs.INSTANCE
														.getHEArgs()
														.getGrmList() + "'.");
					}
					DataInputStream bigEndianDataStream = new DataInputStream(
							fileStream);
					LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(
							bigEndianDataStream, Float.SIZE);
					binList.add(littleEndianDataStream);
				}
			}
		}
		catch (IOException e)
		{
			Logger.handleException(	e,
									"An exception occurred when reading the GRM file '" + CmdArgs.INSTANCE
											.getHEArgs().getGrmList() + "'.");
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmList() || CmdArgs.INSTANCE
				.getHEArgs().isGrmTxtList())
		{
			if (grmList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmBinaryList())
		{
			if (binList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
			isBinGrm = true;
		}
		grmID = idFileList.get(0);
	}
}
