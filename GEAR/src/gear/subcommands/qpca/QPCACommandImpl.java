package gear.subcommands.qpca;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;

import gear.ConstValues;
import gear.data.InputDataSet;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class QPCACommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		qpcaArgs = (QPCACommandArguments)cmdArgs;
		InputDataSet data = new InputDataSet();
		data.readSubjectIDFile(qpcaArgs.getGrmID());
		readFam();
		readGrm(data.getNumberOfSubjects());
		EigenAnalysis();
	}

	private void readGrm(int numSubjects)
	{
		if (qpcaArgs.getGrmBin() != null)
		{
			readGrmBin(qpcaArgs.getGrmBin(), numSubjects);
		}
		else
		{
			BufferedReader reader = qpcaArgs.getGrmText() == null ?
					BufferedReader.openGZipFile(qpcaArgs.getGrmGZ(), "GRM (gzip)") :
					BufferedReader.openTextFile(qpcaArgs.getGrmText(), "GRM");
			readGrm(reader, numSubjects);
		}
	}

	private void readGrmBin(String fileName, int numSubjects)
	{
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /*littleEndian*/true);
		grmMat = new double[numSubjects][numSubjects];
		Logger.printUserLog("Constructing A matrix: a " + numSubjects + " X " + numSubjects + " matrix.");
		for (int i = 0; i < grmMat.length; i++) 
		{
			for (int j = 0; j <= i; j++)
			{
				if (grmBin.available() >= ConstValues.FLOAT_SIZE)
				{
					grmMat[i][j] = grmMat[j][i] = grmBin.readFloat();
				}
			}
		}
	}

	private void readGrm(BufferedReader reader, int numSubjects)
	{
		grmMat = new double[numSubjects][numSubjects];
		String[] tokens = null;
		for (int i = 0; i < grmMat.length; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if ((tokens = reader.readTokens(4)) != null)
				{
					grmMat[i][j] = grmMat[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
	}
	
	private void EigenAnalysis()
	{
		if(famID.size() != grmMat.length)
		{
			Logger.printUserLog("Inconsisitent sample size.\nGEAR quitted");
			System.exit(0);
		}
		DecimalFormat fmt = new DecimalFormat("0.0000");
		DecimalFormat fmtp = new DecimalFormat("0.000E000");
		PrintStream grmWriter = FileUtil.CreatePrintStream(new String(qpcaArgs.getOutRoot() + ".crm"));
		Array2DRowRealMatrix rm = new Array2DRowRealMatrix(grmMat);

		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			for (int j = 0; j < rm.getColumnDimension(); j++)
			{
				grmWriter.print(rm.getEntry(i, j) + " ");
			}
			grmWriter.println();
		}
		grmWriter.close();

		EigenDecompositionImpl ed = new EigenDecompositionImpl(rm.copy(), 1e-6);

		double[][] ev = new double[rm.getRowDimension()][rm.getColumnDimension()];
		for (int i = 0; i < rm.getColumnDimension(); i++)
		{
			ev[i] = ed.getEigenvector(i).toArray();
		}

		Array2DRowRealMatrix evM = new Array2DRowRealMatrix(ev);
		Array2DRowRealMatrix evMat = (Array2DRowRealMatrix) evM.transpose();

		PrintStream evaWriter = FileUtil.CreatePrintStream(new String(qpcaArgs.getOutRoot() + ".eigenval"));
		double[] eR=ed.getRealEigenvalues();

		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			if(Math.abs(eR[i]) >= 0.0001)
			{
				evaWriter.println(fmt.format(eR[i]));
			}
			else
			{
				evaWriter.println(fmtp.format(eR[i]));				
			}
		}
		evaWriter.close();

		PrintStream eveWriter = FileUtil.CreatePrintStream(new String(qpcaArgs.getOutRoot() + ".eigenvec"));

		for (int i = 0; i < evMat.getRowDimension(); i++)
		{
			eveWriter.print(famID.get(i).get(0) + "\t" + famID.get(i).get(1) + "\t");
			for (int j = 0; j < qpcaArgs.getEV(); j++)
			{
				if(Math.abs(evMat.getEntry(i, j)) >= 0.0001)
				{
					eveWriter.print(fmt.format(evMat.getEntry(i, j)) + "\t");					
				}
				else
				{
					eveWriter.print(fmtp.format(evMat.getEntry(i, j)) + "\t");					
				}
			}
			eveWriter.println();
		}
		eveWriter.close();
	}

	private void readFam()
	{
		BufferedReader reader = BufferedReader.openTextFile(qpcaArgs.getGrmID(), "fam");
		
		String[] tokens;
		while ((tokens = reader.readTokensAtLeast(2)) != null)
		{
			ArrayList<String> f = NewIt.newArrayList();
			f.add(tokens[0]);
			f.add(tokens[1]);
			famID.add(f);
		}
		reader.close();
	}

	private double[][] grmMat;
	ArrayList<ArrayList<String>> famID = NewIt.newArrayList();
	private QPCACommandArguments qpcaArgs;
}
