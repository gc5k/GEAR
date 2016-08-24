package gear.subcommands.qpca;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.data.SubjectID;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;

public class QPCACommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		qpcaArgs = (QPCACommandArguments)cmdArgs;
		data = new InputDataSet2();
		data.addFile(qpcaArgs.getGrmID());
		if (qpcaArgs.getKeepFile() != null)
		{
			data.addFile(qpcaArgs.getKeepFile());
		}
		data.LineUpFiles();

		readGrm();
		
		EigenAnalysis();
	}

	private void readGrm()
	{
		if (qpcaArgs.getGrmBin() != null)
		{
			readGrmBin(qpcaArgs.getGrmBin());
		}
		else
		{
			BufferedReader reader = qpcaArgs.getGrmText() == null ?
					BufferedReader.openGZipFile(qpcaArgs.getGrmGZ(), "GRM (gzip)") :
					BufferedReader.openTextFile(qpcaArgs.getGrmText(), "GRM");
			readGrm(reader);
		}
	}

	private void readGrmBin(String fileName)
	{
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /*littleEndian*/true);
		double[][] grmMat = new double[data.getFileSampleSize(0)][data.getFileSampleSize(0)];
		Logger.printUserLog("Constructing A matrix: a " + data.getFileSampleSize(0) + " X " + data.getFileSampleSize(0) + " matrix.");
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
		grmBin.close();
		
		A = lineUpMatrix(grmMat, data.getMatchedSubjectIdx(0));
	}

	private void readGrm(BufferedReader reader)
	{
		double[][] grmMat = new double[data.getFileSampleSize(0)][data.getFileSampleSize(0)];
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
		
		A = lineUpMatrix(grmMat, data.getMatchedSubjectIdx(0));
	}

	private double[][] lineUpMatrix(double[][] B, int[] subIdx)
	{
		double[][] a = new double[subIdx.length][subIdx.length];

		for (int i = 0; i < subIdx.length; i++)
		{
			for (int j = 0; j < subIdx.length; j++)
			{
				a[i][j] = a[j][i] = B[subIdx[i]][subIdx[j]];
			}
		}
		return a;
	}

	private void EigenAnalysis()
	{
		DecimalFormat fmt = new DecimalFormat("0.0000");
		DecimalFormat fmtp = new DecimalFormat("0.000E000");
		PrintStream grmWriter = FileUtil.CreatePrintStream(new String(qpcaArgs.getOutRoot() + ".crm"));
		Array2DRowRealMatrix rm = new Array2DRowRealMatrix(A);

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
			
			if((i+1) <= qpcaArgs.getEV())
			{
				Logger.printUserLog("The " + (i+1) + "th eigenvalue is " + eR[i] + "."); 
			}
		}
		evaWriter.close();

		PrintStream eveWriter = FileUtil.CreatePrintStream(new String(qpcaArgs.getOutRoot() + ".eigenvec"));

		ArrayList<SubjectID> SID = data.getMatchedSubjectID(0);
		
		for (int i = 0; i < evMat.getRowDimension(); i++)
		{
			eveWriter.print(SID.get(i).toString() + "\t");
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

	private double[][] A;
//	ArrayList<ArrayList<String>> famID = NewIt.newArrayList();
	private QPCACommandArguments qpcaArgs;
	private InputDataSet2 data;
}
