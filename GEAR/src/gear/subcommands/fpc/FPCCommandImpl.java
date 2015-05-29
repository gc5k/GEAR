package gear.subcommands.fpc;

import java.io.PrintStream;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;

public class FPCCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		FPCArgs = (FPCCommandArguments) cmdArgs;
//		batch = FPCArgs.getMetaBatch();
		refIdx = FPCArgs.getReference();
		co = FPCArgs.getCoordinates();
		System.out.println(co[0][0] + " " + co[0][1]);
		System.out.println(co[1][0] + " " + co[1][1]);
		System.out.println(co[2][0] + " " + co[2][1]);

		readFst();
		FstCartographer();
		writeOut();
	}

	private void FstCartographer()
	{
		Fpc = new double[fstMat.length][2];
		for(int i = 0; i < fstMat.length; i++)
		{
			
			if(i == refIdx[0])
			{
					Fpc[i][0] = co[0][0]; Fpc[i][1] = co[0][1];
			} 
			else if (i == refIdx[1] )
			{
					Fpc[i][0] = co[1][0]; Fpc[i][1] = co[1][1];
			} 
			else if (i == refIdx[2] )
			{
					Fpc[i][0] = co[2][0]; Fpc[i][1] = co[2][1];
			}
			else
			{
				double f0 = fstMat[i][refIdx[0]];
				double f1 = fstMat[i][refIdx[1]];
				double f2 = fstMat[i][refIdx[2]];

				double A1 = (co[1][0] - co[0][0]) * (f0/(f0+f1)) + co[0][0];
				double A2 = (co[1][1] - co[0][1]) * (f0/(f0+f1)) + co[0][1];

				double B1 = (co[2][0] - co[0][0]) * (f0/(f0+f2)) + co[0][0];
				double B2 = (co[2][1] - co[0][1]) * (f0/(f0+f2)) + co[0][1];

				double C1 = (co[2][0] - co[1][0]) * (f1/(f1+f2)) + co[1][0];
				double C2 = (co[2][1] - co[1][1]) * (f1/(f1+f2)) + co[1][1];
				Fpc[i][0] = (A1 + B1 + C1)/3;
				Fpc[i][1] = (A2 + B2 + C2)/3;
				System.out.println("Cohort " + (1+i) + "'s Fst: " + f0 + " " + f1 + " " + f2 + " " + Fpc[i][0] + " " + Fpc[i][1]);
				System.out.println("Coordinates on edge 1-2: " + A1 + " " + A2);
				System.out.println("Coordinates on edge 1-3: " + B1 + " " + B2);
				System.out.println("Coordinates on edge 2-3: " + C1 + " " + C2);
				System.out.println("Gravity: " + Fpc[i][0] + " " + Fpc[i][1] + "\n");
			}
		}
	}

	private void readFst()
	{
		BufferedReader reader = BufferedReader.openTextFile(FPCArgs.getFstFile(), "Summary Statistic file");
		String[] tokens = reader.readTokens();
		int tokensLen;
		int cnt = 0;

		do
		{
			tokensLen = tokens.length;
			if (cnt == 0)
			{
				fstMat = new double[tokensLen][tokensLen];
			}
			
			for(int i = 0; i < tokensLen; i++)
			{
				fstMat[cnt][i] = Double.parseDouble(tokens[i]);
			}
			cnt++;
		} while( (tokens = reader.readTokens(fstMat.length)) != null && cnt < fstMat.length);

	}

	private void writeOut()
	{
		PrintStream FpcOut = FileUtil.CreatePrintStream(FPCArgs.getOutRoot() + ".fpc");
		FpcOut.println("Fpc1\tFpc2");
		for(int i = 0; i < Fpc.length; i++)
		{
			FpcOut.println(Fpc[i][0] + "\t" + Fpc[i][1]);
		}
		FpcOut.close();
	}

	FPCCommandArguments FPCArgs;
	private double[][] fstMat;
	private double[][] Fpc;
	private double[][] co;
//	private ArrayList<String> batch;
	private int[] refIdx;
}
