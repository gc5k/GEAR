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
		System.out.println("The coordinates for the three referen populations are: ");
		System.out.println("Ref pop1: " + co[0][0] + " " + co[0][1]);
		System.out.println("Ref pop2: " + co[1][0] + " " + co[1][1]);
		System.out.println("Ref pop3: " + co[2][0] + " " + co[2][1]);
		g = new double[2];
		g[0] = (co[0][0] + co[1][0] + co[2][0])/3;
		g[1] = (co[0][1] + co[1][1] + co[2][1])/3;
		System.out.println("Gravity: " + g[0] + " " + g[1]);
		readFst();
		FstCartographer();
		writeOut();
	}

	private void FstCartographer()
	{
		Fpc = new double[fstMat.length][3];
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
				System.out.println("Cohort " + (1+i) + "'s Fst: " + f0 + " " + f1 + " " + f2);
				System.out.println("Coordinates on edge 1-2: " + A1 + " " + A2);
				System.out.println("Coordinates on edge 1-3: " + B1 + " " + B2);
				System.out.println("Coordinates on edge 2-3: " + C1 + " " + C2);
				System.out.println("Gravity: " + Fpc[i][0] + " " + Fpc[i][1]);
				if (contains(co[0][0], co[0][1], co[1][0], co[1][1], g[0], g[1], Fpc[i][0], Fpc[i][1]))
				{
					Fpc[i][2] = 1;
				}
				else if (contains(co[0][0], co[0][1], co[2][0], co[2][1], g[0], g[1], Fpc[i][0], Fpc[i][1]))
				{
					Fpc[i][2] = 2;
				}
				else if (contains(co[1][0], co[1][1], co[2][0], co[2][1], g[0], g[1], Fpc[i][0], Fpc[i][1]))
				{
					Fpc[i][2] = 3;
				}
				System.out.println("Inside subspace: " + subspace[(int) Fpc[i][2]] + "\n");
			}
		}
	}

    public boolean  contains(double x1, double y1, double x2, double y2, double g1, double g2, double x, double y) 
    {

    	int npoints = 3;
    	double[] xpoints = new double[npoints];
    	double[] ypoints = new double[npoints];

    	xpoints[0] = x1;
    	ypoints[0] = y1;

    	xpoints[1] = x2;
    	ypoints[1] = y2;
    	
    	xpoints[2] = g1;
    	ypoints[2] = g2;

        int hits = 0;
        double lastx = xpoints[npoints - 1];
        double lasty = ypoints[npoints - 1];
        double curx, cury;

        // Walk the edges of the polygon
        for (int i = 0; i < npoints; lastx = curx, lasty = cury, i++) {
            curx = xpoints[i];
            cury = ypoints[i];

            if (cury == lasty) {
                continue;
            }

            double leftx;
            if (curx < lastx) {
                if (x >= lastx) {
                    continue;
                }

                leftx = curx;
            } else {
                if (x >= curx) {
                    continue;
                }
                leftx = lastx;
            }

            double test1, test2;
            if (cury < lasty) {
                if (y < cury || y >= lasty) {
                    continue;
                }

                if (x < leftx) {
                    hits++;
                    continue;

                }

                test1 = x - curx;
                test2 = y - cury;

            } else {
                if (y < lasty || y >= cury) {
                    continue;
                }

                if (x < leftx) {
                    hits++;
                    continue;
                }
                test1 = x - lastx;
                test2 = y - lasty;
            }
            if (test1 < (test2 / (lasty - cury) * (lastx - curx))) {

                hits++;
            }

        }
        return ((hits & 1) != 0);
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
		FpcOut.println("Fpc1\tFpc2\tSubspace");
		for(int i = 0; i < Fpc.length; i++)
		{
			FpcOut.println(Fpc[i][0] + "\t" + Fpc[i][1] + "\t" + subspace[(int) Fpc[i][2]]);
		}
		FpcOut.close();
	}

	FPCCommandArguments FPCArgs;
	private double[][] fstMat;
	private double[][] Fpc;
	private double[][] co;
	private double[] g;
	private String[] subspace={"Ref", "Ref:1-2", "Ref:1-3", "Ref:2-3"}; 
//	private ArrayList<String> batch;
	private int[] refIdx;
}
