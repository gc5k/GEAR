package regression;

import im.GenomeScan;
import im.IntervalPriorProbability;
import im.population.IMPopulation;
import publicAccess.ToolKit;

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealMatrixImpl;

public class Likelihood {
    IMPopulation imp;
    GenomeScan gs;
    int[] SNPIdx;
    int[][] ChrInt;

    public Likelihood(IMPopulation im, GenomeScan g, String com) {
    	imp = im;
    	gs = g;
    	SNPIdx = ToolKit.StringToIntArray(com);
    	if(SNPIdx.length > 1) {
    		System.err.println("one dimensional search is allowed only.");
    		System.exit(0);
    	}
    	ChrInt = ChrInt();
    }

    public double LogLikelihoodAlternativeCIM(LinearRegression lr, int interval, double[][] weight) {
    	double LODLikelihood = 0;
    	RealMatrix B = lr.estimate;
    	double[] v = new double[B.getRowDimension()];
    	for(int i = 0; i < v.length; i++) {
    		v[i] = B.getEntry(i, 0);
    	}
    	for (int i = 0; i < weight.length; i++) {
    		v[1+i] = 0;
    	}
    	RealMatrix BA = new RealMatrixImpl(v);
    	RealMatrix X = lr.X();
    	RealMatrix Y = lr.Y();
        RealMatrix fit = X.multiply(BA);
        RealMatrix Res1 = Y.subtract(fit);

    	double mse = lr.mse;
        IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
        for (int i = 0; i < SNPIdx.length; i++) {
            iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
        }

        for (int i = 0; i < Res1.getRowDimension(); i++) {
   			int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[0][0], ChrInt[0][1]);
   			double[] PriorProbability = new double[iip[0].NumQTLtypes()];
   			double LogInd = 0;
   			for (int j = 0; j < PriorProbability.length; j++) {
   				double r = Res1.getEntry(i, 0);
   				PriorProbability[j] = iip[0].PriorProbabilityAt(IIPIdx, interval, j);
   				for (int k = 0; k < weight.length; k++) {
   					r -= weight[k][j] * B.getEntry(1+k, 0);
   				}
   				LogInd += GaussPDF(r, 0, mse) * PriorProbability[j];
   			}
    		LODLikelihood += Math.log10(LogInd);
    	}
    	return LODLikelihood;
    }

    public double LogLikelihoodNullCIM(LinearRegression lr) {
    	double LODLikelihood = 0;
    	RealMatrix B = lr.estimate;
    	RealMatrix X = lr.X();
    	RealMatrix Y = lr.Y();
    	RealMatrix residual = lr.getResidual();
    	double mse = lr.mse;
   		for (int i = 0; i < residual.getRowDimension(); i++) {
   			double res = residual.getEntry(i, 0);
       		double LogInd = GaussPDF(res, 0, mse); 
    		LODLikelihood += Math.log10(LogInd);
    	}
    	return LODLikelihood;
    }

    public double LogLikelihoodAlternativeIM(LinearRegression lr, int interval, double[][] weight) {
    	double LODLikelihood = 0;
    	RealMatrix B = lr.estimate;
    	RealMatrix Y = lr.Y();
    	double mse = lr.mse;
        IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
        for (int i = 0; i < SNPIdx.length; i++) {
            iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
        }

        for (int i = 0; i < Y.getRowDimension(); i++) {
   			int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[0][0], ChrInt[0][1]);
   			double[] PriorProbability = new double[iip[0].NumQTLtypes()];
   			for (int j = 0; j < PriorProbability.length; j++) {
   				PriorProbability[j] = iip[0].PriorProbabilityAt(IIPIdx, interval, j);
   			}
   			double LogInd = 0;
   			for (int j = 0; j < PriorProbability.length; j++) {
   				double r = Y.getEntry(i, 0) - B.getEntry(0, 0);
   				for (int k = 0; k < weight.length; k++) {
   					r -= weight[k][j] * B.getEntry(1+k, 0);
   				}
   				LogInd += GaussPDF(r, 0, mse) * PriorProbability[j];
   			}
    		LODLikelihood += Math.log10(LogInd);
    	}
    	return LODLikelihood;
    }

    public double LogLikelihoodNullIM(LinearRegression lr) {
    	double LODLikelihood = 0;
    	RealMatrix B = lr.estimate;
    	RealMatrix X = lr.X();
    	RealMatrix Y = lr.Y();
    	RealMatrix residual = lr.getResidual();
    	double mse = lr.mse;
   		for (int i = 0; i < residual.getRowDimension(); i++) {
   			double res = residual.getEntry(i, 0);
       		double LogInd = GaussPDF(res, 0, mse); 
    		LODLikelihood += Math.log10(LogInd);
    	}
    	return LODLikelihood;
    }
    
    public static double GaussPDF(double mu1, double mu0, double var) {
    	double pdf = 0;
    	pdf = (1/Math.sqrt(2*Math.PI*var)) * Math.exp(-(mu1-mu0)*(mu1-mu0)/(2*var));
    	return pdf;
    }

    public int[][] ChrInt() {
    	int[][] chrint = new int[SNPIdx.length][2];
        int c = 0;
        int idx = 0;
        for (int i = 0; i < imp.ChromosomeNumber(); i++) {
            for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                if (c == SNPIdx[idx]) {
                    chrint[idx][0] = i;
                    chrint[idx][1] = j;
                    idx++;
                    if (idx == SNPIdx.length) {
                        break;
                    }
                }
                c++;
            }
            if (idx == SNPIdx.length) {
                break;
            }
        }
        return chrint;
    }    
}
