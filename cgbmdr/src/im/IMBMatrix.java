package im;

import java.util.ArrayList;

import im.population.IMPopulation;
import im.IntervalPriorProbability;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class IMBMatrix {

    IMPopulation imp;
    GenomeScan gs;
    int[] SNPIdx;
    int[][] ChrInt;
    ArrayList id;
    int order;

    public IMBMatrix(GenomeScan g, IMPopulation im) {
        gs = g;
        imp = im;
    }

    public double[][] getPPMatrix(String com, double[][] coeff) {
        SNPIdx = ToolKit.StringToIntArray(com);
        double[][] matrix = new double[imp.IndividualNumber()][SNPIdx.length * (coeff.length) + 1];
        ChrInt();
        IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
        for (int i = 0; i < SNPIdx.length; i++) {
            iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
        }
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            matrix[i][0] = 1;
            for (int j = 0; j < SNPIdx.length; j++) {
                int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[j][0], ChrInt[j][1]);
                for (int jj = 0; jj < coeff.length; jj++) {
                    for (int k = 0; k < iip[j].NumQTLtypes(); k++) {
                        matrix[i][1 + j * coeff.length +jj] += iip[j].MidPriorProbability(IIPIdx, k) * coeff[jj][k];
                    }
                }
            }
        }
//        for(int i = 0; i < matrix.length; i ++) {
//            System.out.print(imp.MarkerAt(i, ChrInt[0][0], ChrInt[0][1]) + " " + imp.MarkerAt(i, ChrInt[0][0], ChrInt[0][1]+1) + "\t");
//            for(int j = 0; j < matrix[i].length; j++) {
//                System.out.print(matrix[i][j]+"\t");
//            }
//            System.out.println(imp.PhenotypeAt(i, 0));
//        }
        return matrix;
    }

    public double[][] getPPMatrixAtPoint(String com, int w, double[][] coeff) {
        SNPIdx = ToolKit.StringToIntArray(com);
        double[][] matrix = new double[imp.IndividualNumber()][SNPIdx.length * (coeff.length) + 1];
        ChrInt();
        IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
        for (int i = 0; i < SNPIdx.length; i++) {
            iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
        }
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            matrix[i][0] = 1;
            for (int j = 0; j < SNPIdx.length; j++) {
                int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[j][0], ChrInt[j][1]);
                for (int jj = 0; jj < coeff.length; jj++) {
                    for (int k = 0; k < iip[j].NumQTLtypes(); k++) {
                        matrix[i][1 + j * coeff.length +jj] += iip[j].PriorProbabilityAt(IIPIdx, w, k) * coeff[jj][k];
                    }
                }
            }
        }
//        for(int i = 0; i < matrix.length; i ++) {
//            System.out.print(imp.MarkerAt(i, ChrInt[0][0], ChrInt[0][1]) + " " + imp.MarkerAt(i, ChrInt[0][0], ChrInt[0][1]+1) + "\t");
//            for(int j = 0; j < matrix[i].length; j++) {
//                System.out.print(matrix[i][j]+"\t");
//            }
//            System.out.println(imp.PhenotypeAt(i, 0));
//        }
        return matrix;
    }

    public void setOrder(int o) {
        order = o;
        ChrInt = new int[order][2];
    }

    public String getString() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < ChrInt.length; i++) {
            sb.append(ChrInt[i][0] + "," + ChrInt[i][1] + "\t");
        }
        return sb.toString();
    }

    private void ChrInt() {
        int c = 0;
        int idx = 0;
        for (int i = 0; i < imp.ChromosomeNumber(); i++) {
            for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                if (c == SNPIdx[idx]) {
                    ChrInt[idx][0] = i;
                    ChrInt[idx][1] = j;
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
    }
}
