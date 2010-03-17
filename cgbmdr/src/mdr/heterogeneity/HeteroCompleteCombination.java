/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package mdr.heterogeneity;

/**
 *
 * @author Guo-bo Chen
 */
public class HeteroCompleteCombination {
    private double[][][] statistics;
    private double[][] mean_statistics;
    private String name;                        // combination of loci

    HeteroCompleteCombination(int numTraits, int interval) {
        statistics = new double[numTraits][interval][2];
        mean_statistics = new double[numTraits][2];
    }

    public void SetStatistic(int traitidx, int interval, int statisticidx, double value) {
        statistics[traitidx][interval][statisticidx] = value;
    }

    public String GetModelName() {
        return name;
    }

    public void SetModelName(String n) {
        name = n;
    }

    public double get(int traitidx, int interval, int statisticidx) {
        return statistics[traitidx][interval][statisticidx];
    }

    public void summarise() {
        for(int i =0; i <statistics.length; i++) {
            for(int j = 0; j < statistics[i][0].length; j++) {
                for(int k =0; k < statistics[i].length; k++) {
                    mean_statistics[i][j] += statistics[i][k][j];
                }
                mean_statistics[i][j] /= statistics[i].length;
            }
        }
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(name + "\n");
        for(int i = 0; i < mean_statistics.length; i++) {
            for(int j = 0; j < mean_statistics[i].length; j++) {
                sb.append(mean_statistics[i][j] + "\t");
            }
        }
        sb.append("\n======================\n");
        for(int i = 0; i < statistics.length; i++) {
            for(int j = 0; j < statistics[i].length; j++) {
                for(int k = 0; k < statistics[i][j].length; k++) {
                    sb.append(statistics[i][j][k] + "\t");
                }
                sb.append("\n");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
