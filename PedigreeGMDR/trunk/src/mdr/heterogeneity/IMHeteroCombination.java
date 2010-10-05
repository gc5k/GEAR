/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package mdr.heterogeneity;

/**
 *
 * @author gc5k
 */
public class IMHeteroCombination {
    private double[][][][] statistics;
    private String name;                        // combination of loci

    IMHeteroCombination(int numTraits, int interval, int dimension) {
        statistics = new double[numTraits][interval][dimension][2];
    }

    public void SetStatistic(int traitidx, int interval, int dimidx, int statisticidx, double value) {
        statistics[traitidx][interval][dimidx][statisticidx] = value;
    }

    public String GetModelName() {
        return name;
    }

    public void SetModelName(String n) {
        name = n;
    }

    public double get(int traitidx, int interval, int dimidx, int statisticidx) {
        return statistics[traitidx][interval][dimidx][statisticidx];
    }
}
