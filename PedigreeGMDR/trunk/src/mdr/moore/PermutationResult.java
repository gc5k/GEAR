
package mdr.moore;

import publicAccess.ToolKit;
import publicAccess.PublicData;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
/**
 *
 * @author Guo-Bo Chen
 */

public class PermutationResult {

    private double[][] obs_stat;
    private double[][][] stats;
    private double[][] mean_stat;
    private double[][] sum_stat;
    private int[][] idx_stat;

    private int interval;
    private int replication;
    StringBuffer sb = new StringBuffer();
    public PermutationResult(int numPermutation, int intvl) {
        interval = intvl;
        replication = numPermutation - 1;
        stats = new double[replication][PublicData.NumOfStatistics][interval];
        obs_stat = new double[PublicData.NumOfStatistics][interval];

        sum_stat = new double[PublicData.NumOfStatistics][2];
        mean_stat = new double[PublicData.NumOfStatistics][replication];
        idx_stat = new int[2][replication];
    }

    public void set(int rep, int stat, int intvl, double v) {
        if(rep == 0) {
            obs_stat[stat][intvl] = v;
        } else {
            stats[rep-1][stat][intvl] = v;
            mean_stat[stat][rep-1] += v;
        }
    }

    public void summarise() {
        for (int i = 0; i < mean_stat.length; i++) {
            for (int j = 0; j < mean_stat[i].length; j++) {
                mean_stat[i][j] /= interval;
            }
        }

        for (int i = 0; i < sum_stat.length; i++) {
            sum_stat[i][0] = ToolKit.Mean(mean_stat[i]);
            sum_stat[i][1] = ToolKit.Variance(mean_stat[i]);
            idx_stat[i] = ToolKit.Sort(mean_stat[i]);
        }
    }

    public String toString() {
        sb.append("***********permutation**********" + System.getProperty("line.separator"));
        sb.append("the observed statistics:" + System.getProperty("line.separator"));
        for (int i = 0; i < obs_stat.length; i++) {
        	if ( i == PublicData.TestingAccuIdx ) {
        		sb.append("Mean of Testing Accuracy: ");
        	} else {
        		sb.append("Mean of Training Accuracy: ");
        	}
            for (int j = 0; j < obs_stat[i].length; j++) {
                sb.append(obs_stat[i][j] + " ");
            }
            sb.append(System.getProperty("line.separator"));
        }

//        for (int i = 0; i < stats.length; i++) {
//            sb.append("statistics at permutation " + (i + 1) + System.getProperty("line.separator"));
//            for (int j = 0; j < stats[i].length; j++) {
//                sb.append(j + " Accuracy: ");
//                for (int k = 0; k < stats[i][j].length; k++) {
//                    sb.append(stats[i][j][k] + " ");
//                }
//                sb.append(mean_stat[j][i] +System.getProperty("line.separator"));
//            }
//        }
        sb.append("summary of permutation" + System.getProperty("line.separator"));
        for (int i = 0; i < sum_stat.length; i++) {
            sb.append("Mean is " + sum_stat[i][0] + " +/- " + Math.sqrt(sum_stat[i][1]) + System.getProperty("line.separator"));
        }
        int idx_001 = (new Double(replication * (1-0.01))).intValue();
        int idx_005 = (new Double(replication * (1-0.05))).intValue();
        double pvalue = 0;
        double obs_mean = ToolKit.Mean(obs_stat[0]);
        for (int i = 0; i < replication; i++) {
//            sb.append((i+1) + ": ");
//            for (int j = 0; j < PublicData.NumOfStatistics; j++) {
//                sb.append(mean_stat[j][idx_stat[j][i]] + " ");
//            }
//            sb.append(System.getProperty("line.separator"));
            if(mean_stat[0][idx_stat[0][i]] > obs_mean) {
                pvalue +=1;
            }
        }
 
        sb.append("threshold at 0.05: " + mean_stat[0][idx_stat[0][idx_005]] + System.getProperty("line.separator"));
        sb.append("threshold at 0.01: " + mean_stat[0][idx_stat[0][idx_001]] + System.getProperty("line.separator"));
        sb.append("P value (permutation) is " + (pvalue/replication) + System.getProperty("line.separator"));
        NormalDistribution nd = new NormalDistributionImpl(sum_stat[PublicData.TestingAccuIdx][0], Math.sqrt(sum_stat[PublicData.TestingAccuIdx][1]));
        double mean_tr = 0;
        for (int i = 0; i < obs_stat[PublicData.TestingAccuIdx].length; i++) {
        	mean_tr += obs_stat[PublicData.TestingAccuIdx][i];
        }
        mean_tr /= obs_stat[PublicData.TestingAccuIdx].length;
        try { 
        	sb.append("P value (z-score) is " + (1-nd.cumulativeProbability(mean_tr)));
        } catch (MathException E) {
        	E.printStackTrace(System.err);
        }
        return sb.toString();
    }
}
