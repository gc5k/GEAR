package im;

import java.util.HashMap;
import im.population.IMPopulation;

/**
 *
 * @author Guo-Bo Chen
 */
public class GenomeScan {

    private String phenotypeType;
    private double step;
    private HashMap<String, IntervalPriorProbability> intervalPP;
    private IMPopulation imp;

    public GenomeScan(IMPopulation ip, double s) {
        imp = ip;
        step = s;
        intervalPP = new HashMap();
    }

    public void CalculateIPP() {
        if(!imp.IsPicked()) {
            for (int i = 0; i < imp.ChromosomeNumber(); i++) {
                for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                    IntervalPriorProbability ipp = new IntervalPriorProbability(imp, step, (imp.DistanceAt(i, j + 1) - imp.DistanceAt(i, j)), i, j);
                    String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, i, j);
                    intervalPP.put(key, ipp);
                    for (int k = 0; k < imp.ObservedGenotype(); k++) {
                        int markerScore = (imp.MarkerAt(k, i, j)) * 10 + imp.MarkerAt(k, i, j + 1);
                        if (!ipp.ScanInterval(markerScore)) {
                            break;
                        }
                    }
                }
            }
        } else {
            int[][] pmi = imp.getPickedMarkerIndex();
            for(int i = 0; i < pmi.length; i++) {
                for(int j = 0; j < pmi[i].length-1; j++) {
                    IntervalPriorProbability ipp = new IntervalPriorProbability(imp, step, imp.DistanceAt(i, j+1) - imp.DistanceAt(i, j), i, j);
                    String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, i, j);
                    intervalPP.put(key, ipp);
                    for(int k = 0; k < imp.ObservedGenotype(); k++) {
                        int markerScore = (imp.MarkerAt(k, i, j) *10 + imp.MarkerAt(k, i, j+1));
                        if (!ipp.ScanInterval(markerScore)) {
                            break;
                        }
                    }
                }
            }
        }
    }

    public IntervalPriorProbability getIPPTable(int chr, int mk) {
        return intervalPP.get(IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, chr, mk));
    }

    public void getChrInt(int[] SNPIdx) {
        
    }

    public double getStep() {
        return step;
    }

    public void printIPP(int chr, int mk) {
        String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, chr, mk);
        IntervalPriorProbability ipp = intervalPP.get(key);
        System.out.println(ipp);
        return;
    }

    public void printIPP() {
        for (int i = 0; i < imp.ChromosomeNumber(); i++) {
            for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, i, j);
                IntervalPriorProbability ipp = intervalPP.get(key);
                System.out.println(key);
                System.out.println(ipp);
            }
        }
    }

    public int getIntervalPriorProbability(int chr, int mk, int markerScore) {
        String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, chr, mk);
        IntervalPriorProbability ipp = intervalPP.get(key);
        return ipp.getMarkerScoreIndex(markerScore);
    }

    public int getIPPRowIndexForIndividual(int id, int chr, int mk) {
        int markerScore;
        int midx = id;
        HashMap hsIndGeno = imp.getIndGeno();
        if (hsIndGeno != null) {
            midx = Integer.parseInt((String) hsIndGeno.get(new Integer(id)));
        }
        markerScore = imp.MarkerAt(midx, chr, mk)*10 + imp.MarkerAt(midx, chr, mk+1);
        String key = IMToolKit.MakeKey(IMToolKit.chr_prefix, IMToolKit.interval_prefix, chr, mk);
        IntervalPriorProbability ipp = intervalPP.get(key);
        return ipp.getMarkerScoreIndex(markerScore);
    }
}
