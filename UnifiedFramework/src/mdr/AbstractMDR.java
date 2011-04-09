package mdr;

import java.util.HashMap;
import java.util.HashSet;

import mdr.algorithm.CombinationGenerator;
import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
/**
 *
 * @author Guo-Bo Chen
 */
public abstract class AbstractMDR {
    protected boolean isMooreMDR;
    protected int order;
    protected int numTraits;
    protected double[] offset;
    protected double[][] currBestStatistic;//keeps currently the best statistic given the order of interaction

    protected HashSet<Integer> testedInteraction = new HashSet();//Interaction has already tested
    protected HashMap dataPartitionMap;

    protected CombinationGenerator comGenerator;
    protected DataFile data;
    protected Subdivision subdivision;

    public AbstractMDR(DataFile dr, Subdivision sd, CombinationGenerator cg, int traits, double[] os, boolean ismooremdr) {
        data = dr;
        subdivision = sd;
        comGenerator = cg;
        numTraits = traits;
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);
        isMooreMDR = ismooremdr;
    }

    public boolean isMooreMDR() {
        return isMooreMDR;
    }

    public abstract void search(int or);
}
