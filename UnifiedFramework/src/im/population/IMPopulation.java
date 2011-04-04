package im.population;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.io.PrintStream;
import mdr.TraitStatistic;
import im.reader.IMReaderAbstract;
import im.reader.IMConstant;
import im.IMToolKit;
import regression.LinearRegression;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class IMPopulation {

    protected String[][] markerName;
    String[] MDRselectedMarkerName;
    protected int[][][] marker;
    protected double[][] distance;
    protected double[][] recombination;
    protected double[][] phenotype;
    protected String[][] strPhenotype;
    protected double[][] nuisance;
    protected double[][] score;
    protected String populationType;
    protected int mapping_function;
    protected int Environment;
    protected int Replication;
    protected HashMap IndNestedEnv;
    protected HashMap IndGeno;
    protected ArrayList IDs;
    protected HashMap traitStatistics;
    protected ArrayList permutationIDs;
    protected boolean isPermutation = false;
    protected boolean picked = false;
    protected IMReaderAbstract.CrossParameter crp;
    protected String[][] ancillaryPheno;
    protected int[][] pickedMarkerIndex;
    protected ArrayList pickedMarkerIdx;
    private boolean isSimulation;

    public static class MarkerIndex {

        private int chr;
        private int index;

        public MarkerIndex(int c, int i) {
            chr = c;
            index = i;
        }

        public int getChromosome() {
            return chr;
        }

        public int getIndex() {
            return index;
        }
    }

    public IMPopulation() {
        initial0();
    }

    public IMPopulation(IMReaderAbstract ira) {
        initial0();
        initial(ira);
    }

    private void initial0() {
        isSimulation = true;
        IDs = new ArrayList();
        permutationIDs = new ArrayList();
        IndNestedEnv = new HashMap();
        traitStatistics = new HashMap();
    }

    private void initial(IMReaderAbstract ira) {
        isSimulation = false;
        IndGeno = new HashMap();
        strPhenotype = ira.getPhenotype();
        markerName = ira.getMarkerName();
        distance = ira.getDistance();
        recombination = ira.getRecombination();
        IMReaderAbstract.ChromosomeParameter chrp = ira.getChromosomeParameter();
        mapping_function = chrp.getFunction();
        crp = ira.getCrossParameter();
        populationType = crp.getCross();
        ancillaryPheno = ira.getAncillaryPheno();

        HashSet hsEnv = new HashSet();
        HashSet hsRep = new HashSet();
        HashSet hsInd = new HashSet();
        for (int i = 0; i < ancillaryPheno.length; i++) {
            hsEnv.add(ancillaryPheno[i][0]);
            hsRep.add(ancillaryPheno[i][1]);
            hsInd.add(ancillaryPheno[i][2]);
            int ind = Integer.parseInt(ancillaryPheno[i][2]);
            ind--;
            IndNestedEnv.put(new Integer(i), ancillaryPheno[i][0]);
            IndGeno.put(new Integer(i), (new Integer(ind)).toString());
            IDs.add(i);
            permutationIDs.add(i);
        }

        Replication = hsRep.size();
        Environment = hsEnv.size();
        String[][] markercodingtable = crp.getMarkerCoding();
        String[][][] m = ira.getMarkers();
        marker = new int[m.length][][];
        for (int i = 0; i < m.length; i++) {
            marker[i] = new int[m[i].length][];
            for (int j = 0; j < m[i].length; j++) {
                marker[i][j] = new int[m[i][j].length];
                for (int k = 0; k < m[i][j].length; k++) {
                    for (int v = 0; v < 3; v++) {
                        if (m[i][j][k].compareTo(markercodingtable[IMConstant.AA][2]) == 0) {
                            marker[i][j][k] = IMToolKit.AA;
                        } else if (m[i][j][k].compareTo(markercodingtable[IMConstant.Aa][2]) == 0) {
                            marker[i][j][k] = IMToolKit.Aa;
                        } else {
                            marker[i][j][k] = IMToolKit.aa;
                        }
                    }
                }
            }
        }
        score = new double[strPhenotype.length][];
        for (int i = 0; i < strPhenotype.length; i++) {
            score[i] = new double[strPhenotype[i].length];
        }

        phenotype = new double[strPhenotype.length][];
        for (int i = 0; i < strPhenotype.length; i++) {
            phenotype[i] = new double[strPhenotype[i].length];
            for (int j = 0; j < strPhenotype[i].length; j++) {
                if (strPhenotype[i][j].compareTo(crp.getMissingTrait()) == 0) {
                    phenotype[i][j] = 0;
                } else {
                    phenotype[i][j] = Double.parseDouble(strPhenotype[i][j]);
                }
            }
        }
        BasicTraitStatistics();
    }

    public void PickAllMarkers() {
        pickedMarkerIdx = new ArrayList();
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                MarkerIndex mi = new MarkerIndex(i, j);
                pickedMarkerIdx.add(mi);
            }
        }

        pickedMarkerIndex = new int[distance.length][];
        for (int i = 0; i < distance.length; i++) {
            pickedMarkerIndex[i] = new int[distance[i].length];
            for (int j = 0; j < distance[i].length; j++) {
                pickedMarkerIndex[i][j] = j;
            }
        }
        picked = true;
        MDRSetSelectedMarker(pickedMarkerIndex);
    }

    public void PickMarkers(double cM) {
        if (cM > 0) {
            pickedMarkerIndex = new int[distance.length][];
            for (int i = 0; i < distance.length; i++) {
                ArrayList m = new ArrayList();
                double d = 0;
                int idx = 0;
                m.add(new Integer(0));
                for (int j = 0; j < distance[i].length; j++) {
                    if (Math.abs(distance[i][j] - d) >= cM) {
                        if (idx != j - 1) {
                            idx = j - 1;
                        } else {
                            idx = j;
                        }
                        m.add(new Integer(idx));
                        d = distance[i][idx];
                    }
                }
                if (m.size() < 2 || idx != (distance[i].length - 1)) {
                    m.add(new Integer(distance[i].length - 1));
                }
                pickedMarkerIndex[i] = new int[m.size()];
                for (int j = 0; j < m.size(); j++) {
                    pickedMarkerIndex[i][j] = ((Integer) m.get(j)).intValue();
                }
            }
            picked = true;
            MDRSetSelectedMarker(pickedMarkerIndex);
        } else {
            pickedMarkerIdx = new ArrayList();
            for (int i = 0; i < distance.length; i++) {
                for (int j = 0; j < distance[i].length; j++) {
                    MarkerIndex mi = new MarkerIndex(i, j);
                    pickedMarkerIdx.add(mi);
                }
            }

            pickedMarkerIndex = new int[distance.length][];
            for (int i = 0; i < distance.length; i++) {
                pickedMarkerIndex[i] = new int[distance[i].length];
                for (int j = 0; j < distance[i].length; j++) {
                    pickedMarkerIndex[i][j] = j;
                }
            }
            picked = true;
            MDRSetSelectedMarker(pickedMarkerIndex);
        }
    }

    public void SetSelectedMarker(int[][] loc) {
        pickedMarkerIdx = new ArrayList();
        for (int i = 0; i < loc.length; i++) {
            IMPopulation.MarkerIndex mi = new IMPopulation.MarkerIndex(loc[i][0], loc[i][1]);
            pickedMarkerIdx.add(mi);
        }

        if (picked) {
            MDRselectedMarkerName = new String[pickedMarkerIdx.size()];
            for (int i = 0; i < pickedMarkerIdx.size(); i++) {
                IMPopulation.MarkerIndex mi = (IMPopulation.MarkerIndex) pickedMarkerIdx.get(i);
                MDRselectedMarkerName[i] = markerName[mi.getChromosome()][mi.getIndex()];
            }
        }
    }

    public ArrayList getPickedMarkerIdx() {
        return pickedMarkerIdx;
    }

    public int[][] getPickedMarkerIndex() {
        return pickedMarkerIndex;
    }

    public void Swith2Permutation(boolean ip, long sd) {
        isPermutation = ip;
        if (isPermutation) {
            Random rnd = new Random(sd);
            Collections.shuffle(permutationIDs, rnd);
        }
    }

    /**
     * This one is designed for the my IF2 data
     */
    public void tunePermutation(long sd) {
        int ogi = ObservedGenotype();
        int env = getEnvironment();
        ArrayList PI = new ArrayList();
        for (int i = 0; i < env; i++) {
            for (int j = 0; j < ogi; j++) {
                PI.add(new Integer(ogi * i * Replication + j));
            }
        }
        Random rnd = new Random(sd);
        Collections.shuffle(PI, rnd);
        for (int i = 0; i < env; i++) {
            for (int j = 0; j < ogi; j++) {
                Integer I = (Integer) PI.get(i * ogi + j);
                for (int k = 0; k < Replication; k++) {
                    int idx = i * Replication * ogi + k * ogi + j;
                    int v = I.intValue() + k * ogi;
                    permutationIDs.set(idx, new Integer(v));
                }
            }
        }
    }

    public HashMap getIndGeno() {
        return IndGeno;
    }

    public String getCrossParameter() {
        return populationType;
    }

    public int IndividualNumber() {
        return phenotype.length;
    }

    public int ObservedGenotype() {
        return marker.length;
    }

    public int MapFunction() {
        return mapping_function;
    }

    public int ChromosomeNumber() {
        if (picked) {
            return pickedMarkerIndex.length;
        } else {
            return distance.length;
        }
    }

    public int ChromosomeNumberSelected() {
        if (picked) {
            return pickedMarkerIndex.length;
        } else {
            return 0;
        }
    }

    public int ChromosomeNumberOriginal() {
        return distance.length;
    }

    public int IntervalNumberAtChromosome(int chr) {
        if (picked) {
            return pickedMarkerIndex[chr].length - 1;
        } else {
            return distance[chr].length - 1;
        }
    }

    public int IntervalNumberAtChromosomeSelected(int chr) {
        if (picked) {
            return pickedMarkerIndex[chr].length - 1;
        } else {
            return 0;
        }
    }

    public int IntervalNumberAtChromosomeOriginal(int chr) {
        return distance[chr].length - 1;
    }

    public int SumIntevals() {
        int s = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                s += IntervalNumberAtChromosome(i);
            }
        } else {
            for (int i = 0; i < distance.length; i++) {
                s += IntervalNumberAtChromosome(i);
            }
        }
        return s;
    }

    public int SumIntevalsSelected() {
        int s = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                s += IntervalNumberAtChromosome(i);
            }
        }
        return s;
    }

    public int SumIntevalsOriginal() {
        int s = 0;
        for (int i = 0; i < distance.length; i++) {
            s += IntervalNumberAtChromosome(i);
        }
        return s;
    }

    public int SumMarkers() {
        int s = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                s += pickedMarkerIndex[i].length;
            }
        } else {
            for (int i = 0; i < marker[0].length; i++) {
                s += marker[0][i].length;
            }
        }
        return s;
    }

    public int SumMarkersSelected() {
        int s = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                s += pickedMarkerIndex[i].length;
            }
        }
        return s;
    }

    public int SumMarkersOriginal() {
        int s = 0;
        for (int i = 0; i < marker[0].length; i++) {
            s += marker[0][i].length;
        }
        return s;
    }

    public boolean IsPicked() {
        return picked;
    }

    public double[][] Distance() {
        double[][] d = null;
        if (picked) {
            d = new double[pickedMarkerIndex.length][];
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                d[i] = new double[pickedMarkerIndex[i].length];
                for (int j = 0; j < pickedMarkerIndex[i].length; j++) {
                    d[i][j] = distance[i][pickedMarkerIndex[i][j]];
                }
            }
            return d;
        } else {
            return distance;
        }
    }

    public double[][] DistanceSelected() {
        double[][] d = null;
        if (picked) {
            d = new double[pickedMarkerIndex.length][];
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                d[i] = new double[pickedMarkerIndex[i].length];
                for (int j = 0; j < pickedMarkerIndex[i].length; j++) {
                    d[i][j] = distance[i][pickedMarkerIndex[i][j]];
                }
            }
            return d;
        } else {
            return distance;
        }
    }

    public double[][] DistanceOriginal() {
        return distance;
    }

    public double DistanceAt(int chr, int mk) {
        if (picked) {
            return distance[chr][pickedMarkerIndex[chr][mk]];
        } else {
            return distance[chr][mk];
        }
    }

    public double DistanceAtSelected(int chr, int mk) {
        if (picked) {
            return distance[chr][pickedMarkerIndex[chr][mk]];
        } else {
            return 0.0;
        }
    }

    public double DistanceAtOriginal(int chr, int mk) {
        return distance[chr][mk];
    }

    public int MarkerAt(int idv, int chr, int mk) {
        if (picked) {
            return marker[idv][chr][pickedMarkerIndex[chr][mk]];
        }
        return marker[idv][chr][mk];
    }

    public int MarkerAtSelected(int idv, int chr, int mk) {
        if (picked) {
            return marker[idv][chr][pickedMarkerIndex[chr][mk]];
        } else {
            return -1;
        }
    }

    public int MarkerAtOriginal(int idv, int chr, int mk) {
        return marker[idv][chr][mk];
    }

    public int TraitsNumber() {
        return phenotype[0].length;
    }

    public boolean PhenotypeExist(int idv, int sI) {
        String missingtrait;
        if (isSimulation) {
            missingtrait = IMToolKit.MissingTrait;
        } else {
            missingtrait = crp.getMissingTrait();
        }
        if (strPhenotype[idv][sI].compareTo(missingtrait) == 0) {
            return false;
        } else {
            return true;
        }
    }

    public double PhenotypeAt(int idv, int sI) {
        String missingtrait;
        if (isSimulation) {
            missingtrait = IMToolKit.MissingTrait;
        } else {
            missingtrait = crp.getMissingTrait();
        }
        if (strPhenotype[idv][sI].compareTo(missingtrait) == 0) {
            return 0;
        } else {
            return phenotype[idv][sI];
        }
    }

    public double ScoreAt(int idv, int sI) {
        String missingtrait;
        if (isSimulation) {
            missingtrait = IMToolKit.MissingTrait;
        } else {
            missingtrait = crp.getMissingTrait();
        }
        if (strPhenotype[idv][sI].compareTo(missingtrait) == 0) {
            return 0;
        } else {
            return score[idv][sI];
        }
    }

    public void setScoreAt(int idv, int sI, double sc) {
        score[idv][sI] = sc;
    }

    public ArrayList getIDs() {
        return IDs;
    }

    public String NestedEnvironment(int idx) {
        return (String) IndNestedEnv.get(new Integer(idx));
    }

    public ArrayList getPermutatedIDs() {
        return permutationIDs;
    }

    protected void BasicTraitStatistics() {
        for (int i = 0; i < TraitsNumber(); i++) {
            TraitStatistic ts = new TraitStatistic();
            traitStatistics.put(new Integer(i), ts);
            ts.DefaultScore(this, i);
        }
    }

    protected void BasicTraitStatistics(int i) {
        TraitStatistic ts = new TraitStatistic();
        traitStatistics.put(new Integer(i), ts);
        ts.DefaultScore(this, i);
    }

    public double getMean(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        return ts.CalculateMU(this, idx);
    }

    public int getEnvironment() {
        return Environment;
    }

    public double getTraitRatio(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        double ps = ts.getTraitPositiveScore();
        double ns = ts.getTraitNegativeScore();
        return ps / ((-1) * ns);
    }
    
    public double getPositiveScore(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        return ts.getTraitPositiveScore();    	
    }

    public double getNegativeScore(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        return ts.getTraitNegativeScore();    	
    }

    public int getPositiveSubjects(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        return ts.getTraitPositiveSubjects();    	
    }

    public double getNegativeSubjects(int idx) {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(idx));
        return ts.getTraitNegativeSubjects();    	
    }

    
    public void BuildScore(int sIdx, boolean env) {
        ArrayList<Integer> observed = new ArrayList();
        for (int i = 0; i < IndividualNumber(); i++) {
            if (PhenotypeExist(i, sIdx)) {
                observed.add(new Integer(i));
            }
        }
        double[][] phe_y = new double[observed.size()][1];

        double[][] pre_x;
        if (env) {
            pre_x = new double[observed.size()][2];
            for (int i = 0; i < observed.size(); i++) {
                int idx = observed.get(i).intValue();
                phe_y[i][0] = PhenotypeAt(idx, sIdx);
                pre_x[i][0] = 1;
                pre_x[i][1] = Integer.parseInt(NestedEnvironment(idx));
            }
        } else {
            pre_x = new double[observed.size()][1];
            for (int i = 0; i < observed.size(); i++) {
                int idx = observed.get(i).intValue();
                phe_y[i][0] = PhenotypeAt(idx, sIdx);
                pre_x[i][0] = 1;
            }
        }

        LinearRegression lm = new LinearRegression(pre_x, phe_y);
        lm.MLE();
        double[] s = lm.getResiduals();
        for (int i = 0; i < observed.size(); i++) {
            setScoreAt(observed.get(i).intValue(), sIdx, s[i]);
        }
    }

    public String[][] MarkerName() {
        if (picked) {
            String[][] mn = new String[pickedMarkerIndex.length][];
            for (int i = 0; i < mn.length; i++) {
                mn[i] = new String[pickedMarkerIndex[i].length];
                for (int j = 0; j < mn[i].length; j++) {
                    mn[i][j] = new String(markerName[i][pickedMarkerIndex[i][j]]);
                }
            }
            return mn;
        } else {
            return markerName;
        }
    }

    public String[][] MarkerNameSelected() {
        String[][] mn = null;
        if (picked) {
            mn = new String[pickedMarkerIndex.length][];
            for (int i = 0; i < mn.length; i++) {
                mn[i] = new String[pickedMarkerIndex[i].length];
                for (int j = 0; j < mn[i].length; j++) {
                    mn[i][j] = new String(markerName[i][pickedMarkerIndex[i][j]]);
                }
            }
        }
        return mn;
    }

    public String[][] MarkerNameOriginal() {
        return markerName;
    }

    public MarkerInformation CombinationGeneticInformationOriginal(String com) {
        int[] SNPIndex = ToolKit.StringToIntArray(com);
        int[][] ChrInt = new int[SNPIndex.length][2];
        double[] coordinate = new double[SNPIndex.length];
        double[][] dis = DistanceOriginal();
        for (int i = 0; i < ChrInt.length; i++) {
            int c = 0;
            boolean flag = false;
            for (int j = 0; j < dis.length; j++) {
                for (int k = 0; k < dis[j].length; k++) {
                    if (c == SNPIndex[i]) {
                        flag = true;
                        ChrInt[i][0] = j;
                        ChrInt[i][1] = k;
                        coordinate[i] = dis[j][k];
                    }
                    c++;
                }
                if (flag) {
                    break;
                }
            }
        }
        MarkerInformation mkinfor = new MarkerInformation(ChrInt, coordinate);
        return mkinfor;
    }

    public MarkerInformation CombinationGeneticInformationSelected(String com) {
        int[] SNPIndex = ToolKit.StringToIntArray(com);
        int[][] ChrInt = new int[SNPIndex.length][2];
        double[] coordinate = new double[SNPIndex.length];
        double[][] dis = DistanceSelected();
        for (int i = 0; i < ChrInt.length; i++) {
            int c = 0;
            boolean flag = false;
            for (int j = 0; j < dis.length; j++) {
                for (int k = 0; k < dis[j].length; k++) {
                    if (c == SNPIndex[i]) {
                        flag = true;
                        ChrInt[i][0] = j;
                        ChrInt[i][1] = k;
                        coordinate[i] = dis[j][k];
                    }
                    c++;
                }
                if (flag) {
                    break;
                }

            }
        }
        MarkerInformation mkinfor = new MarkerInformation(ChrInt, coordinate);
        return mkinfor;
    }

    public int[][][] Markers() {
        int[][][] mk = null;
        if (picked) {
            mk = new int[marker.length][][];
            for (int i = 0; i < marker.length; i++) {
                mk[i] = new int[marker[i].length][];
                for (int j = 0; j < pickedMarkerIndex.length; j++) {
                    mk[i][j] = new int[pickedMarkerIndex[j].length];
                    for (int k = 0; k < pickedMarkerIndex[j].length; k++) {
                        mk[i][j][k] = marker[i][j][pickedMarkerIndex[j][k]];
                    }
                }
            }
            return mk;
        } else {
            return marker;
        }
    }

    public int[][][] MarkersSelected() {
        int[][][] mk = null;
        if (picked) {
            mk = new int[marker.length][][];
            for (int i = 0; i < marker.length; i++) {
                mk[i] = new int[marker[i].length][];
                for (int j = 0; j < pickedMarkerIndex.length; j++) {
                    mk[i][j] = new int[pickedMarkerIndex[j].length];
                    for (int k = 0; k < marker[i][j].length; k++) {
                        mk[i][j][k] = marker[i][j][pickedMarkerIndex[j][k]];
                    }
                }
            }
        }
        return mk;
    }

    public int[][][] MarkersOriginal() {
        return marker;
    }

    public int MDRMarkerIndexOriginal(int chr, int idx) {
        int index = 0;
        boolean flag = false;
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                if (i == chr && j == idx) {
                    flag = true;
                    break;
                }
                index++;
            }
            if (flag) {
                break;
            }
        }
        return index;
    }

    public void MDRSelecteMarker(int[][] loc) {
        pickedMarkerIdx = new ArrayList();
        for (int i = 0; i < loc.length; i++) {
            for (int j = 0; j < loc[i].length; j++) {
                IMPopulation.MarkerIndex mi = new IMPopulation.MarkerIndex(i, loc[i][j]);
                pickedMarkerIdx.add(mi);
            }
        }

        if (picked) {
            MDRselectedMarkerName = new String[pickedMarkerIdx.size()];
            for (int i = 0; i < pickedMarkerIdx.size(); i++) {
                IMPopulation.MarkerIndex mi = (IMPopulation.MarkerIndex) pickedMarkerIdx.get(i);
                MDRselectedMarkerName[i] = markerName[mi.getChromosome()][mi.getIndex()];
            }
        }
    }

    private void MDRSetSelectedMarker(int[][] loc) {
        pickedMarkerIdx = new ArrayList();
        for (int i = 0; i < loc.length; i++) {
            for (int j = 0; j < loc[i].length; j++) {
                IMPopulation.MarkerIndex mi = new IMPopulation.MarkerIndex(i, loc[i][j]);
                pickedMarkerIdx.add(mi);
            }
        }

        if (picked) {
            MDRselectedMarkerName = new String[pickedMarkerIdx.size()];
            for (int i = 0; i < pickedMarkerIdx.size(); i++) {
                IMPopulation.MarkerIndex mi = (IMPopulation.MarkerIndex) pickedMarkerIdx.get(i);
                MDRselectedMarkerName[i] = markerName[mi.getChromosome()][mi.getIndex()];
            }
        }
    }

    public int MDRMarkerNumber() {
        return MDRselectedMarkerName.length;
    }

    public String[] MDRMarkerName() {
        return MDRselectedMarkerName;
    }

    public ArrayList MDRPickedMarkerIdx() {
        return pickedMarkerIdx;
    }

    public int MarkerNumber() {
        int c = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                c += pickedMarkerIndex[i].length;
            }
            return c;
        } else {
            for (int i = 0; i < distance.length; i++) {
                c += distance[i].length;
            }
            return c;
        }
    }

    public int MarkerNumberSelected() {
        int c = 0;
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                c += pickedMarkerIndex[i].length;
            }
        }
        return c;
    }

    public int MarkerNumberOriginal() {
        int c = 0;
        for (int i = 0; i < distance.length; i++) {
            c += distance[i].length;
        }
        return c;
    }

    public int MarkerNumber(int chr) {
        if (picked) {
            return pickedMarkerIndex[chr].length;
        } else {
            return distance[chr].length;
        }
    }

    public int MarkerNumberSelected(int chr) {
        if (picked) {
            return pickedMarkerIndex[chr].length;
        } else {
            return 0;
        }
    }

    public int MarkerNumberOriginal(int chr) {
        return distance[chr].length;
    }

    public String MarkerNameAt(int chr, int mk) {
        if (picked) {
            return markerName[chr][pickedMarkerIndex[chr][mk]];
        }
        return markerName[chr][mk];
    }

    public String MarkerNameAtSelected(int chr, int mk) {
        if (picked) {
            return markerName[chr][pickedMarkerIndex[chr][mk]];
        } else {
            return null;
        }
    }

    public String MarkerNameAtOriginal(int chr, int mk) {
        return markerName[chr][mk];
    }

    public double[][] getPhenotype() {
        return phenotype;
    }

    public int getReplication() {
        return Replication;
    }

    public void printDistance() {
        if (picked) {
            for (int i = 0; i < pickedMarkerIndex.length; i++) {
                for (int j = 0; j < pickedMarkerIndex[i].length; j++) {
                    System.out.print(distance[i][j] + " ");
                }
                System.out.println();
            }
        }
    }
    
    public void CIMFormat(PrintStream ps) {
    	ps.println("#FileID 285296362");
    	ps.println("#bychromosome");
    	ps.println("-type position");
    	ps.println("-function 1");
    	ps.println("-Units M");
    	ps.println("-chromosomes 1");
    	ps.println("-maximum " + MarkerNumber());
    	ps.println("-named yes");
    	ps.println("-start");
    	for (int i = 1; i <= ChromosomeNumber(); i++) {
    		ps.println("-Chromosome C" + i);
    		for (int j = 0; j <MarkerNumber(i-1); j++) {
    			ps.println(MarkerNameAt(i-1, j) + "\t" + distance[i-1][j]);
    		}
    	}
    	ps.println("-stop\n---------------------------------------------------");
    	ps.println("#bycross");
    	ps.println("-SampleSize "+ IDs.size());
    	ps.println("-Cross " + populationType);
    	ps.println("-traits 1");
    	ps.println("-missingtrait .");
    	ps.println("-case yes");
    	ps.println("-TranslationTable");
    	ps.println("AA    2     2");
    	ps.println("Aa    1     1");
    	ps.println("aa    0     0");
    	ps.println("A-    12    12");
    	ps.println("a-    10    10");
    	ps.println("--    -1    -1");
    	ps.println("-start individuals markers");
    	for (int i = 0; i < IndividualNumber(); i++) {
    		ps.print("I-" + (i+1));
    		for (int j = 0; j < ChromosomeNumber(); j++) {
    			for (int k = 0; k < MarkerNumber(j); k++) {
    				ps.print(" " + MarkerAt(i, j, k));
    			}
    		}
    		ps.print("\n");
    	}
    	ps.println("-stop individuals markers");
    	ps.println("-start traits");
    	ps.print("t1");
    	for (int i = 0; i < IndividualNumber(); i++) {
    		ps.print(" " + PhenotypeAt(i, 0));
    	}
    	ps.println("\n-stop traits");
    	ps.println("-quit");
    	ps.println("-end");
    	
    }
    
    public String toString() {
    	StringBuffer bs = new StringBuffer();
    	for (int i = 0; i < IndividualNumber(); i++) {
    		for (int j = 0; j < ChromosomeNumberOriginal(); j++) {
    			for (int k = 0; k < MarkerNumber(j); k++) {
    				bs.append(MarkerAt(i,j,k) + " ");
    			}
    			bs.append(PhenotypeAt(i, 0) + "\n");
    		}
    	}
    	return bs.toString();
    }
}
