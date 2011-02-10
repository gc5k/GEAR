package mdr.alternative;

import algorithm.*;
import mdr.TraitStatistic;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen
 */
public class DataSet extends AbstractList implements Cloneable {

    private ArrayList sample = new ArrayList();
    private String markerFile;
    private String phenotypeFile;
    private String[] SNPID;
    private String[] phenotype;
    private int[][] missing;
    private HashMap traitStatistics = new HashMap();

    public static class Subject extends AbstractList {

        int id;
        private List genotypes;
        private Object status;
        private Object mu;
        private ArrayList score;

        public Subject(String[] content) {
            this.status = content[content.length - 1];
            genotypes = new ArrayList(content.length - 1);
            genotypes.addAll(Arrays.asList(content).subList(0, content.length - 1));
        }

        public void addScore(String[] s) {
            score = new ArrayList(s.length);
            score.addAll(Arrays.asList(s).subList(0, s.length));
        }

        public Object getGenotype(int idx) {
            return genotypes.get(idx);
        }

        public Object getScore(int idx) {
            return score.get(idx);
        }

        protected boolean isScoreExisted(int idx) {
            return (((String) score.get(idx)).compareTo(PublicData.MissingPhenotypeValue) == 0) ? false : true;
        }

        public double getDoubleScore(int idx) {
            return Double.parseDouble((String) score.get(idx));
        }
        // cgb begin
        public void iniMU(double mu) {
            this.mu = new Double(mu);
        }

        public void iniMU(Object mu) {
            this.mu = mu;
        }

        public Object getMU() {
            return mu;
        }

        public double getMU_double() {
            return getDouble(mu);
        }

        private static double getDouble(Object o) {
            return ((Double) o).doubleValue();
        }

        // cgb end
        public Object getStatus() {
            return status;
        }

        public Integer getIntegerID() {
            return new Integer(id);
        }

        public int size() {
            return genotypes == null ? 0 : genotypes.size();
        }

        public Object get(int idx) {
            return genotypes == null ? null : genotypes.get(idx);
        }

        public void add(int index, String s) {
            modCount++;
            genotypes.add(index, s);
        }

        public Object remove(int idx) {
            modCount++;
            return genotypes.remove(idx);
        }

        protected void setID(int idx) {
            id = idx;
        }

        public String toString() {
            StringBuffer b = new StringBuffer();
            for (int i = 0; i < size(); ++i) {
                b.append(get(i));
                b.append(PublicData.delim);
            }
            b.append(status);
            b.append(PublicData.delim);
            Iterator e = score.iterator();
            for (; e.hasNext();) {
                b.append(e.next());
                b.append(PublicData.delim);
            }
            return b.toString();
        }
    }

    public DataSet(String fileM, String fileP) {
        markerFile = new String(fileM);
        phenotypeFile = new String(fileP);
    }

    public void initial() throws IOException {
        FileReader fr = new FileReader(markerFile);
        LineNumberReader l = new LineNumberReader(fr);
        BufferedReader b = new BufferedReader(l);
        String line;
        if ((line = b.readLine()) == null) {
            System.out.println("It's empty.");
            System.exit(0);
        }
        String[] s = line.split(PublicData.delim);
        SNPID = new String[s.length - 1];
        System.arraycopy(s, 0, SNPID, 0, s.length - 1);
        TreeMap tempMap = new TreeMap();
        int count = 0;
        while ((line = b.readLine()) != null) {
            String[] content = line.split(PublicData.delim);
            Subject sub = new Subject(content);
            sub.setID(count);
            add(sub);
            for (int i = 0; i < content.length - 1; i++) {
                Integer Ii = new Integer(i);
                if (content[i].compareTo(PublicData.MissingGenotype) == 0) {
                    if (tempMap.containsKey(Ii)) {
                        ArrayList tempArray = (ArrayList) tempMap.get(Ii);
                        tempArray.add(new Integer(count));
                    } else {
                        ArrayList tempArray = new ArrayList();
                        tempArray.add(new Integer(count));
                        tempMap.put(Ii, tempArray);
                    }
                }
            }
            count++;
        }
        Set keys = tempMap.keySet();
        missing = new int[SNPID.length][];
        for (Iterator e = keys.iterator(); e.hasNext();) {
            Integer Ii = (Integer) e.next();
            ArrayList tempArray = (ArrayList) tempMap.get(Ii);
            missing[Ii.intValue()] = new int[tempArray.size()];
            for (int i = 0; i < tempArray.size(); i++) {
                missing[Ii.intValue()][i] = ((Integer) tempArray.get(i)).intValue();
            }
        }
        addScore();
    }

    public String getMarkerFileName() {
        return markerFile;
    }

    public String getPhenotypeFileName() {
        return phenotypeFile;
    }

    public void add(int idx, Object o) {
        modCount++;
        sample.add(idx, o);
    }

    public double[] calculateDefaultMu(int[] sI) throws DataSetException {
        double[] defaultMu = new double[sI.length];
        for (int i = 0; i < sI.length; i++) {
            int scrIdx = sI[i];
            double sum = 0;
            int c = 0;
            for (Iterator e = sample.iterator(); e.hasNext();) {
                DataSet.Subject s = (DataSet.Subject) e.next();
                if (s.isScoreExisted(scrIdx)) {
                    sum += s.getDoubleScore(scrIdx);
                    c++;
                }
            }
            if (c == PublicData.epsilon) {
                throw new DataSetException("too small denominator");
            }
            defaultMu[i] = sum / c;
        }
        return defaultMu;
    }

    /**
     * This method should be invoked when a new trait is selected for running.
     */
    public void TraitSummary(int[] sI, double[] os) throws DataSetException {
        for (int i = 0; i < sI.length; i++) {
            Integer I = new Integer(i);
            TraitStatistic t = (TraitStatistic) traitStatistics.get(I);
            if (t != null && Math.abs(t.getOffset() - os[i]) < PublicData.epsilon) {
                continue;
            }

            TraitStatistic ts = new TraitStatistic(os[i], sI[i], sample);
            traitStatistics.put(I, ts);

            int posSub = ts.getTraitPositiveSubjects();
            int negSub = ts.getTraitNegativeSubjects();
            double posScr = ts.getTraitPositiveScore();
            double negScr = ts.getTraitNegativeScore();
            if (posSub == 0) {
                throw new DataSetException("Score " + phenotype[i] + "is abnormal that no any subject has a positive score.");
            } else if (negSub == 0) {
                throw new DataSetException("Score " + phenotype[i] + "is abnormal that no any subject has a negative score.");
            } else if (posScr < PublicData.epsilon) {
                throw new DataSetException("Abnormal magnitude of the sum of the positive score: " + posScr + " at phenotype " + phenotype[i] + " .");
            } else if (Math.abs(negScr) < PublicData.epsilon) {
                throw new DataSetException("Abnormal magnitude of the sum of the negative score " + negScr + "  at phenotype " + phenotype[i] + " .");
            }
            double r = posScr / Math.abs(negScr) > 1 ? posScr / Math.abs(negScr) : Math.abs(negScr) / posScr;
            if (r > PublicData.tooBigRatio) {
                throw new DataSetException("Abnormal ratio of positive score v.s. negative score.");
            }
            r = posSub / Math.abs(negScr) > 1 ? posSub / Math.abs(negScr) : Math.abs(negScr) / posSub;
            if (r > PublicData.tooBigRatio) {
                throw new DataSetException("Abnormal ration of subjects with positive score v.s. with negative score.");
            }
        }
    }

    public double getTraitRatio(int[] markerIdx, int scrIdx) throws DataSetException {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(scrIdx));
        int posSub = ts.getTraitPositiveSubjects();
        int negSub = ts.getTraitNegativeSubjects();
        double posScr = ts.getTraitPositiveScore();
        double negScr = ts.getTraitNegativeScore();
        HashSet hs = new HashSet();
        for (int i = 0; i < markerIdx.length; i++) {
            if (missing[markerIdx[i]] == null) {
                continue;
            }
            for (int j = 0; j < missing[markerIdx[i]].length; j++) {
                Integer subIdx = new Integer(missing[markerIdx[i]][j]);
                hs.add(subIdx);
            }
        }
        for (Iterator e = hs.iterator(); e.hasNext();) {
            int idx = ((Integer) e.next()).intValue();
            Subject s = (Subject) sample.get(idx);
            double scr = s.getDoubleScore(scrIdx);
            if (scr >= 0) {
                posScr -= scr;
                posSub--;
            } else {
                negScr -= scr;
                negSub--;
            }
        }
        if (Math.abs(negScr) < PublicData.epsilon) {
            throw new DataSetException("too small denominator.");
        }
        return posScr / Math.abs(negScr);
    }

    public double defaultScore(int[] markerIdx, int scrIdx) throws DataSetException {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(scrIdx));
        int posSub = ts.getTraitPositiveSubjects();
        int negSub = ts.getTraitNegativeSubjects();
        double posScr = ts.getTraitPositiveScore();
        double negScr = ts.getTraitNegativeScore();
        HashSet hs = new HashSet();
        for (int i = 0; i < markerIdx.length; i++) {
            if (missing[markerIdx[i]] == null) {
                continue;
            }
            for (int j = 0; j < missing[markerIdx[i]].length; j++) {
                Integer subIdx = new Integer(missing[markerIdx[i]][j]);
                hs.add(subIdx);
            }
        }
        for (Iterator e = hs.iterator(); e.hasNext();) {
            int idx = ((Integer) e.next()).intValue();
            Subject s = (Subject) sample.get(idx);
            double scr = s.getDoubleScore(scrIdx);
            if (scr >= 0) {
                posScr -= scr;
                posSub--;
            } else {
                negScr -= scr;
                negSub--;
            }
        }
        if (posScr + Math.abs(negScr) < PublicData.epsilon) {
            throw new DataSetException("too small denominator.");
        }
        return posScr / (posScr + Math.abs(negScr));
    }

    public int getMarkerNum() {
        return SNPID.length;
    }

    public String[] getMarkers(int[] idx) {
        String[] mkIdx = new String[idx.length];
        for (int i = 0; i < idx.length; i++) {
            mkIdx[i] = SNPID[idx[i]];
        }
        return mkIdx;
    }

    public Object get(int idx) {
        return sample.get(idx);
    }

    protected void addScore() throws IOException {
        FileReader fr = new FileReader(phenotypeFile);
        LineNumberReader l = new LineNumberReader(fr);
        BufferedReader b = new BufferedReader(l);
        String line;
        if ((line = b.readLine()) == null) {
            System.out.println("It's empty.");
            System.exit(0);
        }
        phenotype = line.split(PublicData.delim);
        List pl = new ArrayList();
        while ((line = b.readLine()) != null) {
            pl.add(line);
        }
        if (pl.size() != size()) {
            System.err.println("The size of the sample unmatch the size of the score");
            return;
        }

        Iterator e = pl.iterator();
        int c = 0;
        for (; e.hasNext();) {
            String p = (String) e.next();
            String[] ps = p.split(PublicData.delim);
            if (ps.length != phenotype.length) {
                System.err.println("The " + c + "th phenotype column is not match the total size of the phenotypes");
            }
            Subject sub = (Subject) sample.get(c++);
            sub.addScore(ps);
        }
    }

    public Object remove(int index) {
        modCount++;
        return sample.remove(index);
    }

    public void print() {
        for (Iterator e = sample.iterator(); e.hasNext();) {
            System.out.println(e.next());
        }
    }

    public void printMissing() {
        for (int i = 0; i < missing.length; i++) {
            if (missing[i] != null) {
                for (int j = 0; j < missing[i].length; j++) {
                    System.out.print(missing[i][j] + "\t");
                }
                System.out.println();
            }
        }
    }

    public List getSample() {
        return sample;
    }

    public int size() {
        return sample.size();
    }

    public Object clone() {
        try {
            DataSet o = (DataSet) super.clone();
            o.sample = (ArrayList) sample.clone();
            return o;
        } catch (CloneNotSupportedException c) {
            c.printStackTrace(System.err);
            return null;
        }
    }
}
