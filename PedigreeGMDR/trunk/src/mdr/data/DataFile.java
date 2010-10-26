package mdr.data;

import mdr.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import publicAccess.PublicData;
import im.reader.IMReaderAbstract;
import im.population.IMPopulation;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class DataFile {

    protected ArrayList sample = new ArrayList();
    protected ArrayList shuffledPhenoIDs;
    protected static boolean isPermutation = false;
    protected String markerFile;
    protected String phenotypeFile;
    protected String[] SNPID;
    protected String[] traitName;
    protected int[][] missing;
    protected static int[] scrIndex;
    protected double[] offset;
    protected HashMap traitStatistics = new HashMap();

    public static class Subject extends AbstractList {

        int id;
        private List genotypes;
        private Double status;
        private Object mu;
        private ArrayList<Double> score;
        private ArrayList<Double> shuffledScore;

        public Subject(String[] content) {
            if (content[content.length - 1].equals(PublicData.MissingPhenotypeValue)) {
                this.status = null;
            } else {
                this.status = new Double(Double.parseDouble(content[content.length - 1]));
            }
            genotypes = new ArrayList(content.length - 1);
            genotypes.addAll(Arrays.asList(content).subList(0, content.length - 1));
        }

        public void addScore(String[] s) {
            score = new ArrayList<Double>(s.length);
            for (int i = 0; i < s.length; i++) {
                if (s[i].equals(PublicData.MissingPhenotypeValue)) {
                    score.add(null);
                } else {
                    score.add(Double.parseDouble(s[i]));
                }
            }
        }

        public void addShuffledScore(ArrayList<Double> s) {
            shuffledScore = s;
        }

        public Object getGenotype(int idx) {
            return genotypes.get(idx);
        }

        public ArrayList<Double> getScore() {
            return score;
        }

        public Double getScore(int idx) {
            if (scrIndex[idx] >= 0) {
                return isPermutation ? shuffledScore.get(scrIndex[idx]) : score.get(scrIndex[idx]);
            } else {
                return status;
            }
        }

        /**
         * If the score does not exists, then null is returned.
         * 
         * @param idx
         * @return
         */
        public Double getDoubleScore(int idx) {
            if (scrIndex[idx] >= 0) {
                return isPermutation ? shuffledScore.get(scrIndex[idx]) : score.get(scrIndex[idx]);
            } else {
                return status;
            }
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

    public DataFile(String fileM, String fileP, int[] scrIdx) {
        markerFile = fileM;
        phenotypeFile = fileP;
        try {
            initial(scrIdx);
        } catch (IOException E) {
            E.printStackTrace(System.err);
        }
    }

    public DataFile() {}

    public DataFile(ArrayList<String> mkInformation, ArrayList<ArrayList> marker, ArrayList<Integer> statue, ArrayList<String> traitInformation, ArrayList<ArrayList> phenotype, int[] si) {
    	SNPID = mkInformation.toArray(new String[0]);
    	traitName = traitInformation.toArray(new String[0]);
    	initial1(marker, statue, phenotype, si);
    }

    private void initial1(ArrayList<ArrayList> marker, ArrayList statue, ArrayList<ArrayList> phenotype, int[] si) {
        TreeMap tempMap = new TreeMap();        
    	for (int i= 0; i < marker.size(); i++) {
    		ArrayList<String> geno = (ArrayList) marker.get(i).clone();
    		geno.add((String) statue.get(i).toString());
    		String[] content = geno.toArray(new String[0]);
            Subject sub = new Subject(content);
            sub.setID(i);
            sample.add(sub);
            //count alleles;
            for (int j = 0; j < content.length - 1; j++) {
                Integer Ii = new Integer(j);
                if (content[j].compareTo(PublicData.MissingGenotype) == 0) {
                    if (tempMap.containsKey(Ii)) {
                        ArrayList tempArray = (ArrayList) tempMap.get(Ii);
                        tempArray.add(new Integer(i));
                    } else {
                        ArrayList tempArray = new ArrayList();
                        tempArray.add(new Integer(i));
                        tempMap.put(Ii, tempArray);
                    }
                }
            }
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

        for (int i = 0; i < sample.size(); i++ ) {
            Subject sub = (Subject) sample.get(i);
            ArrayList<Double> traits = phenotype.get(i);
            String[] ps = traits.toArray(new String[0]);
            sub.addScore(ps);
        }
        setPhenotypeIndex(si);
        double[] os = calculateDefaultMu();
        TraitSummary(os);
        setScore(os);
    }
    

    public DataFile(IMPopulation imp, int[] scrIdx) {
        markerFile = "simulationMarker.txt";
        phenotypeFile = "simulationPhenotype.txt";

        //---------
        //read the marker names
        //---------
        SNPID = new String[imp.MDRMarkerNumber()];
        String[] m = imp.MDRMarkerName();
        int c = 0;
        System.arraycopy(m, 0, SNPID, c, m.length);

        GregorianCalendar calendar = new GregorianCalendar();
        Random rnd = new Random();
        rnd.setSeed(calendar.get(GregorianCalendar.SECOND));
        String cs = imp.getCrossParameter();
        TreeMap tempMap = new TreeMap();
        int count = 0;
        int[][][] marker = imp.Markers();
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            int idx = i - marker.length * (i / marker.length);
            String[] content = new String[imp.MarkerNumber() + 1]; //"+1" operation is to coordinate with the construction function of Subject, where the last one will be chopped off and set as status.
            c = 0;
            for (int k = 0; k < marker[idx].length; k++) {
                for (int h = 0; h < marker[idx][k].length; h++) {
                    content[c + h] = Integer.toString(marker[idx][k][h]);
                }
                c += marker[idx][k].length;
            }
            content[content.length - 1] = rnd.nextBoolean() ? new String("1") : new String("0");
            Subject sub = new Subject(content);
            sub.setID(count);
            sample.add(sub);
            for (int k = 0; k < content.length; k++) {
                Integer Ii = new Integer(k);
                if (content[k].compareTo(PublicData.MissingGenotype) == 0) {
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
        c = 0;
        double[][] pheno = imp.getPhenotype();
        for (int i = 0; i < pheno.length; i++) {
            String[] ps = new String[pheno[i].length];
            for (int j = 0; j < pheno[i].length; j++) {
                ps[j] = Double.toString(pheno[i][j]);
            }
            Subject sub = (Subject) sample.get(c++);
            sub.addScore(ps);
        }

        setPhenotypeIndex(scrIdx);
        double[] os = calculateDefaultMu();
        TraitSummary(os);
        setScore(os);
    }

    public DataFile(IMReaderAbstract imra, int[] scrIdx) {
        markerFile = imra.getFile();
        phenotypeFile = null;

        //---------
        //read the marker names
        //---------
        SNPID = new String[imra.MarkerNumber()];
        String[][] m = imra.getMarkerName();
        int c = 0;
        for (int i = 0; i < m.length; i++) {
            System.arraycopy(m[i], 0, SNPID, c, m[i].length);
            c += m[i].length;
        }
        //-----------
        //read markers
        //-----------
        Random rnd = new Random();
        GregorianCalendar calendar = new GregorianCalendar();
        rnd.setSeed(calendar.get(GregorianCalendar.SECOND));
        IMReaderAbstract.CrossParameter cs = imra.getCrossParameter();
        TreeMap tempMap = new TreeMap();
        int count = 0;
        String[][][] marker = imra.getMarkers();
        for (int i = 0; i < imra.IndividualNumber(); i++) {
            int idx = i - marker.length * (i / marker.length);
            String[] content = new String[imra.MarkerNumber() + 1]; //"+1" operation is to coordinate with the construction function of Subject, where the last one will be chopped off and set as status.
            c = 0;
            for (int k = 0; k < marker[idx].length; k++) {
                System.arraycopy(marker[idx][k], 0, content, c, marker[idx][k].length);
                c += marker[idx][k].length;
            }
            content[content.length - 1] = rnd.nextBoolean() ? new String("1") : new String("0");
            Subject sub = new Subject(content);
            sub.setID(count);
            sample.add(sub);
            for (int k = 0; k < content.length; k++) {
                Integer Ii = new Integer(k);
                if (content[k].compareTo(PublicData.MissingGenotype) == 0) {
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
        c = 0;
        String[][] pheno = imra.getPhenotype();
        for (int i = 0; i < pheno.length; i++) {
            String[] ps = pheno[i];
            Subject sub = (Subject) sample.get(c++);
            sub.addScore(ps);
        }

        setPhenotypeIndex(scrIdx);
        double[] os = calculateDefaultMu();
        TraitSummary(os);
        setScore(os);
    }

    public void initial(int[] scrIdx) throws IOException {
        FileReader fr = new FileReader(markerFile);
        LineNumberReader l = new LineNumberReader(fr);
        BufferedReader b = new BufferedReader(l);
        String line;

        // --------------------------------------------------------------------
        // Read the titles.
        // --------------------------------------------------------------------
        if ((line = b.readLine()) == null) {
            System.out.println("It's empty.");
            System.exit(0);
        }
        String[] s = line.split(PublicData.delim);
        SNPID = new String[s.length - 1];
        System.arraycopy(s, 0, SNPID, 0, s.length - 1);
        // -------------------------------------------------------------------

        TreeMap tempMap = new TreeMap();
        int count = 0;
        while ((line = b.readLine()) != null) {
            String[] content = line.split(PublicData.delim);
            Subject sub = new Subject(content);
            sub.setID(count);
            sample.add(sub);
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

        setPhenotypeIndex(scrIdx);
        double[] os = calculateDefaultMu();
        TraitSummary(os);
        setScore(os);
    }

    protected void addScore() throws IOException {
        if (phenotypeFile == null) {
            return;
        }
        FileReader fr = new FileReader(phenotypeFile);
        LineNumberReader l = new LineNumberReader(fr);
        BufferedReader b = new BufferedReader(l);
        String line;
        if ((line = b.readLine()) == null) {
            System.out.println("It's empty.");
            System.exit(0);
        }
        traitName = line.split(PublicData.delim);
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
            if (ps.length != traitName.length) {
                System.err.println("The " + c + "th phenotype row does not match the total size of the phenotypes");
            }
            Subject sub = (Subject) sample.get(c++);
            sub.addScore(ps);
        }
    }

    public int getNumberTraits() {
        return traitStatistics.size();
    }

    public double[] getOffset() {
        return offset;
    }

    public void setPhenotypeIndex(int[] sI) {
        scrIndex = new int[sI.length];
        System.arraycopy(sI, 0, scrIndex, 0, sI.length);
    }

    public void setScore(double[] os) {
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);
    }

    public String getMarkerFileName() {
        return markerFile;
    }

    public String getPhenotypeFileName() {
        return phenotypeFile;
    }

    public void add(int idx, Object o) {
        sample.add(idx, o);
    }

    public double[] calculateDefaultMu() {
        double[] defaultMu = new double[scrIndex.length];
        for (int i = 0; i < scrIndex.length; i++) {
            double sum = 0;
            int c = 0;
            for (Iterator e = sample.iterator(); e.hasNext();) {
                DataFile.Subject s = (DataFile.Subject) e.next();
                Double sscore = s.getDoubleScore(i);
                if (sscore != null) {
                    sum += sscore;
                    c++;
                }
            }
            defaultMu[i] = sum / c;
        }
        return defaultMu;
    }

    /**
     * This method should be invoked when a new trait is selected for running.
     */
    public void TraitSummary(double[] os) {
        for (int i = 0; i < scrIndex.length; i++) {
            Integer I = new Integer(i);
            TraitStatistic t = (TraitStatistic) traitStatistics.get(I);
            if (t != null && Math.abs(t.getOffset() - os[i]) < PublicData.epsilon) {
                continue;
            }

            TraitStatistic ts = new TraitStatistic(os[i], i, sample);
            traitStatistics.put(I, ts);

            int posSub = ts.getTraitPositiveSubjects();
            int negSub = ts.getTraitNegativeSubjects();
            double posScr = ts.getTraitPositiveScore();
            double negScr = ts.getTraitNegativeScore();

//            if (posSub == 0) {
//                throw new DataFileException("Score " + phenotype[i] + "is abnormal that no any subject has a positive score.");
//            } else if (negSub == 0) {
//                throw new DataFileException("Score " + phenotype[i] + "is abnormal that no any subject has a negative score.");
//            } else if (posScr < PublicData.epsilon) {
//                throw new DataFileException("Abnormal magnitude of the sum of the positive score: " + posScr + " at phenotype " + phenotype[i] + " .");
//            } else if (Math.abs(negScr) < PublicData.epsilon) {
//                throw new DataFileException("Abnormal magnitude of the sum of the negative score " + negScr + "  at phenotype " + phenotype[i] + " .");
//            }
            double r = posScr / Math.abs(negScr) > 1 ? posScr / Math.abs(negScr) : Math.abs(negScr) / posScr;
//            if (r > PublicData.tooBigRatio) {
//                throw new DataFileException("Abnormal ratio of positive score v.s. negative score.");
//            }
            r = posSub / Math.abs(negSub) > 1 ? posSub / Math.abs(negSub) : Math.abs(negSub) / posSub;
//            if (r > PublicData.tooBigRatio) {
//                throw new DataFileException("Abnormal ratio of subjects with positive score v.s. with negative score.");
//            }
        }
    }

    /**
     * Get the ratio of total positive scores to total negative scores.
     * 
     * @param scrIdx
     * @param markerIdx
     * @return
     * @throws algorithm.DataFileException
     */
    public double getTraitRatio(int scrIdx, int[] markerIdx) throws DataFileException {
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
            throw new DataFileException("too small denominator.");
        }
        return posScr / Math.abs(negScr);
    }

    public double defaultScore(int[] markerIdx, int[] scrIdx, int sI) throws DataFileException {
        TraitStatistic ts = (TraitStatistic) traitStatistics.get(new Integer(sI));
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
            double scr = s.getDoubleScore(scrIdx[sI]);
            if (scr >= 0) {
                posScr -= scr;
                posSub--;
            } else {
                negScr -= scr;
                negSub--;
            }
        }
        if (posScr + Math.abs(negScr) < PublicData.epsilon) {
            throw new DataFileException("too small denominator.");
        }
        return posScr / (posScr + Math.abs(negScr));
    }

    public int getMarkerNum() {
        return SNPID.length;
    }

    public void PickMarker(double cM) {

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

    public Object remove(int index) {
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

    public void switchPermutation(boolean flag) {
        isPermutation = flag;
    }

    public void setShuffledIDs(ArrayList shuffledids) {
        if (isPermutation) {
            shuffledPhenoIDs = shuffledids;
            for (int i = 0; i < shuffledPhenoIDs.size(); i++) {
                Subject sub1 = (Subject) sample.get(i);
                Subject sub2 = (Subject) sample.get(((Integer) shuffledPhenoIDs.get(i)).intValue());
                ArrayList<Double> s = sub2.getScore();
                sub1.addShuffledScore(s);
            }
        }
    }

    public ArrayList getSample() {
        return sample;
    }

    public int size() {
        return sample.size();
    }
}
