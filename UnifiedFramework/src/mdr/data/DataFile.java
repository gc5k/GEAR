package mdr.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;

import publicAccess.PublicData;

import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class DataFile {

    protected ArrayList<DataFile.Subject> sample = NewIt.newArrayList();
    protected static boolean isPermutation = false;
    protected String markerFile;
    protected String phenotypeFile;
    protected String[] SNPID;
    protected String[] traitName;

    protected static int currScrIdx = -1;

    public static class Subject {
        int id;
        private byte[] genotypes;
        private byte status;
        private double[] score;
        private double defaultScore;
        

        public Subject(String[] content) {
            if (content[content.length - 1].equals(PublicData.MissingPhenotypeValue)) {
                this.status = -1;
            } else {
                status = Byte.parseByte(content[content.length - 1]);
            }
            genotypes = new byte[content.length - 1];
            for(int i = 0; i < genotypes.length; i++) {
            	genotypes[i] = Byte.parseByte(content[i]);
            }
        }

        public Subject(byte[] g, byte s) {
        	genotypes = g;
        	status = s;
        }

        private void setDefaultScore(double m) {
        	defaultScore = status - m;
        }

        public void setScore(double s) {
        	if(currScrIdx<0) {
        		defaultScore = s;
        	} else {
        		score[currScrIdx] = s;
        	}
        }
        
        public void addScore(double[] s) {
        	score = s;
        }

        public void addScore(String[] s) {
            score = new double[s.length];
            for (int i = 0; i < s.length; i++) {
                if (s[i].equals(PublicData.MissingPhenotypeValue)) {
                    score[i] = Double.NaN;
                } else {
                    score[i] = Double.parseDouble(s[i]);
                }
            }
        }


        public byte getGenotype(int idx) {
            return genotypes[idx];
        }

        public double getSelectedScore() {
        	if (currScrIdx < 0) {
        		return defaultScore;
        	} else if (score[currScrIdx] != Double.NaN) {
                return score[currScrIdx];
            } else {
                return status;
            }
        }

        public byte getStatus() {
            return status;
        }

        public Integer getIntegerID() {
            return new Integer(id);
        }

        public int size() {
            return genotypes == null ? 0 : genotypes.length;
        }

        public byte get(int idx) {
            return genotypes == null ? null : genotypes[idx];
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

            for (Double s:score) {
                b.append(s);
                b.append(PublicData.delim);
            }
            return b.toString();
        }
    }

    public DataFile(String fileM, String fileP) {
        markerFile = fileM;
        phenotypeFile = fileP;
        try {
            initial();
        } catch (IOException E) {
            E.printStackTrace(System.err);
        }
    }

    public DataFile() {}

    public DataFile(String[] mkInformation, byte[][] marker, byte[] status, String[] traitInformation, double[][] phenotype) {
    	SNPID = mkInformation;
    	traitName = traitInformation;
    	initial1(marker, status, phenotype);
    }

    private void initial1(byte[][] marker, byte[] status, double[][] score) {
    	double mu = 0;
    	for (int i= 0; i < marker.length; i++) {
    		byte[] g = marker[i];
            Subject sub = new Subject(g, status[i]);
            mu += status[i];
            sub.setID(i);
            sample.add(sub);
    	}
    	mu /= status.length;

        for (int i = 0; i < score.length; i++ ) {
            Subject sub = (Subject) sample.get(i);
            sub.setDefaultScore(mu);
            sub.addScore(score[i]);
        }
    }

    public void initial() throws IOException {
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

        int count = 0;
        while ((line = b.readLine()) != null) {
            String[] content = line.split(PublicData.delim);
            Subject sub = new Subject(content);
            sub.setID(count);
            sample.add(sub);

            count++;
        }
        addScore();
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
        ArrayList<String> pl = NewIt.newArrayList();

        while ((line = b.readLine()) != null) {
            pl.add(line);
        }
        if (pl.size() != size()) {
            System.err.println("The size of the sample unmatch the size of the score");
            return;
        }

        int c = 0;
        for (String p:pl) {
            String[] ps = p.split(PublicData.delim);
            if (ps.length != traitName.length) {
                System.err.println("The " + c + "th phenotype row does not match the total size of the phenotypes");
            }
            Subject sub = (Subject) sample.get(c++);
            sub.addScore(ps);
        }
    }

    public String getMarkerFileName() {
        return markerFile;
    }

    public void setScore(double[] s) {
    	for(int i = 0; i < size(); i++) {
    		Subject sub = get(i);
    		sub.setScore(s[i]);
    	}
    }

    public static void setScoreIndex(int i) {
    	currScrIdx = i;
    }

    public String getPhenotypeFileName() {
        return phenotypeFile;
    }

    public void add(int idx, DataFile.Subject o) {
        sample.add(idx, o);
    }

    public int getNumMarker() {
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

    public Subject get(int idx) {
        return sample.get(idx);
    }

    public Subject remove(int index) {
        return sample.remove(index);
    }

    public void print() {
        for (DataFile.Subject e:sample) {
            System.out.println(e);
        }
    }

    public ArrayList<DataFile.Subject> getSample() {
        return sample;
    }

    public int size() {
        return sample.size();
    }

    public String getTraitName(int idx) {
    	return traitName[idx];
    }
}
