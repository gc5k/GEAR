package mdr.data;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import mdr.data.DataFile;
import im.population.MarkerInformation;
import im.population.IMPopulation;

import publicAccess.PublicData;
import publicAccess.ToolKit;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class ExDataFile extends DataFile {

    protected ArrayList MarkerInfor;
    private IMPopulation imp;
    public ExDataFile(IMPopulation im, int[] scrIdx, double cM) {
        super();
        markerFile = "simulationMarker.txt";
        phenotypeFile = "simulationPhenotype.txt";

        imp = im;
        if (cM > 0) {
            imp.PickMarkers(cM);
        } else {
            imp.PickMarkers(cM);
        }
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
            String[] content = new String[imp.MarkerNumber() + 2]; //"+2" operation is to coordinate with the construction function of Subject, where the last second one is environemnt and the last one is chopped off and set as status.
            c = 0;
            for (int k = 0; k < marker[idx].length; k++) {
                for (int h = 0; h < marker[idx][k].length; h++) {
                    content[c + h] = Integer.toString(marker[idx][k][h]);
                }
                c += marker[idx][k].length;
            }
            content[content.length - 2] = (String) imp.NestedEnvironment(i);
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

    public void PickMarker(double cM) {
        if (cM < 0) {
            imp.PickAllMarkers();
        } else {
            imp.PickMarkers(cM);
        }
    }

    public int getMarkerNum() {
        DataFile.Subject sub = (DataFile.Subject) sample.get(0);
        int numAttri = sub.size();
        return numAttri - 1;
    }

    public IMPopulation getIMPopulation() {
        return imp;
    }

    public MarkerInformation CombinationGeneticInformationOriginal(String com) {
        int[] SNPIndex = ToolKit.StringToIntArray(com);
        int[][] ChrInt = new int[SNPIndex.length][2];
        double[] coordinate = new double[SNPIndex.length];
        double[][] dis = imp.DistanceOriginal();
        for (int i = 0; i < ChrInt.length; i++) {
            int c = 0;
            boolean flag = false;
            for (int j = 0; j < dis.length; j++) {
                for (int k = 0; k < dis[j].length; k++ ) {
                    if (c == SNPIndex[i]) {
                        flag = true;
                        ChrInt[i][0] = j;
                        ChrInt[i][1] = k;
                        coordinate[i] = dis[j][k];
                    }
                    c++;
                }
                if (flag ) {
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
        double[][] dis = imp.DistanceSelected();
        for (int i = 0; i < ChrInt.length; i++) {
            int c = 0;
            boolean flag = false;
            for (int j = 0; j < dis.length; j++) {
                for (int k = 0; k < dis[j].length; k++ ) {
                    if (c == SNPIndex[i]) {
                        flag = true;
                        ChrInt[i][0] = j;
                        ChrInt[i][1] = k;
                        coordinate[i] = dis[j][k];
                    }
                    c++;
                }
                if (flag ) {
                    break;
                }
            }
        }
        MarkerInformation mkinfor = new MarkerInformation(ChrInt, coordinate);
        return mkinfor;
    }

    public double[][] DistanceOriginal() {
        return imp.DistanceOriginal();
    }

    public double[][] DistanceSelected() {
        return imp.DistanceSelected();
    }

    public int MDRMarkerIndexOriginal(int chr, int idx) {
        int index = 0;
        boolean flag = false;
        double[][] distance = imp.DistanceOriginal();
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                if( i == chr && j == idx ) {
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
}
