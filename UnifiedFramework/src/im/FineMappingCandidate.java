package im;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import im.IMSet;
import im.IMToolKit;
import im.MDR2IM;
import im.population.IMPopulation;
import im.population.MarkerInformation;

import mdr.data.ExDataFile;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class FineMappingCandidate {

    ArrayList sigPoints;
    ArrayList<ArrayList> sigCombinations;
    ExDataFile edf;

    MDR2IM mdr2im;
    ArrayList<ArrayList> fineMappingCombinations;
    double windowSize;

    public FineMappingCandidate(ArrayList sigpoints, ExDataFile ed, double ws) {
        sigPoints = sigpoints;
        edf = ed;
        windowSize = ws;
        IMPopulation imp = edf.getIMPopulation();
        mdr2im = new MDR2IM(sigPoints, imp);
    }

    public ArrayList<ArrayList> candidateCombinations() {
        sigCombinations = mdr2im.summary();
        if (sigCombinations == null) {
            return null;
        }
        fineMappingCombinations = new ArrayList();
        TreeSet candidates = new TreeSet();
        for (Iterator e = sigCombinations.iterator(); e.hasNext();) {
            ArrayList sigCom = (ArrayList) e.next();
            TreeSet hs = new TreeSet();
            for (Iterator e1 = sigCom.iterator(); e1.hasNext();) {
                String com = (String) e1.next();
                MarkerInformation mi = mdr2im.CombinationGeneticInformation(com);
//                System.out.println(com + "***");
                String[] combinations = searchCandidateCombination(mi);
                for (int i = 0; i < combinations.length; i++) {
                    hs.add(combinations[i]);
                }
            }

            ArrayList Coms = new ArrayList();
            for (Iterator e1 = hs.iterator(); e1.hasNext(); ) {
                Coms.add((String) e1.next());
            }
            candidates.addAll(hs);
        }

        int[][][] location = new int[candidates.size()][][];

        HashMap coms = new HashMap();
        int idx = 0;
        for (Iterator e = candidates.iterator(); e.hasNext(); ) {
            String com = (String) e.next();
            MarkerInformation mi = edf.CombinationGeneticInformationOriginal(com);
            int[][] c = mi.atLocation();
            location[idx] = new int[c.length][2];
            StringBuffer loc = new StringBuffer();
            for (int i = 0; i < c.length; i++) {
                location[idx][i][0] = c[i][0];
                location[idx][i][1] = c[i][1];
                loc.append(c[i][0] + IMToolKit.separator + c[i][1] + IMToolKit.separator);
            }
            coms.put(loc.toString(), com);
//            System.out.println(com);
            idx++;
        }

        IMSet imSet = new IMSet();
        HashSet<TreeSet> locSet = imSet.findConnectedComponents(location);
        ArrayList<ArrayList> SigComs = new ArrayList();
        for (Iterator e = locSet.iterator(); e.hasNext();) {
            TreeSet ele = (TreeSet) e.next();
            ArrayList sigCom = new ArrayList();
            for (Iterator f = ele.iterator(); f.hasNext();) {
                int i = ((Integer) f.next()).intValue();
                StringBuffer sb = new StringBuffer();
                for (int j = 0; j < location[i].length; j++) {
                    sb.append(String.valueOf(location[i][j][0]) + IMToolKit.separator + location[i][j][1] + IMToolKit.separator);
                }
                String key = (String) coms.get(sb.toString());
                sigCom.add(key);
            }
            SigComs.add(sigCom);
        }
        return SigComs;
    }

    private String[] searchCandidateCombination(MarkerInformation mi) {
        int[][] chr = mi.atLocation();
        double[] co = mi.atCoordinate();
        int[][] s = new int[co.length][];
        double[][] distance = edf.DistanceOriginal();
        int[] ai = new int[co.length];
        for (int i = 0; i < chr.length; i++) {
            ai[i] = edf.MDRMarkerIndexOriginal(chr[i][0], chr[i][1]);
        }
        ArrayList sms = new ArrayList();
        for (int i = 0; i < chr.length; i++) {
            ArrayList sm = new ArrayList();
            for (int j = 0; j < distance[chr[i][0]].length; j++) {
                if (Math.abs(distance[chr[i][0]][j] - distance[chr[i][0]][chr[i][1]]) < windowSize) {
                    Integer index = edf.MDRMarkerIndexOriginal(chr[i][0], j);
                    if (i == 0) { // first marker in the combintion
                        for (int k = i + 1; k < chr.length; k++) {
                            if (index > ai[k]) {
                                break;
                            }
                        }
                    } else if (i == (chr.length - 1)) {// last marker in the combination
                        for (int k = 0; k < chr.length - 1; k++) {
                            if (index < ai[k]) {
                                break;
                            }
                        }
                    } else { // markers in the middle of the combination
                        for (int k = 0; k < i ; k++) {
                            if (index < ai[k]) {
                                break;
                            }
                        }
                        for (int k = i + 1; k < chr.length; k++) {
                            if (index > ai[k]) {
                                break;
                            }
                        }
                    }
                    sm.add(index);
                }
            }
            sms.add(sm);
        }

        int[][] comElements = new int[sms.size()][];
        for (int i = 0; i < sms.size(); i++) {
            ArrayList sm = (ArrayList) sms.get(i);
            comElements[i] = new int[sm.size()];
            for (int j = 0; j < sm.size(); j++) {
                Integer index = (Integer) sm.get(j);
                comElements[i][j] = index.intValue();
            }
        }

        String[] Xcom = IMToolKit.XCombination(comElements);
        return Xcom;
    }
}
