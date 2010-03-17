package im;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;
import im.population.IMPopulation;
import im.population.MarkerInformation;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class MDR2IM {
    ArrayList sigPoints;
    HashMap coms;
    IMPopulation imp;
    int[][] ChrInt;
    double[] coordinate;
    int[][][] location;//length of the first dimension = number of significant points; length of the second dimension = order of interaction; length of the third dimension = 2, chr and interval, respectively.
    public MDR2IM(ArrayList sp, IMPopulation im) {
        sigPoints = sp;
        imp = im;
        coms = new HashMap();
        Convert2GeneticInformation();
    }

    public MarkerInformation CombinationGeneticInformation(String com) {
        int[] SNPIndex = ToolKit.StringToIntArray(com);
        ChrInt = new int[SNPIndex.length][2];
        ChrInt(SNPIndex);
        coordinate = new double[SNPIndex.length];
        double[][] distance = imp.DistanceOriginal();
        for (int i = 0; i < ChrInt.length; i++) {
            coordinate[i] = distance[ChrInt[i][0]][ChrInt[i][1]];
        }
        MarkerInformation mkinfor = new MarkerInformation(ChrInt, coordinate);
        return mkinfor;
    }

    public ArrayList<ArrayList> summary() {
        IMSet imSet = new IMSet();
        HashSet locSet = imSet.findConnectedComponents(location);
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

    private void Convert2GeneticInformation() {
        location = new int[sigPoints.size()][][];
        for(int i = 0; i < sigPoints.size(); i++) {
            String com = (String) sigPoints.get(i);
            int[] SNPIndex = ToolKit.StringToIntArray(com);
            ChrInt = new int[SNPIndex.length][2];
            ChrInt(SNPIndex);

            StringBuffer loc = new StringBuffer();
            location[i] = new int[ChrInt.length][2];
            for(int j = 0; j < ChrInt.length; j++) {
                location[i][j][0] = ChrInt[j][0];
                location[i][j][1] = ChrInt[j][1];
                loc.append(ChrInt[j][0] + IMToolKit.separator + ChrInt[j][1] + IMToolKit.separator);
            }
            coms.put(loc.toString(), com);
        }
        return;
    }

    private void ChrInt(int[] snpindex) {
        ArrayList mdrSelectedMarkerIndex = imp.MDRPickedMarkerIdx();
        for (int i = 0; i < snpindex.length; i++) {
            for (int j = 0; j < mdrSelectedMarkerIndex.size(); j++ ) {
                IMPopulation.MarkerIndex mi = (IMPopulation.MarkerIndex) mdrSelectedMarkerIndex.get(j);
                if (snpindex[i] == j) {
                    ChrInt[i][0] = mi.getChromosome();
                    ChrInt[i][1] = mi.getIndex();
                    break;
                }
            }
        }
    }
}
