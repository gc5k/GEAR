package im;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Iterator;
import java.util.Set;
import im.IMSet;
import mdr.heterogeneity.IMHeteroLinearMergeSearch;
import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Summary {

    double[] peakStat;
    int[][][] location;//length of the first dimension = number of significant points; length of the second dimension = order of interaction; length of the third dimension = 2, chr and interval, respectively.
    HashMap coms;
    HashMap InvertComs;
    IMHeteroLinearMergeSearch ihm;
    ArrayList<ArrayList> SigComs;

    public Summary(IMHeteroLinearMergeSearch imhlms, double[] ps) {
        ihm = imhlms;
        if (ps != null) {
            peakStat = new double[ps.length];
            System.arraycopy(ps, 0, peakStat, 0, ps.length);
            Arrays.sort(peakStat);
        }
    }

    public void significantpoint(Double t) {
        double threshold = 0;
        if (peakStat != null) {
            int idxPS = (new Double(peakStat.length * 0.95)).intValue();
            threshold = peakStat[idxPS];
        }
        if (t != null) {
            threshold = t.doubleValue();
        }

        ArrayList sp = ihm.SignificantPoints(threshold);
        int order = ihm.getOrder();
        location = new int[sp.size()][order][2];
        coms = new HashMap();
        InvertComs = new HashMap();
        int idx = 0;
        for (Iterator e = sp.iterator(); e.hasNext();) {
            String str = (String) e.next();
            String s[] = str.split(PublicData.seperator);
            StringBuffer sb = new StringBuffer();
            StringBuffer sk = new StringBuffer();
            for (int i = 0; i < order; i++) {
                location[idx][i][0] = Integer.parseInt(s[i * 2]);
                sb.append(s[i * 2]);
                sb.append(IMToolKit.separator);
                location[idx][i][1] = Integer.parseInt(s[i * 2 + 1]);
                sb.append(s[i * 2 + 1]);
                sb.append(IMToolKit.separator);
            }
            for (int i = 0; i < order; i++) {
                sk.append(s[s.length - order + i]);
                if (i < (order - 1)) {
                    sk.append(IMToolKit.separator);
                }
            }
            coms.put(sb.toString(), sk.toString());
            InvertComs.put(sk.toString(), sb.toString());
            idx++;
        }
        IMSet imSet = new IMSet();
        HashSet locSet = imSet.findConnectedComponents(location);
        SigComs = new ArrayList();
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
    }

    public ArrayList<ArrayList> getSignificantCombinations() {
        return SigComs;
    }

    public HashMap getSignificantCombinationMap() {
        return InvertComs;
    }
}
