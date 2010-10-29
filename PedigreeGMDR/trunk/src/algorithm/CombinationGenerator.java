package algorithm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Iterator;
import java.util.Set;

import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen
 */
public class CombinationGenerator {

    private HashMap<Integer, List> combination = new HashMap();

    private int start;
    private int end;
    private int N;
    private int[] a;

    public CombinationGenerator(int s, int e, int n) {
        start = s;
        end = e;
        N = n;
        a = new int[N];
        for (int i = 0; i < N; i++) {
            a[i] = i;
        }
    }

    public void fill(Integer order, List com) {
    	combination.put(order, com);
    }

    public List get(Object o) {
        return combination.get(o);
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int size() {
        return combination.size();
    }

    public void generateCombination() {
        for (int i = start; i <= end; i++) {
            combination.put(new Integer(i), combine(a, N, i));
        }
    }

    public void print() {
        Set keys = combination.keySet();
        for (Iterator i = keys.iterator(); i.hasNext();) {
            List c = (List) i.next();
            for (Iterator j = c.iterator(); j.hasNext();) {
                System.out.println((String) j.next());
            }
            System.out.println();
        }
    }

    public List combine(int[] a, int n, int m) {
        List Comb = new ArrayList();
        int[] order = new int[m + 1];
        for (int i = 0; i <= m; i++) {
            order[i] = i - 1;
        }
        int count = 0;
        int k = m;
        boolean flag = true;
        String combi;
        while (order[0] == -1) {
            if (flag) {
                combi = new String();
                for (int i = 1; i <= m; i++) {
                    if (i != m) {
                        combi += (new Integer(a[order[i]])).toString() + PublicData.seperator;
                    } else {
                        combi += (new Integer(a[order[i]])).toString();
                    }
                }
                Comb.add(combi);
                count++;
                flag = false;
            }
            order[k]++;
            if (order[k] == n) {
                order[k--] = 0;
                continue;
            }
            if (k < m) {
                order[k + 1] = order[k];
                k++;
                continue;
            }
            if (k == m) {
                flag = true;
            }
        }
        return Comb;
    }
}
