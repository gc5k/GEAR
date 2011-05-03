package mdr.algorithm;

import java.util.ArrayList;

import publicAccess.PublicData;
import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class CombinationGenerator {

	ArrayList<String> Comb;
    private int start;
    private int N;
    private int[] a;

    public CombinationGenerator(int s, int n) {
        start = s;
        N = n;
        a = new int[N];
        for (int i = 0; i < N; i++) {
            a[i] = i;
        }
    }

    public int getStart() {
        return start;
    }

    public ArrayList<String> getCombination() {
    	return Comb;
    }

    public void generateCombination() {
    	combine(a, N, start);
    }

    private ArrayList<String> combine(int[] a, int n, int m) {
        Comb = NewIt.newArrayList();
        int[] order = new int[m + 1];
        for (int i = 0; i <= m; i++) {
            order[i] = i - 1;
        }
        int count = 0;
        int k = m;
        boolean flag = true;
        StringBuilder combi;
        while (order[0] == -1) {
            if (flag) {
                combi = new StringBuilder();
                for (int i = 1; i <= m; i++) {
                    combi.append(Integer.toString(a[order[i]]));
                    if (i != m) {
                        combi.append(PublicData.seperator);
                    }
                }
                Comb.add(combi.toString());
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
