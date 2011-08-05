package mdr.algorithm;

import java.util.ArrayList;

import publicAccess.PublicData;
import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class CombinationGeneratorII {

	ArrayList<String> Comb;
    private int start;
    private int N;
    private int[] a;

    private int[] exclude = null;
    private int[] include = null;
    public CombinationGeneratorII(int s, int n, int[] in, int[] ex) {
        start = s;
        N = n;
        int len = n;
        if(ex != null) {
        	exclude = new int[ex.length];
        	System.arraycopy(ex, 0, exclude, 0, ex.length);
        	len -= exclude.length;
        }
        if(in != null) {
        	include = new int[in.length];
        	System.arraycopy(in, 0, include, 0, in.length);
        	len -= include.length;
        	start -= include.length;
        }
        a = new int[len];
        int idx = 0;
        for (int i = 0; i < N; i++) {
        	boolean flag = true;
        	if( include != null) {
        		for(int j = 0; j < include.length; j++) {
        			if(include[j] == i) {
        				flag = false;
        				break;
        			}
        		}
        	}
        	if( exclude != null) {
        		for(int j = 0; j < exclude.length; j++) {
        			if(exclude[j] == i) {
        				flag = false;
        				break;
        			}
        		}
        	}
        	if(flag) {
        		a[idx++] = i;
        	}
        }
    }

    public int getStart() {
        return start;
    }

    public ArrayList<String> getCombination() {
    	return Comb;
    }

    public void generateCombination() {
    	combine(a, a.length, start);
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
                if (include != null) {
                	for (int i = 0; i < include.length; i++) {
                		combi.append(include[i]);
                		combi.append(PublicData.seperator);
                	}
                }
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

    public static void main(String[] args) {
    	int[] in = {1, 2};
    	int[] ex = {8, 9};
    	CombinationGeneratorII cg = new CombinationGeneratorII(2, 10, in, ex);
    	cg.generateCombination();
    	System.out.println("done");
    }
}
