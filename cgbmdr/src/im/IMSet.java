package im;

import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

/**
 *
 * @author Zhi-Xiang Zhu
 */
public class IMSet {

    /**
     * This method performs the algorithm of finding connected components of a graph.
     * The judgement of whether two vetices are connected are defined in the helper function
     * isConnected(). A component has to be constructed as a TreeSet because the input
     * are passed in in some sorted order. This order of points should be kept in each
     * component.
     *
     * @param points    a set of input points. Each point is a two-dimension int array.
     * @return          a HashSet with each element is type of java.util.HashSet<Integer>
     */
    public static HashSet<TreeSet> findConnectedComponents(int points[][][]) {
        assert points != null;

        TreeSet [] componentIndex = new TreeSet [points.length];
        for (int i = 0; i < points.length; i++) {
            componentIndex[i] = new TreeSet<Integer>();
            componentIndex[i].add(i);
            for (int j = 0; j < i; j++) {
                if (componentIndex[i] != componentIndex[j] && isConnected(points[i], points[j])) {
                    // union operation
                    for (Iterator iter = componentIndex[i].iterator(); iter.hasNext();) {
                        int k = ((Integer) iter.next()).intValue();
                        componentIndex[k] = componentIndex[j];
                        componentIndex[j].add(k);
                    }
                }
            }
        }

        HashSet<TreeSet> disjointSet = new HashSet<TreeSet>();
        for (int i = 0; i < points.length; i++) {
            if (!disjointSet.contains(componentIndex[i]))
                disjointSet.add(componentIndex[i]);
        }
        return disjointSet;
    }

    /**
     * helper method, determining whether two points are connected.
     *
     * @param p1
     * @param p2
     * @return
     */
    private static boolean isConnected(int p1[][], int p2[][]) {
        assert p1.length == p2.length;
//        for (int i = 0; i < p1.length; ++i) {
//            assert p1[i].length == 2;
//            for (int j = 0; j < p2.length; ++j) {
//                assert p2[j].length == 2;
//                if (p1[i][0] == p2[j][0] &&
//                        Math.abs(p1[i][1] - p2[j][1]) <= 1)
//                    return true;
//            }
//        }
//
        for (int i = 0; i < p1.length; i++) {
            assert p1[i].length == 2;
            assert p2[i].length == 2;
            if (p1[i][0] != p2[i][0] || Math.abs(p1[i][1] - p2[i][1]) > 1) {
                return false;
            }
        }
        return true;
    }
}
