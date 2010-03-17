package im;

/**
 *
 * @author Guo-Bo Chen
 */
public class IMToolKit {

    public static String separator = ",";
    public static String chr_prefix = "chr_";
    public static String interval_prefix = "interval_";
    public static String MissingGenotype = ".";
    public static String MissingTrait = ".";
    public static int AA = 2;
    public static int Aa = 1;
    public static int aa = 0;
    public static int B1 = 0;
    public static int B2 = 1;
    public static int F2 = 2;
    public static int IF2 = 3;
    public static int DH = 4;
    public static int[] populationindicator = {B1, B2, F2, IF2, DH};

    public static double Felsenstein(double d, int mf) {//k=1,haldane; k=0, kosambi; Genetics 91:769-775
        return 0.5 * (1 - Math.exp(2 * (mf - 2) * d)) / (1 - (mf - 1) * Math.exp(2 * (mf - 2) * d));
    }

    public static String MakeKey(String chr, String mk, int c, int m) {
        String key = new String(chr + Integer.toString(c) + mk +
                Integer.toString(m));
        return key;
    }

    //XCombination was developed by Zhi-Xiang Zhu
    public static int[][] XCombination(int[] orders) {
        int len = 1;
        for (int i = 0; i < orders.length; ++i) {
            len *= orders[i];
        }

        int[][] indices = new int[len][];

        indices[0] = new int[orders.length];
        for (int i = 0; i < indices[0].length; ++i) {
            indices[0][i] = 0;
        }

        for (int i = 1; i < indices.length; ++i) {
            indices[i] = new int[orders.length];
            int j = indices[i].length;
            do {
                indices[i][--j] = 0;
            } while (j >= 0 && indices[i - 1][j] == orders[j] - 1);
            if (j >= 0) {
                indices[i][j] = indices[i - 1][j] + 1;
            }
            while (--j >= 0) {
                indices[i][j] = indices[i - 1][j];
            }
        }
        return indices;
    }

    //XCombination was developed by Zhi-Xiang Zhu
    public static String[] XCombination(int[][] groups) {
        int[] orders = new int[groups.length];
        for (int i = 0; i < orders.length; ++i) {
            orders[i] = groups[i].length;
        }
        int[][] indices = XCombination(orders);
        String[] combinations = new String[indices.length];
        for (int i = 0; i < indices.length; ++i) {
            combinations[i] = String.valueOf(groups[0][indices[i][0]]);
            for (int j = 1; j < indices[i].length; ++j) {
                combinations[i] += ",";
                combinations[i] += String.valueOf(groups[j][indices[i][j]]);
            }
        }
        return combinations;
    }
}

