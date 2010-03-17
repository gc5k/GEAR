package im.population;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class MarkerInformation {
    int[][] chr;
    double[] coordinate;
    public MarkerInformation(int[][] c, double[] co) {
        chr = c;
        coordinate = co;
    }

    public int[][] atLocation() {
        return chr;
    }

    public double[] atCoordinate() {
        return coordinate;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < chr.length; i++) {
           sb.append(chr[i][0] + "," + chr[i][1] + "," + coordinate[i]+";");
        }
        return sb.toString();
    }
}
