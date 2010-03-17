package im.reader;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Random;

import im.IMToolKit;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Imputation {

    IMReaderAbstract imReader;
    String[][] markercoding;
    String[][][] marker;
    double[][] distance;

    Random rnd;

    public Imputation(IMReaderAbstract imr) {
        imReader = imr;
        markercoding = imReader.crp.getMarkerCoding();
        marker = imReader.getMarkers();
        distance = imReader.getDistance();
        rnd = new Random();
        rnd.setSeed(imReader.crp.getSampleSize());
    }

    public void Impute() {
        for (int i = 0; i < marker.length; i++) {
            for (int j = 0; j < marker[i].length; j++) {
                for (int k = 0; k < marker[i][j].length; k++) {
                    if (   marker[i][j][k].compareTo(markercoding[IMConstant.AA][2]) != 0
                        && marker[i][j][k].compareTo(markercoding[IMConstant.Aa][2]) != 0
                        && marker[i][j][k].compareTo(markercoding[IMConstant.aa][2]) != 0
                        ) {
                        impute(i, j, k);
                    }
                }
            }
        }
    }

    public String impute(int idx, int cI, int mI) {
        Pattern pf2 = Pattern.compile("^f2", Pattern.CASE_INSENSITIVE);
        Pattern pif2 = Pattern.compile("^if2", Pattern.CASE_INSENSITIVE);
        Matcher match;

        String mk = null;
        String[] fmk = {null, null};
        double dis_M1Q = 0;
        double dis_QM2 = 0;
        double dis_M1M2 = 0;

        int missing_type;
        if (marker[idx][cI][mI].compareTo(markercoding[IMConstant.A_][2]) == 0) {
            missing_type = IMConstant.A_;
        } else if (marker[idx][cI][mI].compareTo(markercoding[IMConstant.a_][2]) == 0) {
            missing_type = IMConstant.a_;
        } else {
            missing_type = IMConstant.Missing;
        }

        int i1 = mI;
        while (--i1 >= 0) {
            if ((marker[idx][cI][i1].compareTo(markercoding[IMConstant.AA][2]) == 0) ||
                    (marker[idx][cI][i1].compareTo(markercoding[IMConstant.Aa][2]) == 0) ||
                    (marker[idx][cI][i1].compareTo(markercoding[IMConstant.aa][2]) == 0)) {
                fmk[0] = marker[idx][cI][i1];
                break;
            }
        }
        if (i1 >= 0) {
            dis_M1Q = distance[cI][mI] - distance[cI][i1];
        }

        int i2 = mI;
        while (++i2 < marker[idx][cI].length) {
            if ((marker[idx][cI][i2].compareTo(markercoding[IMConstant.AA][2]) == 0) ||
                    (marker[idx][cI][i2].compareTo(markercoding[IMConstant.Aa][2]) == 0) ||
                    (marker[idx][cI][i2].compareTo(markercoding[IMConstant.aa][2]) == 0)) {
                fmk[1] = marker[idx][cI][i2];
                break;
            }
        }
        if (i2 < marker[idx][cI].length) {
            dis_QM2 = distance[cI][i2] - distance[cI][mI];
        }

        if (i1 >= 0 && i2 < marker[idx][cI].length) {
            dis_M1M2 = distance[cI][i2] - distance[cI][mI];
        }

        match = pf2.matcher(imReader.crp.getCross());
        String imputed_marker = null;
        if (match.matches()) {
            imputed_marker = F2Imputation(missing_type, fmk, dis_M1M2, dis_M1Q, dis_QM2);
        }

        match = pif2.matcher(imReader.crp.getCross());
        if (match.matches()) {
            imputed_marker = F2Imputation(missing_type, fmk, dis_M1M2, dis_M1Q, dis_QM2);
        }
        marker[idx][cI][mI] = imputed_marker;
        return mk;
    }

    public String F2Imputation(int missing_type, String[] fmk, double dis_M1M2, double dis_M1Q, double dis_QM2) {
        String mk = null;
        String[] M = {markercoding[IMConstant.AA][2], markercoding[IMConstant.Aa][2], markercoding[IMConstant.aa][2]};
        double rM1M2 = IMToolKit.Felsenstein(dis_M1M2, imReader.chrp.getFunction());
        double rM1Q = IMToolKit.Felsenstein(dis_M1Q, imReader.chrp.getFunction());
        double rQM2 = IMToolKit.Felsenstein(dis_QM2, imReader.chrp.getFunction());

        double[] pM = {0, 0, 1};
        double p0 = 0, p1 = 0, p2 = 0;
        if (fmk[0] == null && fmk[1] == null) {
            pM[0] = 0.25;
            pM[1] = 0.75;
            pM[2] = 1;
        } else if (fmk[0] == null) {//only available the flanking marker at the right side
            double p = IMToolKit.Felsenstein(dis_QM2, imReader.chrp.getFunction());
            if (fmk[1].compareTo(markercoding[IMConstant.AA][2]) == 0) {//AA
                p0 = (1 - rQM2) * (1 - rQM2);
                p1 = 2 * (1 - rQM2) * rQM2;
                p2 = rQM2 * rQM2;
            } else if (fmk[1].compareTo(markercoding[IMConstant.Aa][2]) == 0) {//Aa
                p0 = (1 - rQM2) * rQM2;
                p1 = (1 - rQM2) * (1 - rQM2) + rQM2 * rQM2;
                p2 = (1 - rQM2) * rQM2;
            } else {//aa
                p0 = rQM2 * rQM2;
                p1 = 2 * (1 - rQM2) * rQM2;
                p2 = (1 - rQM2) * (1 - rQM2);
            }
            if (missing_type == IMConstant.A_) {
                pM[0] = p0 / (p0 + p1);
                pM[1] = 1;
                pM[2] = 1;
            } else if (missing_type == IMConstant.a_) {
                pM[0] = 0;
                pM[1] = p1 / (p1 + p2);
                pM[2] = 1;
            } else {
                pM[0] = p0;
                pM[1] = p0 + p1;
                pM[2] = 1;
            }
        } else if (fmk[1] == null) {// only avialable the flanking marker at the left side
            double p = IMToolKit.Felsenstein(dis_M1Q, imReader.chrp.getFunction());
            if (fmk[0].compareTo(markercoding[IMConstant.AA][2]) == 0) {//AA
                p0 = (1 - rM1Q) * (1 - rM1Q);
                p1 = 2 * (1 - rM1Q) * rM1Q;
                p2 = rM1Q * rM1Q;
            } else if (fmk[0].compareTo(markercoding[IMConstant.Aa][2]) == 0) {//Aa
                p0 = (1 - rM1Q) * rM1Q;
                p1 = (1 - rM1Q) * (1 - rM1Q) + rM1Q * rM1Q;
                p2 = (1 - rM1Q) * rM1Q;
            } else {//aa
                p0 = rM1Q * rM1Q;
                p1 = 2 * (1 - rM1Q) * rM1Q;
                p2 = (1 - rM1Q) * (1 - rM1Q);
            }
            if (missing_type == IMConstant.A_) {
                pM[0] = p0 / (p0 + p1);
                pM[1] = 1;
                pM[2] = 1;
            } else if (missing_type == IMConstant.a_) {
                pM[0] = 0;
                pM[1] = p1 / (p1 + p2);
                pM[2] = 1;
            } else {
                pM[0] = p0;
                pM[1] = p0 + p1;
                pM[2] = 1;
            }
        } else {//both available
            int markerscore = Integer.parseInt(fmk[0]) * 10 + Integer.parseInt(fmk[1]);
            double theta = rM1Q / rM1M2;
            double eta = (rM1M2 * rM1M2) / ((1 - rM1M2) * (1 - rM1M2) + rM1M2 * rM1M2);
            switch (markerscore) {
                case 22:
                    p0 = 1;
                    p1 = 0;
                    p2 = 0;
                    break;
                case 21:
                    p0 = 1 - theta;
                    p1 = theta;
                    p2 = 0;
                    break;
                case 20:
                    p0 = (1 - theta) * (1 - theta);
                    p1 = 2 * theta * (1 - theta);
                    p2 = theta * theta;
                    break;
                case 12:
                    p0 = theta;
                    p1 = 1 - theta;
                    p2 = 0;
                    break;
                case 11:
                    p0 = eta * theta * (1 - theta);
                    p1 = 1 - 2 * eta * theta * (1 - theta);
                    p2 = eta * theta * (1 - theta);
                    break;
                case 10:
                    p0 = 0;
                    p1 = 1 - theta;
                    p2 = theta;
                    break;
                case 2:
                    p0 = theta * theta;
                    p1 = 2 * theta * (1 - theta);
                    p2 = (1 - theta) * (1 - theta);
                    break;
                case 1:
                    p0 = 0;
                    p1 = theta;
                    p2 = 1 - theta;
                    break;
                case 0:
                    p0 = 0;
                    p1 = 0;
                    p2 = 1;
            }
            if (missing_type == IMConstant.A_) {
                pM[0] = p0 / (p0 + p1);
                pM[1] = 1;
                pM[2] = 1;
            } else if (missing_type == IMConstant.a_) {
                pM[0] = 0;
                pM[1] = p1 / (p1 + p2);
                pM[2] = 1;
            } else {
                pM[0] = p0;
                pM[1] = p0 + p1;
                pM[2] = 1;
            }
        }
        double r = rnd.nextDouble();
        int i = 0;
        while (i < pM.length ) {
            if (r < pM[i]) {
                mk = M[i];
                break;
            }
            i++;
        }
        return mk;
    }
}
