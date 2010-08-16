package im;

import java.util.ArrayList;
import java.util.HashMap;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import im.population.IMPopulation;
import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class IntervalPriorProbability {

    int chromosome;
    int interval;
    int walks;
    String populationType;
    int populationIndicator;
    double step;
    double length;
    int map_function;
    HashMap<Integer, Integer> genoPPMap;
    int[] countFlankingMarkerType;
    ArrayList QTLtype;
    Pattern[] population;
    double[][][] interval_pp;   // [NumGenotype][Interval/step][NumQTL]
    boolean isLastInterval;
    
    public IntervalPriorProbability(IMPopulation imp, double s, double l, int chr, int inter, boolean isLast) {
        chromosome = chr;
        interval = inter;
        length = l;
        populationType = imp.getCrossParameter();
        step = s;
        map_function = imp.MapFunction();
        walks = CalculateWalks(l, s);
        isLastInterval = isLast;
        if (isLastInterval) {
        	walks++;
        }
        initial();
    }
    
    public IntervalPriorProbability(IMPopulation imp, double s, double l, int chr, int inter) {
        chromosome = chr;
        interval = inter;
        length = l;
        populationType = imp.getCrossParameter();
        step = s;
        map_function = imp.MapFunction();
        walks = CalculateWalks(l, s);
        initial();
    }

    public int getWalks() {
        return walks;
    }

    private int CalculateWalks(double l, double s) {
        assert l >= 0 && s >= 0;
        if (l < publicAccess.PublicData.epsilon) {
            return 1;
        }
        int len = ((int) Math.ceil(l / s)) > 0 ? ((int) Math.round(l / s)) : 1;

        return len;
    }

    /**
     * See "Statistical Methods for Mapping Quantitative Trait Loci", Zhao-Bang Zeng, August 21,
     * 2000, Section 6.1.
     */
    private void initial() {

        population = new Pattern[5];
        population[0] = Pattern.compile("b1", Pattern.CASE_INSENSITIVE);
        population[1] = Pattern.compile("b2", Pattern.CASE_INSENSITIVE);
        population[2] = Pattern.compile("f2", Pattern.CASE_INSENSITIVE);
        population[3] = Pattern.compile("if2", Pattern.CASE_INSENSITIVE);
        population[4] = Pattern.compile("dh", Pattern.CASE_INSENSITIVE);
        Matcher match;
        for (int i = 0; i < population.length; i++) {
            match = population[i].matcher(populationType);
            if (match.matches()) {
                populationIndicator = IMToolKit.populationindicator[i];
                break;
            }
        }

        genoPPMap = new HashMap();
        QTLtype = new ArrayList();
        if (populationIndicator == IMToolKit.B1) {
            QTLtype.add(new String("QQ"));
            QTLtype.add(new String("Qq"));
            countFlankingMarkerType = new int[4];
            interval_pp = new double[4][][];
            for (int i = 0; i < interval_pp.length; i++) {
                interval_pp[i] = new double[walks][];
                for (int j = 0; j < interval_pp[i].length; j++) {
                    interval_pp[i][j] = new double[2];
                }
            }
            genoPPMap.put(new Integer(22), new Integer(0));
            genoPPMap.put(new Integer(21), new Integer(1));
            genoPPMap.put(new Integer(12), new Integer(2));
            genoPPMap.put(new Integer(11), new Integer(3));
        } else if (populationIndicator == IMToolKit.B2) {
            QTLtype.add(new String("qq"));
            QTLtype.add(new String("Qq"));
            countFlankingMarkerType = new int[4];
            interval_pp = new double[4][][];
            for (int i = 0; i < interval_pp.length; i++) {
                interval_pp[i] = new double[walks][];
                for (int j = 0; j < interval_pp[i].length; j++) {
                    interval_pp[i][j] = new double[2];
                }
            }
            genoPPMap.put(new Integer(0), new Integer(0));
            genoPPMap.put(new Integer(1), new Integer(1));
            genoPPMap.put(new Integer(10), new Integer(2));
            genoPPMap.put(new Integer(11), new Integer(3));
        } else if (populationIndicator == IMToolKit.F2 || populationIndicator == IMToolKit.IF2) {
            QTLtype.add(new String("QQ"));
            QTLtype.add(new String("Qq"));
            QTLtype.add(new String("qq"));
            countFlankingMarkerType = new int[9];
            interval_pp = new double[9][][];
            for (int i = 0; i < interval_pp.length; i++) {
                interval_pp[i] = new double[walks][];
                for (int j = 0; j < interval_pp[i].length; j++) {
                    interval_pp[i][j] = new double[3];
                }
            }
            genoPPMap.put(new Integer(22), new Integer(0));
            genoPPMap.put(new Integer(21), new Integer(1));
            genoPPMap.put(new Integer(20), new Integer(2));
            genoPPMap.put(new Integer(12), new Integer(3));
            genoPPMap.put(new Integer(11), new Integer(4));
            genoPPMap.put(new Integer(10), new Integer(5));
            genoPPMap.put(new Integer(2), new Integer(6));
            genoPPMap.put(new Integer(1), new Integer(7));
            genoPPMap.put(new Integer(0), new Integer(8));
        } else if (populationIndicator == IMToolKit.DH) {
            countFlankingMarkerType = new int[4];
            QTLtype.add(new String("QQ"));
            QTLtype.add(new String("qq"));
            interval_pp = new double[4][][];
            for (int i = 0; i < interval_pp.length; i++) {
                interval_pp[i] = new double[walks][];
                for (int j = 0; j < interval_pp[i].length; j++) {
                    interval_pp[i][j] = new double[2];
                }
            }
            genoPPMap.put(new Integer(22), new Integer(0));
            genoPPMap.put(new Integer(20), new Integer(1));
            genoPPMap.put(new Integer(02), new Integer(2));
            genoPPMap.put(new Integer(0), new Integer(3));
        }
    }

    public int getMarkerScoreIndex(int markerScore) {
        return ((Integer) genoPPMap.get(new Integer(markerScore))).intValue();
    }

    public int getQTLScoreIndex(String qtl) {
        int idx = 0;
        for (int i = 0; i < QTLtype.size(); i++) {
            if (qtl.compareTo((String) QTLtype.get(i)) == 0) {
                idx = i;
                break;
            }
        }
        return idx;
    }
    
    public boolean ScanInterval(int markerScore) {
        boolean flag = false;
        int idx = ((Integer) genoPPMap.get(new Integer(markerScore))).intValue();
        if (countFlankingMarkerType[idx] == 0) {
            if (populationIndicator == IMToolKit.B1) {
                B1PriorProbability(markerScore, idx);
            } else if (populationIndicator == IMToolKit.B2) {
                B2PriorProbability(markerScore, idx);
            } else if (populationIndicator == IMToolKit.F2 || populationIndicator == IMToolKit.IF2) {
                F2PriorProbability(markerScore, idx);
            } else if (populationIndicator == IMToolKit.DH) {
                DHPriorProbability(markerScore, idx);
            }
            countFlankingMarkerType[idx] = 1;
        }
        for (int i = 0; i < countFlankingMarkerType.length; i++) {
            if (countFlankingMarkerType[i] == 0) {
                flag = true;
                break;
            }
        }
        return flag;
    }

    public void B1PriorProbability(int markerScore, int idx) {
        //          QQ    Qq
        //===================
        //M1M1M2M2
        //M1M1M2m2
        //M1m1M2M2
        //M1m1M2m2

        /*
         * This function writes its computational result in interval_pp[idx].
         */

        double[][] pp = new double[walks][2];
        double rM1M2 = IMToolKit.Felsenstein(length, map_function);
        double s = 0;
        for (int i = 0; i < walks; i++) {
            double rM1Q = IMToolKit.Felsenstein(s, map_function);
            double rQM2 = IMToolKit.Felsenstein(length - s, map_function);
            switch (markerScore) {
                case 22:
                    pp[i][0] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
                    pp[i][1] = rM1Q * rQM2 / (1 - rM1M2);
                    break;
                case 21:
                    pp[i][0] = (1 - rM1Q) * rQM2 / rM1M2;
                    pp[i][1] = rM1Q * (1 - rQM2) / rM1M2;
                    break;
                case 12:
                    pp[i][0] = rM1Q * (1 - rQM2) / rM1M2;
                    pp[i][1] = (1 - rM1Q) * rQM2 / rM1M2;
                    break;
                case 11:
                    pp[i][0] = rM1Q * rQM2 / (1 - rM1M2);
                    pp[i][1] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
            }
            if ((walks - i)==2) {
            	if (isLastInterval) {
            		s  = length;
            	} else {
            		s += step;
            	}
            } else  {
            	s += step;
            }
        }
        BayesProbability(pp);
        for (int i = 0; i < walks; i++) {
            System.arraycopy(pp[i], 0, interval_pp[idx][i], 0, pp[i].length);
        }
        return;
    }

    public void B2PriorProbability(int markerScore, int idx) {
        //          qq    Qq
        //===================
        //m1m1m2m2
        //m1m1M2m2
        //M1m1m2m2
        //M1m1M2m2

        /*
         * This function writes its computational result in interval_pp[idx].
         */

        double[][] pp = new double[walks][2];
        double rM1M2 = IMToolKit.Felsenstein(length, map_function);
        double s = 0;
        for (int i = 0; i < pp.length; i++) {
            double rM1Q = IMToolKit.Felsenstein(s, map_function);
            double rQM2 = IMToolKit.Felsenstein(length - s, map_function);
            switch (markerScore) {
                case 0:
                    pp[i][0] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
                    pp[i][1] = rM1Q * rQM2 / (1 - rM1M2);
                    break;
                case 1:
                    pp[i][0] = (1 - rM1Q) * rQM2 / rM1M2;
                    pp[i][1] = rM1Q * (1 - rQM2) / rM1M2;
                    break;
                case 10:
                    pp[i][0] = rM1Q * (1 - rQM2) / rM1M2;
                    pp[i][1] = (1 - rM1Q) * rQM2 / rM1M2;
                    break;
                case 11:
                    pp[i][0] = rM1Q * rQM2 / (1 - rM1M2);
                    pp[i][1] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
            }
            if ((walks - i)==2) {
            	if (isLastInterval) {
            		s  = length;
            	} else {
            		s += step;
            	}
            } else  {
            	s += step;
            }
        }
        BayesProbability(pp);
        for (int i = 0; i < walks; i++) {
            System.arraycopy(pp[i], 0, interval_pp[idx][i], 0, pp[i].length);
        }
        return;
    }

    public void F2PriorProbability(int markerScore, int idx) {//it can be applied to IF2
        //           QQ  Qq qq
        //====================
        //M1M1M2M2
        //M1M1M2m2
        //M1M1m2m2
        //M1m1M2M2
        //M1m1M2m2
        //M1m1m2m2
        //m1m1M2M2
        //m1m1M2m2
        //m1m1m2m2

        /*
         * This function writes its computational result in interval_pp[idx].
         */

        double[][] pp = new double[walks][3];
        double rM1M2 = IMToolKit.Felsenstein(length, map_function);
        double s = 0;
        for (int i = 0; i < pp.length; i++) {
            double rM1Q = IMToolKit.Felsenstein(s, map_function);
            double theta = rM1Q / rM1M2;
            double eta = (rM1M2 * rM1M2) / ((1 - rM1M2) * (1 - rM1M2) + rM1M2 * rM1M2);
            switch (markerScore) {
                case 22:
                    pp[i][0] = 1;
                    pp[i][1] = 0;
                    pp[i][2] = 0;
                    break;
                case 21:
                    pp[i][0] = 1 - theta;
                    pp[i][1] = theta;
                    pp[i][2] = 0;
                    break;
                case 20:
                    pp[i][0] = (1 - theta) * (1 - theta);
                    pp[i][1] = 2 * theta * (1 - theta);
                    pp[i][2] = theta * theta;
                    break;
                case 12:
                    pp[i][0] = theta;
                    pp[i][1] = 1 - theta;
                    pp[i][2] = 0;
                    break;
                case 11:
                    pp[i][0] = eta * theta * (1 - theta);
                    pp[i][1] = 1 - 2 * eta * theta * (1 - theta);
                    pp[i][2] = eta * theta * (1 - theta);
                    break;
                case 10:
                    pp[i][0] = 0;
                    pp[i][1] = 1 - theta;
                    pp[i][2] = theta;
                    break;
                case 2:
                    pp[i][0] = theta * theta;
                    pp[i][1] = 2 * theta * (1 - theta);
                    pp[i][2] = (1 - theta) * (1 - theta);
                    break;
                case 1:
                    pp[i][0] = 0;
                    pp[i][1] = theta;
                    pp[i][2] = 1 - theta;
                    break;
                case 0:
                    pp[i][0] = 0;
                    pp[i][1] = 0;
                    pp[i][2] = 1;
            }
            if ((walks - i)==2) {
            	if (isLastInterval) {
            		s  = length;
            	} else {
            		s += step;
            	}
            } else  {
            	s += step;
            }
        }
        BayesProbability(pp);
        for (int i = 0; i < walks; i++) {
            System.arraycopy(pp[i], 0, interval_pp[idx][i], 0, pp[i].length);
        }
        return;
    }

    public void DHPriorProbability(int markerScore, int idx) {
        //          QQ    qq
        //===================
        //M1M1M2M2
        //M1M1m2m2
        //m1m1M2M2
        //m1m1m2m2

        /*
         * This function writes its computational result in interval_pp[idx].
         */

        double[][] pp = new double[walks][2];
        double rM1M2 = IMToolKit.Felsenstein(length, map_function);
        double s = 0;
        for (int i = 0; i < pp.length; i++) {
            double rM1Q = IMToolKit.Felsenstein(s, map_function);
            double rQM2 = IMToolKit.Felsenstein(length - s, map_function);
            switch (markerScore) {
                case 22:
                    pp[i][0] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
                    pp[i][1] = rM1Q * rQM2 / (1 - rM1M2);
                    break;
                case 21:
                    pp[i][0] = (1 - rM1Q) * rQM2 / rM1M2;
                    pp[i][1] = rM1Q * (1 - rQM2) / rM1M2;
                    break;
                case 12:
                    pp[i][0] = rM1Q * (1 - rQM2) / rM1M2;
                    pp[i][1] = (1 - rM1Q) * rQM2 / rM1M2;
                    break;
                case 11:
                    pp[i][0] = rM1Q * rQM2 / (1 - rM1M2);
                    pp[i][1] = (1 - rM1Q) * (1 - rQM2) / (1 - rM1M2);
            }
            if ((walks - i)==2) {
            	if (isLastInterval) {
            		s  = length;
            	} else {
            		s += step;
            	}
            } else  {
            	s += step;
            }
        }
        BayesProbability(pp);
        for (int i = 0; i < walks; i++) {
            System.arraycopy(pp[i], 0, interval_pp[idx][i], 0, pp[i].length);
        }
        return;
    }

    public int NumQTLtypes() {
        return QTLtype.size();
    }

    public ArrayList QTLGenoType() {
        return QTLtype;
    }

    public double PriorProbabilityAt(int FlankingMarker, int interval, int qtlType) {
        return interval_pp[FlankingMarker][interval][qtlType];
    }

    public double MidPriorProbability(int FlankingMarker, int qtlType) {
        int i = (new Double(walks / 2)).intValue();
        return interval_pp[FlankingMarker][i][qtlType];
    }

    private void BayesProbability(double[][] pp) {
        for (int i = 0; i < pp.length; i++) {
            double sum = 0;
            for (int j = 0; j < pp[i].length; j++) {
                sum += pp[i][j];
            }
            for (int j = 0; j < pp[i].length; j++) {
                pp[i][j] /= sum;
            }
        }
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < interval_pp.length; i++) {
            for (int j = 0; j < interval_pp[i].length; j++) {
                for (int k = 0; k < interval_pp[i][j].length; k++) {
                    sb.append(interval_pp[i][j][k]);
                    sb.append(" ");
                }
                sb.append("\n");
            }
        }
        return sb.toString();
    }
}
