
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.math.linear.RealMatrix;

import algorithm.*;
import im.population.simulation.AbstractLoci;
import im.population.IMPopulation;
import im.IMBMatrix;
import im.IntervalPriorProbability;
import im.GenomeScan;
import im.population.simulation.*;
import publicAccess.PublicData;
import publicAccess.ToolKit;
import regression.LinearRegression;

/**
 *
 * @author Guo-Bo Chen
 */
public class RegPopulation {

    public static void main(String[] args) throws IOException {

        File config = null;
        if (args.length > 0) {
            config = new File(args[0]);
        }
        BufferedReader Reader = null;
        if (config != null) {
            Reader = new BufferedReader(new FileReader(config));
        }
        ArrayList<String> param = new ArrayList();
        String Line;

        while ((Reader != null) && (Line = Reader.readLine()) != null) {
            if (Line.length() == 0 || Line.startsWith("#")) {
                continue;
            }
            param.add(Line);
        }

// population size
        int ps = 500;
        if (param.size() > 0) {
            ps = Integer.parseInt(param.get(0));
        }
// phenotype number
        int[] pheNum = {0};
        double[] os = {0};
        if (param.size() > 1) {
            pheNum[0] = Integer.parseInt(param.get(1));
            os[0] = 0;
        }
// population type
        String pt = new String("F2");
        if (param.size() > 2) {
            pt = param.get(2);
        }
// population mu
        double mu = 10;
        if (param.size() > 3) {
            mu = Double.parseDouble(param.get(3));
        }
// residual
        double sd = 1;
        if (param.size() > 4) {
            sd = Double.parseDouble(param.get(4));
        }
// mapping function
        int mf = 1;
        if (param.size() > 5) {
            mf = Integer.parseInt(param.get(5));
        }
// step
        double step = 0.01;
        if (param.size() > 6) {
            step = Double.parseDouble(param.get(6));
        }
// search the same chromosome?
        boolean ssci = true;
        if (param.size() > 7) {
            ssci = Boolean.parseBoolean(param.get(7));
        }
// interval
        int interval = 5;
        if (param.size() > 8) {
            interval = Integer.parseInt(param.get(8));
        }
// seed
        long seed = 8;
        if (param.size() > 9) {
            seed += Long.parseLong(param.get(9));
        }
// selective mu
        double MU = 0;
        if (param.size() > 10) {
            MU = Double.parseDouble(param.get(10));
        }
// selective Threshold
        double T = 0;
        if (param.size() > 11) {
            T = Double.parseDouble(param.get(11));
        }
// search start
        int search_start = 1;
        if (param.size() > 12) {
            search_start = Integer.parseInt(param.get(12));
        }
// search end
        int search_end = 1;
        if (param.size() > 13) {
            search_end = Integer.parseInt(param.get(13));
        }
// switch to permutation
        boolean switch2permutation = false;
        if (param.size() > 14) {
            switch2permutation = Boolean.parseBoolean(param.get(14));
        }
// replication
        int rep = 1;
        if (param.size() > 15) {
            rep = Integer.parseInt(param.get(15));
        }
// permutation
        int permutation = 0;
        if (param.size() > 16) {
            permutation = Integer.parseInt(param.get(16));
        }

////////////////////////////////////// the second file;//////////////////////////////////////////////////////
        File config2 = null;
        if (args.length > 1) {
            config2 = new File(args[1]);
        }
        BufferedReader Reader2 = null;
        if (config != null) {
            Reader2 = new BufferedReader(new FileReader(config2));
        }
        ArrayList<String> param2 = new ArrayList();
        String Line2;

        while ((Reader2 != null) && (Line2 = Reader2.readLine()) != null) {
            if (Line2.length() == 0 || Line2.startsWith("#")) {
                continue;
            }
            param2.add(Line2);
        }
//linkage map
        double d[][] = {{0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}};

        if (param2.size() > 1) {
            d = new double[Integer.parseInt(param2.get(0))][];
            for (int k = 0; k < Integer.parseInt(param2.get(0)); k++) {
                String[] distance = param2.get(1+k).split("[,\\s]++");
                d[k] = new double[distance.length];
                for (int kk = 0; kk < distance.length; kk++) {
                    d[k][kk] = Double.parseDouble(distance[kk]);
                }
            }
        }
//QTL
        int[] chr1 = {0};
        int[] loci1 = {3};
        int[] genotype1 = {1};
        double[] effect1 = {0.5};
        int environment1 = 0;
        AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1, effect1, environment1);
        ArrayList QTL = new ArrayList();
        QTL.add(al);

        int pl = 0;
        if (param2.size() > 0) {
            pl = Integer.parseInt(param2.get(0));
        }
        if (param2.size() > (pl + 1)) {
            QTL.clear();
            int qtlnumber = Integer.parseInt(param2.get(pl + 1));
            for (int k = 0; k < qtlnumber; k++) {
                String[] chr = param2.get((pl + 2 + k * 5)).split("[,\\s]++");
                chr1 = new int[chr.length];
                for (int kk = 0; kk < chr.length; kk++) {
                    chr1[kk] = Integer.parseInt(chr[kk]);
                }

                String[] loc = param2.get((pl + 2 + k * 5) + 1).split("[,\\s]++");
                loci1 = new int[loc.length];
                for (int kk = 0; kk < loc.length; kk++) {
                    loci1[kk] = Integer.parseInt(loc[kk]);
                }

                String[] fun = param2.get((pl + 2 + k * 5) + 2).split("[,\\s]++");
                genotype1 = new int[fun.length];
                for (int kk = 0; kk < fun.length; kk++) {
                    genotype1[kk] = Integer.parseInt(fun[kk]);
                }

                String[] eff = param2.get((pl + 2 + k * 5) + 3).split("[,\\s]++");
                effect1 = new double[eff.length];
                for (int kk = 0; kk < eff.length; kk++) {
                    effect1[kk] = Double.parseDouble(eff[kk]);
                }
                int envi1 = new Integer((param2.get(pl +2 + k * 5) + 4));
                al = new AbstractLoci(chr1, loci1, genotype1, effect1, envi1);
                QTL.add(al);
            }
        }

        int Env = 2;
        double[] env={0.0, 0.0};
        if ((MU - 0) > PublicData.epsilon) {
            AbstractPopulation ap;
            if (pt.compareTo("F2") == 0) {
                ap = new F2Population(ps, pheNum.length, pt, d, seed + rep, mu, env, sd, QTL, mf);
            } else if (pt.compareTo("DH") == 0) {
                ap = new DHPopulation(ps, pheNum.length, pt, d, seed + rep, mu, env, sd, QTL, mf);
            } else {
                ap = new BackCrossPopulation(ps, pheNum.length, pt, d, seed + rep, mu, env, sd, QTL, mf);
            }
            ap.ProducePopulation();

            for (int i = 0; i < pheNum.length; i++) {
                ap.ProducePhenotype(i, MU, T);
            }
            MU = ap.getMean(0);
        }

        for (int i_rep = 0; i_rep < rep; i_rep++) {
            AbstractPopulation ap;
            if (pt.compareTo("F2") == 0) {
                ap = new F2Population(ps, pheNum.length, pt, d, seed + i_rep, mu, env, sd, QTL, mf);
            } else if (pt.compareTo("DH") == 0) {
                ap = new DHPopulation(ps, pheNum.length, pt, d, seed + i_rep, mu, env, sd, QTL, mf);
            } else {
                ap = new BackCrossPopulation(ps, pheNum.length, pt, d, seed + i_rep, mu, env, sd, QTL, mf);
            }
            ap.ProducePopulation();

            for (int i = 0; i < pheNum.length; i++) {
                ap.ProducePhenotype(i, MU, T);
            }

            GenomeScan gs = new GenomeScan(ap, step);
            gs.CalculateIPP();
            double[][] Y = new double[ap.IndividualNumber()][1];
            double[][] X;

            if (permutation == 0) {
                ArrayList ids = ap.getIDs();
                for (int i = 0; i < ids.size(); i++) {
                    Integer id = (Integer) ids.get(i);
                    Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
                }
                IMBMatrix imb = new IMBMatrix(gs, ap);
                for (int i = search_start; i <= search_end; i++) {
                    imb.setOrder(i);
                    CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumIntevals());
                    cg.generateCombination();
                    List com = cg.get(i);
                    double[][] Coeff={{1,0,-1},{0.5,-0.5,0.5}};

                    for (Iterator e = com.iterator(); e.hasNext(); ) {
                        String s = (String) e.next();
                        for(int w = 0; w < 10; w++) {
//                        X = imb.getPPMatrix(s, Coeff);
                            X = imb.getPPMatrixAtPoint(s, w, Coeff);
                            System.out.print(imb.getString());
                            LinearRegression lm = new LinearRegression(X, Y);
                            lm.MLE();
                            System.out.println(w + " " + lm.getEstimate()+" "+lm.getSSTO()+" "+lm.getSSR()+" "+lm.getSSE() + " " +lm.get_F_Statistic() + " " + lm.getP_F());
                        }
                    }
                }
            } else {
                for (int i_permu = 0; i_permu < permutation; i_permu++) {
                    ap.Swith2Permutation(switch2permutation, seed * (i_rep * 100) + i_permu);
                    ArrayList ids = ap.getIDs();
                    for (int i = 0; i < ids.size(); i++) {
                        Integer id = (Integer) ids.get(i);
                        Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
                    }
                    for (int i = search_start; i <= search_end; i++) {
                        CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumIntevals());
                        cg.generateCombination();
                    }
                }
            }
        }
    }
}
