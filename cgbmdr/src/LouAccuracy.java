
import java.io.*;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;
import im.GenomeScan;
import im.Summary;
import im.population.simulation.AbstractPopulation;
import im.population.simulation.AbstractLoci;
import im.population.simulation.BackCrossPopulation;
import im.population.simulation.F2Population;
import im.population.simulation.DHPopulation;

import mdr.heterogeneity.IMHeteroLinearMergeSearch;
import mdr.heterogeneity.LouIMHeteroLinearCompleteMergeSearch;
import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class LouAccuracy {

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

// env mu
        double[] E = {0};
        if (param.size() > 0) {
            String[] e = param.get(0).split("[,\\s]++");
            E = new double[e.length];
            for (int i = 0; i < e.length; i++) {
                E[i] = Double.parseDouble(e[i]);
            }
        }

// population size
        int ps = 1000;
        if (param.size() > 1) {
            ps = Integer.parseInt(param.get(1));
        }

// phenotype number
        int[] pheNum = {0};
        double[] os = {0};
        if (param.size() > 2) {
            pheNum[0] = Integer.parseInt(param.get(2));
            os[0] = 0;
        }

// population type
        String pt = new String("B1");
        if (param.size() > 3) {
            pt = param.get(3);
        }
// population mu
        double mu = 10;
        if (param.size() > 4) {
            mu = Double.parseDouble(param.get(4));
        }
// residual
        double sd = 0;
        if (param.size() > 5) {
            sd = Double.parseDouble(param.get(5));
        }
// mapping function
        int mf = 1;
        if (param.size() > 6) {
            mf = Integer.parseInt(param.get(6));
        }
// step
        double step = 0.01;
        if (param.size() > 7) {
            step = Double.parseDouble(param.get(7));
        }
// search the same chromosome?
        boolean ssci = true;
        if (param.size() > 8) {
            ssci = Boolean.parseBoolean(param.get(8));
        }
// interval
        int interval = 5;
        if (param.size() > 9) {
            interval = Integer.parseInt(param.get(9));
        }
// seed
        long seed = 0;
        if (param.size() > 10) {
            seed = Long.parseLong(param.get(10));
        }
// selective mu
        double MU = 0;
        if (param.size() > 11) {
            MU = Double.parseDouble(param.get(11));
        }
// selective Threshold
        double T = 0.0;
        if (param.size() > 12) {
            T = Double.parseDouble(param.get(12));
        }
// search start
        int search_start = 1;
        if (param.size() > 13) {
            search_start = Integer.parseInt(param.get(13));
        }
// search end
        int search_end = 1;
        if (param.size() > 14) {
            search_end = Integer.parseInt(param.get(14));
        }
// switch to permutation
        boolean switch2permutation = true;
        if (param.size() > 15) {
            switch2permutation = Boolean.parseBoolean(param.get(15));
        }
// replication
        int rep = 1;
        if (param.size() > 16) {
            rep = Integer.parseInt(param.get(16));
        }
// permutation
        int permutation = 0;
        if (param.size() > 17) {
            permutation = Integer.parseInt(param.get(17));
        }
// threshold
        double threshold = 0.6;
        if (param.size() > 18) {
            threshold = Double.parseDouble(param.get(18));
        }

////////////////////////////////////// the second file;
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
        double d[][] = {
            {0, 0.1, 0.2, 0.299, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}
        };

        if (param2.size() > 1) {
            d = new double[Integer.parseInt(param2.get(0))][];
            for (int k = 0; k < Integer.parseInt(param2.get(0)); k++) {
                String[] distance = param2.get(1 + k).split("[,\\s]++");
                d[k] = new double[distance.length];
                for (int kk = 0; kk < distance.length; kk++) {
                    d[k][kk] = Double.parseDouble(distance[kk]);
                }
            }
        }
//QTL
        ArrayList QTL = new ArrayList();
        int[] chr1 = {0};
        int[] loci1 = {3};
        int[] genotype1 = {1};
        double[] effect1 = {2};
        int environment1 = 0;
        AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1, effect1, environment1);
        QTL.add(al);

//        ArrayList QTL = new ArrayList();
//        int[] chr1 = {0,1};
//        int[] loci1 = {3,5};
//        int[] genotype1 = {1,10,12,21};
//        double[] effect1 = {0.5, 0.5, 0.5, 0.5};
//        int environment1 = 0;
//        AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1, effect1, environment1);
//        QTL.add(al);
//
//        int[] chr2 = {3,4};
//        int[] loci2 = {8,5};
//        int[] genotype2 = {0,11,22};
//        double[] effect2 = {1,0.5,1};
//        int environment2 = 0;
//        AbstractLoci al2 = new AbstractLoci(chr2, loci2, genotype2, effect2, environment2);
//        QTL.add(al2);

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

                int envi = (new Integer(param2.get((pl + 2 + k * 5) + 4))).intValue();
                al = new AbstractLoci(chr1, loci1, genotype1, effect1, envi);
                QTL.add(al);
            }
        }

        if ((MU - 0) > PublicData.epsilon) {//when phenotype should be selected, calculate the population mean first.
            AbstractPopulation ap;
            if (pt.compareTo("F2") == 0) {
                ap = new F2Population(ps, pheNum.length, pt, d, seed + rep, mu, E, sd, QTL, mf);
            } else if (pt.compareTo("DH") == 0) {
                ap = new DHPopulation(ps, pheNum.length, pt, d, seed + rep, mu, E, sd, QTL, mf);
            } else {
                ap = new BackCrossPopulation(ps, pheNum.length, pt, d, seed + rep, mu, E, sd, QTL, mf);
            }
            ap.ProducePopulation();

            for (int i = 0; i < pheNum.length; i++) {
                ap.ProducePhenotype(i, MU, T);
            }
            MU = ap.getMean(0);
        }

//        StringBuffer OUT = new StringBuffer();
//        OUT.append("result.txt");
        PrintStream POUT = new PrintStream(new BufferedOutputStream(new FileOutputStream("cim.txt")));
        PrintStream LouAccu = new PrintStream(new BufferedOutputStream(new FileOutputStream("LouAccu.txt")));
        for (int i_rep = 0; i_rep < rep; i_rep++) {
            AbstractPopulation ap;
            if (pt.compareTo("F2") == 0) {
                ap = new F2Population(ps, pheNum.length, pt, d, seed + i_rep, mu, E, sd, QTL, mf);
            } else if (pt.compareTo("DH") == 0) {
                ap = new DHPopulation(ps, pheNum.length, pt, d, seed + i_rep, mu, E, sd, QTL, mf);
            } else {
                ap = new BackCrossPopulation(ps, pheNum.length, pt, d, seed + i_rep, mu, E, sd, QTL, mf);
            }
            ap.ProducePopulation();

            for (int i = 0; i < pheNum.length; i++) {
                ap.ProducePhenotype(i, MU, T);
            }
            ap.CIMFormat(POUT);
            boolean env = ap.getEnvironment() > 1 ? true : false;
            ap.BuildScore(0, env);
            Subdivision subdivision = new Subdivision(interval, seed + i_rep, ap);
            subdivision.RandomPartition();
            GenomeScan gs = new GenomeScan(ap, step); //set up genome-wide scan
            gs.CalculateIPP(); //calculate conditional probability of putative QTL for given loci
            StringBuffer out = new StringBuffer();
            out.append(search_start);
            out.append("loci_result");
            out.append(i_rep);
            out.append(".txt");
            PrintStream Pout = new PrintStream(new BufferedOutputStream(new FileOutputStream(out.toString())));
            if (i_rep == 0) {
                StringBuffer pout = new StringBuffer();
                pout.append("permutation");
                pout.append(i_rep);
                pout.append(".txt");
                PrintStream POut = new PrintStream(new BufferedOutputStream(new FileOutputStream(pout.toString())));
//                System.setOut(POut);

                for (int i_permu = 0; i_permu < permutation; i_permu++) {
                    ap.Swith2Permutation(switch2permutation, seed * (i_rep * 100) + i_permu);
//                    ap.tunePermutation(seed + i_permu);
                    for (int i = search_start; i <= search_end; i++) {
                        CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumMarkers());
                        cg.generateCombination();
                        IMHeteroLinearMergeSearch imh = new IMHeteroLinearMergeSearch(gs, ap, subdivision, cg, pheNum, os, false, ap.getPermutatedIDs());
                        imh.setSearch(ssci, env);
                        imh.search(i, pheNum[0]);
//                        System.out.println(imh.peakStatistic());
                    }
                }
                POut.close();
            }

            for (int i = search_start; i <= search_end; i++) {
                CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumMarkers());
                cg.generateCombination();  // combinations for intervals;
                LouIMHeteroLinearCompleteMergeSearch imh = new LouIMHeteroLinearCompleteMergeSearch(gs, ap, subdivision, cg, pheNum, os, false);
                imh.setSearch(ssci, env);
                imh.search(i, pheNum[0]);
                imh.printWholeAccuracyLou(LouAccu, 0);
            }
            Pout.close();
        }
        LouAccu.close();
    }
}
