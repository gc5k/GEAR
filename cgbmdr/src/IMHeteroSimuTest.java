
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import algorithm.CombinationGenerator;
import algorithm.Subdivision;
import mdr.data.ExDataFile;
import mdr.SchematicPermutation;

import mdr.heterogeneity.IMHeteroMarkerLinearMergeSearch;

import im.FineMappingCandidate;
import im.population.MarkerInformation;

import im.population.simulation.AbstractLoci;
import im.population.simulation.AbstractPopulation;
import im.population.simulation.BackCrossPopulation;
import im.population.simulation.DHPopulation;
import im.population.simulation.F2Population;
import im.population.simulation.ParameterReader;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class IMHeteroSimuTest {

    public static void main(String[] args) throws IOException {
        ParameterReader pr = null;
        if (args.length > 0) {
            pr = new ParameterReader();
            pr.ConfigPopulation(args[0]);
            pr.ConfigGeneticParameter(args[1]);
        }

        int rep = 1;
        int permutation = 100;
        int ps = 500;
        int[] pheNum = {0};
        double mu = 10;
        double[] E = {0};
        double sd = 1;
        long seed = 12;
        int mf = 1;
        double step = 0.02;
        double windowSize = 0.05;
        int interval = 5;
        boolean ssci = true;
        double MU_sel = 0;
        double T_sel = 0;

        String pt = new String("F2");
        int search_start = 1;
        int search_end = 1;
        boolean isMooreMDR = false;
        boolean switch2permutation = true;
        double threshold = 0.6;

        double d[][] = {
            {0.0,0.0999,0.1998,0.2997,0.3996,0.45,0.4995,0.5994,0.6993,0.7992,0.8991, 0.99},
        };

        ArrayList QTL = new ArrayList();
        int[] chr1 = {0};
        int[] loci1 = {5};
        int[] genotype1 = {1};
        double[] effect1 = {0.5};
        int environment1 = 0;
        AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1, effect1, environment1);
        QTL.add(al);

        if (pr != null) {
            rep = pr.SimulationReplication();
            permutation = pr.PermutationReplication();
            ps = pr.PopulationSize();
            pt = pr.PopulationType();
            pheNum = pr.PhenotypeNumber();
            mu = pr.Mu();
            E = pr.EMu();
            sd = pr.StandardDeviation();
            seed = pr.Seed();
            mf = pr.MappingFunction();
            step = pr.Step();
            windowSize = pr.WindowSize();
            interval = pr.Interval();
            ssci = pr.SearchSameChromosomeInteraction();
            MU_sel = pr.MuSelective();
            T_sel = pr.ThresholdSelective();
            search_start = pr.SearchStart();
            search_end = pr.SearchEnd();
            switch2permutation = true;
            threshold = pr.Threhold();

            d = pr.Distance();
            QTL = pr.QTLInformation();
        }

        int[] scrIdx = {0};

        for (int i_rep = 0; i_rep <= rep; i_rep++) {

            StringBuffer out = new StringBuffer();
            out.append(search_start);
            out.append("loci_result");
            out.append(i_rep);
            out.append(".txt");
            PrintStream Out = new PrintStream(new BufferedOutputStream(new FileOutputStream(out.toString())));
            System.setOut(Out);
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
                ap.ProducePhenotype(i, MU_sel, T_sel);
            }

            ExDataFile exdr = new ExDataFile(ap, scrIdx, 0.1);
            exdr.PickMarker(0.1);
            Subdivision sud = new Subdivision(interval, seed, exdr);
            sud.RandomPartition();

            CombinationGenerator cg = new CombinationGenerator(search_start, search_end, exdr.getMarkerNum());
            cg.generateCombination();
            if (permutation > 0 && i_rep == 0) {
                exdr.switchPermutation(true);
                double[] ts = new double[permutation];
                for (int i = 0; i < permutation; i++) {
                    IMHeteroMarkerLinearMergeSearch as = new IMHeteroMarkerLinearMergeSearch(exdr, sud, cg, scrIdx.length, exdr.getOffset(), isMooreMDR);
                    SchematicPermutation sp = new SchematicPermutation(ap);
                    exdr.setShuffledIDs(sp.getShuffledIDs(i));
                    for (int j = search_start; j <= search_end; j++) {
                        as.search(j, scrIdx[0]);
                        ts[i] = as.getTopStatistic(scrIdx[0], i);
                    }
                }
                Arrays.sort(ts);
                String pout = "permutation.txt";
                PrintStream POut = new PrintStream(new BufferedOutputStream(new FileOutputStream(pout.toString())));
                for (int i = 0; i < ts.length; i++) {
                    POut.println(ts[i]);
                }
                POut.close();
                int thresholdIndex = (int) (ts.length * 0.95);
                threshold = ts[thresholdIndex];
            }

            System.out.println(threshold);
            exdr.switchPermutation(false);
            IMHeteroMarkerLinearMergeSearch as = new IMHeteroMarkerLinearMergeSearch(exdr, sud, cg, scrIdx.length, exdr.getOffset(), isMooreMDR);
            ArrayList sigPoints = null;

            for (int i = search_start; i <= search_end; i++) {
                as.search(i, scrIdx[0]);
                sigPoints = as.getSignificantPoints(i, scrIdx[0], threshold);
                ArrayList<ArrayList> candidateComs = null;
                if (sigPoints != null) {
                    System.out.println("selected intervals");
                    for (Iterator e = sigPoints.iterator(); e.hasNext();) {
                        MarkerInformation mi = exdr.CombinationGeneticInformationSelected((String) e.next());
                        System.out.println(mi);
                    }
                    FineMappingCandidate fmc = new FineMappingCandidate(sigPoints, exdr, windowSize);
                    candidateComs = fmc.candidateCombinations();
                    System.out.println("There are " + candidateComs.size() + " candidates");
                    int c = 0;
                    for (Iterator e = candidateComs.iterator(); e.hasNext();) {
                        ArrayList candidateCom = (ArrayList) e.next();
                        System.out.println("Candidate Interval " + ++c);
                        for (Iterator e1 = candidateCom.iterator(); e1.hasNext();) {
                            String com = (String) e1.next();
                            MarkerInformation mi = exdr.CombinationGeneticInformationOriginal(com);
                            System.out.println(mi);
                        }
                    }
                    ExDataFile cdr = new ExDataFile(ap, scrIdx, -1);
                    Subdivision csud = new Subdivision(interval, seed, cdr);
                    csud.RandomPartition();

                    CombinationGenerator ccg = new CombinationGenerator(search_start, search_end, cdr.getMarkerNum());
                    ccg.generateCombination();
                    IMHeteroMarkerLinearMergeSearch cas = new IMHeteroMarkerLinearMergeSearch(cdr, csud, ccg, scrIdx.length, exdr.getOffset(), isMooreMDR);
                    for (Iterator e = candidateComs.iterator(); e.hasNext();) {
                        ArrayList com = (ArrayList) e.next();
                        cas.search(i, scrIdx[0], com);
                        MarkerInformation mi = ap.CombinationGeneticInformationOriginal(cas.getTopCombination());
                        System.out.println("Best estimate is " + mi + " " + cas.getTopStatistic(i, i));
                    }
                }
            }
            Out.close();
        }
    }
}
