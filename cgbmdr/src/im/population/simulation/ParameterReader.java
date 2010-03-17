package im.population.simulation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class ParameterReader {

    String fileConfig;
    String fileParameter;

    int rep = 1;                        //1-0
    int permutation = 100;              //1-1
    String pt = new String("F2");       //1-2
    int ps = 1000;                      //1-3
    int[] pheNum = {0};                 //1-4
    double mu = 10;                     //1-5
    double[] E = {0};                   //1-6
    double sd = 1;                      //1-7
    long seed = 0;                      //1-8
    int mf = 1;                         //1-9
    double step = 0.01;                 //1-10
    double windowSize = 0.05;           //1-11
    int interval = 5;                   //1-12
    boolean ssci = true;                //1-13
    double mu_selective = 0;            //1-14
    double threshold_selective = 0.0;   //1-15
    int search_start = 1;               //1-16
    int search_end = 1;                 //1-17
    boolean switch2permutation = true;  //1-18
    double threshold = 0.6;             //1-19

    double distance[][] = {             //2-0
            {0.0, 0.02, 0.04, 0.06, 0.08,
                0.0999, 0.12, 0.14, 0.16, 0.18,
                0.1998, 0.22, 0.24, 0.26, 0.28,
                0.2997, 0.32, 0.34, 0.35, 0.36, 0.38,
                0.3996, 0.42, 0.44, 0.46, 0.48,
                0.4995, 0.52, 0.54, 0.56, 0.58,
                0.5994, 0.62, 0.64, 0.66, 0.68,
                0.6993, 0.72, 0.74, 0.76, 0.78,
                0.7992, 0.82, 0.84, 0.86, 0.88,
                0.8991, 0.94, 0.96, 0.98, 0.99},
            {0.0, 0.02, 0.04, 0.06, 0.08,
                0.999, 0.12, 0.14, 0.16, 0.18,
                0.1998, 0.22, 0.24, 0.26, 0.28,
                0.2997, 0.32, 0.34, 0.35, 0.36, 0.38,
                0.3996, 0.42, 0.44, 0.46, 0.48,
                0.4995, 0.52, 0.54, 0.56, 0.58,
                0.5994, 0.62, 0.64, 0.66, 0.68,
                0.6993, 0.72, 0.74, 0.76, 0.78,
                0.7992, 0.82, 0.84, 0.86, 0.88,
                0.8991, 0.94, 0.96, 0.98, 0.99}
        };

    ArrayList QTL;                      //2-1

    public ParameterReader() {
    }

    public void ConfigPopulation(String file) {
        fileConfig = file;
        BufferedReader Reader = null;
        if (fileConfig != null) {
            try {
                Reader = new BufferedReader(new FileReader(fileConfig));
            } catch (IOException E) {
                E.printStackTrace(System.err);
            }
        }
        ArrayList<String> param = new ArrayList();
        String Line;
        try {
            while ((Reader != null) && (Line = Reader.readLine()) != null) {
                if (Line.length() == 0 || Line.startsWith("#")) {
                    continue;
                }
                param.add(Line);
            }
        } catch (IOException E) {
            E.printStackTrace(System.err);
        }



// replication
        if (param.size() > 0) {
            rep = Integer.parseInt(param.get(0));
        }

// permutation
        if (param.size() > 1) {
            permutation = Integer.parseInt(param.get(1));
        }

// population type
        if (param.size() > 2) {
            pt = param.get(2);
        }

// population size
        if (param.size() > 3) {
            ps = Integer.parseInt(param.get(3));
        }

// phenotype number
        if (param.size() > 4) {
            pheNum[0] = Integer.parseInt(param.get(4));
        }

// population mu
        if (param.size() > 5) {
            mu = Double.parseDouble(param.get(5));
        }

// env mu
        if (param.size() > 6) {
            String[] e = param.get(6).split("[,\\s]++");
            E = new double[e.length];
            for (int i = 0; i < e.length; i++) {
                E[i] = Double.parseDouble(e[i]);
            }
        }

// residual
        if (param.size() > 7) {
            sd = Double.parseDouble(param.get(7));
        }


// seed
        if (param.size() > 8) {
            seed = Long.parseLong(param.get(8));
        }

// mapping function
        if (param.size() > 9) {
            mf = Integer.parseInt(param.get(9));
        }

// step
        if (param.size() > 10) {
            step = Double.parseDouble(param.get(10));
        }

// window size
        if (param.size() > 11 ) {
            windowSize = Double.parseDouble(param.get(11));
        }

// interval
        if (param.size() > 12) {
            interval = Integer.parseInt(param.get(12));
        }

// search the same chromosome?
        if (param.size() > 13) {
            ssci = Boolean.parseBoolean(param.get(13));
        }

// selective mu
        if (param.size() > 14) {
            mu_selective = Double.parseDouble(param.get(14));
        }

// selective Threshold
        if (param.size() > 15) {
            threshold_selective = Double.parseDouble(param.get(15));
        }

// search start
        if (param.size() > 16) {
            search_start = Integer.parseInt(param.get(16));
        }

// search end
        if (param.size() > 17) {
            search_end = Integer.parseInt(param.get(17));
        }

// switch to permutation
        if (param.size() > 18) {
            switch2permutation = Boolean.parseBoolean(param.get(18));
        }

// threshold
        if (param.size() > 19) {
            threshold = Double.parseDouble(param.get(19));
        }
    }

    public void ConfigGeneticParameter(String file) {
        fileParameter = file;
        BufferedReader Reader2 = null;
        if (fileParameter != null) {
            try {
                Reader2 = new BufferedReader(new FileReader(fileParameter));
            } catch (IOException E) {
                E.printStackTrace(System.err);
            }
        }
        ArrayList<String> param2 = new ArrayList();
        String Line2;

        try {
            while ((Reader2 != null) && (Line2 = Reader2.readLine()) != null) {
                if (Line2.length() == 0 || Line2.startsWith("#")) {
                    continue;
                }
                param2.add(Line2);
            }
        } catch (IOException E) {
            E.printStackTrace(System.err);
        }

//linkage map
        if (param2.size() > 1) {
            distance = new double[Integer.parseInt(param2.get(0))][];
            for (int k = 0; k < Integer.parseInt(param2.get(0)); k++) {
                String[] dis = param2.get(1 + k).split("[,\\s]++");
                distance[k] = new double[dis.length];
                for (int kk = 0; kk < dis.length; kk++) {
                    distance[k][kk] = Double.parseDouble(dis[kk]);
                }
            }
        }

//QTL
        QTL = new ArrayList();
        int[] chr1 = {0};
        int[] loci1 = {3};
        int[] genotype1 = {1};
        double[] effect1 = {0.5};
        int environment1 = 0;
        AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1, effect1, environment1);
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

                int envi = (new Integer(param2.get((pl + 2 + k * 5) + 4))).intValue();
                al = new AbstractLoci(chr1, loci1, genotype1, effect1, envi);
                QTL.add(al);
            }
        }
    }

    public int SimulationReplication() {
        return rep;
    }

    public int PermutationReplication() {
        return permutation;
    }

    public String PopulationType() {
        return pt;
    }

    public int PopulationSize() {
        return ps;
    }

    public int[] PhenotypeNumber() {
        return pheNum;
    }

    public double Mu() {
        return mu;
    }

    public double[] EMu() {
        return E;
    }

    public double StandardDeviation() {
        return sd;
    }

    public long Seed() {
        return seed;
    }

    public int MappingFunction() {
        return mf;
    }

    public double Step() {
        return step;
    }

    public double WindowSize() {
        return windowSize;
    }

    public int Interval() {
        return interval;
    }

    public boolean SearchSameChromosomeInteraction() {
        return ssci;
    }

    public double MuSelective() {
        return mu_selective;
    }

    public double ThresholdSelective() {
        return threshold_selective;
    }

    public int SearchStart() {
        return search_start;
    }

    public int SearchEnd() {
        return search_end;
    }

    public boolean Switch2Permutation() {
        return switch2permutation;
    }

    public double Threhold() {
        return threshold;
    }

    public double[][] Distance() {
        return distance;
    }

    public ArrayList QTLInformation() {
        return QTL;
    }
}
