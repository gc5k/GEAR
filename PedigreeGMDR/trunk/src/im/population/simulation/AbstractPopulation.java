package im.population.simulation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;

import im.population.IMPopulation;
import im.IMToolKit;
import regression.LinearRegression;

/**
 *
 * @author Guo-Bo Chen
 */
public abstract class AbstractPopulation extends IMPopulation {

    // The size of the 1st dimension is the number of chromosomes in an individual.  -- by Zhixiang Zhu
    double[][] _distance;
    double[][] _recombination;
    
    int populationSize;
    int phenotypeSize;
    long seed;
    double mean;
    double stand_deviation;
    double[] Env_mean;
    Random rnd;
    int[][][] parent;
    int[][][] marker_QTL;
    ArrayList QTL;
    HashMap<Integer, TreeSet> QTLLoci;  // <chr, loci>

    public AbstractPopulation(int ps, int pheNum, String pt, double[][] d, long s, double mu, double[] env, double sd, ArrayList Q, int mf) {
        super();
        populationSize = ps;
        populationType = pt;
        phenotypeSize = ps * env.length;
        _distance = new double[d.length][];
        for (int i = 0; i < _distance.length; i++) {
            _distance[i] = new double[d[i].length];
            System.arraycopy(d[i], 0, _distance[i], 0, d[i].length);
        }
        seed = s;
        mean = mu;
        stand_deviation = sd;
        QTL = Q;
        rnd = new Random();
        rnd.setSeed(seed);
        mapping_function = mf;
        Environment = env.length;
        phenotype = new double[phenotypeSize][pheNum];
        strPhenotype = new String[phenotypeSize][pheNum];
        score = new double[phenotypeSize][pheNum];
        QTLLoci = new HashMap();
        Env_mean = new double[env.length];
        System.arraycopy(env, 0, Env_mean, 0, env.length);

        initial();
    }

    private void initial() {

        for (int i = 0; i < phenotypeSize; i++) {
            Integer id = new Integer(i);
            IDs.add(id);
            Integer pid = new Integer(i);
            permutationIDs.add(pid);
            int curr_env = i / populationSize;
            IndNestedEnv.put(new Integer(i), (new Integer(curr_env)).toString());
        }

        parent = new int[2][][];
        for (int i = 0; i < parent.length; i++) {
            parent[i] = new int[_distance.length][];
            for (int j = 0; j < parent[i].length; j++) {
                parent[i][j] = new int[_distance[j].length];
                for (int k = 0; k < parent[i][j].length; k++) {
                    if (i == 0) {
                        parent[i][j][k] = 1;
                    } else {
                        parent[i][j][k] = 0;
                    }
                }
            }
        }

        marker_QTL = new int[phenotypeSize][][];
        for (int i = 0; i < phenotypeSize; i++) {
            marker_QTL[i] = new int[_distance.length][];
            for (int j = 0; j < marker_QTL[i].length; j++) {
                marker_QTL[i][j] = new int[_distance[j].length];
            }
        }

        _recombination = new double[_distance.length][];
        for (int i = 0; i < _recombination.length; i++) {
            _recombination[i] = new double[_distance[i].length];
            _recombination[i][0] = 0.5;
            for (int j = 1; j < _recombination[i].length; j++) {
                _recombination[i][j] = IMToolKit.Felsenstein(_distance[i][j] - _distance[i][j - 1], mapping_function);
            }
        }

        // Transform the pairs of <chromosom-id, locus-id> to pairs of <chromosome-id, loci-set>.  -- by Zhixiang Zhu
        for (Iterator it = QTL.iterator(); it.hasNext();) {
            AbstractLoci al = (AbstractLoci) it.next();
            int[] chr = al.getChr();
            int[] loci = al.getLoci();
            for (int i = 0; i < chr.length; i++) {
                Integer c = new Integer(chr[i]);
                TreeSet ts;
                if (QTLLoci.containsKey(c)) {
                    ts = (TreeSet) QTLLoci.get(c);
                } else {
                    ts = new TreeSet();
                    QTLLoci.put(c, ts);
                }
                ts.add(new Integer(loci[i]));
            }
        }
    }

    protected int[][] ProduceHaplotype() {
        int[][] haplotype;
        haplotype = new int[_distance.length][];
        for (int i = 0; i < haplotype.length; i++) {
            haplotype[i] = new int[_distance[i].length];
        }
        int PI = 0;
        for (int j = 0; j < _distance.length; j++) {
            for (int k = 0; k < _distance[j].length; k++) {
                double r = rnd.nextFloat();
                PI = (r < _recombination[j][k]) ? (1 - PI) : PI;
                haplotype[j][k] = parent[PI][j][k];
            }
        }
        return haplotype;
    }

    protected void FilterOutQTL() {
        marker = new int[marker_QTL.length][][];
        for (int i = 0; i < marker_QTL.length; i++) {   // for each instance
            int dim_chr = marker_QTL[i].length;
            marker[i] = new int[dim_chr][];
            for (int j = 0; j < marker_QTL[i].length; j++) {    // for each chromosome of an instance
                int dim_marker = marker_QTL[i][j].length;
                TreeSet ts = null;
                if (QTLLoci.containsKey(new Integer(j))) {      // this chromosome harbours QTLs
                    ts = (TreeSet) QTLLoci.get(new Integer(j));
                    dim_marker -= ts.size();
                }
                marker[i][j] = new int[dim_marker];
                int idx = 0;
                for (int k = 0; k < marker_QTL[i][j].length; k++) { // for each marker on a chromosome
                    if (ts != null && ts.contains(new Integer(k))) {
                        continue;
                    }
                    marker[i][j][idx++] = marker_QTL[i][j][k];
                }
            }
        }

        distance = new double[_distance.length][];
        recombination = new double[_recombination.length][];
        for (int i = 0; i < _distance.length; i++) {
            int dim_marker = _distance[i].length;
            TreeSet ts = null;
            if (QTLLoci.containsKey(new Integer(i))) {  // this chromosome harbours QTLs
                ts = (TreeSet) QTLLoci.get(new Integer(i));
                dim_marker -= ts.size();
            }
            distance[i] = new double[dim_marker];
            recombination[i] = new double[dim_marker];
            int idx = 0;
            for (int j = 0; j < _distance[i].length; j++) {
                if (ts != null && ts.contains(new Integer(j))) {
                    continue;
                }
                distance[i][idx] = _distance[i][j];
                if (idx == 0) {
                    recombination[i][idx] = 0.5;
                } else {
                    recombination[i][idx] = IMToolKit.Felsenstein(distance[i][idx] - distance[i][idx - 1], mapping_function);
                }
                idx++;
            }
        }
    }

    public abstract void ProducePopulation();

    public abstract void print();

    /**
     *
     * @param sI    score index
     * @param MU
     * @param T
     */
    public void ProducePhenotype(int sI, double MU, double T) {
        double mu = 0;
        int i = 0;
        while (i < phenotypeSize) {
            int curr_env = (new Double(Math.floor(i / populationSize))).intValue();
            phenotype[i][sI] = mean + Env_mean[curr_env] + rnd.nextGaussian() * stand_deviation;
            strPhenotype[i][sI] = new Double (phenotype[i][sI]).toString();
            int[] chr;
            int[] loci;
            int[] genotype;
            double[] effect;
            int nested_env;
            for (Iterator it = QTL.iterator(); it.hasNext();) {
                AbstractLoci al = (AbstractLoci) it.next();
                chr = al.getChr();
                loci = al.getLoci();
                genotype = al.getGenotype();
                effect = al.getEffect();
                nested_env = al.getEnvironment();
                int genoscore = 0;
                if (curr_env != nested_env) {
                    continue;
                }
                for (int j = 0; j < chr.length; j++) {
                    genoscore *= 10;
                    genoscore += marker_QTL[i][chr[j]][loci[j]];
                }
                for (int j = 0; j < genotype.length; j++) {
                    if (genoscore == genotype[j]) {
                        phenotype[i][sI] += effect[j];
                        break;
                    }
                }
            }
            if (Math.abs(phenotype[i][sI] - MU) > T) {
                mu += phenotype[i][sI];
                i++;
            }
        }
        BasicTraitStatistics(sI);
    }

    protected void AutomatedName() {
        markerName = new String[marker[0].length][];
        for(int i = 0; i < marker[0].length; i++) {
            markerName[i] = new String[marker[0][i].length];
            for(int j = 0; j < marker[0][i].length; j++) {
                markerName[i][j] = new String("Chr_" + Integer.toString(i) + "_Marker_" + Integer.toString(j));
            }
        }
    }

    public void printPhenotype() {
        for (int i = 0; i < phenotypeSize; i++) {
            for (int j = 0; j < phenotype[i].length; j++) {
                System.out.print(phenotype[i][j] + "\t");
            }
            System.out.println();
        }
    }
}
