package im.population.simulation;

import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen
 */
public class F2Population extends AbstractPopulation {

    public F2Population(int ps, int pheNum, String pt, double[][] d, long s, double mu, double[] env, double sd, ArrayList Q, int mf) {
        super(ps, pheNum, pt, d, s, mu, env, sd, Q, mf);
    }

    public void ProducePopulation() {
        for (int h = 0; h < Environment; h++) {
            for (int i = 0; i < populationSize; i++) {
                for (int f = 0; f < 2; f++) {
                    int[][] haplotype = ProduceHaplotype();
                    for (int j = 0; j < marker_QTL[h * populationSize + i].length; j++) {
                        for (int k = 0; k < marker_QTL[h * populationSize + i][j].length; k++) {
                            marker_QTL[h * populationSize + i][j][k] += haplotype[j][k];
                        }
                    }
                }
            }
        }
        FilterOutQTL();
        AutomatedName();
    }

    public void print() {
        for (int i = 0; i < marker.length; i++) {
            for (int j = 0; j < marker[i].length; j++) {
                for (int k = 0; k < marker[i][j].length; k++) {
                    System.out.print(marker[i][j][k] + " ");
                }
            }
            for (int j = 0; j < phenotype[i].length; j++) {
                System.out.print(phenotype[i][j] + " ");
            }
            for (int j = 0; j < score[i].length; j++) {
                System.out.print(score[i][j] + " ");
            }
            System.out.println();
        }
    }
}
