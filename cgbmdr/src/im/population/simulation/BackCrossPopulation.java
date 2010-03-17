package im.population.simulation;

import java.util.ArrayList;
/**
 *
 * @author Guo-Bo Chen
 */
public class BackCrossPopulation extends AbstractPopulation {

    public BackCrossPopulation(int ps, int pheNum, String pt, double[][] d, long s, double mu, double[] env, double sd, ArrayList Q, int mf) {
        super(ps, pheNum, pt, d, s, mu, env, sd, Q, mf);
    }

    public void ProducePopulation() {
        for (int i = 0; i < IndividualNumber(); i++) {
            int[][] haplotype = ProduceHaplotype();
            int PI = 0;
            PI = (populationType.compareTo("B1") == 0) ? 0:1;
            for (int j = 0; j < marker_QTL[i].length; j++) {
                for (int k = 0; k < marker_QTL[i][j].length; k++) {
                    marker_QTL[i][j][k] = parent[PI][j][k] + haplotype[j][k];
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
            System.out.println();
        }
    }
}
