
import java.io.*;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.heterogeneity.AbstractHeteroMergeSearch;
import mdr.heterogeneity.HeteroLinearMergeSearch;
import im.reader.IMReaderWinQTLCart;
import im.reader.Imputation;
import im.population.IMPopulation;

/**
 *
 * @author Guo-Bo Chen
 */
public class IMHeteroTest {

    public static void main(String[] args) throws IOException {
        IMReaderWinQTLCart wq = new IMReaderWinQTLCart("IF2Rice.mcd");
        try {
            wq.readData();
        } catch (IOException E) {
            E.printStackTrace(System.err);
        }
        int search_start = 1;
        int search_end = 1;
        int seed = 10;
        int interval = 5;
        double step = 0.02;
        int scIdx = 0;
        double window = 0.2;
        boolean ssci = true;
        boolean env = true;
        double[] os = {0, 0, 0, 0};
        int[] pheNum = {0, 1, 2, 3};
        Imputation im = new Imputation(wq);
        im.Impute();
        IMPopulation imp = new IMPopulation(wq);

        imp.BuildScore(pheNum[0], true);

        int fromC = 1;
        int toC = 1;
        int scrIdx[] = {0};
        double offset[] = {0.5};
        int partitionMethod = 0;
        boolean isMooreMDR = false;

        DataFile dr = new DataFile(wq, scrIdx);

        Subdivision sd = new Subdivision(interval, seed, dr);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(fromC, toC, dr.getMarkerNum());
        cg.generateCombination();
        AbstractHeteroMergeSearch as = new HeteroLinearMergeSearch(dr, sd, cg, scrIdx.length, offset, isMooreMDR);
        for (int i = search_start; i <= search_end; i++) {
            as.search(i, scrIdx[0]);
        }
        as.print();
    }
}
