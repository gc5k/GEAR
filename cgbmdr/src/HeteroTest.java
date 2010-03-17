
import java.io.*;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.heterogeneity.AbstractHeteroMergeSearch;
import mdr.heterogeneity.HeteroLinearMergeSearch;

/**
 *
 * @author Guo-Bo Chen
 */
public class HeteroTest {

    public static void main(String[] args) throws IOException {
        String mFile = "ped0.txt";
        String pFile = "phe0.txt";
        int fromC = 1;
        int toC = 1;
        int scrIdx[] = {0};
        int interval = 10;
        long seed = 7;
        int partitionMethod = 0;
        boolean isMooreMDR = false;

//        pFile = null;
        DataFile dr = new DataFile(mFile, pFile, scrIdx);
        Subdivision sd = new Subdivision(interval, seed, dr);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(fromC, toC, dr.getMarkerNum());
        cg.generateCombination();
        AbstractHeteroMergeSearch as = new HeteroLinearMergeSearch(dr, sd, cg, scrIdx.length, dr.getOffset(), isMooreMDR);
        for (int i = fromC; i <= toC; i++) {
            as.search(i, scrIdx[0]);
        }
        as.print();
    }
}
