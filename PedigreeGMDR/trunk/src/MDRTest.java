
import java.io.*;
import java.util.ArrayList;
import java.util.StringTokenizer;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Partition;
import algorithm.Subdivision;

import mdr.moore.AbstractMergeSearch;
import mdr.moore.LinearMergeSearch;

import publicAccess.PublicData;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class MDRTest {

    public static void main(String[] args) throws IOException {
        File Config = new File(args[0]);
        ArrayList<String> Parameters = new ArrayList();
        BufferedReader Reader = new BufferedReader(new FileReader(Config));
        String Line;
        while ((Line = Reader.readLine()) != null) {
            if (Line.length() == 0 || Line.startsWith("#")) {
                continue;
            }
            Parameters.add(Line);
        }
        int ParamNum = Parameters.size();
        if (ParamNum != 9) {
            throw new IOException("Problems in reading paremeters.");
        }

        String mFile = Parameters.get(0);
        String pFile = Parameters.get(1);
        if (pFile.compareTo(PublicData.IgnorFile) == 0) {
            pFile = null;
        }
        int fromC = Integer.parseInt(Parameters.get(2));
        int toC = Integer.parseInt(Parameters.get(3));
        StringTokenizer Tokenizer = new StringTokenizer((String) Parameters.get(4),
                ",\t ");
        String[] scr_idx = new String[Tokenizer.countTokens()];
        int ii = 0;
        while (Tokenizer.hasMoreTokens()) {
            scr_idx[ii++] = Tokenizer.nextToken();
        }
        int scrIdx[] = ToolKit.StringArrayTOIntArray(scr_idx);
        int interval = Integer.parseInt(Parameters.get(5));
        long seed = Integer.parseInt(Parameters.get(6));
        int partitionMethod = Integer.parseInt(Parameters.get(7));
        int searchMethod = Integer.parseInt(Parameters.get(8));
        boolean isMooreMDR = true;
        DataFile dr = new DataFile(mFile, pFile, scrIdx);
        Subdivision sd = new Subdivision(interval, seed, dr);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(fromC, toC, dr.getMarkerNum());
        cg.generateCombination();
        AbstractMergeSearch as = new LinearMergeSearch(dr, sd, cg, scrIdx.length, dr.getOffset(), isMooreMDR);
        for (int i = fromC; i <= toC; i++) {
            as.search(i);
            as.summarise();
        }
    }
}
