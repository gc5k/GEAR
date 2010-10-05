package mdrAlternative;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Partition;
import algorithm.Subdivision;

import mdr.Combination;
import mdr.moore.LinearMergeSearch;
import mdr.moore.AbstractMergeSearch;

import mdrAlternative.AbstractCrossValidation;
import mdrAlternative.AbstractSearch;
import mdrAlternative.CV;
import mdrAlternative.CVHolder;
import mdrAlternative.DynamicSearch;
import mdrAlternative.GenotypeSearch;
import mdrAlternative.LinearSearch;
import mdrAlternative.TreeSearch;

import publicAccess.PublicData;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class main {

    public static void main(String[] args) throws IOException {
        String mFile = "sample1.txt";
        String pFile = "phenotype1.txt";
        int fromC = 1;
        int toC = 3;
        int scrIdx[] = {0};
        int interval = 5;
        long seed = 56;
        int partitionMethod = PublicData.PairedPartition;
        boolean isMooreMDR = true;
        int searchMethod = 0;
        DataFile dr = new DataFile(mFile, pFile, scrIdx);
        Subdivision sd = new Subdivision(interval, seed, dr);
        sd.RandomPartition();

        Partition partition = new Partition(interval, seed, dr, scrIdx.length, dr.getOffset());
        partition.partition(0, partitionMethod);
        CombinationGenerator cg = new CombinationGenerator(fromC, toC, dr.getMarkerNum());
        cg.generateCombination();
        AbstractSearch as;
        if (searchMethod == PublicData.LinearSearch) {
            as = new LinearSearch(dr, cg, scrIdx);
        } else if (searchMethod == PublicData.TreeSearch) {
            as = new TreeSearch(dr, cg, scrIdx);
        } else if (searchMethod == PublicData.DynamicSearch) {
            as = new DynamicSearch(dr, cg, scrIdx);
        } else {
            as = new GenotypeSearch(dr, cg, scrIdx);
        }
        ArrayList cvResultList = new ArrayList();
        for (int i = fromC; i <= toC; i++) {
            as.search(i);
            CVHolder cvHolder = new CVHolder(i, scrIdx);
            cvResultList.add(cvHolder);
        }
        
        HashMap modelMap = (HashMap) as.getModelMap();
        Set orderedSet = new TreeSet(modelMap.keySet());
        int idx = 0;
        for (Iterator e = orderedSet.iterator(); e.hasNext();) {
            CVHolder cvHolder = (CVHolder) cvResultList.get(idx++);
            Integer order = (Integer) e.next();
            HashMap models = (HashMap) modelMap.get(order);
            Set keys = new TreeSet(models.keySet());
            for (Iterator e1 = keys.iterator(); e1.hasNext();) {
                String key = (String) e1.next();
                int[] c = ToolKit.StringToIntArray(key);
                Combination model = (Combination) models.get(key);
                CV cv = new CV(dr, sd, c, scrIdx, dr.getOffset(), model, true);
//                AbstractCrossValidation cv = new AbstractCrossValidation(dr, partition, c, scrIdx, offset, model, true);
                cv.kFolderSubdivision();
                cv.kFolder();
                cv.calculate();
                cvHolder.add(cv);
            }
        }
       
        for (Iterator e = cvResultList.iterator(); e.hasNext();) {
            CVHolder cvHolder = (CVHolder) e.next();
            cvHolder.printTopModels();
        }
    }

}
