package mdr.alternative;

import mdr.alternative.CV;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author Guo-Bo Chen
 */
public class CVHolder extends HashMap {

    private int order;
    private int[] scrIdx;
    private double topMeanTA;
    private double topMeanTrA;

    private ArrayList topModels;

    public CVHolder(int or, int[] sI) {
        order = or;
        topModels = new ArrayList();
        scrIdx = new int[sI.length];
        System.arraycopy(sI, 0, scrIdx, 0, sI.length);
        for( int i=0; i<scrIdx.length; i++) {
            topModels.add(new ArrayList());
        }
    }

    public int getOrder() {
        return order;
    }

    public void add(CV cv) {
        String modelKey = cv.getCVKey();
        put(modelKey, cv);
        int i = 0;
        for (Iterator e = cv.iterator(); e.hasNext(); ) {
            CV.CVOneTraitResult cvTrait = (CV.CVOneTraitResult) e.next();
            ArrayList topList = (ArrayList) topModels.get(i);
            int j = 0;
            for (Iterator e1 = cvTrait.iterator(); e1.hasNext();) {
                CV.CVPair cvPair = (CV.CVPair) e1.next();
                double trAccu = cvPair.getTrainingAccuracy();
                if (size() > 1) {
                    CV topCV = (CV) topList.get(j);
                    CV.CVOneTraitResult topResult = (CV.CVOneTraitResult) topCV.get(i);
                    CV.CVPair topCVPair = (CV.CVPair) topResult.get(j);
                    double topTrAccu = topCVPair.getTrainingAccuracy();
                    if (topTrAccu < trAccu) {
                        topList.set(j, cv);
                    }
                } else {
                    topList.add(cv);
                }
                j++;
            }
            i++;
        }
    }
    
    public void printTopModels() {
        int c=0;
        System.out.println("Order of Interaction is " + order);
        for (Iterator e = topModels.iterator(); e.hasNext();) {
            ArrayList topResult = (ArrayList) e.next();
            int count = 0;
            for (Iterator e1 = topResult.iterator(); e1.hasNext();) {
                CV cv = (CV) e1.next();
                if(cv.isMooreMDR()) {
                    cv.testingCalculate(c, count);
                }
                CV.CVOneTraitResult cvResult = (CV.CVOneTraitResult) cv.get(c);
                CV.CVPair cvPair = (CV.CVPair) cvResult.get(count);
                System.out.println(cv.getCVKey() + ", " + cvPair.getTestingAccuracy() + ", " + cvPair.getTrainingAccuracy());
                count++;
                topMeanTA  = topMeanTA + cvPair.getTestingAccuracy();
                topMeanTrA = topMeanTrA + cvPair.getTrainingAccuracy();
            }
            topMeanTA /= topResult.size();
            topMeanTrA /= topResult.size();
            System.out.println("mean of testing accuracy: "+ topMeanTA);
            System.out.println("mean of training accuracy: "+ topMeanTrA);
            c++;
        }
    }
}
