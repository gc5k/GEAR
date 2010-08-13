package publicAccess;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import mdr.Cell;
/**
 *
 * @author Guo-Bo Chen
 */
public class ToolKit {

    public static int[] StringToIntArray(String s) {
        String[] unit = s.split(PublicData.seperator);
        int[] value = new int[unit.length];
        for (int i = 0; i < unit.length; i++) {
            value[i] = Integer.parseInt(unit[i]);
        }
        return value;
    }

    public static int[] StringArrayTOIntArray(String[] s) {
        if (s == null) {
            return null;
        }
        int[] d = new int[s.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = Integer.parseInt(s[i]);
        }
        return d;
    }
    
    public static String IntArrayToString(int[] array) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < array.length; i++) {
            Integer I = new Integer(array[i]);
            sb.append(I);
            if (i == (array.length - 1)) {
                continue;
            }
            sb.append(PublicData.seperator);
        }
        return sb.toString();
    }
    
    public static int Acertainment (double posScr, double negScr, double threshold) {
        int status = 0;
        if (negScr == 0) {
            if (posScr == 0) {
                status = -1;
            } else {
                status = 1;
            }
        } else {
            if ((posScr / Math.abs(negScr)) == threshold) {
                status = PublicData.tieValue;
            } else {
                status = (posScr / Math.abs(negScr)) > threshold ? 1 : 0;
            }
        }
        return status;
    }
    
    public static double Mean(double[] vector) {
        double sum = 0;
        if (vector.length == 0) {
            return 0;
        }
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        return sum / (double) vector.length;
    }

    public static double Variance(double[] vector) {
        double sum = 0, sumSquared = 0;
        if (vector.length <= 1) {
            return 0;
        }
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
            sumSquared += (vector[i] * vector[i]);
        }
        double result = (sumSquared - (sum * sum / (double) vector.length)) /
                (double) (vector.length - 1);
        // We don't like negative variance
        if (result < 0) {
            return 0;
        } else {
            return result;
        }
    }

    public static double Accuracy(HashMap model) throws ToolKitException {
        if (model == null) {
            throw new ToolKitException("It is an empty model.");
        }
        double accuracy = 0;
        double truePos = 0;
        double trueNeg = 0;
        double falsePos = 0;
        double falseNeg = 0;
        Set keys = model.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Cell cell = (Cell) model.get(key);
            if (cell.getStatus() == 1) {
                truePos += cell.getPositiveScore();
                falseNeg += Math.abs(cell.getNegativeScore());
            } else if (cell.getStatus() == 0) {
                trueNeg += Math.abs(cell.getNegativeScore());
                falsePos += cell.getPositiveScore();
            } else {
                //what's going on if for a negative status;
            }
        }
        double denominator = truePos + trueNeg + falsePos + falseNeg;
        if(denominator == 0) {
            throw new ToolKitException("denominator is zero.");
        }
        accuracy = (truePos + trueNeg) / (denominator);
        return accuracy;
    }

    public static HashMap AccuracyLou(HashMap model, HashMap QTLProbabilityMatrix) throws ToolKitException {
    	HashMap posteriorProbability = new HashMap();
    	double accuracy = 0;
        if (model == null) {
            throw new ToolKitException("It is an empty model.");
        }
        double weightH = 0;
        double weightL = 0;
        double truePos = 0;
        double trueNeg = 0;
        double falsePos = 0;
        double falseNeg = 0;
        Set keys = model.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Cell cell = (Cell) model.get(key);
            if (cell.getStatus() == 1) {
            	weightH += ((Double) QTLProbabilityMatrix.get(key)).doubleValue();
                truePos += cell.getPositiveScore();
                falseNeg += Math.abs(cell.getNegativeScore());
            } else if (cell.getStatus() == 0) {
            	weightL += ((Double) QTLProbabilityMatrix.get(key)).doubleValue();
                trueNeg += Math.abs(cell.getNegativeScore());
                falsePos += cell.getPositiveScore();
            } else {
                //what's going on if for a negative status;
            }
        }
        for (Iterator e = keys.iterator(); e.hasNext();) {
        	String key = (String) e.next();
        	Cell cell = (Cell) model.get(key);
        	if (cell.getStatus() == 1) {
        		posteriorProbability.put(key, new Double(weightH));
        	} else if (cell.getStatus() == 0) {
        		posteriorProbability.put(key, new Double(weightL));
        	} else {
        		
        	}
        }
        return posteriorProbability;
    }

    public static double CGBStatistic(HashMap model) throws ToolKitException {
        double cgbStatistic = 0;

        return cgbStatistic;
    }

    public static int[] Sort(/*@non_null@*/double[] array) {
        int[] index = new int[array.length];
        array = (double[]) array.clone();
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
            if (Double.isNaN(array[i])) {
                array[i] = Double.MAX_VALUE;
            }
        }
        QuickSort(array, index, 0, array.length - 1);
        return index;
    }

    private static void QuickSort(double[] array,  int[] index, int left, int right) {
        if (left < right) {
            int middle = Partition(array, index, left, right);
            QuickSort(array, index, left, middle);
            QuickSort(array, index, middle + 1, right);
        }
    }

    
    private static int Partition(double[] array, int[] index, int l, int r) {
        double pivot = array[index[(l + r) / 2]];
        int help;
        while (l < r) {
            while ((array[index[l]] < pivot) && (l < r)) {
                l++;
            }
            while ((array[index[r]] > pivot) && (l < r)) {
                r--;
            }
            if (l < r) {
                help = index[l];
                index[l] = index[r];
                index[r] = help;
                l++;
                r--;
            }
        }
        if ((l == r) && (array[index[r]] > pivot)) {
            r--;
        }
        return r;
    }
}
