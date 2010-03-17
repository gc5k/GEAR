package im;

import im.population.*;
import java.util.HashMap;
import java.util.ArrayList;
import im.population.IMPopulation;
import regression.LinearRegression;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BuildScore {

    protected IMPopulation imp;
    public BuildScore(IMPopulation im) {
        imp = im;
    }

    public void BuildScore(int sIdx, boolean env) {
        ArrayList<Integer> observed = new ArrayList();
        for(int i = 0; i < imp.IndividualNumber(); i++) {
            if(imp.PhenotypeExist(i, sIdx)) {
                observed.add(new Integer(i));
            }
        }
        double[][] phe_y = new double[observed.size()][1];

        double[][] pre_x;
        if(env) {
            pre_x = new double[observed.size()][2];
            for (int i = 0; i < observed.size(); i++) {
                int idx = observed.get(i).intValue();
                phe_y[i][0] = imp.PhenotypeAt(idx, sIdx);
                pre_x[i][0] = 1;
                pre_x[i][1] = Integer.parseInt(imp.NestedEnvironment(idx));
            }
        } else {
            pre_x = new double[observed.size()][1];
            for (int i = 0; i < observed.size(); i++) {
                int idx = observed.get(i).intValue();
                phe_y[i][0] = imp.PhenotypeAt(idx, sIdx);
                pre_x[i][0] = 1;
            }
        }

        LinearRegression lm = new LinearRegression(pre_x, phe_y);
        lm.MLE();
        double[] s = lm.getResiduals();
        for (int i = 0; i < observed.size(); i++) {
            imp.setScoreAt(observed.get(i).intValue(), sIdx, s[i]);
        }
    }
}
