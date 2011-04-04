package mdr;

import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen
 */
public class BestKFoldCVResult extends ArrayList {

    int order;

    public BestKFoldCVResult(int or, int traits, int interval) {
        order = or;
        for (int i = 0; i < traits; i++) {
            OneTraitResult i_thResult = new OneTraitResult(i);
            for (int j = 0; j < interval; j++) {
                OneCVSet cvPair = new OneCVSet(j, null);
                i_thResult.add(cvPair);
            }
            add(i_thResult);
        }
    }

    public String getKeyAt(int i, int j) {
        OneCVSet cvSet = (OneCVSet) ((OneTraitResult) get(i)).get(j);
        return cvSet.getCombination();
    }

    public OneCVSet get(int i, int j) {
        return (OneCVSet) ((OneTraitResult) get(i)).get(j);
    }

    public void set(int i, int j, OneCVSet cvSet) {
        ((OneTraitResult) get(i)).set(j, cvSet);
    }

    public int getOrder() {
        return order;
    }

    public void summarize() {
        for (int i = 0; i < size(); i++) {
            OneTraitResult i_thResult = (OneTraitResult) get(i);
            i_thResult.summarise();
        }
    }
}
