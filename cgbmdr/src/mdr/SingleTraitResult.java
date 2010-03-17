package mdr;

import java.util.AbstractList;
import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen
 */
public class SingleTraitResult extends AbstractList{
    ArrayList cvResult = new ArrayList();
    int interval;
    public SingleTraitResult (int itvl) {
        interval = itvl;
    }

    public void add(int idx, Object o) {
        modCount++;
        cvResult.add(idx, o);
    }

    public Object get(int idx) {
        return cvResult.get(idx);
    }

    public Object remove(int idx) {
        modCount++;
        return cvResult.remove(idx);
    }

    public int size() {
        return cvResult.size();
    }
}
