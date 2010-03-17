package mdr.moore;


import mdr.Suite;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;

import mdr.data.DataFile;

/**
 *
 * @author Guo-Bo Chen
 */
public class PSuite extends Suite {

    private int numTraits;
    private ArrayList refs = new ArrayList();

    public PSuite(Suite s, int traits) {
        super(s);
        initial(traits);
    }

    public PSuite(int traits) {
        super(traits);
        initial(traits);
    }

    public void initial(int traits) {
        numTraits = traits;
        for (int i = 0; i < numTraits; i++) {
            refs.add(new ArrayList());
        }        
    }

    public void addSubject(DataFile.Subject sub, ArrayList traitsIndex) {
        add(sub);
        Integer id = sub.getIntegerID();
        for (int i = 0; i < numTraits; i++) {
            ArrayList ref = (ArrayList) refs.get(i);
            HashMap traitindex = (HashMap) traitsIndex.get(i);
            ref.add( traitindex.get(id));
        }
    }

    public void summarize(int idx, double offset, ArrayList Trait) {
        posSubjects[idx] = 0;
        negSubjects[idx] = 0;
        posScore[idx] = 0;
        negScore[idx] = 0;
        ArrayList trait = (ArrayList) Trait.get(idx);
        ArrayList ref = (ArrayList) refs.get(idx);
        for (Iterator e = ref.iterator(); e.hasNext();) {
            int i = ((Integer) e.next()).intValue();
            double scr = ((Double) (trait.get(i))).doubleValue() - offset;
            if ((scr) >= 0) {
                posSubjects[idx]++;
                posScore[idx] += scr;
            } else {
                negSubjects[idx]++;
                negScore[idx] += scr;
            }
        }
    }
    
    public void print(int idx) {
        System.out.println("Trait " + idx + " " + posSubjects[idx] + " " + negSubjects[idx] + " " + posScore[idx] + " " + negScore[idx]);     
    }
}
