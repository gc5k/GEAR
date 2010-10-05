package family.RabinowitzLairdAlgorithm;

import java.util.TreeMap;

public class Rabinowitz1 extends AbstractGenoDistribution {

    TreeMap parentGenoMap;
    String parentgeno1;

    public Rabinowitz1(TreeMap child, TreeMap parent) {
        super(child);
        parentGenoMap = new TreeMap(parent);
        this.parentgeno1 = (String) parentGenoMap.firstKey();
        countChildrenAllele(childrenGenoMap);
        countAllele(childrenGenoMap);
        countParentAllele(parentGenoMap);
        countAllele(parentGenoMap);
    }

    protected void genotypeParents() {
    }

    public String[] getNontransmitted(final String transmitted) {
        return null;
    }

    public String[] getNontransmitted() {
        String control[] = new String[getChildrenNum()];
        TreeMap controlMap = new TreeMap();
        // how to deal with missing data???

        if (childrenGenoMap.size() == 1) {// situation 1
            // System.out.println("Table1 s1");

            String[] genopool = new String[1];
            genopool[0] = (String) childrenGenoMap.firstKey();
            double[] freq = new double[1];
            freq[0] = 1.0;
            Produce(control, controlMap, genopool, freq);
        } else if (childrenGenoMap.size() == 2) {// situation 2, 3
            // System.out.println("Table1 s2,3");

            String[] genopool = new String[2];
            genopool[0] = (String) childrenGenoMap.firstKey();
            genopool[1] = (String) childrenGenoMap.lastKey();
            double[] freq = new double[2];
            freq[0] = 0.5;
            freq[1] = 1.0;
            do {
                controlMap.clear();
                Produce(control, controlMap, genopool, freq);
            } while (!(controlMap.size() > 1));
        } else {
            System.err.println("Wrecked in Rabinowitz table 1");
        }
        return control;
    }
}
