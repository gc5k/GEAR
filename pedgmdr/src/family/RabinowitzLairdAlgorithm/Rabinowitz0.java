package family.RabinowitzLairdAlgorithm;

import java.util.TreeMap;
//both parents' genotype were observed
public class Rabinowitz0 extends AbstractGenoDistribution {

    TreeMap parentGenoMap;
    public Rabinowitz0(TreeMap children, TreeMap parents) {
        super(children);
        this.parentGenoMap = new TreeMap(parents);
        genotypeParents();
    }

    public void genotypeParents() {
        parentGeno.add((String) parentGenoMap.firstKey());
        parentGeno.add((String) parentGenoMap.lastKey());
    }

    public String[] getNontransmitted(final String Transmited) {
        return null;
    }

    public String[] getNontransmitted() {
        String[] control = new String[getChildrenNum()];
        TreeMap controlMap = new TreeMap();
        for (int i = 0; i < control.length; i++) {
            control[i] = Transmit(controlMap);
        }
        return control;
    }

    protected String Transmit(TreeMap cM) {
        String geno = new String(RandomAssign());

        if (cM.containsKey(geno)) {
            Integer c = ((Integer) cM.get(geno));
            int v = (c.intValue());
            v++;
            c = new Integer(v);
            cM.put(geno, c);
        } else {
            Integer c = new Integer(1);
            cM.put(new String(geno), c);
        }
        return geno;
    }
}
