package family.RabinowitzLairdAlgorithm.rabinowitz;

import java.util.TreeMap;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;

import util.NewIt;
//both parents' genotype were observed
public class Rabinowitz0 extends AbstractGenoDistribution {

    TreeMap<String, Integer> parentGenoMap;
    public Rabinowitz0(TreeMap<String, Integer> children, TreeMap<String, Integer> parents) {
        super(children);
        this.parentGenoMap = new TreeMap<String, Integer>(parents);
        genotypeParents();
    }

    public void genotypeParents() {
        parentGeno.add(parentGenoMap.firstKey());
        parentGeno.add(parentGenoMap.lastKey());
    }

    public String[] getNontransmitted(final String Transmited) {
        return null;
    }

    public String[] getNontransmitted() {
        String[] control = new String[getChildrenNum()];
        TreeMap<String, Integer> controlMap = NewIt.newTreeMap();
        for (int i = 0; i < control.length; i++) {
            control[i] = Transmit(controlMap);
        }
        return control;
    }

    protected String Transmit(TreeMap<String, Integer> cM) {
        String geno = new String(RandomAssign());

        if (cM.containsKey(geno)) {
            Integer c = ((Integer) cM.get(geno));
            c++;
            cM.put(geno, c);
        } else {
            Integer c = new Integer(1);
            cM.put(new String(geno), c);
        }
        return geno;
    }
}
