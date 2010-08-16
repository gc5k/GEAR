package im.population.simulation;

/**
 *
 * @author Guo-Bo Chen
 */
public class AbstractLoci {
    int[] chr;
    int[] loci;
    double[] location;
    int[] genotype;
    double[] effect;

    int env;
    public AbstractLoci(int[] c, int[] l, int[] g, double[] e, int environment) {
        chr = new int[c.length];
        System.arraycopy(c, 0, chr, 0, c.length);
        loci = new int[l.length];
        System.arraycopy(l, 0, loci, 0, l.length);
        genotype = new int[g.length];
        System.arraycopy(g, 0, genotype, 0, g.length);
        effect = new double[e.length];
        System.arraycopy(e, 0, effect, 0, e.length);
        env = environment;
    }

    public AbstractLoci(int[] c, double[] lc, int[] g, double[] e, int environment) {
        chr = new int[c.length];
        System.arraycopy(c, 0, chr, 0, c.length);
        loci = new int[c.length];
        location = new double[lc.length];
        System.arraycopy(lc, 0, location, 0, lc.length);
        genotype = new int[g.length];
        System.arraycopy(g, 0, genotype, 0, g.length);
        effect = new double[e.length];
        System.arraycopy(e, 0, effect, 0, e.length);
        env = environment;
    }
    
    public void setLoci(int lc) {
    	loci[0] = lc;
    }

    public int[] getChr() {
        return chr;
    }

    public int[] getLoci() {
        return loci;
    }

    public double[] getLocation() {
    	return location;
    }

    public int[] getGenotype() {
        return genotype;
    }

    public double[] getEffect() {
        return effect;
    }

    public int getEnvironment() {
        return env;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("chr"+ " " + "loci"+"\n");
        for(int i = 0; i < chr.length; i++) {
            sb.append(chr[i] + " " + loci[i] + "\n");
        }
        sb.append("genotype effect\n");
        for(int i = 0; i < genotype.length; i++) {
            sb.append(genotype[i]+ " " + effect[i]);
        }
        sb.append("environment " + env + "\n");
        return sb.toString();
    }
}
