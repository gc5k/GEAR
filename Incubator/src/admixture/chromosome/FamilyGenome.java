package admixture.chromosome;

import java.util.ArrayList;
import java.util.Iterator;

public class FamilyGenome extends ArrayList<FamilySingleChromosome> {
	private static final long serialVersionUID = 1L;
	private int FamID;
	private int numKid;
//	private ArrayList<FamilySingleChromosome> FG = new ArrayList<FamilySingleChromosome>();
	private double[][] ancestry_overall_p;//[2][ancestry]
	private double[][] ancestry_overall_o;//[number of kids][ancestry]
	public FamilyGenome(int fid, int num_kid) {
		FamID = fid;
		numKid = num_kid;
	}
	
	public void AscertainAncestry(double[][][] post_snp_ancestry) {
		ancestry_overall_p = new double[2][post_snp_ancestry[0][0].length];
		ancestry_overall_o = new double[numKid][post_snp_ancestry[0][0].length];
		for (Iterator<FamilySingleChromosome> i = this.iterator(); i.hasNext(); ) {
			FamilySingleChromosome fsc = i.next();
			fsc.AscertainParentSingleChromosomeAncestry(post_snp_ancestry);
			fsc.AscertainOffspringSingleChromosomeAncestry(post_snp_ancestry);
			double[][] a_p = fsc.getParentSingleChromosomeAncestry();
			double[][] a_o = fsc.getOffspringSingleChromosomeAncestry();
			
			for(int j = 0; j < a_p[0].length; j++) {
				ancestry_overall_p[0][j] += a_p[0][j];
				ancestry_overall_p[1][j] += a_p[1][j];
			}
			
			for(int j = 0; j < a_o.length; j++) {
				for (int k = 0; k < a_o[j].length; k++) {
					ancestry_overall_p[j][k] += a_o[j][k];
				}
			}
		}
	}
}
