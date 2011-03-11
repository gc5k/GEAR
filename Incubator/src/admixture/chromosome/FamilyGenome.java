package admixture.chromosome;

import java.util.ArrayList;
import java.util.Iterator;

public class FamilyGenome extends ArrayList<FamilySingleChromosome> {
	private static final long serialVersionUID = 1L;
	private int FamID;
	private int numKid;
//	private ArrayList<FamilySingleChromosome> FG = new ArrayList<FamilySingleChromosome>();
	private double[][][] ancestry_chr_p;//[chr][2; 0 for father, 1 for mother][ancestry]
	private double[][][] ancestry_chr_o;//[chr][number of kids][ancestry]
	public FamilyGenome(int fid, int num_kid) {
		FamID = fid;
		numKid = num_kid;
	}

	public void AscertainGenomeAncestry() {
		ancestry_chr_p = new double[NumChromosome()][][];
		ancestry_chr_o = new double[NumChromosome()][][];
		int chr = 0;
		for (Iterator<FamilySingleChromosome> i = this.iterator(); i.hasNext(); ) {
			FamilySingleChromosome fsc = i.next();
			double[][] a_p = fsc.getParentChromosomeAncestry();
			double[][] a_o = fsc.getOffspringChromosomeAncestry();
			ancestry_chr_p[chr] = a_p;
			ancestry_chr_o[chr] = a_o;
			chr++;
		}
	}

	public int NumChromosome() {
		return this.size();
	}

	public double[] FatherAncestry() {
		return cal_ancestry(ancestry_chr_p, 0);
	}

	public double[] MotherAncestry() {
		return cal_ancestry(ancestry_chr_p, 1);
	}

	public double[][] OffspringAncestry() {
		double[][] a_o = new double[numKid][];
		for(int i = 0; i < numKid; i++) {
			a_o[i] = cal_ancestry(ancestry_chr_o, i);
		}
		return a_o;
	}

	private double[] cal_ancestry(double[][][] an, int idx) {
		double fAnc[] = new double[an[0][idx].length];

		for(int i = 0; i < an.length; i++) {
			for(int j = 0; j < an[i][idx].length; j++) {
				fAnc[j] += an[i][idx][j];
			}
		}

		for(int i = 0; i < fAnc.length; i++) {
			fAnc[i] /= NumChromosome();
		}
		return fAnc;
	}
	
	public String ParentGenotype(int pi, int[] chr, int[] loc) {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < chr.length; i++) {
			FamilySingleChromosome fsc = get(chr[i]);
			int[] g = fsc.ParentGenotype(pi, loc[i]);
			sb.append(g[0]);
			sb.append(g[1]);
		}
		return sb.toString();
	}

	public String OffspringGenotype(int pi, int[] chr, int[] loc) {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < chr.length; i++) {
			FamilySingleChromosome fsc = get(chr[i]);
			int[] g = fsc.OffspringGenotype(pi, loc[i]);
			sb.append(g[0]);
			sb.append(g[1]);
		}
		return sb.toString();
	}

	public int getNumberOffspring() {
		return numKid;
	}

	public void printAncestry() {
		double[] a_f = FatherAncestry();
		double[] a_m = MotherAncestry();
		double[][] a_o = OffspringAncestry();
		System.out.println("FamID: " + FamID);
		System.out.println("Father:");
		for(int i = 0; i < a_f.length; i++) {
			System.out.print(a_f[i] + " "); 
		}
		System.out.println();

		System.out.println("Mother: ");		
		for(int i = 0; i < a_m.length; i++) {
			System.out.print(a_m[i] + " ");
		}
		System.out.println();

		for(int i = 0; i < a_o.length; i++) {
			System.out.println("Kid: " + i);
			for(int j = 0; j < a_o[i].length; j++) {
				System.out.print(a_o[i][j] + " ");
			}
			System.out.println();
		}
	}
}
