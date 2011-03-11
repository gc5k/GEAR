package admixture.chromosome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import org.apache.commons.lang3.ArrayUtils;

public class FamilyGenome extends ArrayList<FamilySingleChromosome> {
	private static final long serialVersionUID = 1L;
	private int FamID;
	private int numKid;
	// private ArrayList<FamilySingleChromosome> FG = new
	// ArrayList<FamilySingleChromosome>();
	private double[][][] ancestry_chr_p;// [control_chr][2; 0 for father, 1 for
										// mother][ancestry]
	private double[][][] ancestry_chr_o;// [control_chr][number of
										// kids][ancestry]
	private int[] disease_linked_chr;

	public FamilyGenome(int fid, int num_kid) {
		FamID = fid;
		numKid = num_kid;
		disease_linked_chr = null;
	}

	public void AscertainGenomeAncestry() {
		ancestry_chr_p = new double[NumChromosome()][][];
		ancestry_chr_o = new double[NumChromosome()][][];
		int chr = 0;
		for (Iterator<FamilySingleChromosome> i = this.iterator(); i.hasNext();) {
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

	public double[] ParentAncestry(int idx) {
		return cal_ancestry(ancestry_chr_p, idx);
	}

	public double[][] OffspringAncestry() {
		double[][] a_o = new double[numKid][];
		for (int i = 0; i < numKid; i++) {
			a_o[i] = cal_ancestry(ancestry_chr_o, i);
		}
		return a_o;
	}

	public double[] OffspringAncestry(int idx) {
		return cal_ancestry(ancestry_chr_o, idx);
	}

	private double[] cal_ancestry(double[][][] an, int idx) {
		double fAnc[] = new double[an[0][idx].length];
		int dlc = disease_linked_chr == null ? 0 : disease_linked_chr.length;
		for (int i = 0; i < an.length; i++) {
			if (dlc != 0 && Arrays.binarySearch(disease_linked_chr, i) < 0)
				continue;
			for (int j = 0; j < an[i][idx].length; j++) {
				fAnc[j] += an[i][idx][j];
			}
		}

		for (int i = 0; i < fAnc.length; i++) {
			fAnc[i] /= NumChromosome() - dlc;
		}
		return fAnc;
	}

	public String ParentGenotype(int pi, int[] chr, int[] loc) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < chr.length; i++) {
			FamilySingleChromosome fsc = get(chr[i]);
			int[] g = fsc.ParentGenotype(pi, loc[i]);
			if (g[0] > g[1]) {
				sb.append(g[0]);
				sb.append(g[1]);
			} else {
				sb.append(g[1]);
				sb.append(g[0]);
			}
		}
		return sb.toString();
	}

	public String OffspringGenotype(int pi, int[] chr, int[] loc) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < chr.length; i++) {
			FamilySingleChromosome fsc = get(chr[i]);
			int[] g = fsc.OffspringGenotype(pi, loc[i]);
			if (g[0] > g[1]) {
				sb.append(g[0]);
				sb.append(g[1]);
			} else {
				sb.append(g[1]);
				sb.append(g[0]);
			}
		}
		return sb.toString();
	}

	public void addFamilyChromosome(FamilySingleChromosome fsc) {
		if (fsc.isDiseaseLinked()) {
			disease_linked_chr = ArrayUtils.add(disease_linked_chr, fsc.getChrID());
		}
		add(fsc);
	}

	public int getNumberOffspring() {
		return numKid;
	}

	public int getFamilyID() {
		return FamID;
	}

	public void printGenome() {
		System.out.println("FamID " + FamID);
		for (Iterator<FamilySingleChromosome> i = this.iterator(); i.hasNext();) {
			FamilySingleChromosome fsc = i.next();
			int[][][] g_p = fsc.getParentChromosome();
			int[][][] g_o = fsc.getOffspingChromosome();
			for (int j = 0; j < g_p.length; j++) {
				System.out.println("Parent: " + j);
				for (int k = 0; k < g_p[j].length; k++) {
					for (int l = 0; l < g_p[j][k].length; l++) {
						System.out.print(g_p[j][k][l] + " ");
					}
					System.out.println();
				}
			}
			for (int j = 0; j < g_o.length; j++) {
				System.out.println("Kid: " + j);
				for (int k = 0; k < g_o[j].length; k++) {
					for (int l = 0; l < g_o[j][k].length; l++) {
						System.out.print(g_o[j][k][l] + " ");
					}
					System.out.println();
				}
			}
		}
	}

	public void printAncestry() {
		double[] a_f = ParentAncestry(0);
		double[] a_m = ParentAncestry(1);
		double[][] a_o = OffspringAncestry();
		System.out.println("FamID: " + FamID);
		System.out.println("Father:");
		for (int i = 0; i < a_f.length; i++) {
			System.out.print(a_f[i] + " ");
		}
		System.out.println();

		System.out.println("Mother: ");
		for (int i = 0; i < a_m.length; i++) {
			System.out.print(a_m[i] + " ");
		}
		System.out.println();

		for (int i = 0; i < a_o.length; i++) {
			System.out.println("Kid: " + i);
			for (int j = 0; j < a_o[i].length; j++) {
				System.out.print(a_o[i][j] + " ");
			}
			System.out.println();
		}
	}
}
