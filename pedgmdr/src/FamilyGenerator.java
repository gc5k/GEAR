import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FamilyGenerator {

	int samplesize;
	int numLocus;
	int dimension;
	String[] cell;
	int[] AffLoci;
	int[] Kid_Diff_Family;
	int[] AffKid_Diff_Family;
	String[][] Allele;
	double[][] AlleleFreq;
	double[] recombination;
	double[] LD;
	double[] DPrime;
	double[] corMarkers;// it is calculated based on LD and allele frequencies.
	int[] FamNum;
	double[][] ParentGenotypeMissingRate;
	double KidGenotypeMissingRate;

	ArrayList Parents;
	ArrayList Children;
	ArrayList PCovariate;
	ArrayList PPhenotype;
	ArrayList CCovariate;
	ArrayList CPhenotype;
	ArrayList Traits;
	ArrayList PTraits;
	double error;
	double deviation;
	double intercept;
	double covariable;
	double geneAffect;
	double threshold;
	String model;
	Random randomData;

	FamilyGenerator(long sd) {
		randomData = new Random(sd);
		Parents = new ArrayList();
		Children = new ArrayList();
		CCovariate = new ArrayList();
		CPhenotype = new ArrayList();
		PCovariate = new ArrayList();
		PPhenotype = new ArrayList();
		Traits = new ArrayList();
		PTraits = new ArrayList();
	}

	public static class Parameter {
		protected String filename;
		protected BufferedReader buffer;
		protected ArrayList<String> lines;

		protected long seed; // 0
		protected int simu_replication;// 1
		protected String model; // 2
		protected int[] AffLoci; // 3
		protected double[] AlleleFreq; // 4
		protected double[] Recombination; // 5
		protected double[] LD; //6
		protected String[] FunctionalGenotype; //7 
		
		protected double intercept; // 8
		protected double gene; // 9
		protected double cov; // 10
		protected double dev; // 11
		protected double err; // 12
		protected double threshold; // 13
		protected int family_size; // 14
		protected int[] FamNum; // 15	
		protected int[] Kid; // 16
		protected int[] AffKid; // 17
		protected double[][] ParentMissingRate; // 18
		protected double KidGenotypeMissingRate;// 19

		private double[] corMarkers;
		private double[] DPrime;
		

		public Parameter() {
			lines = new ArrayList();
		}

		public void read(String file) throws IOException {
			filename = file;
			buffer = new BufferedReader(new FileReader(new File(file)));
			sweepComments();
			parseValue();
		}

		public void sweepComments() throws IOException {
			boolean flag = true;
			String line;
			while ((line = buffer.readLine()) != null) {
				if (Pattern.matches("\\s*", line)) {// empty line
					continue;
				} else if (Pattern.matches("^\\s*//*.*", line)) {// multiple line
					// comments
					flag = false;
				} else {
					lines.add(line);
				}
			}
		}

		public void parseValue() {

			seed = Long.parseLong(lines.get(0));
			simu_replication = Integer.parseInt(lines.get(1));
			model = lines.get(2);
			String[] AL = lines.get(3).split(",");
			AffLoci = new int[AL.length];
			for (int i = 0; i < AL.length; i++) {
				AffLoci[i] = Integer.parseInt(AL[i]);
			}
			String[] AF = lines.get(4).split(",");
			AlleleFreq = new double[AF.length];
			for (int i = 0; i < AlleleFreq.length; i++) {
				AlleleFreq[i] = Double.parseDouble(AF[i]);
			}
			String[] Rec = lines.get(5).split(",");
			Recombination = new double[Rec.length];
			for (int i = 0; i < AlleleFreq.length; i++) {
				Recombination[i] = Double.parseDouble(Rec[i]);
			}
			String[] ld = lines.get(6).split(",");
			LD = new double[ld.length];
			for (int i = 0; i < LD.length; i++) {
				LD[i] = Double.parseDouble(ld[i]);
			}
			String[] fl = lines.get(7).split(",");
			FunctionalGenotype = new String[fl.length];
			System.arraycopy(fl, 0, FunctionalGenotype, 0, fl.length);			

			intercept = Double.parseDouble(lines.get(8));
			gene = Double.parseDouble(lines.get(9));
			cov = Double.parseDouble(lines.get(10));			
			dev = Double.parseDouble(lines.get(11));
			err = Double.parseDouble(lines.get(12));
			threshold = Double.parseDouble(lines.get(13));			

			family_size = Integer.parseInt(lines.get(14));
			String[] FM = lines.get(15).split(",");
			FamNum = new int[FM.length];
			for (int i = 0; i < FM.length; i++) {
				FamNum[i] = Integer.parseInt(FM[i]);
			}

			String[] K = lines.get(16).split(",");
			Kid = new int[K.length];
			for (int i = 0; i < K.length; i++) {
				Kid[i] = Integer.parseInt(K[i]);
			}
			String[] AK = lines.get(17).split(",");
			AffKid = new int[AK.length];
			for (int i = 0; i < AK.length; i++) {
				AffKid[i] = Integer.parseInt(AK[i]);
			}
			String[] MS = lines.get(18).split(",");
			ParentMissingRate = new double[MS.length / 2][2];
			for (int i = 0; i < MS.length / 2; i++) {
				ParentMissingRate[i][0] = Double.parseDouble(MS[i * 2 + 0]);
				ParentMissingRate[i][1] = Double.parseDouble(MS[i * 2 + 1]);
			}
			KidGenotypeMissingRate = Double.parseDouble(lines.get(19));
			calculateCorrelation_for_Markers();
			calculateDPrime();
		}
		
		public void recipe() throws IOException{
			PrintWriter out = new PrintWriter(new File("recipe.txt"));
			out.println(this.toString());
			out.close();
		}
		
		private void calculateCorrelation_for_Markers() {
			corMarkers = new double[LD.length];
			for (int i = 0; i < corMarkers.length; i++) {
				corMarkers[i] = LD[i]/
				Math.sqrt(AlleleFreq[i]*(1 - AlleleFreq[i]) * AlleleFreq[i+1] * (1 - AlleleFreq[i+1]));
				if (corMarkers[i] >= 1) {
					System.err.println("LD between marker " + i + " and marker " + (i+1) + " is irrational, and has been set to 0.");
					LD[i] = 0;
					corMarkers[i] = 0;
				}
			}
		}
		
		private void calculateDPrime() {
			DPrime = new double[LD.length];
			for (int i = 0; i < DPrime.length; i++) {
				if (LD[i] > 0) {
					double numerator1 = AlleleFreq[i] * (1 - AlleleFreq[i+1]);
					double numerator2 = (1 - AlleleFreq[i]) * AlleleFreq[i+1];
					if (numerator1 < numerator2) {
						DPrime[i] = LD[i] / numerator1;
					}
				} else {
					double numerator1 = AlleleFreq[i] * AlleleFreq[i+1];
					double numerator2 = (1 - AlleleFreq[i]) * (1 - AlleleFreq[i+1]);
					if (numerator1 < numerator2) {
						DPrime[i] = LD[i] / numerator2;
					}
				}
			}
		}
		
		public String toString() {
			StringBuffer sb = new StringBuffer();
			sb.append("seed : " + seed);
			sb.append("\nreplication of simulation : " + simu_replication);
			sb.append("\nmodel : " + model );
			sb.append("\nFunctional loci : ");
			for (int i = 0; i < AffLoci.length; i++) {
				sb.append(AffLoci[i] + ",");
			}
			sb.append("\nAllele Frequency :");
			for (int i = 0; i < AlleleFreq.length; i++) {
				sb.append(AlleleFreq[i] + ",");
			}			
			sb.append("\nRecombination : ");
			for (int i = 0; i < Recombination.length; i++) {
				sb.append(Recombination[i] + ",");
			}
			sb.append("\nLD(D'): ");
			for (int i = 0; i < LD.length; i++) {
				sb.append(LD[i] + "(" + DPrime[i] + "),");
			}
			sb.append("\nFunctional Genotype : ");
			for (int i = 0; i < FunctionalGenotype.length; i++) {
				sb.append(FunctionalGenotype[i] + ",");
			}
			sb.append("\nintercept : " + intercept);
			sb.append("\ngene : " + gene);
			sb.append("\nCovariable : " + cov);
			sb.append("\nSD Cov : " + dev);
			sb.append("\nS.D. Residual error : " + err);
			sb.append("\nThreshold : " + threshold);
			sb.append("\nfamily size : " + family_size);
			sb.append("\nFamily number in different categories : ");
			for (int i = 0; i < FamNum.length; i++) {
				sb.append(FamNum[i] + ",");
			}
			sb.append("\nKids number (affected) in different categories : ");
			for (int i = 0; i < Kid.length; i++) {
				sb.append(Kid[i] + "(" + AffKid[i] + ")" + "," );
			}
			sb.append("\nParent genotype Missing Rate : ");
			for (int i = 0; i < ParentMissingRate.length; i++) {
				sb.append(ParentMissingRate[i][0] + " " + ParentMissingRate[i][1] + ",");
			}
			sb.append("\ngenotype missing rate for kids : " + KidGenotypeMissingRate);			
			return sb.toString();
		}
	}

	public void SettingUpParameter(Parameter p) {
		model = p.model;
		numLocus = p.AlleleFreq.length;
		samplesize = p.family_size;
		geneAffect = p.gene;
		error = p.err;
		covariable = p.cov;
		deviation = p.dev;
		intercept = p.intercept;
		AffLoci = new int[p.AffLoci.length];
		System.arraycopy(p.AffLoci, 0, AffLoci, 0, p.AffLoci.length);
		KidGenotypeMissingRate = p.KidGenotypeMissingRate;
		FamNum = new int[p.FamNum.length];
		System.arraycopy(p.FamNum, 0, FamNum, 0, p.FamNum.length);
		ParentGenotypeMissingRate = new double[p.ParentMissingRate.length][2];
		for (int i = 0; i < ParentGenotypeMissingRate.length; i++) {
			System.arraycopy(p.ParentMissingRate[i], 0,
					ParentGenotypeMissingRate[i], 0,
					p.ParentMissingRate[i].length);
		}
		threshold = p.threshold;

		AlleleFreq = new double[p.AlleleFreq.length][2];
		for (int i = 0; i < AlleleFreq.length; i++) {
			AlleleFreq[i][0] = p.AlleleFreq[i];
			AlleleFreq[i][1] = 1 - p.AlleleFreq[i];
		}
		Allele = new String[AlleleFreq.length][2];
		for (int i = 0; i < Allele.length; i++) {
			Allele[i][0] = "1";
			Allele[i][1] = "2";
		}
		Kid_Diff_Family = new int[p.Kid.length];
		System.arraycopy(p.Kid, 0, Kid_Diff_Family, 0, p.Kid.length);
		AffKid_Diff_Family = new int[p.AffKid.length];
		System.arraycopy(p.AffKid, 0, AffKid_Diff_Family, 0, p.AffKid.length);
		recombination = new double[p.Recombination.length];
		System.arraycopy(p.Recombination, 0, recombination, 0,
				p.Recombination.length);
		LD = new double[p.LD.length];
		System.arraycopy(p.LD, 0, LD, 0, p.LD.length);
		cell = new String[p.FunctionalGenotype.length];
		System.arraycopy(p.FunctionalGenotype, 0, cell, 0, p.FunctionalGenotype.length);
		
		corMarkers = new double[p.corMarkers.length];
		System.arraycopy(p.corMarkers, 0, corMarkers, 0, p.corMarkers.length);
		
		DPrime = new double[p.DPrime.length];
		System.arraycopy(p.DPrime, 0, DPrime, 0, p.DPrime.length);
	}

	public int affection(double obs) {
		if (model.compareTo("B") == 0) {
			double prob = Math.exp(obs) / (1D + Math.exp(obs));
			if (randomData.nextFloat() < prob) {
				return 1;
			}
			return 0;
		}

		if (obs > threshold) {
			return 1;
		}
		return 0;
	}

	private void generateFounderGenotype(String[][] pair) {
		for (int i = 0; i < 2; i++) {
			int chrIdx = 0;
			for (int k = 0; k < numLocus; k++) {
				if(k == 0) {
					if (randomData.nextFloat() < 0.5) {
						chrIdx = 0;
					} else {
						chrIdx = 1;
					}
					pair[k][i] = Allele[k][chrIdx];
				} else {
					double r = (1 - corMarkers[k-1])/2;
					if (randomData.nextFloat() < r) {
						chrIdx = 1 - chrIdx;
					}
					pair[k][i] = Allele[k][chrIdx];
				}
			}
		}
	}

	public void create() {
		int NumAffectedKid = 0;
		int i = 0;
		while (i < samplesize) {
			int sib;
			ArrayList p_temp = new ArrayList();
			ArrayList p_trait = new ArrayList();
			ArrayList p_covariate = new ArrayList();
			ArrayList p_phenotype = new ArrayList();
			for (int j = 0; j < 2; ++j) {
				double covariate;
				double obs;
				String[][] pair = new String[numLocus][2];
				generateFounderGenotype(pair);
//				for (int l = 0; l < numLocus; ++l) {
//					for (int k = 0; k < 2; ++k) {
//						double rd = randomData.nextFloat();
//						int index = 0;
//						double fq = AlleleFreq[l][index];
//						while (rd > fq) {
//							fq += AlleleFreq[l][(++index)];
//						}
//						pair[l][k] = Allele[l][index];
//					}
//				}
				p_temp.add(pair);
				Integer genestatus = getStatus(pair);

				if (model.compareTo("B") == 0) {
					covariate = randomData.nextGaussian() * deviation;
					obs = genestatus.doubleValue() * geneAffect + covariate
							* covariable + intercept;
				} else {
					covariate = randomData.nextGaussian() * deviation;
					obs = genestatus.doubleValue() * geneAffect + covariate
							* covariable + intercept
							+ randomData.nextGaussian() * error;
					p_phenotype.add(new Double(obs));
				}
				Integer status = new Integer(affection(obs));
				p_covariate.add(new Double(covariate));
				p_trait.add(status);
			}

			ArrayList c_temp = new ArrayList();
			ArrayList c_trait = new ArrayList();
			ArrayList c_covariate = new ArrayList();
			ArrayList c_phenotype = new ArrayList();
			int affect = 0;

			if (i < FamNum[0]) {
				NumAffectedKid = AffKid_Diff_Family[0];
				sib = Kid_Diff_Family[0];
			} else if (i < FamNum[1]) {
				NumAffectedKid = AffKid_Diff_Family[1];
				sib = Kid_Diff_Family[1];
			} else {
				NumAffectedKid = AffKid_Diff_Family[2];
				sib = Kid_Diff_Family[2];
			}

			for (int ii = 0; ii < sib; ++ii) {
				double covariate;
				double obs;
				String[][] child = new String[numLocus][2];
				for (int j = 0; j < 2; ++j) {
					String[][] P = (String[][]) (String[][]) p_temp.get(j);
					int chr = 0;
					for (int k = 0; k < numLocus; ++k) {
						double rd = randomData.nextFloat();
						if (rd > recombination[k]) {
							chr = 1 - chr;
						}
						child[k][j] = P[k][chr];
					}
				}

				Integer genestatus = getStatus(child);

				if (model.compareTo("B") == 0) {
					covariate = randomData.nextGaussian() * deviation;
					obs = genestatus.doubleValue() * geneAffect + covariate
							* covariable + intercept;
				} else {
					covariate = randomData.nextGaussian() * deviation;
					obs = genestatus.doubleValue() * geneAffect + covariate
							* covariable + intercept
							+ randomData.nextGaussian() * error;
					c_phenotype.add(new Double(obs));
				}

				Integer status = new Integer(affection(obs));
				affect += status.intValue();
				c_covariate.add(new Double(covariate));
				c_trait.add(status);
				c_temp.add(child);
			}

			if (affect >= NumAffectedKid) {
				Children.add(c_temp);
				Parents.add(p_temp);
				PCovariate.add(p_covariate);
				PPhenotype.add(p_phenotype);
				CCovariate.add(c_covariate);
				CPhenotype.add(c_phenotype);
				PTraits.add(p_trait);
				Traits.add(c_trait);
				++i;
			}
		}
	}

	public Integer getStatus(String[][] geno) {
		String genostr = new String();
		for (int i = 0; i < AffLoci.length; ++i) {
			if (geno[AffLoci[i]][0].compareTo(geno[AffLoci[i]][1]) > 0) {
				genostr = genostr + geno[AffLoci[i]][1] + geno[AffLoci[i]][0];
			} else {
				genostr = genostr + geno[AffLoci[i]][0] + geno[AffLoci[i]][1];
			}
		}

		for (int i = 0; i < cell.length; ++i) {
			if (genostr.compareTo(cell[i]) == 0) {
				return new Integer(1);
			}
		}
		return new Integer(0);
	}

	public void print(String ped, String phe) throws IOException {
		PrintWriter pedout = new PrintWriter(new File(ped));
		PrintWriter pheout = new PrintWriter(new File(phe));

		pedout.print("FID\tID\tFA\tM0\tSex\tAffection\t");
		for (int i = 0; i < Allele.length; ++i) {
			if (i != Allele.length - 1) {
				pedout.print("M" + (i + 1) + "\t");
			} else {
				pedout.print("M" + (i + 1));
			}
		}
		pedout.println();

		if (model.compareTo("B") == 0) {
			pheout.println("FIP\tID\taffection\tcov");
		} else {
			pheout.println("FIP\tID\taffection\tphenotype\tcov");
		}

		for (int i = 0; i < samplesize; ++i) {
			double[] MR = null;
			for (int j = 0; j < FamNum.length; j++) {
				if (i < FamNum[j]) {
					MR = ParentGenotypeMissingRate[j];
					break;
				}
			}
			int FID = 20000 + i;
			int ID = FID * 10;
			int FA = ID;
			int MO = ID + 1;
			ArrayList temp = (ArrayList) Parents.get(i);
			ArrayList pt_temp = (ArrayList) PTraits.get(i);
			ArrayList pc_temp = (ArrayList) PCovariate.get(i);
			ArrayList pp_temp = null;
			if (model.compareTo("B") != 0) {
				pp_temp = (ArrayList) PPhenotype.get(i);
			}
			String[][] pair1 = (String[][]) (String[][]) temp.get(0);
			String[][] pair2 = (String[][]) (String[][]) temp.get(1);

			pedout.print(FID + "\t" + ID + "\t0\t0\t0\t"
					+ (((Integer) pt_temp.get(0)).intValue() + 1) + "\t");

			for (int j = 0; j < Allele.length; ++j) {
				double rd = randomData.nextFloat();
				if (rd > MR[0]) {
					pedout.print(pair1[j][0] + " " + pair1[j][1]);
				} else {
					pedout.print("0 0");
				}
				if (j != Allele.length - 1) {
					pedout.print("\t");
				}
			}

			pedout.println();
			if (model.compareTo("B") == 0) {
				pheout.println(FID + "\t" + ID + "\t" + pt_temp.get(0) + "\t"
						+ pc_temp.get(0));
			} else {
				pheout.println(FID + "\t" + ID + "\t" + pt_temp.get(0) + "\t"
						+ pp_temp.get(0) + "\t" + pc_temp.get(0));
			}
			++ID;

			pedout.print(FID + "\t" + ID + "\t0\t0\t1\t"
					+ (((Integer) pt_temp.get(1)).intValue() + 1) + "\t");
			for (int j = 0; j < Allele.length; ++j) {
				double rd = randomData.nextFloat();
				if (rd > MR[1]) {
					pedout.print(pair2[j][0] + " " + pair2[j][1]);
				} else {
					pedout.print("0 0");
				}
				if (j != Allele.length - 1) {
					pedout.print("\t");
				}
			}

			pedout.println();
			if (model.compareTo("B") == 0) {
				pheout.println(FID + "\t" + ID + "\t" + pt_temp.get(1) + "\t"
						+ pc_temp.get(1));
			} else {
				pheout.println(FID + "\t" + ID + "\t" + pt_temp.get(1) + "\t"
						+ pp_temp.get(1) + "\t" + pc_temp.get(1));
			}

			++ID;

			ArrayList children = (ArrayList) Children.get(i);
			ArrayList t_temp = (ArrayList) Traits.get(i);
			ArrayList cc_temp = (ArrayList) CCovariate.get(i);
			ArrayList cp_temp = null;
			if (model.compareTo("B") != 0) {
				cp_temp = (ArrayList) CPhenotype.get(i);
			}
			for (int j = 0; j < children.size(); ++j) {
				String[][] child = (String[][]) (String[][]) children.get(j);
				int sex = (randomData.nextBoolean()) ? 1 : 0;
				pedout.print(FID + "\t" + ID + "\t" + FA + "\t" + MO + "\t"
						+ sex + "\t"
						+ (((Integer) t_temp.get(j)).intValue() + 1) + "\t");

				for (int k = 0; k < child.length; ++k) {
					double rd = randomData.nextFloat();
					if (rd > KidGenotypeMissingRate) {
						pedout.print(child[k][0] + " " + child[k][1]);
					} else {
						pedout.print("0 0");
					}
					if (k != child.length - 1) {
						pedout.print("\t");
					}
				}

				pedout.println();
				if (model.compareTo("B") == 0) {
					pheout.println(FID + "\t" + ID + "\t" + t_temp.get(j)
							+ "\t" + cc_temp.get(j));
				} else {
					pheout.println(FID + "\t" + ID + "\t" + t_temp.get(j)
							+ "\t" + cp_temp.get(j) + "\t" + cc_temp.get(j));
				}
				++ID;
			}
		}
		pedout.close();
		pheout.close();
	}

	public void setAffect(double gene) {
		geneAffect = gene;
	}

	public void setParameters(double intc, double cov, double dev, double err) {
		error = err;
		covariable = cov;
		deviation = dev;
		intercept = intc;
	}

	public void setSeed(long seed) {
		randomData.setSeed(seed);
	}

	public void setThreshold(double t) {
		threshold = t;
	}

	public static void main(String[] args) {
		Parameter pr = new Parameter();
		if (args.length > 0) {
			try {
				pr.read(args[0]);
			} catch (IOException E) {
				E.printStackTrace(System.err);
			}
		} else {
			pr.model = "B";
			if (pr.model.compareTo("B") == 0) {
				pr.intercept = -5.29;
				pr.dev = 3.16;
				pr.cov = 1;
				pr.err = 1;
				pr.gene = 1.09;
			} else {
				pr.intercept = 0;
				pr.dev = 1;
				pr.cov = 1;
				pr.err = 1;
				pr.gene = 2.05;
			}
			pr.seed = 2008;
			pr.simu_replication = 2;
			pr.family_size = 600;
			pr.KidGenotypeMissingRate = 0;
			pr.FamNum = new int[3];
			pr.FamNum[0] = 200;
			pr.FamNum[1] = 400;
			pr.FamNum[2] = pr.family_size;
			pr.Kid = new int[3];
			pr.Kid[0] = 2;
			pr.Kid[1] = 3;
			pr.Kid[2] = 4;
			pr.AffKid = new int[3];
			pr.AffKid[0] = 2;
			pr.AffKid[1] = 2;
			pr.AffKid[2] = 2;
			pr.ParentMissingRate = new double[3][2];
			pr.ParentMissingRate[0][0] = 0;
			pr.ParentMissingRate[0][1] = 0;
			pr.ParentMissingRate[1][0] = 0;
			pr.ParentMissingRate[1][1] = 1;
			pr.ParentMissingRate[2][0] = 1;
			pr.ParentMissingRate[2][1] = 1;

			pr.threshold = 2.05;
		}

		// means AABb, Aabb, AaBB, AABb

		for (int i = 0; i < pr.simu_replication; ++i) {
			System.out.println("Simulation:" + i);
			String Ped = new Integer(i).toString() + ".ped";
			String Phe = new Integer(i).toString() + ".phe";
			FamilyGenerator simu = new FamilyGenerator(pr.seed * 1000 + i);
			simu.SettingUpParameter(pr);
			simu.create();
			try {
				simu.print(Ped, Phe);
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
		try {
			pr.recipe();
		} catch (IOException E) {
			E.printStackTrace(System.err);
		}
	}
}