import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.regex.Pattern;

public class FamilyGenerator {

	boolean isNullHypothesis;
	String model;	
	int[] AffLoci;
	String[] FunctionalGenotype;	
	double[] recombination;	

	double intercept;
	double geneAffect;
	double covariable;	
	double deviation;	
	double error;
	double threshold;
	int FamilySize;
	int[] FamNum;
	int[] Kid_Diff_Family;
	int[] AffKid_Diff_Family;
	double[][] ParentGenotypeMissingRate;
	double KidGenotypeMissingRate;
	double[] DPrime;
	double[][] corMarkers;// it is calculated based on LD and allele frequencies.

	String[][] Allele;
	int numLocus;
	ArrayList Parents;
	ArrayList Children;
	ArrayList PCovariate;
	ArrayList PPhenotype;
	ArrayList CCovariate;
	ArrayList CPhenotype;
	ArrayList Traits;
	ArrayList PTraits;
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

		protected boolean isNullHypothesis; // 0
		protected long seed; // 1
		protected int simu_replication;// 2
		protected String model; // 3
		protected int[] AffLoci; // 4
		protected String[] FunctionalGenotype; // 5
		protected double[] Recombination; // 6		
		protected double[][] AlleleFreq; // 7
		protected double[] DPrime; // 8

		protected double intercept; // 9
		protected double gene; // 10
		protected double cov; // 11
		protected double dev; // 12
		protected double err; // 13
		protected double threshold; // 14
		protected int family_size; // 15
		protected int[] FamNum; // 16
		protected int[] Kid; // 17
		protected int[] AffKid; // 18
		protected double[][] ParentMissingRate; // 19
		protected double KidGenotypeMissingRate;// 20

		private double[][] corMarkers;

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
				} else if (Pattern.matches("^//.*", line)) {// multiple line
					// comments
					flag = false;
				} else {
					lines.add(line);
				}
			}
		}

		public void parseValue() {
			isNullHypothesis = Boolean.parseBoolean(lines.get(0));
			seed = Long.parseLong(lines.get(1));
			simu_replication = Integer.parseInt(lines.get(2));
			model = lines.get(3);
			String[] AL = lines.get(4).split(",");
			AffLoci = new int[AL.length];
			for (int i = 0; i < AL.length; i++) {
				AffLoci[i] = Integer.parseInt(AL[i]);
			}
			String[] fl = lines.get(5).split(",");
			FunctionalGenotype = new String[fl.length];
			System.arraycopy(fl, 0, FunctionalGenotype, 0, fl.length);	
			String[] Rec = lines.get(6).split(",");
			Recombination = new double[Rec.length];
			for (int i = 0; i < Recombination.length; i++) {
				Recombination[i] = Double.parseDouble(Rec[i]);
			}
			String[] AF = lines.get(7).split(",");
			AlleleFreq = new double[AF.length/Recombination.length][Recombination.length];
			for (int i = 0; i < AlleleFreq.length; i++) {
				for (int j = 0; j < AlleleFreq[i].length; j++) {
					AlleleFreq[i][j] = Double.parseDouble(AF[i*AlleleFreq[i].length + j]);
				}
			}
			String[] ld = lines.get(8).split(",");
			DPrime = new double[ld.length];
			for (int i = 0; i < DPrime.length; i++) {
				DPrime[i] = Double.parseDouble(ld[i]);
			}

			intercept = Double.parseDouble(lines.get(9));
			gene = Double.parseDouble(lines.get(10));
			cov = Double.parseDouble(lines.get(11));
			dev = Double.parseDouble(lines.get(12));
			err = Double.parseDouble(lines.get(13));
			threshold = Double.parseDouble(lines.get(14));

			family_size = Integer.parseInt(lines.get(15));
			String[] FM = lines.get(16).split(",");
			FamNum = new int[FM.length];
			for (int i = 0; i < FM.length; i++) {
				FamNum[i] = Integer.parseInt(FM[i]);
			}

			String[] K = lines.get(17).split(",");
			Kid = new int[K.length];
			for (int i = 0; i < K.length; i++) {
				Kid[i] = Integer.parseInt(K[i]);
			}
			String[] AK = lines.get(18).split(",");
			AffKid = new int[AK.length];
			for (int i = 0; i < AK.length; i++) {
				AffKid[i] = Integer.parseInt(AK[i]);
			}
			String[] MS = lines.get(19).split(",");
			ParentMissingRate = new double[MS.length / 2][2];
			for (int i = 0; i < MS.length / 2; i++) {
				ParentMissingRate[i][0] = Double.parseDouble(MS[i * 2 + 0]);
				ParentMissingRate[i][1] = Double.parseDouble(MS[i * 2 + 1]);
			}
			KidGenotypeMissingRate = Double.parseDouble(lines.get(20));
			calculateCorrelation_for_Markers_with_known_DPrime();
		}

		public void recipe() throws IOException{
			PrintWriter out = new PrintWriter(new File("recipe.txt"));
			out.println(this.toString());
			out.close();
		}

		protected void calculateCorrelation_for_Markers_with_known_DPrime() {
			corMarkers = new double[AlleleFreq.length][DPrime.length];
			for (int i = 0; i < corMarkers.length; i++) {
				for (int j = 0; j < corMarkers[i].length; j++) {
					double d = 0;
					if(DPrime[j] > 0) {
						d = DPrime[j] * Math.min(AlleleFreq[i][j]*(1 - AlleleFreq[i][j+1]), AlleleFreq[i][j+1]*(1-AlleleFreq[i][j]));
					} else {
						d = DPrime[j] * Math.min(AlleleFreq[i][j]* AlleleFreq[i][j+1], (1-AlleleFreq[i][j+1])*(1-AlleleFreq[i][j]));
					}
					corMarkers[i][j] = d/Math.sqrt(AlleleFreq[i][j]*(1 - AlleleFreq[i][j]) * AlleleFreq[i][j+1] * (1 - AlleleFreq[i][j+1]));					
				}
			}
		}

		public String toString() {
			StringBuffer sb = new StringBuffer();
			sb.append("was the data simulated under the null hypothesis: " + isNullHypothesis);
			sb.append("seed : " + seed);
			sb.append("\nreplication of simulation : " + simu_replication);
			sb.append("\nmodel : " + model );
			sb.append("\nFunctional loci : ");
			for (int i = 0; i < AffLoci.length; i++) {
				if (i == AffLoci.length -1 ) {
					sb.append(AffLoci[i]);
				} else {
					sb.append(AffLoci[i] + ",");
				}
			}
			sb.append("\nFunctional Genotype :");
			for (int i = 0; i < FunctionalGenotype.length; i++) {
				sb.append(FunctionalGenotype[i] + ",");
			}
			sb.append("\nAllele Frequency :");
			for (int i = 0; i < AlleleFreq.length; i++) {
				for (int j = 0; j < AlleleFreq[i].length; j++) {
					sb.append(AlleleFreq[i][j] + ",");
				}
				sb.append("\n");
			}			
			sb.append("\nRecombination : ");
			for (int i = 0; i < Recombination.length; i++) {
				sb.append(Recombination[i] + ",");
			}
			sb.append("\nD': ");
			for (int i = 0; i < DPrime.length; i++) {
				sb.append(DPrime[i] + ",");
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
		isNullHypothesis = p.isNullHypothesis;
		model = p.model;
		numLocus = p.Recombination.length;
		AffLoci = new int[p.AffLoci.length];
		System.arraycopy(p.AffLoci, 0, AffLoci, 0, p.AffLoci.length);

		Allele = new String[numLocus][2];
		for (int i = 0; i < Allele.length; i++) {
			Allele[i][0] = "1";
			Allele[i][1] = "2";
		}
		recombination = new double[p.Recombination.length];
		System.arraycopy(p.Recombination, 0, recombination, 0, p.Recombination.length);
		DPrime = new double[p.DPrime.length];
		System.arraycopy(p.DPrime, 0, DPrime, 0, p.DPrime.length);
		FunctionalGenotype = new String[p.FunctionalGenotype.length];
		System.arraycopy(p.FunctionalGenotype, 0, FunctionalGenotype, 0, p.FunctionalGenotype.length);

		intercept = p.intercept;
		geneAffect = p.gene;
		covariable = p.cov;
		deviation = p.dev;		
		error = p.err;
		threshold = p.threshold;
		
		FamilySize = p.family_size;
		Kid_Diff_Family = new int[p.Kid.length];
		System.arraycopy(p.Kid, 0, Kid_Diff_Family, 0, p.Kid.length);
		AffKid_Diff_Family = new int[p.AffKid.length];
		System.arraycopy(p.AffKid, 0, AffKid_Diff_Family, 0, p.AffKid.length);
		ParentGenotypeMissingRate = new double[p.ParentMissingRate.length][2];
		for (int i = 0; i < ParentGenotypeMissingRate.length; i++) {
			System.arraycopy(p.ParentMissingRate[i], 0,
					ParentGenotypeMissingRate[i], 0,
					p.ParentMissingRate[i].length);
		}

		KidGenotypeMissingRate = p.KidGenotypeMissingRate;
		FamNum = new int[p.FamNum.length];
		System.arraycopy(p.FamNum, 0, FamNum, 0, p.FamNum.length);

		corMarkers = new double[p.corMarkers.length][p.corMarkers[0].length];
		for (int i = 0; i < corMarkers.length; i++) {
			System.arraycopy(p.corMarkers[i], 0, corMarkers[i], 0, p.corMarkers[i].length);
		}
		
		DPrime = new double[p.DPrime.length];
		System.arraycopy(p.DPrime, 0, DPrime, 0, p.DPrime.length);
	}

	public int affection(double obs) {
		if (!isNullHypothesis ) {
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
		} else {
			return randomData.nextFloat() < 0.5 ? 0:1;
		}
	}

	private void generateFounderGenotype(String[][] pair, int FamCategory) {
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
					double r = (1 - corMarkers[FamCategory][k-1])/2;
					if (randomData.nextFloat() < r) {
						chrIdx = 1 - chrIdx;
					}
					pair[k][i] = Allele[k][chrIdx];
				}
			}
		}
	}

	public void create() {
		int i = 0;
		while (i < FamilySize) {
			ArrayList p_temp = new ArrayList();
			ArrayList p_trait = new ArrayList();
			ArrayList p_covariate = new ArrayList();
			ArrayList p_phenotype = new ArrayList();
			int FamCategory = 0;

			for (int j = 0; j < FamNum.length; j++) {
				if (i < FamNum[j]) {
					FamCategory = j;
					break;
				}
			}

			for (int j = 0; j < 2; ++j) {
				double covariate;
				double obs;
				String[][] pair = new String[numLocus][2];
				generateFounderGenotype(pair, FamCategory);
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


			for (int ii = 0; ii < Kid_Diff_Family[FamCategory]; ++ii) {
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
			if (AffKid_Diff_Family[FamCategory] > 0) {
				if (affect >= AffKid_Diff_Family[FamCategory]) {
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
			} else {
				if (affect == Math.abs(AffKid_Diff_Family[FamCategory])) {
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

		for (int i = 0; i < FunctionalGenotype.length; ++i) {
			if (genostr.compareTo(FunctionalGenotype[i]) == 0) {
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

		for (int i = 0; i < FamilySize; ++i) {
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
			pr.isNullHypothesis = false;
			pr.seed = 2008;
			pr.model = "B";
			pr.AffLoci = new int[2];
			pr.AffLoci[0] = 0;
			pr.AffLoci[1] = 1;
			pr.FunctionalGenotype = new String[4];
			pr.FunctionalGenotype[0] = "1112";
			pr.FunctionalGenotype[1] = "1211";
			pr.FunctionalGenotype[2] = "1222";
			pr.FunctionalGenotype[3] = "2212";
			pr.Recombination = new double[10];
			Arrays.fill(pr.Recombination, 0.5);
			pr.AlleleFreq = new double[3][10];
			for (int i = 0; i < pr.AlleleFreq.length; i++) {
				Arrays.fill(pr.AlleleFreq[i], 0.5);
			}
			pr.DPrime = new double[9];
			Arrays.fill(pr.DPrime, 0.7);
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
			pr.threshold = 2.05;
			pr.simu_replication = 2;
			pr.family_size = 200;
			pr.FamNum = new int[1];
			pr.FamNum[0] = 200;
			pr.Kid = new int[1];
			pr.Kid[0] = 2;
			pr.AffKid = new int[1];
			pr.AffKid[0] = -2;
			pr.ParentMissingRate = new double[1][2];
			pr.ParentMissingRate[0][0] = 0;
			pr.ParentMissingRate[0][1] = 0;
			pr.KidGenotypeMissingRate = 0;
			pr.calculateCorrelation_for_Markers_with_known_DPrime();
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