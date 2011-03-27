package admixture.population;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import admixture.population.genome.DNAStirrer;
import admixture.population.genome.GeneFlow;
import admixture.population.genome.HotSpot;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.genome.chromosome.FamilyGenome;
import admixture.population.genome.chromosome.FamilySingleChromosome;
import admixture.population.phenotype.FamilyPhenotype;
import admixture.population.phenotype.PhenotypeGenerator;
import admixture.population.phenotype.QualityControl;

public class GenerateColony {

	protected int n_chr;
	protected int n_phe;
	protected long seed;
	protected PhenotypeGenerator pheGenerator;
	protected HotSpot hs;
	protected Habitat FamHab;
	protected Habitat CaseControlHab;
	protected ArrayList<DNAStirrer> DNAPool;
	protected ArrayList<ChromosomeGenerator> ChrGenerator;
	protected int control_chr;
	protected double[] disease_rate;
	protected boolean recombinationFree;
	protected boolean isNullHypothesis;
	protected static int CurrFam = 0;

	public static class Builder {
		private int n_phe = 1;
		private long seed = 2011;
		private int control_chr = 0;
		private double[] disease_rate;
		private HotSpot hotspot;
		private ArrayList<DNAStirrer> DNAPool;
		private ArrayList<ChromosomeGenerator> ChrGenerator;
		private PhenotypeGenerator pheGenerator;
		private boolean recombinationFree = true;
		private boolean isNullHypothesis = true;
		private ArrayList<GeneFlow> GF;
		public Builder(double[] dr, HotSpot h, ArrayList<DNAStirrer> dp, ArrayList<ChromosomeGenerator> cg, PhenotypeGenerator p) {
			disease_rate = new double[dr.length];
			System.arraycopy(dr, 0, disease_rate, 0, dr.length);
			hotspot = h;
			DNAPool = dp;
			ChrGenerator = cg;
			pheGenerator = p;
		}
		
		public Builder numPhenotype(int np) { n_phe = np; return this;}
		
		public Builder seed(long s) { seed = s; return this;}
		
		public Builder diseaseChr(int dc) { control_chr = dc; return this;}
		
		public Builder recombinationFree(boolean rf) { recombinationFree = rf; return this;}
		
		public Builder isNullHypothesis(boolean isNull) { isNullHypothesis = isNull; return this;}

		public GenerateColony build() { return new GenerateColony(this); }
	}

	public GenerateColony (Builder builder) {
		disease_rate = builder.disease_rate;
		hs = builder.hotspot;
		DNAPool = builder.DNAPool;
		ChrGenerator = builder.ChrGenerator;
		pheGenerator = builder.pheGenerator;

		n_chr = DNAPool.size();
		
		n_phe = builder.n_phe;
		seed = builder.seed;
		control_chr = builder.control_chr;

		recombinationFree = builder.recombinationFree;
		isNullHypothesis = builder.isNullHypothesis;
	}

	public GenerateColony (int np, long s, int dc, double[] dr, HotSpot h, ArrayList<DNAStirrer> dp, ArrayList<ChromosomeGenerator> cg,
			PhenotypeGenerator p, boolean rf, boolean isNull) {
		n_phe = np;
		n_chr = dp.size();
		seed = s;
		control_chr = dc;
		disease_rate = dr;
		hs = h;
		DNAPool = dp;
		ChrGenerator = cg;
		pheGenerator = p;
		recombinationFree = rf;
		isNullHypothesis = isNull;

		hs.setSeed(seed);
		pheGenerator.setSeed(seed);
	}

	public static void setCurrFamilyID(int ci) {
		CurrFam = ci;
	}

	public void GenerateFamHab(int N_Fam, int N_Kid, QualityControl qc) {
		FamHab = new Habitat();
		generateFamilies(FamHab, N_Fam, N_Kid, qc);
	}
	
	public void GenerateCCHab(int N_Fam, int N_Kid, QualityControl qc) {
		CaseControlHab = new Habitat();
		generateFamilies(CaseControlHab, N_Fam, N_Kid, qc);
	}

	private void generateFamilies(Habitat hab, int N_Fam, int N_Kid, QualityControl qc) {
		for (int i = 0; i < N_Fam; i++) {
			FamilyGenome fg = new FamilyGenome(CurrFam + i + 1, N_Kid);
			FamilyPhenotype fp;
			int r = 0;
			do {
				for (int j = 0; j < DNAPool.size(); j++) {
					if (r > 0 && control_chr == j) {
						continue;
					}
					int chrID = j;
					DNAStirrer ds = DNAPool.get(j);
					ChromosomeGenerator cg = ChrGenerator.get(j);
					hs.rev(ds.NumberOfSNP());
//					hs.GenerateRecombination(AdmixtureConstant.free_recombination);
//					int[] f_hotspot = hs.getHotSpot();
//					hs.GenerateRecombination(AdmixtureConstant.free_recombination);
//					int[] m_hotspot = hs.getHotSpot();
					if (r == 0) {
						fg.addFamilyChromosome(cg.generateFamilySingleChromosome(chrID, N_Kid, hs,
								control_chr != j));
					} else {
						fg.setFamilyChromosome(j, cg.generateFamilySingleChromosome(chrID, N_Kid, hs, control_chr != j));
					}
				}
				if(!isNullHypothesis) {
					fp = pheGenerator.getGeneratePhenotypeAdmixtureLogistic(fg, disease_rate);
				} else {
					fp = pheGenerator.getGeneratePhenotypeAncestry(fg, disease_rate);
				}
				r++;
			} while (!qc.Accept(fp));
			hab.AddFamilyGenome(fg);
			hab.AddFamilyPhenotype(fp);
		}
		CurrFam += N_Fam;
	}

	public Habitat getFamHab() {
		return FamHab;
	}

	public Habitat getCCHab() {
		return CaseControlHab;
	}
	
	public int getNumberChromosome() {
		return n_chr;
	}
	
	public int getNumberPhenotype() {
		return n_phe;
	}

	public int getControlChr() {
		return control_chr;
	}

	public ArrayList<DNAStirrer> getDNAPool() {
		return DNAPool;
	}

	public void printFamilyGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();
		
		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(!printLinked) {
				if(i != control_chr) continue;
			} else {
				if(i == control_chr) continue;
			}
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();
		
		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}
		pheout.println();
		for(int f = 0; f < FamP.size(); f++) {
			//print phenotype
			FamilyPhenotype fp = FamP.get(f);
			StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
			}
			sp[0].append(fp.getStringParentPhenotype(0));
			sp[1].append(fp.getStringParentPhenotype(1));

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i + 2].append(fp.getStringOffspringPhenotype(i));
			}
			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = FamG.get(f);
			StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
			}
			sb[0].append(0 + " " + 0 + " " + 1 + " " + fp.getParentStatus(0) + " ");
			sb[1].append(0 + " " + 0 + " " + 2 + " " + fp.getParentStatus(1) + " ");

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sb[2+i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if(!printLinked) {
					if(fsc.isDiseaseLinked()) continue;
				} else {
					if(!fsc.isDiseaseLinked()) continue;
				}
				if(isAllele) {
					sb[0].append(fsc.getStringParentChromosome(0));
					sb[1].append(fsc.getStringParentChromosome(1));					
				} else {
					sb[0].append(fsc.getGenotypeStringParentChromosome(0));
					sb[1].append(fsc.getGenotypeStringParentChromosome(1));
				}
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i + 2].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();	
	}	

	public void printCCGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		
		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(!printLinked) {
				if(i != control_chr) continue;
			} else {
				if(i == control_chr) continue;
			}
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();
		
		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}
		pheout.println();
		ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
		for(int f = 0; f < CCP.size(); f++) {
			FamilyPhenotype fp = CCP.get(f);
			StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
				sp[i].append(fp.getStringOffspringPhenotype(i));
			}

			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			FamilyGenome fg = CCG.get(f);
			StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

			for(int i = 0; i < fg.getNumberOffspring(); i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				sb[i].append(0 + " " + 0 + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if(!printLinked) {
					if(fsc.isDiseaseLinked()) continue;
				} else {
					if(!fsc.isDiseaseLinked()) continue;
				}
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();	
	}

	public void printGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();
		
		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(!printLinked) {
				if(i != control_chr) continue;
			} else {
				if(i == control_chr) continue;
			}
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();
		
		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}
		pheout.println();
		for(int f = 0; f < FamP.size(); f++) {
			//print phenotype
			FamilyPhenotype fp = FamP.get(f);
			StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
			}
			sp[0].append(fp.getStringParentPhenotype(0));
			sp[1].append(fp.getStringParentPhenotype(1));

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i + 2].append(fp.getStringOffspringPhenotype(i));
			}
			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = FamG.get(f);
			StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
			}
			sb[0].append(0 + " " + 0 + " " + 1 + " " + fp.getParentStatus(0) + " ");
			sb[1].append(0 + " " + 0 + " " + 2 + " " + fp.getParentStatus(1) + " ");

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sb[2+i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if(!printLinked) {
					if(fsc.isDiseaseLinked()) continue;
				} else {
					if(!fsc.isDiseaseLinked()) continue;
				}
				if(isAllele) {
					sb[0].append(fsc.getStringParentChromosome(0));
					sb[1].append(fsc.getStringParentChromosome(1));					
				} else {
					sb[0].append(fsc.getGenotypeStringParentChromosome(0));
					sb[1].append(fsc.getGenotypeStringParentChromosome(1));
				}
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i + 2].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}

		//print case-control population 
		ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
		for(int f = 0; f < CCP.size(); f++) {
			FamilyPhenotype fp = CCP.get(f);
			StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
				sp[i].append(fp.getStringOffspringPhenotype(i));
			}

			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			FamilyGenome fg = CCG.get(f);
			StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

			for(int i = 0; i < fg.getNumberOffspring(); i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				sb[i].append(0 + " " + 0 + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if(!printLinked) {
					if(fsc.isDiseaseLinked()) continue;
				} else {
					if(!fsc.isDiseaseLinked()) continue;
				}
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}					
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();	
	}	

	public void printOffspringCC(String ped, String phe, boolean isAllele) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(i != control_chr) continue;
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();

		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}

		pheout.println();
		for(int f = 0; f < FamP.size(); f++) {
			//print phenotype
			FamilyPhenotype fp = FamP.get(f);
			StringBuffer[] sp = new StringBuffer[fp.getNumberAffectedOffspring()];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getOffspringID(i) + " ");
				sp[i].append(fp.getStringOffspringPhenotype(i));
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = FamG.get(f);
			StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getOffspringID(i) + " " 
						+ fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				for(int i = 0; i < sb.length; i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
					pedout.println(sb[i].toString());
				}
			}
		}

		//print case-control population 
		ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
		for(int f = 0; f < CCP.size(); f++) {
			FamilyPhenotype fp = CCP.get(f);
			StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getOffspringID(i) + " ");
				sp[i].append(fp.getStringOffspringPhenotype(i));
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = CCG.get(f);
			StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getOffspringID(i) + " ");
				sb[i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				for(int i = 0; i < sb.length; i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
					pedout.println(sb[i].toString());
				}
			}
		}
		pedout.close();
		pheout.close();
	}

	public void printFounder(String ped, String phe, boolean isAllele, boolean unrelatedOnly) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(i != control_chr) continue;
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();

		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}

		pheout.println();
		for(int f = 0; f < FamP.size(); f++) {
			//print phenotype
			FamilyPhenotype fp = FamP.get(f);
			StringBuffer[] sp = new StringBuffer[2];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
			}
			sp[0].append(fp.getStringParentPhenotype(0));
			sp[1].append(fp.getStringParentPhenotype(1));

			for(int i = 0; i < sp.length; i++) {
				if(unrelatedOnly && i >= 2) {
					continue;
				}
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = FamG.get(f);
			StringBuffer[] sb = new StringBuffer[2];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
			}
			sb[0].append(0 + " " + 0 + " " + 1 + " " + fp.getParentStatus(0) + " ");
			sb[1].append(0 + " " + 0 + " " + 2 + " " + fp.getParentStatus(1) + " ");

			int cnt = 0;
			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				if(isAllele) {
					sb[0].append(fsc.getStringParentChromosome(0));
					sb[1].append(fsc.getStringParentChromosome(1));					
				} else {
					sb[0].append(fsc.getGenotypeStringParentChromosome(0));
					sb[1].append(fsc.getGenotypeStringParentChromosome(1));
				}
				if(unrelatedOnly && cnt >= 2) {
					continue;
				}
				cnt++;
			}

			for(int i = 0; i < sb.length; i++) {
				if(unrelatedOnly && i >= 2) {
					continue;
				}
				pedout.println(sb[i].toString());
			}
		}

		//print case-control population 
		ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
		for(int f = 0; f < CCP.size(); f++) {
			FamilyPhenotype fp = CCP.get(f);
			StringBuffer[] sp = new StringBuffer[2];

			for(int i = 0; i < 2; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
				sp[i].append(fp.getStringParentPhenotype(i));
			}

			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = CCG.get(f);
			StringBuffer[] sb = new StringBuffer[2];

			for(int i = 0; i < 2; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				sb[i].append(0 + " " + 0 + " " + 1 + " " + fp.getParentStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				for(int i = 0; i < 2; i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringParentChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringParentChromosome(i));
					}
				}
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();
	}
	public void printUnrelatedIndividual(String ped, String phe, boolean isAllele, boolean unrelatedOnly) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

		pedout.print("FID ID FA MO SEX Affection ");
		for(int i = 0; i < n_chr; i++) {
			if(i != control_chr) continue;
			DNAStirrer ds = DNAPool.get(i);
			String[] SN = ds.getSNPNames();
			for(int j = 0; j < SN.length; j++) {
				pedout.print(SN[j] + " ");
			}
		}
		pedout.println();

		pheout.print("FID ID ");
		for(int i = 0; i < n_phe; i++) {
			pheout.print("phe" + i + " ");
		}

		pheout.println();
		for(int f = 0; f < FamP.size(); f++) {
			//print phenotype
			FamilyPhenotype fp = FamP.get(f);
			StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

			for(int i = 0; i < sp.length; i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
			}
			sp[0].append(fp.getStringParentPhenotype(0));
			sp[1].append(fp.getStringParentPhenotype(1));

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i + 2].append(fp.getStringOffspringPhenotype(i));
			}

			for(int i = 0; i < sp.length; i++) {
				if(unrelatedOnly && i >= 2) {
					continue;
				}
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = FamG.get(f);
			StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
			}
			sb[0].append(0 + " " + 0 + " " + 1 + " " + fp.getParentStatus(0) + " ");
			sb[1].append(0 + " " + 0 + " " + 2 + " " + fp.getParentStatus(1) + " ");

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sb[2+i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			int cnt = 0;
			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				if(isAllele) {
					sb[0].append(fsc.getStringParentChromosome(0));
					sb[1].append(fsc.getStringParentChromosome(1));					
				} else {
					sb[0].append(fsc.getGenotypeStringParentChromosome(0));
					sb[1].append(fsc.getGenotypeStringParentChromosome(1));
				}
				if(unrelatedOnly && cnt >= 2) {
					continue;
				}
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i + 2].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}
				cnt++;
			}

			for(int i = 0; i < sb.length; i++) {
				if(unrelatedOnly && i >= 2) {
					continue;
				}
				pedout.println(sb[i].toString());
			}
		}

		//print case-control population 
		ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
		for(int f = 0; f < CCP.size(); f++) {
			FamilyPhenotype fp = CCP.get(f);
			StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sp[i] = new StringBuffer();
				sp[i].append(fp.getFamilyID() + " " + fp.getFamilyID() * 10000 + i + " ");
				sp[i].append(fp.getStringOffspringPhenotype(i));
			}

			for(int i = 0; i < sp.length; i++) {
				pheout.println(sp[i].toString());
			}

			//print genotype
			FamilyGenome fg = CCG.get(f);
			StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

			for(int i = 0; i < fg.getNumberOffspring(); i++) {
				sb[i] = new StringBuffer();
				sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				sb[i].append(0 + " " + 0 + " " + 1 + " " + fp.getOffspringStatus(i) + " ");
			}

			for(FamilySingleChromosome fsc:fg) {
				if (fsc.isDiseaseLinked()) continue;
				for(int i = 0; i < fg.getNumberOffspring(); i++) {
					if(isAllele) {
						sb[i].append(fsc.getStringOffspringChromosome(i));
					} else {
						sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
					}
				}
			}

			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();
	}
}
