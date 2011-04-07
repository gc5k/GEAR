package admixture.population;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import admixture.population.genome.DNAStirrer;
import admixture.population.genome.chromosome.FamilyGenome;
import admixture.population.genome.chromosome.FamilySingleChromosome;
import admixture.population.phenotype.FamilyPhenotype;
import arsenal.NewIt;

public class Print2File {
	ArrayList<GenerateColony> colony;

	public Print2File() {
		colony = NewIt.newArrayList();
	}

	public void addColony(GenerateColony gc) {
		colony.add(gc);
	}

	public void printGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));

		int colony_count = 0;
		for (GenerateColony gc : colony) {
			Habitat FamHab = gc.getFamHab();
			ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

			ArrayList<DNAStirrer> DNAPool = gc.getDNAPool();
			if (colony_count == 0) {
				pedout.print("FID ID FA MO SEX Affection ");
				for (int i = 0; i < gc.getNumberChromosome(); i++) {
					if(!printLinked) {
						if (i != gc.getControlChr()) continue;
					} else {
						if (i == gc.getControlChr()) continue;
					}
					DNAStirrer ds = DNAPool.get(i);
					String[] SN = ds.getSNPNames();
					for (int j = 0; j < SN.length; j++) {
						pedout.print(SN[j] + " ");
					}
				}
				pedout.println();

				pheout.print("FID ID ");
				for (int i = 0; i < gc.getNumberPhenotype(); i++) {
					pheout.print("phe" + i + " ");
				}
				pheout.println();
			}

			for (int f = 0; f < FamP.size(); f++) {
				// print phenotype
				FamilyPhenotype fp = FamP.get(f);
				StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

				for (int i = 0; i < sp.length; i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getIndividualID(i) + " ");
				}
				sp[0].append(fp.getStringParentPhenotype(0));
				sp[1].append(fp.getStringParentPhenotype(1));

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i + 2].append(fp.getStringOffspringPhenotype(i));
				}
				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				// print genotype
				FamilyGenome fg = FamG.get(f);
				StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
				for (int i = 0; i < sb.length; i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				}
				sb[0].append(0 + " " + 0 + " " + 1 + " " + (fp.getParentStatus(0) + 1) + " ");
				sb[1].append(0 + " " + 0 + " " + 2 + " " + (fp.getParentStatus(1) + 1) + " ");

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sb[2 + i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " "
							+ (fp.getOffspringStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (!printLinked) {
						if (fsc.isDiseaseLinked()) continue;
					} else {
						if (!fsc.isDiseaseLinked()) continue;
					}
					if (isAllele) {
						sb[0].append(fsc.getStringParentChromosome(0));
						sb[1].append(fsc.getStringParentChromosome(1));
					} else {
						sb[0].append(fsc.getGenotypeStringParentChromosome(0));
						sb[1].append(fsc.getGenotypeStringParentChromosome(1));
					}
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i + 2].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}

			// print case-control population
			Habitat CaseControlHab = gc.getCCHab();
			ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
			for (int f = 0; f < CCP.size(); f++) {
				FamilyPhenotype fp = CCP.get(f);
				StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getOffspringID(i) + " ");
					sp[i].append(fp.getStringOffspringPhenotype(i));
				}

				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				FamilyGenome fg = CCG.get(f);
				StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

				for (int i = 0; i < fg.getNumberOffspring(); i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getOffspringID(i) + " ");
					sb[i].append(0 + " " + 0 + " " + 1 + " " + (fp.getOffspringStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (!printLinked) {
						if (fsc.isDiseaseLinked()) continue;
					} else {
						if (!fsc.isDiseaseLinked()) continue;
					}
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}
			colony_count++;
		}
		pedout.close();
		pheout.close();
	}
	public void printCCGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));

		int colony_count = 0;
		for (GenerateColony gc : colony) {

			ArrayList<DNAStirrer> DNAPool = gc.getDNAPool();
			if (colony_count == 0) {
				pedout.print("FID ID FA MO SEX Affection ");
				for (int i = 0; i < gc.getNumberChromosome(); i++) {
					if(!printLinked) {
						if (i != gc.getControlChr()) continue;
					} else {
						if (i == gc.getControlChr()) continue;
					}
					DNAStirrer ds = DNAPool.get(i);
					String[] SN = ds.getSNPNames();
					for (int j = 0; j < SN.length; j++) {
						pedout.print(SN[j] + " ");
					}
				}
				pedout.println();

				pheout.print("FID ID ");
				for (int i = 0; i < gc.getNumberPhenotype(); i++) {
					pheout.print("phe" + i + " ");
				}
				pheout.println();
			}

			// print case-control population
			Habitat CaseControlHab = gc.getCCHab();
			ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
			for (int f = 0; f < CCP.size(); f++) {
				FamilyPhenotype fp = CCP.get(f);
				StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getOffspringID(i) + " ");
					sp[i].append(fp.getStringOffspringPhenotype(i));
				}

				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				FamilyGenome fg = CCG.get(f);
				StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

				for (int i = 0; i < fg.getNumberOffspring(); i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getOffspringID(i) + " ");
					sb[i].append(0 + " " + 0 + " " + 1 + " " + (fp.getOffspringStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (!printLinked) {
						if (fsc.isDiseaseLinked()) continue;
					} else {
						if (!fsc.isDiseaseLinked()) continue;
					}
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}
			colony_count++;
		}
		pedout.close();
		pheout.close();
	}

	public void printFamilyGenotype2file(String ped, String phe, boolean isAllele, boolean printLinked) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));

		int colony_count = 0;
		for (GenerateColony gc : colony) {
			Habitat FamHab = gc.getFamHab();
			ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

			ArrayList<DNAStirrer> DNAPool = gc.getDNAPool();
			if (colony_count == 0) {
				pedout.print("FID ID FA MO SEX Affection ");
				for (int i = 0; i < gc.getNumberChromosome(); i++) {
					if(!printLinked) {
						if (i != gc.getControlChr()) continue;
					} else {
						if (i == gc.getControlChr()) continue;
					}
					DNAStirrer ds = DNAPool.get(i);
					String[] SN = ds.getSNPNames();
					for (int j = 0; j < SN.length; j++) {
						pedout.print(SN[j] + " ");
					}
				}
				pedout.println();

				pheout.print("FID ID ");
				for (int i = 0; i < gc.getNumberPhenotype(); i++) {
					pheout.print("phe" + i + " ");
				}
				pheout.println();
			}

			for (int f = 0; f < FamP.size(); f++) {
				// print phenotype
				FamilyPhenotype fp = FamP.get(f);
				StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

				for (int i = 0; i < sp.length; i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getIndividualID(i) + " ");
				}
				sp[0].append(fp.getStringParentPhenotype(0));
				sp[1].append(fp.getStringParentPhenotype(1));

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i + 2].append(fp.getStringOffspringPhenotype(i));
				}
				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				// print genotype
				FamilyGenome fg = FamG.get(f);
				StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
				for (int i = 0; i < sb.length; i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				}
				sb[0].append(0 + " " + 0 + " " + 1 + " " + (fp.getParentStatus(0) + 1) + " ");
				sb[1].append(0 + " " + 0 + " " + 2 + " " + (fp.getParentStatus(1) + 1) + " ");

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sb[2 + i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " "
							+ (fp.getOffspringStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (!printLinked) {
						if (fsc.isDiseaseLinked()) continue;
					} else {
						if (!fsc.isDiseaseLinked()) continue;
					}
					if (isAllele) {
						sb[0].append(fsc.getStringParentChromosome(0));
						sb[1].append(fsc.getStringParentChromosome(1));
					} else {
						sb[0].append(fsc.getGenotypeStringParentChromosome(0));
						sb[1].append(fsc.getGenotypeStringParentChromosome(1));
					}
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i + 2].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}
			colony_count++;
		}
		pedout.close();
		pheout.close();
	}

	
	public void printUnrelatedIndividual(String ped, String phe, boolean isAllele, boolean unrelatedOnly)
			throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));

		int colony_count = 0;
		for (GenerateColony gc : colony) {

			ArrayList<DNAStirrer> DNAPool = gc.getDNAPool();
			if (colony_count == 0) {
				pedout.print("FID ID FA MO SEX Affection ");
				for (int i = 0; i < gc.getNumberChromosome(); i++) {
					if (i == gc.getControlChr())
						continue;
					DNAStirrer ds = DNAPool.get(i);
					String[] SN = ds.getSNPNames();
					for (int j = 0; j < SN.length; j++) {
						pedout.print(SN[j] + " ");
					}
				}
				pedout.println();

				pheout.print("FID ID ");
				for (int i = 0; i < gc.getNumberPhenotype(); i++) {
					pheout.print("phe" + i + " ");
				}
				pheout.println();
			}
			Habitat FamHab = gc.getFamHab();
			ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

			for (int f = 0; f < FamP.size(); f++) {
				// print phenotype
				FamilyPhenotype fp = FamP.get(f);
				StringBuffer[] sp = new StringBuffer[2 + fp.getNumberOffspring()];

				for (int i = 0; i < sp.length; i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getIndividualID(i) + " ");
				}
				sp[0].append(fp.getStringParentPhenotype(0));
				sp[1].append(fp.getStringParentPhenotype(1));

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i + 2].append(fp.getStringOffspringPhenotype(i));
				}

				for (int i = 0; i < sp.length; i++) {
					if (unrelatedOnly && i >= 2) {
						continue;
					}
					pheout.println(sp[i].toString());
				}

				// print genotype

				FamilyGenome fg = FamG.get(f);
				StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
				for (int i = 0; i < sb.length; i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				}
				sb[0].append(0 + " " + 0 + " " + 1 + " " + (fp.getParentStatus(0) + 1) + " ");
				sb[1].append(0 + " " + 0 + " " + 2 + " " + (fp.getParentStatus(1) + 1) + " ");

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sb[2 + i].append(fg.getFatherID() + " " + fg.getMotherID() + " " + 1 + " "
							+ (fp.getOffspringStatus(i) + 1) + " ");
				}

				int cnt = 0;
				for (FamilySingleChromosome fsc : fg) {
					if (fsc.isDiseaseLinked())
						continue;
					if (isAllele) {
						sb[0].append(fsc.getStringParentChromosome(0));
						sb[1].append(fsc.getStringParentChromosome(1));
					} else {
						sb[0].append(fsc.getGenotypeStringParentChromosome(0));
						sb[1].append(fsc.getGenotypeStringParentChromosome(1));
					}
					if (unrelatedOnly && cnt >= 2) {
						continue;
					}
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i + 2].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i + 2].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
					cnt++;
				}

				for (int i = 0; i < sb.length; i++) {
					if (unrelatedOnly && i >= 2) {
						continue;
					}
					pedout.println(sb[i].toString());
				}
			}

			// print case-control population
			Habitat CaseControlHab = gc.getCCHab();
			ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
			for (int f = 0; f < CCP.size(); f++) {
				FamilyPhenotype fp = CCP.get(f);
				StringBuffer[] sp = new StringBuffer[fp.getNumberOffspring()];

				for (int i = 0; i < fp.getNumberOffspring(); i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getOffspringID(i) + " ");
					sp[i].append(fp.getStringOffspringPhenotype(i));
				}

				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				// print genotype
				FamilyGenome fg = CCG.get(f);
				StringBuffer[] sb = new StringBuffer[fg.getNumberOffspring()];

				for (int i = 0; i < fg.getNumberOffspring(); i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getOffspringID(i) + " ");
					sb[i].append(0 + " " + 0 + " " + 1 + " " + (fp.getOffspringStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (fsc.isDiseaseLinked())
						continue;
					for (int i = 0; i < fg.getNumberOffspring(); i++) {
						if (isAllele) {
							sb[i].append(fsc.getStringOffspringChromosome(i));
						} else {
							sb[i].append(fsc.getGenotypeStringOffspringChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}
			colony_count++;
		}
		pedout.close();
		pheout.close();
	}

	public void printFounder(String ped, String phe, boolean isAllele, boolean unrelatedOnly) throws IOException {
		PrintWriter pedout = new PrintWriter(new BufferedWriter(new FileWriter(ped)));
		PrintWriter pheout = new PrintWriter(new BufferedWriter(new FileWriter(phe)));

		int colony_count = 0;
		for (GenerateColony gc : colony) {

			Habitat FamHab = gc.getFamHab();
			ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();

			ArrayList<DNAStirrer> DNAPool = gc.getDNAPool();
			if (colony_count == 0) {
				pedout.print("FID ID FA MO SEX Affection ");
				for (int i = 0; i < gc.getNumberChromosome(); i++) {
					if (i == gc.getControlChr())
						continue;
					DNAStirrer ds = DNAPool.get(i);
					String[] SN = ds.getSNPNames();
					for (int j = 0; j < SN.length; j++) {
						pedout.print(SN[j] + " ");
					}
				}
				pedout.println();

				pheout.print("FID ID ");
				for (int i = 0; i < gc.getNumberPhenotype(); i++) {
					pheout.print("phe" + i + " ");
				}
				pheout.println();
			}
			for (int f = 0; f < FamP.size(); f++) {
				// print phenotype
				FamilyPhenotype fp = FamP.get(f);
				StringBuffer[] sp = new StringBuffer[2];

				for (int i = 0; i < sp.length; i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getIndividualID(i) + " ");
				}
				sp[0].append(fp.getStringParentPhenotype(0));
				sp[1].append(fp.getStringParentPhenotype(1));

				for (int i = 0; i < sp.length; i++) {
					if (unrelatedOnly && i >= 2) {
						continue;
					}
					pheout.println(sp[i].toString());
				}

				// print genotype
				FamilyGenome fg = FamG.get(f);
				StringBuffer[] sb = new StringBuffer[2];
				for (int i = 0; i < sb.length; i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
				}
				sb[0].append(0 + " " + 0 + " " + 1 + " " + (fp.getParentStatus(0) + 1) + " ");
				sb[1].append(0 + " " + 0 + " " + 2 + " " + (fp.getParentStatus(1) + 1) + " ");

				int cnt = 0;
				for (FamilySingleChromosome fsc : fg) {
					if (fsc.isDiseaseLinked())
						continue;
					if (isAllele) {
						sb[0].append(fsc.getStringParentChromosome(0));
						sb[1].append(fsc.getStringParentChromosome(1));
					} else {
						sb[0].append(fsc.getGenotypeStringParentChromosome(0));
						sb[1].append(fsc.getGenotypeStringParentChromosome(1));
					}
					if (unrelatedOnly && cnt >= 2) {
						continue;
					}
					cnt++;
				}

				for (int i = 0; i < sb.length; i++) {
					if (unrelatedOnly && i >= 2) {
						continue;
					}
					pedout.println(sb[i].toString());
				}
			}

			// print case-control population
			Habitat CaseControlHab = gc.getCCHab();
			ArrayList<FamilyGenome> CCG = CaseControlHab.getFamilyGenome();
			ArrayList<FamilyPhenotype> CCP = CaseControlHab.getFamilyPhenotype();
			for (int f = 0; f < CCP.size(); f++) {
				FamilyPhenotype fp = CCP.get(f);
				StringBuffer[] sp = new StringBuffer[2];

				for (int i = 0; i < 2; i++) {
					sp[i] = new StringBuffer();
					sp[i].append(fp.getFamilyID() + " " + fp.getIndividualID(i) + " ");
					sp[i].append(0 + " " + 0 + " " + 1 + (fp.getParentStatus(i) + 1)+ " ");
				}

				for (int i = 0; i < sp.length; i++) {
					pheout.println(sp[i].toString());
				}

				// print genotype
				FamilyGenome fg = CCG.get(f);
				StringBuffer[] sb = new StringBuffer[2];

				for (int i = 0; i < 2; i++) {
					sb[i] = new StringBuffer();
					sb[i].append(fg.getFamilyID() + " " + fg.getIndividualID(i) + " ");
					sb[i].append(0 + " " + 0 + " " + 1 + " " + (fp.getParentStatus(i) + 1) + " ");
				}

				for (FamilySingleChromosome fsc : fg) {
					if (fsc.isDiseaseLinked())
						continue;
					for (int i = 0; i < 2; i++) {
						if (isAllele) {
							sb[i].append(fsc.getStringParentChromosome(i));
						} else {
							sb[i].append(fsc.getGenotypeStringParentChromosome(i));
						}
					}
				}

				for (int i = 0; i < sb.length; i++) {
					pedout.println(sb[i].toString());
				}
			}
			pedout.close();
			pheout.close();
			colony_count++;
		}
	}
}
