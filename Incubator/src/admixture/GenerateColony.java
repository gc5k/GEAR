package admixture;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import admixture.chromosome.FamilyGenome;
import admixture.chromosome.FamilySingleChromosome;
import admixture.phenotype.FamilyPhenotype;
import admixture.phenotype.PhenotypeGenerator;
import admixture.phenotype.QualityControl;

public class GenerateColony {

	private long seed;
	private PhenotypeGenerator pg;
	private HotSpot hs;
	private Habitat FamHab;
	private Habitat CaseControlHab;
	private ArrayList<DNAStirrer> DNAPool;
	private ArrayList<ChromosomeGenerator> ChrGenerator;
	private int disease_chr;
	private double[] disease_rate;
	private boolean recombination_free;
	private int CurrFam;
	public GenerateColony (long s, int dc, double[] dr, HotSpot h, ArrayList<DNAStirrer> dp, ArrayList<ChromosomeGenerator> cg,
			PhenotypeGenerator p, boolean rf) {
		seed = s;
		disease_chr = dc;
		disease_rate = dr;
		hs = h;
		DNAPool = dp;
		ChrGenerator = cg;
		pg = p;
		recombination_free = rf;

		hs.setSeed(seed);
		pg.setSeed(seed);
	}

	public void GenerateNewFamHab(int N_Fam, int N_Kid, QualityControl qc) {
		FamHab = new Habitat();
		generateFamilies(FamHab, N_Fam, N_Kid, qc);
	}
	
	public void GenerateCCHab(int N_Fam, int N_Kid, QualityControl qc) {
		CaseControlHab = new Habitat();
		generateFamilies(CaseControlHab, N_Fam, N_Kid, qc);
	}
	
	private void generateFamilies(Habitat hab, int N_Fam, int N_Kid, QualityControl qc) {
		for (int i = 0; i < N_Fam; i++) {
			FamilyGenome fg = new FamilyGenome(i, N_Kid);
			FamilyPhenotype fp;
			int r = 0;
			do {
				for (int j = 0; j < DNAPool.size(); j++) {
					if (r > 0 && disease_chr != j) {
						continue;
					}
					int chrID = j;
					DNAStirrer ds = DNAPool.get(j);
					ChromosomeGenerator cg = ChrGenerator.get(j);
					hs.rev(ds.NumberOfSNP());
					hs.GenerateRecombination(AdmixtureConstant.free_recombination);
					int[] f_hotspot = hs.getHotSpot();
					hs.GenerateRecombination(AdmixtureConstant.free_recombination);
					int[] m_hotspot = hs.getHotSpot();
					if (r == 0) {
						fg.addFamilyChromosome(cg.generateFamilySingleChromosome(chrID, N_Kid, f_hotspot, m_hotspot,
								ds.PostSNPAncestralProb(), disease_chr == j));
					} else {
						fg.setFamilyChromosome(j, cg.generateFamilySingleChromosome(chrID, N_Kid, f_hotspot,
								m_hotspot, ds.PostSNPAncestralProb(), disease_chr == j));
					}
				}
				fp = pg.getGeneratePhenotypeAncestry(fg, disease_rate);
				r++;
			} while (!qc.Accept(fp));
			fp.print();
			fg.printGenome();
			hab.AddFamilyGenome(fg);
			hab.AddFamilyPhenotype(fp);
		}
		CurrFam += N_Fam;
	}

	public void print2file(String ped, String phe) throws IOException {
		PrintWriter pedout = new PrintWriter(new File(ped));
		PrintWriter pheout = new PrintWriter(new File(phe));
		ArrayList<FamilyGenome> FamG = FamHab.getFamilyGenome();
		ArrayList<FamilyPhenotype> FamP = FamHab.getFamilyPhenotype();
		
		for(FamilyGenome fg:FamG) {
			StringBuffer[] sb = new StringBuffer[2 + fg.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
			}
			for(FamilySingleChromosome fsc:fg) {
				sb[0].append(fsc.getStringParentChromosome(0));
				sb[1].append(fsc.getStringParentChromosome(1));
				for(int i = 0; i < fg.getNumberOffspring(); i++)
				sb[i+2].append(fsc.getStringOffspringChromosome(i));
			}
			
			for(int i = 0; i < sb.length; i++) {
				pedout.println(sb[i].toString());
			}
		}
		
		for(FamilyPhenotype fp:FamP) {
			StringBuffer[] sb = new StringBuffer[2 + fp.getNumberOffspring()];
			for(int i = 0; i < sb.length; i++) {
				sb[i] = new StringBuffer();
			}
			sb[0].append(fp.getStringParentPhenotype(0));
			sb[1].append(fp.getStringParentPhenotype(1));
			for(int i = 0; i < fp.getNumberOffspring(); i++) {
				sb[i + 2].append(fp.getStringOffspringPhenotype(i));
			}
			for(int i = 0; i < sb.length; i++) {
				pheout.println(sb[i].toString());
			}
		}
		pedout.close();
		pheout.close();
	}
}
