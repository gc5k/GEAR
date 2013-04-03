package family.qc.colqc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;

import parameter.Parameter;
import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;
import gear.util.NewIt;

public class SNPFilter implements SNPFilterInterface {

	protected MapFile mapData;
	protected int[] WSNP;

	protected int[] wseq = null;
	protected int[] bgseq = null;
	protected boolean filterFlag = false;
	ArrayList<SNP> snpList;
	protected HashSet<Integer> selectedSNPSet;
	protected HashSet<Integer> excludedSNPSet;
	protected HashSet<Integer> bgSNPSet;

	protected int[][] in_range = null;
	protected int[][] ex_range = null;

	public SNPFilter(MapFile mapData) {
		this.mapData = mapData;
		snpList = mapData.getMarkerList();
		selectedSNPSet = NewIt.newHashSet();
		excludedSNPSet = NewIt.newHashSet();
		bgSNPSet = NewIt.newHashSet();
	}

	public void Select() {
/*
		if (Parameter.bgsnpFlag) {
			selectBackgroundSNP();
			filterFlag = true;
		}
*/
		if (Parameter.INSTANCE.inchrFlag || Parameter.INSTANCE.exchrFlag) {
			selectChromosome();
			filterFlag = true;
		}
/*
//		if (Parameter.snpwindowFlag) {
//			selectSNPWindow();
//			filterFlag = true;
//		}

		if (Parameter.snpFlag) {
			selectSNPs();
			filterFlag = true;
		}


		if (Parameter.regionFlag) {
			selectSNPRegion();
			filterFlag = true;
		}

		if (Parameter.geneFlag) {
			selectGene();
			filterFlag = true;
		}
		
		if (filterFlag && selectedSNPSet.size() == 0
				&& excludedSNPSet.size() == 0 && bgSNPSet.size() == 0) {
			System.err.println("No snps selected. GMDR quit.");
			Test.LOG.append("No snps selected. GMDR quit.\n");
			Test.printLog();
			System.exit(0);
		}
		*/
		makeWSNPList();
		return;
	}
/*
	private void selectBackgroundSNP() {

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getName();
			int idx = ArrayUtils.indexOf(Parameter.bgsnp, chr);
			if (idx >= 0) {
				addBackgroundSNP(i);
			}
		}

	}
*/

	private void selectChromosome() {

		if (Parameter.INSTANCE.inchrFlag) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int idx = ArrayUtils.indexOf(Parameter.INSTANCE.inchr, chr);
				if (idx >= 0) {
					includeSNP(i);
				}
			}
		}

		if (Parameter.INSTANCE.exchrFlag) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int idx = ArrayUtils.indexOf(Parameter.INSTANCE.exchr, chr);
				if (idx >= 0) {
					excludeSNP(i);
				}
			}
		}
	}

	/*
	private void selectSNPWindow() {

		ArrayList<SNP> snpwindowList = NewIt.newArrayList();
		ArrayList<Integer> snpwindowIndex = NewIt.newArrayList();

		int c = 0;
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String rs = snp.getName();
			int idx = ArrayUtils.indexOf(Parameter.snpwindow, rs);
			if (idx >= 0) {
				snpwindowList.add(snp);
				snpwindowIndex.add(new Integer(idx));
				c++;
			}
			if (c == Parameter.snpwindow.length) {
				break;
			}
		}
		if (snpwindowList.size() == 0) {
			System.err.println("could not identify the window size.");
		}
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			for (int j = 0; j < snpwindowList.size(); j++) {
				int d = snpwindowList.get(j).compareTo(snp);
				if (d != Integer.MAX_VALUE && d != Integer.MIN_VALUE) {
					int idx = snpwindowIndex.get(j).intValue();
					if (d > Parameter.snp_window[idx][0]
							&& d < Parameter.snp_window[idx][1])
						includeSNP(i);
				}
			}
		}
	}
	*/
    /*
	private void selectGene() {

		int count = 0;
		ArrayList<ArrayList<Integer>> xsnps = NewIt.newArrayList();

		if(!Parameter.snp2genefileFlag && !Parameter.snp2genefilesFlag) {
			for (int i = 0; i < Parameter.gene.length; i++) {
				ArrayList<Integer> s = NewIt.newArrayList();
				xsnps.add(s);
			}

			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int pos = snp.getPosition();
				for (int j = 0; j < Parameter.gene_chr.length; j++) {

					if (chr.compareTo(Parameter.gene_chr[j]) == 0) {
						if (pos >= (Parameter.gene_begin[j] - Parameter.genewindow) * 1000
								&& pos <= (Parameter.gene_end[j] + Parameter.genewindow) * 1000) {
							xsnps.get(j).add(i);
							includeSNP(i);
							count++;
						}
					}
				}
			}
		} else {
			for (int i = 0; i < Parameter.gene.length; i++) {
				ArrayList<Integer> s = NewIt.newArrayList();
				xsnps.add(s);
			}

			for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext(); ) {
				int i = e.next().intValue();
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int pos = snp.getPosition();
				for (int j = 0; j < Parameter.gene_chr.length; j++) {

					if (chr.compareTo(Parameter.gene_chr[j]) == 0) {
						if (pos >= (Parameter.gene_begin[j] - Parameter.genewindow) * 1000
								&& pos <= (Parameter.gene_end[j] + Parameter.genewindow) * 1000) {
							xsnps.get(j).add(i);

							count++;
						}
					}
				}
			}

		}

		int countgene = 0;
		if (Parameter.snp2genefilesFlag) {
			for (int i = 0; i < xsnps.size(); i++) {
				StringBuffer sbsnp = new StringBuffer(Parameter.gene[i] + ".snp");
				StringBuffer sb = new StringBuffer(Parameter.gene[i] + ".gene");
				PrintStream PW = null;
				PrintStream PWsnp = null;
				try {
					PW = new PrintStream(sb.toString());
					PWsnp = new PrintStream(sbsnp.toString());
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}

				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				countgene++;
				System.err.print(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);
				System.err.println(" [writing snps into " + Parameter.gene[i]
						+ ".gene]");
				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);
				Test.LOG.append(" [writing snps into " + Parameter.gene[i]
						+ ".gene]\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PWsnp.println(snp.getName());
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() + " " + Parameter.gene[i]);
				}
				PWsnp.close();
				PW.close();
			}
			System.err.println(count + " snps selected for " + countgene + " selected genes in total.");
			Test.LOG.append(count + " snps selected for " + countgene + " selected genes in total.\n");
			Test.printLog();
			System.exit(1);
		} else {
			StringBuffer sbsnp = new StringBuffer(Parameter.out + ".snp");
			StringBuffer sb = new StringBuffer(Parameter.out + ".gene");
			PrintStream PW = null;
			PrintStream PWsnp = null;
			try {
				PW = new PrintStream(sb.toString());
				PWsnp = new PrintStream(sbsnp.toString());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < xsnps.size(); i++) {
				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				countgene++;
				System.err.println(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);

				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i] + "\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PWsnp.println(snp.getName());
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() + " " + Parameter.gene[i]);
				}
			}
			System.err.println("writing snps into " + Parameter.out + ".gene");
			Test.LOG.append("writing snps into " + Parameter.out + ".gene\n");
			PW.close();
			PWsnp.close();
			System.err.println(count + " snps selected for " + countgene + " selected genes in total.");
			Test.LOG.append(count + " snps selected for " + countgene + " selected genes in total.\n");
			if (Parameter.snp2genefileFlag) {
				Test.printLog();
				System.exit(1);
			}
		}

	}


	private void selectSNPRegion() {
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int pos = snp.getPosition();
			int idx = ArrayUtils.indexOf(Parameter.chr_reg, chr);
			if (idx >= 0) {
				if (pos >= Parameter.begin[idx] * 1000
						&& pos <= Parameter.end[idx] * 1000) {
					includeSNP(i);
				}
			}
		}
		return;
	}

	private void selectSNPs() {

		if (Parameter.includesnp != null) {
			HashSet<String> SS = NewIt.newHashSet();
			for (int i = 0; i < Parameter.includesnp.length; i++) {
				SS.add(Parameter.includesnp[i]);
			}
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				if(SS.contains(rs)) {
					includeSNP(i);
				}
			}
		}

		if (Parameter.excludesnp != null) {
			HashSet<String> SS = NewIt.newHashSet();
			for (int i = 0; i < Parameter.includesnp.length; i++) {
				SS.add(Parameter.includesnp[i]);
			}
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				if(!SS.contains(rs)) {
					includeSNP(i);
				}
			}
		}

		return;
	}
*/

	private void makeWSNPList() {

		if (selectedSNPSet.size() > 0) {
			WSNP = new int[selectedSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = selectedSNPSet.iterator(); e
						.hasNext();) {
				Integer V = e.next();
				WSNP[c++] = V.intValue();
			}
			Arrays.sort(WSNP);
		} else if (excludedSNPSet.size() > 0) {
			WSNP = new int[snpList.size() - excludedSNPSet.size()];
			int c = 0;
			for (int i = 0; i < snpList.size(); i++) {
				if (!excludedSNPSet.contains(i)) {
					WSNP[c++] = i;
				}
			}
		} else {
			WSNP = new int[snpList.size()];
			for (int i = 0; i < snpList.size(); i++) {
				WSNP[i] = i;
			}
		}

		wseq = new int[WSNP.length];
		for (int i = 0; i < WSNP.length; i++) {
			wseq[i] = i;
		}
	}

	private boolean includeSNP(int i) {
		if (selectedSNPSet.contains(i)) {
			return true;
		} else {
			selectedSNPSet.add(i);
			return false;
		}
	}

	private boolean excludeSNP(int i) {

		if (excludedSNPSet.contains(i)) {
			return true;
		} else {
			excludedSNPSet.add(i);
			return false;
		}
	}

	public int[] getWorkingSNP() {
		return WSNP;
	}

	@Override
	public int[] getWSeq() {
		return wseq;
	}

	@Override
	public int[] getBgSeq() {
		return bgseq;
	}

	@Override
	public int[][] getWSeq2() {
		return null;
	}
}
