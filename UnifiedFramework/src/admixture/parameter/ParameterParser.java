package admixture.parameter;

import java.util.ArrayList;
import java.util.Arrays;

import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;

public class ParameterParser {

	public static int[] selectedSNP(MapFile mf, String[] selectedSNP) {
		if(selectedSNP == null) {
			return null;
		}
		boolean[] flag = new boolean[selectedSNP.length];
		Arrays.fill(flag, false);
		int[] index = new int[selectedSNP.length];
		ArrayList<SNP> snpList = mf.getMarkerList();
		for (int i = 0; i < selectedSNP.length; i++) {
			for (int j = 0; j < snpList.size(); j++) {
				SNP snp = snpList.get(j);
				if(snp.getName().compareTo(selectedSNP[i]) == 0) {
					flag[i] = true;
					index[i] = j;
				}
			}
		}
		boolean f = true;
		for (int i = 0; i < flag.length; i++) {
			if (!flag[i]) {
				System.err.println("selected snp: " + selectedSNP[i] + " is unavailable.");
				f = false;
			}
		}
		if(!f) {
			System.exit(0);
			return null;
		} else {
			return index;
		}
	}
}
