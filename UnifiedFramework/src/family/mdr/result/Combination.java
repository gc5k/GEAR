
package family.mdr.result;

import java.util.HashMap;
import java.util.Map.Entry;

import publicAccess.PublicData;

import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;


/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Combination extends HashMap<String, Suite> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	
	public String printModel(int[] idx, MapFile mf) {
		StringBuffer sb = new StringBuffer();
		for(Entry<String, Suite> entry: this.entrySet()) {
			String key = entry.getKey();
			String[] g = key.split(PublicData.seperator);
			sb.append("Geno:");
			for(int i = 0; i < idx.length; i++) {
				SNP snp = mf.getSNP(idx[i]);
				sb.append(snp.getPolymorphism(g[i])+", ");
			}
			sb.append(" " + entry.getValue() + ", ");
		}
		return sb.toString();
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		for(Entry<String, Suite> entry: this.entrySet()) {
			String key = entry.getKey();
			sb.append("Genotype-" + key + " " + entry.getValue() + "; ");
		}
		return sb.toString();
	}
}
