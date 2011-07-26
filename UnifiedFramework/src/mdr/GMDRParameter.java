package mdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class GMDRParameter {

	protected ArrayList<String> lines;
	protected BufferedReader buffer;
	protected String filename;
	protected String marker_file;
	protected String phe_file;
	protected int interaction_from;
	protected int interaction_end;
	protected int[] scr_idx;
	protected int interval;
	protected long seed;
	protected int partition_method;
	protected boolean isMooreMDR;
	protected int replication_permutation;
	protected int search_method;

	public GMDRParameter() {

	}

	public void read(String file) throws IOException {
		filename = file;
		lines = new ArrayList();
		buffer = new BufferedReader(new FileReader(new File(file)));
		sweepComments();
		parseValue();
	}

	public void sweepComments() throws IOException {
		boolean flag = true;
		String line;
		while ((line = buffer.readLine()) != null) {
			if (Pattern.matches("^\\s*//*.*", line)) {// empty line
				continue;
			} else {
				lines.add(line);
			}
		}
	}

	public void parseValue() {
		marker_file = lines.get(0);
		phe_file = lines.get(1);

		interaction_from = Integer.parseInt(lines.get(2));
		interaction_end = Integer.parseInt(lines.get(3));
		String[] c = lines.get(4).split(",");
		scr_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			scr_idx[i] = Integer.parseInt(c[i]) - 1;
		}
		interval = Integer.parseInt(lines.get(5));
		seed = Long.parseLong(lines.get(6));
		partition_method = Integer.parseInt(lines.get(7));
		isMooreMDR = Boolean.parseBoolean(lines.get(8));
		replication_permutation = Integer.parseInt(lines.get(9));
		search_method = Integer.parseInt(lines.get(10));
	}
	
	public void setMarkerFile(String mf) {
		marker_file = mf;
	}
	
	public String getMarkerFile() {
		return marker_file;
	}
	
	public void setPhenotypeFile(String pf) {
		phe_file = pf;
	}
	
	public String getPhenotypeFile() {
		return phe_file;
	}
	
	public void setInteractionFrom(int i) {
		interaction_from = i;
	}
	
	public int getInterctionFrom() {
		return interaction_from;
	}
	
	public void setInteractionEnd(int i) {
		interaction_end = i;
	}
	
	public int getInteractionEnd() {
		return interaction_end;
	}
	
	public void setScoreIndex(String sc) {
		String[] c = sc.split(",");
		scr_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			scr_idx[i] = Integer.parseInt(c[i]) - 1;
		}
	}

	public void setScoreIndex(int[] scrIdx) {		
		scr_idx = scrIdx;
	}
	
	public int[] getScoreIndex() {
		return scr_idx;
	}
	
	public void setInterval(int i) {
		interval = i;
	}
	
	public int getInterval() {
		return interval;
	}
	
	public void setSeed(long s) {
		seed = s;
	}
	
	public long getSeed() {
		return seed;
	}
	
	public void setPartitionMethod(int i) {
		partition_method = i;
	}
	
	public int getPartitionMethod() {
		return partition_method;
	}
	
	public void setMooreMDR(boolean m) {
		isMooreMDR = m;
	}
	
	public boolean isMooreMDR() {
		return isMooreMDR;
	}
	
	public void setReplicationPermutation(int i) {
		replication_permutation = i;
	}
	
	public int getReplicationPermutation() {
		return replication_permutation;
	}
	
	public void setSearchMethod(int i) {
		search_method = i;
	}
	
	public int getSearchMethod() {
		return search_method;
	}
}
