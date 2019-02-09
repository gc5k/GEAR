package gear.family.pedigree.file;

import gear.util.Logger;
import gear.util.NewIt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class MapFile
{
	protected static final String DELIMITER = "\\s+";
	protected ArrayList<SNP> snpList = NewIt.newArrayList();
	protected HashMap<String, Integer> chrSNPCount = NewIt.newHashMap();
	protected ArrayList<Integer> badline;
	protected int[] workingSnpIndexes;
	protected int numMarkerOriginal;
	protected String filename = null;

	public MapFile(String filename) {
		this.filename = filename;
	}
	
	public String getFilename() {
		return filename;
	}

	public void setWorkingSNPs(int[] workingSnpIndexes) {
		this.workingSnpIndexes = workingSnpIndexes;
		ArrayList<SNP> filteredSNPs = NewIt.newArrayList();
		filteredSNPs.ensureCapacity(workingSnpIndexes.length);
		for (int i = 0; i < workingSnpIndexes.length; i++)
		{
			SNP snp = snpList.get(workingSnpIndexes[i]);
			filteredSNPs.add(snp);
		}
		snpList = filteredSNPs;
	}

	public void parseMap() {
		File mapfile = new File(filename);

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(mapfile));
		} catch (IOException e) {
			Logger.handleException(e, "Cannot open the map file '" + mapfile + "'.");
		}

		String line = null;
		try
		{
			int c = 0;
			while ((line = reader.readLine()) != null)
			{
				c++;
				String[] tokens = line.split(DELIMITER);
				if (tokens.length < 4)
				{
					if (badline == null)
						badline = NewIt.newArrayList();
					badline.add(new Integer(c));
				}
				String chr = tokens[0];
				String name = tokens[1];
				// Skip genetic distance field at tokens[2].
				float dis = Float.parseFloat(tokens[2]);
				int pos = Integer.parseInt(tokens[3]);
				addSNP(chr, name, dis, pos);
			}
			reader.close();
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occured when reading the map file '" + filename + "'.");
		}

		if (badline != null)
		{
			Logger.printUserError("The following line(s) have error:");
			String badlines = "";
			for (Integer i : badline)
			{
				badlines += i + ",";
			}
			Logger.printUserError(badlines);
			System.exit(1);
		}
		numMarkerOriginal = snpList.size();
	}

	public HashMap<String, Integer> getChrSNPCount()
	{
		return chrSNPCount;
	}

	private void count(String chr)
	{
		if (chrSNPCount.containsKey(chr))
		{
			Integer c = chrSNPCount.get(chr);
			c++;
			chrSNPCount.put(chr, c);
		} else
		{
			chrSNPCount.put(chr, new Integer(1));
		}
	}

	public void addSNP(String chr, String name, float dis, int pos) {
		count(chr);
		snpList.add(new SNP(chr, name, dis, pos));
	}

	public void addSNP(String chr, String name, float dis, int pos, char a1, char a2) {
		count(chr);
		snpList.add(new SNP(chr, name, dis, pos, a1, a2));
	}

	public void setMarker(int l)
	{
		snpList.ensureCapacity(l);
		for (int i = 0; i < l; i++)
		{
			StringBuilder sb = new StringBuilder();
			sb.append("snp");
			sb.append(i);
			addSNP(sb.toString(), "-1", -1f, -1);
		}
	}

	public ArrayList<SNP> getMarkerList()
	{
		return snpList;
	}

	public String getMarkerName(int i)
	{
		return snpList.get(i).getName();
	}

	public int getMarkerNumberOriginal()
	{
		return numMarkerOriginal;
	}

	public int getMarkerNumber()
	{
		return snpList.size();
	}

	public SNP getSNP(int i)
	{
		return snpList.get(i);
	}

	public void setPolymorphism(char[][] p, short[][] freq)
	{
		if (p.length != snpList.size())
		{
			Logger.printUserError("The map file and the pedigree file do not match.");
			Logger.getDevLogger().severe("p.length != snpList.size()");
			System.exit(1);
		} else
		{
			for (int i = 0; i < p.length; i++)
			{
				SNP snp = snpList.get(i);
				snp.setAllele(p[i], freq[i]);
			}
		}
	}

	public void setAlleleFrequency(double[][] freq)
	{
		for (int i = 0; i < freq.length; i++)
		{
			SNP snp = snpList.get(i);
			snp.setAllele(freq[i]);
		}
	}

	public void setPolymorphismMarker(char[][] p)
	{
		for (int i = 0; i < p.length; i++)
		{
			SNP snp = snpList.get(i);
			snp.setAllelePolymorphism(p[i]);
		}
	}
}
