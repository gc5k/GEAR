package gear.family.pedigree.file;

import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.lang3.ArrayUtils;

import gear.ConstValues;
import gear.CmdArgs;
import gear.data.Person;
import gear.data.Family;
import gear.data.UniqueRecordSet;
import gear.family.pedigree.Hukou;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

/**
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PedigreeFile
{

	protected char[][] alleleSet;
	protected short[][] alleleFreq;
	protected ArrayList<Hukou> hukouBook;
	protected UniqueRecordSet<Family> families = new UniqueRecordSet<Family>();
	protected HashSet<String> sixthCol = NewIt.newHashSet();
	protected boolean isSixthColBinary = true;
	protected int numMarkers;
	protected String titleLine = null;
	protected String filename;
	protected boolean header = true;
	
	public PedigreeFile(String filename) {
		this.filename = filename;
	}
	
	public String getFilename() {
		return filename;
	}

	public Family getFamily(String familystrID)
	{
		return this.families.get(familystrID);
	}

	/**
	 * this method iterates through each family in Hashtable families and adds
	 * up the number of individuals in total across all families
	 * 
	 * @return the total number of individuals in all the family objects in the
	 *         families hashtable
	 */
	public int getNumIndividuals()
	{
		int total = 0;
		for (int ii = 0; ii < families.size(); ++ii)
		{
			Family fam = families.get(ii);
			total += fam.size();
		}
		return total;
	}

	public UniqueRecordSet<Family> getFamilies()
	{
		return families;
	}

	/**
	 * Taking in a pedigree file in the form of a vector of strings and parses
	 * it. The data parsed is stored in families in the member hashtable
	 * families. Note that the "Linkage" here is the relationship between
	 * relatives in a pedigree, but is not that term of genetics.
	 */
	public void parseLinkage(int[] workingSnpIndexes)
	{
		BufferedReader reader = BufferedReader.openTextFile(filename, "ped");
		Person per;
		int k = 0;
		Hukou hukou;
		hukouBook = NewIt.newArrayList();
		String[] tokens = reader.readTokensAtLeast(6);
		
		if (tokens.length % 2 != 0)
		{
			reader.errorPreviousLine("The number of columns is not an even number.");
		}
		
		int numCols = tokens.length;
		int numMarkers = (numCols - 6) / 2;

		this.numMarkers = workingSnpIndexes == null ? numMarkers : workingSnpIndexes.length;
		alleleSet = new char[this.numMarkers][2];
		alleleFreq = new short[this.numMarkers][2];
		for (int i = 0; i < this.numMarkers; i++)
		{
			alleleSet[i][0] = alleleSet[i][1] = ConstValues.MISSING_ALLELE_CHAR;
		}
		
		do
		{
			if (numCols != tokens.length)
			{
				String msg = "";
				msg += "The first row has " + numCols + " columns, ";
				msg += "but this row has " + tokens.length + " columns.";
				reader.errorPreviousLine(msg);
			}

			per = new Person(numMarkers);
			per.setFamilyID(tokens[0]);
			per.setPersonID(tokens[1]);
			per.setDadID(tokens[2]);
			per.setMomID(tokens[3]);

			int Gender = Integer.parseInt(tokens[4]);
			sixthCol.add(tokens[5]);

			per.setGender(Gender);
			per.setAffectedStatus(tokens[5]);

			hukou = new Hukou(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]);

			int c = 0;
			for (int j = 0; j < numMarkers; j++)
			{
				int idx = workingSnpIndexes == null ? j : ArrayUtils.indexOf(workingSnpIndexes, j);
				if (idx < 0)
				{
					continue;
				}
				
				try
				{
					String[] allele = { tokens[6 + j * 2],
							tokens[6 + j * 2 + 1] };
					boolean flag = !allele[0].equals(ConstValues.MISSING_ALLELE_STRING) &&
							       !allele[1].equals(ConstValues.MISSING_ALLELE_STRING);
					if (flag)
					{
						int[] code = recode(c, allele);
						per.addMarker(flag, code[0], code[1], c);
						alleleFreq[c][code[0]]++;
						alleleFreq[c][code[1]]++;
					}
					else
					{
						per.addMarker(flag, 0, 0, c);
					}
				}
				catch (NumberFormatException e)
				{
					Logger.handleException(e,
							"An invalid genotype is found in the ped file at line "
									+ (k + 1) + " for marker " + (c + 1)
									+ ".");
				}
				c++;
			}
			// check if the family exists already in the Hashtable
			Family family = families.get(per.getFamilyID());
			if (family == null)
			{
				// it doesn't exist, so create a new FamilyStruct object
				family = new Family(per.getFamilyID());
				families.put(family);
			}

			if (family.hasPerson(per.getPersonID()))
			{
				String msg = "";
				msg += "Person " + per.getPersonID() + " in family ";
				msg += per.getFamilyID() + " appears more than once.";
				reader.errorPreviousLine(msg);
			}
			hukouBook.add(hukou);
			family.addPerson(per);
			k++;
		} while ((tokens = reader.readTokens(numCols)) != null);
		checkIs6ColBinary();
		Logger.printUserLog("Read " + hukouBook.size() + " individuals from '" + filename + "'.");

	}

	protected void checkIs6ColBinary()
	{
		for (String c : sixthCol) {
			if (CmdArgs.INSTANCE.status_shiftFlag) {
				if (c.compareTo("1") != 0
						&& c.compareTo("0") != 0
						&& c.compareTo(ConstValues.MISSING_PHENOTYPE) != 0) {
					isSixthColBinary = false;
					break;
				}
			} else {
				if (c.compareTo("2") != 0
						&& c.compareTo("1") != 0
						&& c.compareTo("0") != 0
						&& c.compareTo(ConstValues.MISSING_PHENOTYPE) != 0) {
					isSixthColBinary = false;
					break;
				}
			}
		}
	}

	public boolean IsSixthColBinary()
	{
		return isSixthColBinary;
	}

	protected int[] recode(int idx, String[] allele) {
		int[] code = { -1, -1 };
		char[] ref = alleleSet[idx];
		if (ref[1] != ConstValues.MISSING_ALLELE_CHAR)
		{
			// two detected alleles
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (allele[i].charAt(0) == ref[j])
					{
						code[i] = j;
						break;
					}
				}
			}
			if (code[0] == -1 || code[1] == -1)
			{
				Logger.printUserError("There're more than 2 alleles in the marker column "
						+ (idx + 1));
			}
		} else {
			// less than two detected alleles
			if (allele[0].compareTo(allele[1]) == 0)
			{
				// 1 when both alleles are same
				if (ref[0] == ConstValues.MISSING_ALLELE_CHAR)
				{
					// zero detected alleles
					ref[0] = allele[0].charAt(0);
					code[0] = code[1] = 0;
				}
				else
				{
					// one detected alleles
					if (allele[0].charAt(0) == ref[0])
					{
						code[0] = code[1] = 0;
					} else
					{
						code[0] = code[1] = 1;
						ref[1] = allele[1].charAt(0);
					}
				}
			} else
			{
				// 2 when both alleles are different
				if (ref[0] == ConstValues.MISSING_ALLELE_CHAR)
				{
					// zero detected alleles
					ref[0] = allele[0].charAt(0);
					ref[1] = allele[1].charAt(0);
					code[0] = 0;
					code[1] = 1;
				}
				else
				{
					// one detected alleles
					if (ref[0] == allele[0].charAt(0))
					{
						ref[1] = allele[1].charAt(0);
						code[0] = 0;
						code[1] = 1;
					}
					else if (ref[0] == allele[1].charAt(0))
					{
						ref[1] = allele[0].charAt(0);
						code[0] = 1;
						code[1] = 0;
					}
					else
					{
						Logger.printUserError("There're more than 3 alleles in marker column "
								+ (idx + 1));
					}
				}
			}
		}
		return code;
	}

	public char[][] getPolymorphism() { return alleleSet; }

	public short[][] getAlleleFrequency() { return alleleFreq; }

	public int getNumMarkers() { return numMarkers; }

	public void setHeader(boolean flag) { header = flag; }

	public ArrayList<Hukou> getHukouBook() { return hukouBook; }
}