package gear.family.pedigree.file;

import java.util.ArrayList;

import gear.ConstValues;
import gear.data.Person;
import gear.data.Family;
import gear.util.BufferedReader;
import gear.util.NewIt;

public class TransposePedigreeReader extends PedigreeFile
{
	public String famFile;
	private ArrayList<String> familyIDs;
	private ArrayList<String> Individualid;
	private MapFile mapData;

	public TransposePedigreeReader(String tped, String famF, MapFile md)
	{
		super(tped);
		famFile = famF;
		mapData = md;
	}

	public void initial()
	{
		familyIDs = NewIt.newArrayList();
		Individualid = NewIt.newArrayList();

		BufferedReader famReader = BufferedReader.openTextFile(famFile, "fam");
		BufferedReader tpedReader = BufferedReader.openTextFile(filename, "tped");

		String line = tpedReader.readLine();
		while (line != null)
		{
			numMarkers++;
			line = tpedReader.readLine();
		}
		tpedReader.close();
		alleleSet = new char[numMarkers][2];
		for (int i = 0; i < numMarkers; i++)
		{
			alleleSet[i][0] = alleleSet[i][1] = ConstValues.MISSING_ALLELE_CHAR;
		}
		alleleFreq = new short[numMarkers][2];
		
		String[] tokens;
		while ((tokens = famReader.readTokens(6)) != null)
		{
			Person per = new Person(numMarkers);
			familyIDs.add(tokens[0]);
			Individualid.add(tokens[1]);
			per.setFamilyID(tokens[0]);
			per.setPersonID(tokens[1]);
			per.setDadID(tokens[2]);
			per.setMomID(tokens[3]);
			per.setGender(Integer.parseInt(tokens[4]));
			per.setAffectedStatus(tokens[5]);

			Family fam = families.get(tokens[0]);
			if (fam == null)
			{
				fam = new Family(tokens[0]);
				families.put(fam);
			}
			fam.addPerson(per);
		}
	}

	public void parseLinkage(int numMarker)
	{
		initial();
		int numCols = Individualid.size() * 2 + 4;
		BufferedReader reader = BufferedReader.openTextFile(filename, "tped");
		Person per;
		int k = 0;
		String[] tokens;
		while ((tokens = reader.readTokens(numCols)) != null)
		{
			mapData.addSNP(tokens[1], tokens[0], Float.parseFloat(tokens[2]), Integer.parseInt(tokens[3]));

			if (tokens.length > 4)
			{
				for (int i = 0; i < familyIDs.size(); i++)
				{
					Family bf = families.get(familyIDs.get(i));
					per = bf.getPerson(Individualid.get(i));
					String[] allele = { tokens[4 + i * 2], tokens[4 + i * 2 + 1] };
					boolean flag = !allele[0].equals(ConstValues.MISSING_ALLELE_STRING)	&&
							       !allele[1].equals(ConstValues.MISSING_ALLELE_STRING);
					if (flag)
					{
						int[] code = recode(k, allele);
						per.addMarker(flag, code[0], code[1], k);
						alleleFreq[k][code[0]]++;
						alleleFreq[k][code[1]]++;
					}
					else
					{
						per.addMarker(flag, 0, 0, k);
					}
				}
			}
			k++;
		}
		reader.close();
	}
}
