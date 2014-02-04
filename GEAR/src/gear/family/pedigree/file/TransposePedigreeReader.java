package gear.family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import gear.ConstValues;
import gear.family.pedigree.genotype.BFamilyStruct;
import gear.family.pedigree.genotype.BPerson;
import gear.util.Logger;
import gear.util.NewIt;

public class TransposePedigreeReader extends PedigreeFile
{

	public String FamFile;
	public String tPedFile;
	private ArrayList<String> Famid;
	private ArrayList<String> Individualid;
	private MapFile mapData;

	public TransposePedigreeReader(String tped, String famF, MapFile md)
	{
		super();
		tPedFile = tped;
		FamFile = famF;
		mapData = md;
	}

	public void initial() throws IOException
	{
		Famid = NewIt.newArrayList();
		Individualid = NewIt.newArrayList();

		BufferedReader reader = new BufferedReader(new FileReader(new File(
				FamFile)));
		BufferedReader treader = new BufferedReader(new FileReader(new File(
				tPedFile)));

		String line = treader.readLine();
		while (line != null)
		{
			num_marker++;
			line = treader.readLine();
		}
		treader.close();
		AlleleSet = new char[num_marker][2];
		for (int i = 0; i < num_marker; i++)
		{
			AlleleSet[i][0] = AlleleSet[i][1] = ConstValues.MISSING_ALLELE_CHAR;
		}
		AlleleFreq = new short[num_marker][2];

		while ((line = reader.readLine()) != null)
		{
			String[] tokens = line.split("\\s+");

			BPerson per = new BPerson(num_marker);
			Famid.add(tokens[0]);
			Individualid.add(tokens[1]);
			per.setFamilyID(tokens[0]);
			per.setPersonID(tokens[1]);
			per.setDadID(tokens[2]);
			per.setMomID(tokens[3]);
			per.setGender(Integer.parseInt(tokens[4]));
			per.setAffectedStatus(tokens[5]);

			BFamilyStruct fam = familySet.getFamily(tokens[0]);
			if (fam == null)
			{
				fam = new BFamilyStruct(tokens[0]);
				familySet.putFamily(fam);
			}
			fam.addPerson(per);
		}
	}

	public void parseLinkage(String infile, int numMarker) throws IOException
	{
		initial();
		int colNum = Individualid.size() * 2 + 4;
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				tPedFile)));
		String line;
		BPerson per;
		int k = 0;
		while ((line = reader.readLine()) != null)
		{
			String[] tokenizer = line.split("\\s+");
			int numTokens = tokenizer.length;
			if (colNum != numTokens)
			{
				Logger.printUserError("Column number mismatch in tpedfile in line "
						+ (k + 1)
						+ ", expect "
						+ colNum
						+ " column(s), but read " + numTokens);
				System.exit(1);
			}
			mapData.addSNP(tokenizer[1], tokenizer[0],
					Float.parseFloat(tokenizer[2]),
					Integer.parseInt(tokenizer[3]));

			if (tokenizer.length > 4)
			{
				for (int i = 0; i < Famid.size(); i++)
				{
					BFamilyStruct bf = familySet.getFamily(Famid.get(i));
					per = bf.getPerson(Individualid.get(i));
					String[] allele = { tokenizer[4 + i * 2],
							tokenizer[4 + i * 2 + 1] };
					boolean flag = !allele[0].equals(ConstValues.MISSING_ALLELE_STRING)	&&
							       !allele[1].equals(ConstValues.MISSING_ALLELE_STRING);
					if (flag)
					{
						int[] code = recode(k, allele);
						per.addMarker(flag, code[0], code[1], k);
						AlleleFreq[k][code[0]]++;
						AlleleFreq[k][code[1]]++;
					} else
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
