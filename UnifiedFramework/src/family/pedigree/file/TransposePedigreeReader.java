package family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import admixture.parameter.Parameter;

import util.NewIt;
import family.mdr.arsenal.MDRConstant;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;

public class TransposePedigreeReader extends PedigreeFile {

	public String FamFile;
	public String tPedFile;
	private int n_individual = 0;
	private ArrayList<String> Famid;
	private ArrayList<String> Individualid;
	private MapFile mapData;

	public TransposePedigreeReader(String tped, String famF, MapFile md) {
		super();
		tPedFile = tped;
		FamFile = famF;
		mapData = md;
	}

	public void initial() throws IOException {
		Famid = NewIt.newArrayList();
		Individualid = NewIt.newArrayList();

		BufferedReader reader = new BufferedReader(new FileReader(new File(FamFile)));
		BufferedReader treader = new BufferedReader(new FileReader(new File(tPedFile)));

		String line = treader.readLine();
		while(line != null) {
			num_marker++;
			line = treader.readLine();
		}
		treader.close();
		AlleleSet = new char[num_marker][2];
		for(int i = 0; i < num_marker; i++) {
			AlleleSet[i][0] = AlleleSet[i][1] = Parameter.missing_allele.charAt(0);
		}
		AlleleFreq = new short[num_marker][2];

		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split(MDRConstant.delim);

			BPerson per = new BPerson(num_marker);
			Famid.add(tokens[0]);
			Individualid.add(tokens[1]);
			per.setFamilyID(tokens[0]);
			per.setPersonID(tokens[1]);
			per.setDadID(tokens[2]);
			per.setMomID(tokens[3]);
			per.setGender(Integer.parseInt(tokens[4]));
			per.setAffectedStatus(Integer.parseInt(tokens[5]));

			BFamilyStruct fam = familystructure.get(tokens[0]);
			if (fam == null) {
				fam = new BFamilyStruct(tokens[0]);
				familystructure.put(tokens[0], fam);
			}
			fam.addPerson(per);
			n_individual++;
		}
	}

	public void parseLinkage(String infile, int numMarker) throws IOException {
		initial();
		int colNum = Individualid.size() * 2 + 4;
		BufferedReader reader = new BufferedReader(new FileReader(new File(tPedFile)));
		String line;
		BPerson per;
		int k = 0;
		while ((line = reader.readLine()) != null) {
			String[] tokenizer = line.split(MDRConstant.delim);
			int numTokens = tokenizer.length;
			if (colNum != numTokens) {
				System.err.println("Column number mismatch in tpedfile in line " + (k + 1));
				System.exit(0);
			}
			mapData.addSNP(tokenizer[1], tokenizer[0], Float.parseFloat(tokenizer[2]), Integer.parseInt(tokenizer[3]));

			if (tokenizer.length > 4) {
				for (int i = 0; i < Famid.size(); i++) {
					BFamilyStruct bf = familystructure.get(Famid.get(i));
					per = bf.getPerson(Individualid.get(i));
					String[] allele = { tokenizer[4 + i * 2], tokenizer[4 + i * 2 + 1] };
					boolean flag = (allele[0].compareTo(Parameter.missing_allele) != 0) && (allele[1].compareTo(Parameter.missing_allele) != 0);
					if (flag) {
						int[] code = recode(k, allele);
						per.addMarker(flag, code[0], code[1], k);
						AlleleFreq[k][code[0]]++; AlleleFreq[k][code[1]]++;
					} else {
						per.addMarker(flag, 0, 0, k);
					}
				}
			}
			k++;
		}
		reader.close();
	}
}
