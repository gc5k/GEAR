package family.pedigree.file;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;

public class BEDReader extends PedigreeFile {
	public String FamFile;
	public BEDReader (String famF, int numMark) {
		super();
		FamFile = famF;
		this.num_marker = numMark;
	}

	public void intial(File infile) throws IOException {
		
		BufferedReader reader = new BufferedReader(new FileReader(FamFile));
		String line;
		while ((line = reader.readLine()) != null) {
			String[] tokens = line.split("\\s+");

			BFamilyStruct fam = familystructure.get(tokens[0]);
			BPerson per = new BPerson(num_marker);
			per.setFamilyID(tokens[0]);
			per.setPersonID(tokens[1]);
			per.setDadID(tokens[2]);
			per.setMomID(tokens[3]);
			per.setGender(Integer.parseInt(tokens[4]));
			per.setAffectedStatus(Integer.parseInt(tokens[5]));

			if(fam == null) {
				fam = new BFamilyStruct(tokens[0]);
				familystructure.put(tokens[0], fam);
			}
			fam.addPerson(per);	
		}
		
	}
	
	public void parseLinkage() throws IOException {
		BufferedInputStream in = null;
		try {
			in = new BufferedInputStream(new FileInputStream(pedfile));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		byte[] magic = new byte[3];
		int n = in.read(magic, 0, 3);
		if(magic[2] == 1) {
			snp_major(in);
		} else {
			individual_major(in);
		}
	}
	
	private void snp_major(BufferedInputStream in) {
		
	}
	
	private void individual_major(BufferedInputStream in) {
		
	}

}
