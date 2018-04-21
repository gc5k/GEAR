package gear.family.pedigree.file;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;

import gear.data.Person;
import gear.data.Family;
import gear.family.pedigree.Hukou;
import gear.family.plink.PLINKBinaryParser;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

public class BEDReader extends PedigreeFile {
	public String famFile;
	private int n_individual = 0;
	private ArrayList<String> famIDs;
	private ArrayList<Person> persons;
	private MapFile mapData;

	public BEDReader(String famF, int numMark, MapFile mapdata) {
		super();
		famFile = famF;
		num_marker = numMark;
		mapData = mapdata;
	}

	@Override
	public void initial() {
		famIDs = NewIt.newArrayList();
		persons = NewIt.newArrayList();

		BufferedReader reader = BufferedReader.openTextFile(famFile, "fam");
		AlleleSet = new char[num_marker][];
		AlleleFreq = new short[num_marker][2];
		for (int i = 0; i < mapData.snpList.size(); i++) {
			SNP snp = mapData.snpList.get(i);
			AlleleSet[i] = snp.getSNP();
		}

		HukouBook = NewIt.newArrayList();
		Hukou hukou;
		String[] tokens;
		while ((tokens = reader.readTokens(6)) != null) {
			Person person = new Person(num_marker);
			famIDs.add(tokens[0]);

			person.setFamilyID(tokens[0]);
			person.setPersonID(tokens[1]);
			person.setDadID(tokens[2]);
			person.setMomID(tokens[3]);
			person.setGender(Integer.parseInt(tokens[4]));
			person.setAffectedStatus(tokens[5]);
			SixthCol.add(tokens[5]);

			hukou = new Hukou(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]);
			Family family = families.get(tokens[0]);
			if (family == null) {
				family = new Family(tokens[0]);
				families.put(family);
			}
			if (family.hasPerson(person.getPersonID())) {
				String msg = "";
				msg += "Person " + person.getPersonID() + " in family ";
				msg += person.getFamilyID() + " appears more than once.";
				reader.errorPreviousLine(msg);
			}
			HukouBook.add(hukou);
			family.addPerson(person);
			persons.add(person);
			n_individual++;
		}
		Is6ColBinary();
		Logger.printUserLog("Reading " + HukouBook.size() + " individuals from '" + famFile + "'.");
	}

	@Override
	public void parseLinkage(String infile, int numMarkerInFile, int[] WSNP) {
		initial();
		pedfile = infile;
		BufferedInputStream in = null;

		try {
			in = new BufferedInputStream(new FileInputStream(new File(pedfile)));
		} catch (FileNotFoundException e) {
			Logger.handleException(e, "Cannot open the pedigree file '" + pedfile + "'.");
		}

		try {
			byte[] magic = new byte[3];
			in.read(magic, 0, 3);
			if (magic[2] == 1) {
				Logger.printUserLog("Reading data in PLINK SNP-major mode.");
				snp_major(in, numMarkerInFile, WSNP);
			} else {
				Logger.printUserLog("Reading data in PLINK individual-major mode.");
				individual_major(in, numMarkerInFile, WSNP);
			}
			in.close();
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when reading the bed file.");
		}
	}

	private void individual_major(BufferedInputStream in, int numMarkerInFile, int[] WSNP) throws IOException {
		int L = 0;
		if (numMarkerInFile % 4 == 0) {
			L = numMarkerInFile / 4;
		} else {
			L = numMarkerInFile / 4 + 1;
		}
		int exL = 0;
		if (WSNP.length % 4 == 0) {
			exL = WSNP.length / 4;
		} else {
			exL = WSNP.length / 4 + 1;
		}
		byte[] geno = new byte[L];
		byte[] extract_geno = new byte[exL];
		for (int i = 0; i < n_individual; i++) {
			in.read(geno, 0, L);
			extract_geno = extractGenotype(geno, numMarkerInFile, WSNP);
			persons.get(i).addAllMarker(extract_geno);
		}
		famIDs = null;
		persons = null;
	}

	private byte[] extractGenotype(byte[] g, int numMarkerInFile, int[] WSNP) {

		int exL = 0;
		if (WSNP.length % 4 == 0) {
			exL = WSNP.length / 4;
		} else {
			exL = WSNP.length / 4 + 1;
		}
		byte[] Exg = new byte[exL];
		int c = 0;
		for (int i = 0; i < numMarkerInFile; i++) {
			int idx = ArrayUtils.indexOf(WSNP, i);
			if (idx < 0)
				continue;
			int posByte = i >> 2;
			int posBit = (i & 0x3) << 1;
			int g1 = (g[posByte] >> posBit) & 3;

			int ExposByte = c >> 2;
			int ExposBit = (c & 0x3) << 1;
			Exg[ExposByte] |= g1 << ExposBit;
			if (g1 == 0) {
				AlleleFreq[c][0] += 2;
			} else if (g1 == 2) {
				AlleleFreq[c][0]++;
				AlleleFreq[c][1]++;
			} else if (g1 == 3) {
				AlleleFreq[c][1] += 2;
			}
			c++;
		}
		return Exg;
	}

	private static int[][] constructSnpMajorGenotypeByteConvertTable() {
		int[][] table = new int[0x100][4];
		for (int byteValue = 0; byteValue <= 0xff; ++byteValue) {
			for (int indIdx = 0; indIdx < 4; ++indIdx) {
				table[byteValue][indIdx] = PLINKBinaryParser.convertToGearGenotype((byteValue >> (indIdx << 1)) & 0x3);
			}
		}
		return table;
	}

	private void snp_major(BufferedInputStream in, int numMarkerInFile, int[] WSNP) throws IOException {
		byte[] g = new byte[(n_individual + 3) / 4];
		int[][] genoByteCvtTable = constructSnpMajorGenotypeByteConvertTable();
		int snpIdx = 0;
		for (int i = 0; i < numMarkerInFile; i++) {
			in.read(g, 0, g.length);
			if (ArrayUtils.indexOf(WSNP, i) >= 0) {
				int indIdx = 0;
				int posByte = snpIdx >> Person.shift;
				int posBit = (i & 0xf) << 1;
				for (int byteIdx = 0; byteIdx < g.length; ++byteIdx) {
					int[] genoValues = genoByteCvtTable[g[byteIdx] & 0xff]; // 0xff
																			// is
																			// necessary
																			// here,
																			// otherwise
																			// Java
																			// will
																			// sign
																			// extend
																			// the
																			// byte
					for (int j = 0; j < 4 && indIdx < n_individual; ++j, ++indIdx) {
						persons.get(indIdx).addByteGenotype(genoValues[j], posByte, posBit);
					}
				}
				snpIdx++;
			}
		}
	}
}
