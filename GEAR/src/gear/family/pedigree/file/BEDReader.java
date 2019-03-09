package gear.family.pedigree.file;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

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
	private ArrayList<String> famIDs;
	private ArrayList<Person> persons;
	private MapFile mapData;
	private BufferedInputStream inStream;
	private boolean isSnpMajor;
	private long bytesPerRow;

	public BEDReader(String bedFilename, String famFilename, MapFile mapData) {
		super(bedFilename);
		famFile = famFilename;
		this.mapData = mapData;
		try {
			inStream = new BufferedInputStream(new FileInputStream(new File(bedFilename)));
			byte[] magic = new byte[3];
			inStream.read(magic, 0, 3);
			isSnpMajor = magic[2] == 1;
		} catch (FileNotFoundException e) {
			Logger.handleException(e, "Cannot open the bed file '" + bedFilename + "'.");
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when reading the bed file.");
		}
	}
	
	public boolean IsSnpMajor() { return isSnpMajor; }
	
	private void parseFam() {
		famIDs = NewIt.newArrayList();
		persons = NewIt.newArrayList();
		numMarkers = mapData.getMarkerNumberOriginal();
		BufferedReader reader = BufferedReader.openTextFile(famFile, "fam");
		alleleSet = new char[numMarkers][];
		alleleFreq = new short[numMarkers][2];
		for (int i = 0; i < mapData.snpList.size(); i++) {
			SNP snp = mapData.snpList.get(i);
			alleleSet[i] = snp.getSNP();
		}

		hukouBook = NewIt.newArrayList();
		Hukou hukou;
		String[] tokens;
		while ((tokens = reader.readTokens(6)) != null) {
			Person person = new Person(numMarkers);
			famIDs.add(tokens[0]);

			person.setFamilyID(tokens[0]);
			person.setPersonID(tokens[1]);
			person.setDadID(tokens[2]);
			person.setMomID(tokens[3]);
			person.setGender(Integer.parseInt(tokens[4]));
			person.setAffectedStatus(tokens[5]);
			sixthCol.add(tokens[5]);

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
			hukouBook.add(hukou);
			family.addPerson(person);
			persons.add(person);
		}
		checkIs6ColBinary();
		Logger.printUserLog("Read " + hukouBook.size() + " individuals from '" + famFile + "'.");
	}
	
	private void calculateBytesPerRow() {
		if (IsSnpMajor()) {
			bytesPerRow = (persons.size() + 3) / 4;
		} else {
			bytesPerRow = (mapData.getMarkerNumberOriginal()) / 4;
		}
	}
	
	public void prepareToParseGenotypes() {
		parseFam();
		calculateBytesPerRow();
	}
	
	public long getBytesPerRow() { return bytesPerRow; }
	
	public int readNextByte() {
		int nextByte = -1;
		try {
			nextByte = inStream.read();
		} catch (IOException e) {
			Logger.handleException(e, "I/O exception occurred when reading the bed file.");
		}
		if (nextByte == -1) {
			Logger.printUserError("Read past the end of the bed file.");
			System.exit(-1);
		}
		return nextByte;
	}

	@Override
	public void parseLinkage(int[] WSNP) {
		long startTime = System.nanoTime();
		try {
			if (isSnpMajor) {
				Logger.printUserLog("Reading data in PLINK SNP-major mode.");
				parseWithSnpMajor(inStream, mapData.getMarkerNumberOriginal(), WSNP);
			} else {
				Logger.printUserLog("Reading data in PLINK individual-major mode.");
				parseWithIndividualMajor(inStream, mapData.getMarkerNumberOriginal(), WSNP);
			}
			inStream.close();
			long endTime = System.nanoTime();
			Logger.printUserLog(String.format("It takes %.1fs to read the data.", (endTime - startTime) / 1e9));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when reading the bed file.");
		}
	}

	private void parseWithIndividualMajor(BufferedInputStream in, int numMarkersInFile, int[] workingSNPs) throws IOException {
		int L = 0;
		if (numMarkersInFile % 4 == 0) {
			L = numMarkersInFile / 4;
		} else {
			L = numMarkersInFile / 4 + 1;
		}
		int exL = 0;
		int numMarkers = workingSNPs == null ? numMarkersInFile : workingSNPs.length;
		if (numMarkers % 4 == 0) {
			exL = numMarkers / 4;
		} else {
			exL = numMarkers / 4 + 1;
		}
		byte[] geno = new byte[L];
		byte[] extract_geno = new byte[exL];
		for (int i = 0; i < persons.size(); i++) {
			in.read(geno, 0, L);
			extract_geno = extractGenotype(geno, numMarkersInFile, workingSNPs);
			persons.get(i).addAllMarker(extract_geno);
		}
		famIDs = null;
		persons = null;
	}

	private byte[] extractGenotype(byte[] g, int numMarkersInFile, int[] workingSNPs) {
		int exL = 0;
		int numMarkers = workingSNPs == null ? numMarkersInFile : workingSNPs.length;
		if (numMarkers % 4 == 0) {
			exL = numMarkers / 4;
		} else {
			exL = numMarkers / 4 + 1;
		}
		byte[] Exg = new byte[exL];
		int c = 0;
		for (int i = 0; i < numMarkersInFile; i++) {
			int idx = workingSNPs == null ? i : ArrayUtils.indexOf(workingSNPs, i);
			if (idx < 0)
				continue;
			int posByte = i >> 2;
			int posBit = (i & 0x3) << 1;
			int g1 = (g[posByte] >> posBit) & 3;

			int ExposByte = c >> 2;
			int ExposBit = (c & 0x3) << 1;
			Exg[ExposByte] |= g1 << ExposBit;
			if (g1 == 0) {
				alleleFreq[c][0] += 2;
			} else if (g1 == 2) {
				alleleFreq[c][0]++;
				alleleFreq[c][1]++;
			} else if (g1 == 3) {
				alleleFreq[c][1] += 2;
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

	private void parseWithSnpMajor(BufferedInputStream in, int numMarkersInFile, int[] wsnp) throws IOException {
		byte[] g = new byte[(persons.size() + 3) / 4];
		int[][] genoByteCvtTable = constructSnpMajorGenotypeByteConvertTable();
		
		boolean[] workingSnpFlags = new boolean[numMarkersInFile];
		if (wsnp == null) {
			Arrays.fill(workingSnpFlags, true);
		} else {
			Arrays.fill(workingSnpFlags, false);
			for (int idx : wsnp)
				workingSnpFlags[idx] = true;
		}
		
		int snpIdx = 0;
		for (int i = 0; i < numMarkersInFile; i++) {
			in.read(g, 0, g.length);
			if (workingSnpFlags[i]) {
				int indIdx = 0;
				int posByte = snpIdx >> Person.shift;
				int posBit = (i & 0xf) << 1;
				for (int byteIdx = 0; byteIdx < g.length; ++byteIdx) {
					// 0xff is necessary here, otherwise Java will sign extend the byte
					int[] genoValues = genoByteCvtTable[g[byteIdx] & 0xff];
					for (int j = 0; j < 4 && indIdx < persons.size(); ++j, ++indIdx) {
						persons.get(indIdx).addByteGenotype(genoValues[j], posByte, posBit);
					}
				}
				snpIdx++;
			}
		}
	}
	
	public void skipOneRow() {
		try {
			inStream.skip(bytesPerRow);
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when reading the bed file.");
		}
	}
}
