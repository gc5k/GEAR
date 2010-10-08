package family.pedigree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;

public class GMDRPhenoFile {

    private String titleLine;
    private ArrayList traitInfor;
    private ArrayList traitFileStrings;
    private ArrayList traitStrings;
    private File phenoFile;
    private ArrayList allSubjects;
    private Hashtable families;
    private boolean bogusParents = false;

    public GMDRPhenoFile() {
        families = new Hashtable();
    }

    public void Initial(File infile) throws GMDRPhenoFileException, IOException {
        traitFileStrings = new ArrayList();
        BufferedReader reader = new BufferedReader(new FileReader(infile));
        phenoFile = infile;
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.length() == 0) {
                // skip blank lines
                continue;
            }
            if (line.startsWith("#")) {
                // skip comments
                continue;
            }
            traitFileStrings.add(line);
        }
        traitStrings = (ArrayList) traitFileStrings.clone();
        titleLine = (String) traitStrings.get(0);
        traitStrings.remove(0);
        int numLines = traitStrings.size();

        if (numLines < 2) {
            throw new GMDRPhenoFileException(
                    "Phenotype data format error: empty file");
        }
        StringTokenizer tokenizer = new StringTokenizer(titleLine, "\n\t\" \"");

        int numTokens = tokenizer.countTokens();

        if (numTokens < 2) {
            throw new GMDRPhenoFileException(
                    "Phenotype data format error: the title line is incorrect");
        }

        // reading the title line:get the marker number
        int c = 0;
        ArrayList temp = new ArrayList();
        while (tokenizer.hasMoreTokens()) {
            if (c++ < 2) {
                tokenizer.nextToken();
            } else {
                String trait = (String) tokenizer.nextToken();
                temp.add(trait);
            }
        }
        traitInfor = (ArrayList) temp.clone();
    }

    public File getPhenoFile() {
        return phenoFile;
    }

    public void parsePhenotype() throws GMDRPhenoFileException {
        int colNum = -1;
        int numLines = traitStrings.size();
        int numTraits = 0;
        if (numLines < 2) {
            throw new GMDRPhenoFileException(
                    "Pheno Data format error: empty file");
        }
        Subject sub;
        this.allSubjects = new ArrayList();

        for (int k = 0; k < numLines; k++) {
            StringTokenizer tokenizer = new StringTokenizer(
                    (String) traitStrings.get(k), "\n\t\" \"");
            int numTokens = tokenizer.countTokens();

            sub = new Subject(traitInfor.size());
            if (colNum < 1) {
                // only check column number count for the first nonblank line
                colNum = numTokens;
                numTraits = (numTokens - 2);
            }
            if (numTokens < 2) {
                throw new GMDRPhenoFileException(
                        "Incorrect number of fields in phefile. line " + (numLines + 2));
            }

            if (colNum != numTokens || colNum != traitInfor.size() + 2) {
                // this line has a different number of columns
                // should send some sort of error message
                throw new GMDRPhenoFileException(
                        "Column number mismatch in phefile. line " + (numLines + 2));
            }
            if (tokenizer.hasMoreTokens()) {
                sub.setFamilyID(new String(tokenizer.nextToken().trim()));
                sub.setSubjectID(new String(tokenizer.nextToken().trim()));
                while (tokenizer.hasMoreTokens()) {
                    try {
                        String trait = tokenizer.nextToken();
                        sub.AddTrait(trait);
                    } catch (NumberFormatException nfe) {
                        throw new GMDRPhenoFileException(
                                "Phenotype file input error: invalid Phenotype on line " + (numLines + 2));
                    }
                }

                // check if the family exists already in the Hashtable
                FamilyUnit fam = (FamilyUnit) this.families.get(sub.getFamilyID());
                if (fam == null) {
                    // it doesnt exist, so create a new Family object
                    fam = new FamilyUnit(sub.getFamilyID());
                }

                if (fam.getSubjects().containsKey(sub.getSubjectID())) {
                    throw new GMDRPhenoFileException("Individual " + sub.getSubjectID() + " in family " + sub.getFamilyID() + " appears more than once.");
                }

                fam.addSubject(sub);
                this.families.put(sub.getFamilyID(), fam);
                this.allSubjects.add(sub);
            }
        }
    }

    public FamilyUnit getFamilyUnit(String familyUnitID) {
        return (FamilyUnit) this.families.get(familyUnitID);
    }

    public Enumeration getFamUnitList() {
        return this.families.keys();
    }

    public int getNumFamilyUnits() {
        return families.size();
    }

    public int getNumSubjects() {
        return allSubjects.size();
    }

    public int getNumTraits() {
        return traitInfor.size();
    }

    public ArrayList getTraitName() {
    	return traitInfor;
    }
    
    public String getTraitAtI(int index)  {
        if (index < 0 || index > traitInfor.size()) {
            System.err.println("Could not find the phenotype at index " + index);
            System.exit(0);
        }
        return (String) traitInfor.get(index);
    }
}