package family.pedigree;

import java.util.ArrayList;

/**
 * stores the phenotypic outcomes of each individual. this class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen
 */
public class Subject {

    private String familyID;
    private String subjectID;
    private ArrayList traits;
    private ArrayList missing;
    private int numTraits;
    private String missingData = ".";

    public Subject(int num) {
        traits = new ArrayList();
        missing = new ArrayList();
        numTraits = num;
    }

    public void AddTrait(final String t) {
        traits.add(t);
        if (t.equals(missingData)) {
            missing.add(new Boolean(true));
        }
    }

    public int getNumberofTraits() {
        return numTraits;
    }

    public ArrayList getTraits() {
        return traits;
    }

    public ArrayList getmissing() {
        return missing;
    }

    public void setFamilyID(String FID) {
        familyID = FID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public void setSubjectID(String ID) {
        subjectID = ID;
    }

    public String getSubjectID() {
        return subjectID;
    }
}
