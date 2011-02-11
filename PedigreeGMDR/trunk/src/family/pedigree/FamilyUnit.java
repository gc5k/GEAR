package family.pedigree;

import java.util.Enumeration;
import java.util.Hashtable;

import util.NewIt;

/**
 * stores the familyName and the members of a family from a pedigree file this class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class FamilyUnit {
//save phenotype information;
    private Hashtable<String, Subject> subjects;
    private String familyUnitName;
    private int mendErrors;

    public FamilyUnit() {
        this.subjects = NewIt.newHashtable();
    }

    public FamilyUnit(String familyName) {
        this.subjects = NewIt.newHashtable();
        this.familyUnitName = familyName;
    }

    /**
     * returns the name of this family (familyName)
     * 
     * @return family name
     */
    public String getFamilyUnitName() {
        return familyUnitName;
    }

    /**
     * sets the family name (familyName)
     * 
     * @param familyName
     */
    public void setFamilyName(String familyName) {
        this.familyUnitName = familyName;
    }

    /**
     * returns the subjects Hashtable (containing individuals)
     * 
     * @return subjects Hashtable
     */
    public Hashtable<String, Subject> getSubjects() {
        return subjects;
    }

    /**
     * sets the subjects hashtable
     * 
     * @param subjects
     */
    public void setSubjects(Hashtable<String, Subject> subjects) {
        this.subjects = subjects;
    }

    public int getMendErrs() {
        return mendErrors;
    }

    public void setMendErrs(int mendErrors) {
        this.mendErrors = mendErrors;
    }

    /**
     * returns the number of subjects of this family
     * 
     * @return number of subjects in this family
     */
    public int getNumsubjects() {
        return this.subjects.size();
    }

    /**
     * returns a list of individualIDs that are subjects of this family in the form of an enumeration
     * 
     * @return enumeration memberlist
     */
    public Enumeration<String> getSubjectsList() {
        return this.subjects.keys();
    }

    /**
     * adds a member to this family (adds to subjects Vector)
     * 
     * @param ind
     *            Individual to add to subjects Vector
     */
    public void addSubject(Subject ind) {
        this.subjects.put(ind.getSubjectID(), ind);
    }

    public void removeSubject(String id) {
        subjects.remove(id);
    }

    public boolean containsSubject(String id) {
        if (subjects.containsKey(id)) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * get the Individual with individualID
     * 
     * @param individualID
     *            the individualID of the individual we want
     * @return the individual with matching individualID
     */
    public Subject getSubject(String subjectID) {
        return (Subject) this.subjects.get(subjectID);
    }
}
