package gear.data;

interface SubjectOrder
{
	public String getSubjectOrderName();
	
	public int getNumberOfSubjects();
	
	/** 
	 * @return -1 if the subject doesn't exist
	 */
	public int getSubjectIndex(SubjectID subjectID);
	
	public SubjectID getSubjectID(int subjectIdx);
	
	public void swapSubjects(int subjectIdx1, int subjectIdx2);
}
