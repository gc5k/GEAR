package admixture;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public interface AdmixtureConstant {
	static public int Marry_In = 1;
	static public int Mirgrate_In = 2;
	static public double epsilon = 10E-6;
	static public boolean With_Genetic_Drift = true;
	static public boolean Without_Genetic_Drift = false;

	static public boolean free_recombination = true;
	static public boolean With_replacement = true;
	static public boolean Without_replacement = false;
	static public boolean With_LD = true;
	static public boolean Without_LD = false;
	
	static public int FamilyExactAffected = 0;
	static public int FamilyMoreThanAffected = 1;
	static public int CaseControl = 2;
	static public int No_selection = 3;
	
	static public boolean printAllele = true;
}
