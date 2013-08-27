package gear.profile;

import gear.family.pedigree.file.SNP;

public abstract class Data
{
	public abstract SNP[] getSNPs();
	
	public abstract class Iterator
	{
		public abstract boolean next();
		
		public abstract int getIndividualIndex();
		
		public abstract String getFamilyID();
		
		public abstract String getIndividualID();
		
		public abstract int getLocusIndex();
		
		public abstract float getAllele1Fraction();
		
		public String getPhenotype()
		{
			return null;
		}
	}
	
	public abstract Iterator iterator();
}
