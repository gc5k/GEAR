package gear.profile;

public class Score
{
	private char allele;
	private float value;

	public Score(char allele, float value)
	{
		this.allele = allele;
		this.value = value;
	}

	public char getAllele()
	{
		return allele;
	}

	public float getValue()
	{
		return value;
	}
	
	public void setValue(float value)
	{
		this.value = value;
	}
}
