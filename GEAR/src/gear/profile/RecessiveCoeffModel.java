package gear.profile;

public class RecessiveCoeffModel extends CoeffModel
{
	@Override
	public float compute(float scoreAlleleFrac)
	{
		return Math.signum(scoreAlleleFrac - 1.0f);
	}
}
