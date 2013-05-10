package gear.profile;

public class DominanceCoeffModel extends CoeffModel
{
	@Override
	public float compute(float scoreAlleleFrac)
	{
		return Math.signum(scoreAlleleFrac);
	}
}
