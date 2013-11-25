package gear.subcommands.profile;

public class AdditiveCoeffModel extends CoeffModel
{
	@Override
	public float compute(float scoreAlleleFrac)
	{
		return scoreAlleleFrac;
	}
}
