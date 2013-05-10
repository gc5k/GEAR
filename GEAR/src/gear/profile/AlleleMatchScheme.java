package gear.profile;

public enum AlleleMatchScheme
{
	MATCH_NONE,
	MATCH_ALLELE1,			// The score allele of a locus matches the first allele in the data file
	MATCH_ALLELE2,			// The score allele of a locus matches the second allele in the data file
	MATCH_ALLELE1_FLIPPED,	// The score allele of a locus matches the flipped first allele in the data file
	MATCH_ALLELE2_FLIPPED	// The score allele of a locus matches the flipped second allele in the data file
}
