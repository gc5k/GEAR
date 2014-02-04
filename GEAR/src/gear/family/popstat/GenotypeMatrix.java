package gear.family.popstat;

import java.util.ArrayList;

import gear.CmdArgs;
import gear.data.Person;
import gear.family.pedigree.PersonIndex;

public class GenotypeMatrix
{
	protected int[][] genotypeMat;
	protected int lenMat;
	protected int lenMatF;
	protected final int shift = 4;
	protected int numMarker = 0;
	protected ArrayList<PersonIndex> pidx;

	public GenotypeMatrix(ArrayList<PersonIndex> pi)
	{
		pidx = pi;
		initial();
	}

	protected void initial()
	{
		genotypeMat = new int[pidx.size()][];

		int c1 = 0;
		for (PersonIndex pi : pidx)
		{
			if (!pi.isPseudo())
			{
				genotypeMat[c1++] = pi.getPerson().getAlleleArray();
			}
		}
		numMarker = pidx.get(0).getPerson().getNumMarkers();

	}

	public int getNumIndivdial() 
	{
		return getGRow();
	}

	public int getNumMarker()
	{
		return numMarker;
	}

	public void setAdditiveScore(int idx, int i, int v)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		genotypeMat[idx][posByte] &= ~(3 << posBit);
		genotypeMat[idx][posByte] |= (v & 3) << posBit;
	}

	public int getAdditiveScore(int idx, int i)
	{
		// 0 homozygote 1/1
		// 1 heterozygosity 1/2
		// 2 homozygote 2/2
		// 3 missing
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		return g;
		// if (g == 1) {// 01
		// return 3;
		// } else if (g == 0) {
		// return 0;
		// } else if (g == 2) {
		// return 1;
		// } else {
		// return 2;
		// }
	}

	public int getAdditiveScoreOnFirstAllele(int idx, int i)
	{
		// 0 homozygote 1/1
		// 1 heterozygosity 1/2
		// 2 homozygote 2/2
		// 3 missing
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		return g == 3 ? g : (2-g);
		// if (g == 1) {// 01
		// return 3;
		// } else if (g == 0) {
		// return 0;
		// } else if (g == 2) {
		// return 1;
		// } else {
		// return 2;
		// }
	}

	public byte getOriginalGenotypeScore(int ind, int i)
	{// this only for write back to bed file purpose
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (genotypeMat[ind][posByte] >> (posBit)) & 3;
		switch (g)
		{
		case 0:
			g = 0;
			break;
		case 1:
			g = 2;
			break;
		case 2:
			g = 3;
			break;
		default:
			g = 1;
			break; // missing
		}
		return (byte) g;
	}

	public int[] getBiAlleleGenotype(int idx, int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> posBit) & 3;
		int[] b = { 0, 0 };
		switch (g)
		{
		case 0:
			b[0] = 0;
			b[1] = 0;
			break;
		case 1:
			b[0] = 0;
			b[1] = 1;
			break;
		case 2:
			b[0] = 1;
			b[1] = 1;
			break;
		default:
			b[0] = Person.MissingAlleleCode;
			b[1] = Person.MissingAlleleCode;
			break;
		}
		return b;
	}

	public String getGenotypeScoreString(int idx, int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		if (g == 3)
		{// 01
			return CmdArgs.INSTANCE.missingGenotype;
		} else
		{
			return Integer.toString(g);
		}
	}

	public int getGRow()
	{
		return genotypeMat.length;
	}

	public int getGCol()
	{
		return genotypeMat[0].length;
	}

	public int[][] getG()
	{
		return genotypeMat;
	}
}
