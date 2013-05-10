package gear.profile;

import family.pedigree.file.SNP;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.CmdArgs;
import gear.ConstValues;

import java.util.ArrayList;

public class PlinkData extends Data
{
	public static PlinkData create()
	{
		PLINKParser parser = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			parser = new PLINKParser(CmdArgs.INSTANCE.getFileArgs().getPed(),
					                 CmdArgs.INSTANCE.getFileArgs().getMap());
		}
		else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			parser = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0).getBed(),
					                       CmdArgs.INSTANCE.getBFileArgs(0).getBim(),
					                       CmdArgs.INSTANCE.getBFileArgs(0).getFam());
		}
		else
		{
			assert false;
		}
		return new PlinkData(parser);
	}
	
	public PlinkData(PLINKParser parser)
	{
		parser.Parse();
		sampleFilter = new SampleFilter(parser.getPedigreeData(), parser.getMapData());
		genoMatrix = new GenotypeMatrix(sampleFilter.getSample());
		snpList = sampleFilter.getMapFile().getMarkerList();
		calcAllele1Frequencies();
	}
	
	private void calcAllele1Frequencies()
	{
		allele1Freqs = new float[snpList.size()];
		int[] numAlleles = new int[snpList.size()];
		
		for (int indIdx = 0; indIdx < genoMatrix.getGRow(); ++indIdx)
		{
			for (int locusIdx = 0; locusIdx < genoMatrix.getGCol(); ++locusIdx)
			{
				int genoValue = genoMatrix.getAdditiveScore(indIdx, locusIdx);
				if (genoValue != ConstValues.BINARY_MISSING_GENOTYPE)
				{
					allele1Freqs[locusIdx] += 2.0f - genoValue;
					numAlleles[locusIdx] += 2;
				}
			}
		}
		
		for (int locusIdx = 0; locusIdx < genoMatrix.getGCol(); ++locusIdx)
		{
			if (numAlleles[locusIdx] != 0)
			{
				allele1Freqs[locusIdx] /= numAlleles[locusIdx];
			}
		}
	}
	
	public class Iterator extends Data.Iterator
	{
		@Override
		public boolean next()
		{
			if (++locIdx >= snpList.size())
			{
				if (++indIdx >= genoMatrix.getGRow())
				{
					return false;
				}
				locIdx = 0;
			}
			return true;
		}

		@Override
		public int getIndividualIndex()
		{
			return indIdx;
		}

		@Override
		public String getFamilyID()
		{
			return sampleFilter.getSample().get(indIdx).getFamilyID();
		}

		@Override
		public String getIndividualID()
		{
			return sampleFilter.getSample().get(indIdx).getIndividualID();
		}

		@Override
		public int getLocusIndex()
		{
			return locIdx;
		}

		@Override
		public float getAllele1Fraction()
		{
			int genoValue = genoMatrix.getAdditiveScore(indIdx, locIdx);
			return genoValue == ConstValues.BINARY_MISSING_GENOTYPE ? 2 * allele1Freqs[locIdx] : 2 - genoValue;
		}

		@Override
		public String getPhenotype()
		{
			return sampleFilter.getHukouBook().get(indIdx).getCol6();
		}
		
		private int indIdx = 0, locIdx = -1;
	}

	@Override
	public SNP[] getSNPs()
	{
		return snpList.toArray(new SNP[0]);
	}

	@Override
	public Iterator iterator()
	{
		return new Iterator();
	}
	
	private GenotypeMatrix genoMatrix;
	private float[] allele1Freqs;
	private ArrayList<SNP> snpList;
	private SampleFilter sampleFilter;
}
