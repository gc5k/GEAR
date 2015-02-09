package gear.subcommands.profile;

import gear.ConstValues;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public final class ProfileCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		profCmdArgs = (ProfileCommandArguments)cmdArgs;
		
		boolean isKeepSC = true;
		if(profCmdArgs.getExtractFile() != null)
		{
			extractSCsnp = NewIt.newHashSet();
			BufferedReader reader = BufferedReader.openTextFile(profCmdArgs.getExtractFile(), "Profile Extract score");
			String[] tokens = null;

			while((tokens = reader.readTokens(1)) != null)
			{
				for(int i = 0; i < tokens.length; i++)
				{
					extractSCsnp.add(tokens[i]);
				}
			}
			Logger.printUserLog("Read " + extractSCsnp.size() + " SNPs in " + profCmdArgs.getExtractFile());
		}

		if(profCmdArgs.getRemoveFile() != null)
		{
			extractSCsnp = NewIt.newHashSet();
			BufferedReader reader = BufferedReader.openTextFile(profCmdArgs.getRemoveFile(), "Profile Remove score");
			String[] tokens = null;

			while((tokens = reader.readTokens(1)) != null)
			{
				for(int i = 0; i < tokens.length; i++)
				{
					extractSCsnp.add(tokens[i]);
				}
			}
			Logger.printUserLog("Read " + extractSCsnp.size() + " SNPs in " + profCmdArgs.getRemoveFile());
			isKeepSC = false;
		}
		
		ScoreFile scoreFile;
		if (profCmdArgs.getScoreFile() != null)
		{
			scoreFile = ScoreFile.readTextFile(profCmdArgs.getScoreFile(), profCmdArgs.getHasScoreHeader(), extractSCsnp, isKeepSC);
		}
		else
		{
			scoreFile = ScoreFile.readTextFileGZ(profCmdArgs.getScoreFileGZ(), profCmdArgs.getHasScoreHeader(), extractSCsnp, isKeepSC);
		}

		HashMap<String, Float> qScoreMap = readQScores();  // LocusName-to-QScore map
		QRange[] qRanges = readQRanges();
		
		Data genoData = initData();

		SNP[] snps = genoData.getSNPs();
		CoeffModel[] coeffModels = initCoeffModels(snps);

		FilteredSNPs filteredSNPs = new FilteredSNPs(snps, scoreFile, qScoreMap, qRanges, profCmdArgs);

		printSNPFilterResult(scoreFile, qRanges, filteredSNPs);

		ArrayList<String> famIDs = new ArrayList<String>();
		ArrayList<String> indIDs = new ArrayList<String>();
		ArrayList<float[][]> riskProfiles = new ArrayList<float[][]>();
		ArrayList<int[][]> numLociUsed = new ArrayList<int[][]>();
		
		for (int traitIdx = 0; traitIdx < scoreFile.getNumberOfTraits(); ++traitIdx)
		{
			Data.Iterator iter = genoData.iterator();
			while (iter.next())
			{
				int indIdx = iter.getIndividualIndex();
				int locIdx = iter.getLocusIndex();
				
				while (indIdx >= famIDs.size())
				{
					famIDs.add(null);
				}
				famIDs.set(indIdx, iter.getFamilyID());
				
				while (indIdx >= indIDs.size())
				{
					indIDs.add(null);
				}
				indIDs.set(indIdx, iter.getIndividualID());
				
				while (indIdx >= riskProfiles.size())
				{
					riskProfiles.add(new float [scoreFile.getNumberOfTraits()][filteredSNPs.getNumLocusGroups()]);
				}
				
				while (indIdx >= numLociUsed.size())
				{
					numLociUsed.add(new int [scoreFile.getNumberOfTraits()][filteredSNPs.getNumLocusGroups()]);
				}
				
				if (filteredSNPs.getScore(traitIdx, locIdx) == null)
				{
					continue;
				}
				
				float scoreAlleleFrac = 0.0f;
				switch (filteredSNPs.getAlleleMatchScheme(locIdx))
				{
				case MATCH_ALLELE1:
				case MATCH_ALLELE1_FLIPPED:
					scoreAlleleFrac = iter.getAllele1Fraction();
					if (profCmdArgs.isScale())
					{
						if(freq[locIdx][0] != Double.NaN && freq[locIdx][0] != 1.0)
						{
							float p = (float) freq[locIdx][0];
							scoreAlleleFrac = (float) ((scoreAlleleFrac - 2 * p) / Math.sqrt(2 * p * (1-p)));							
						}
						else
						{
							scoreAlleleFrac = 0.0f;
						}
					}
					break;
				case MATCH_ALLELE2:
				case MATCH_ALLELE2_FLIPPED:
					scoreAlleleFrac = 2.0f - iter.getAllele1Fraction();
					if (profCmdArgs.isScale())
					{
						if (freq[locIdx][1] != Double.NaN && freq[locIdx][1] != 1.0)
						{
							float p = (float) freq[locIdx][1];
							scoreAlleleFrac = (float) ((scoreAlleleFrac - 2*p) / Math.sqrt(2 * p * (1-p)));
						}
						else
						{
							scoreAlleleFrac = 0.0f;
						}
					}
					break;
				default:
					continue;
				}

				float riskValue = coeffModels[locIdx].compute(scoreAlleleFrac) * filteredSNPs.getScore(traitIdx, locIdx);
				
				for (int locGrpIdx = 0; locGrpIdx < filteredSNPs.getNumLocusGroups(); ++locGrpIdx)
				{
					if (filteredSNPs.isInLocusGroup(locIdx, locGrpIdx))
					{
						riskProfiles.get(indIdx)[traitIdx][locGrpIdx] += riskValue;
						++numLociUsed.get(indIdx)[traitIdx][locGrpIdx];
					}
				}
			}  // while data
		}  // for each trait

		if (profCmdArgs.getIsWeighted())
		{
			for (int indIdx = 0; indIdx < riskProfiles.size(); ++indIdx)
			{
				for (int traitIdx = 0; traitIdx < scoreFile.getNumberOfTraits(); ++traitIdx)
				{
					for (int locGrpIdx = 0; locGrpIdx < filteredSNPs.getNumLocusGroups(); ++locGrpIdx)
					{
						int denom = numLociUsed.get(indIdx)[traitIdx][locGrpIdx];
						if (denom != 0)
						{
							riskProfiles.get(indIdx)[traitIdx][locGrpIdx] /= profCmdArgs.getIsSameAsPlink() ? (denom << 1) : denom;
						}
					}
				}
			}
		}
		
		// Write result to file
		for (int locGrpIdx = 0; locGrpIdx < filteredSNPs.getNumLocusGroups(); ++locGrpIdx)
		{
			String fileName = profCmdArgs.getResultFile();
			if (qRanges != null)
			{
				fileName += ".q." + qRanges[locGrpIdx].getName();
			}
			fileName += ".profile";
			PrintStream predictorFile = FileUtil.CreatePrintStream(fileName);
			
			// Title Line
			predictorFile.print("FID\tIID");
			for (int traitIdx = 0; traitIdx < scoreFile.getNumberOfTraits(); ++traitIdx)
			{
				predictorFile.print("\tSCORE." + scoreFile.getTrait(traitIdx));
			}
			predictorFile.println();
			
			for (int indIdx = 0; indIdx < riskProfiles.size(); indIdx++)
			{
				predictorFile.print(famIDs.get(indIdx) + "\t" + indIDs.get(indIdx));
				for (int traitIdx = 0; traitIdx < scoreFile.getNumberOfTraits(); ++traitIdx)
				{
					predictorFile.print("\t" + riskProfiles.get(indIdx)[traitIdx][locGrpIdx]);
				}
				predictorFile.println();
			}
			predictorFile.close();
		}
	}

	private CoeffModel[] initCoeffModels(SNP[] snps)
	{
		CoeffModel addModel = new AdditiveCoeffModel();
		CoeffModel domModel = new DominanceCoeffModel();
		CoeffModel recModel = new RecessiveCoeffModel();
		
		CoeffModel[] coeffModels = new CoeffModel[snps.length];
		
		switch (profCmdArgs.getCoeffModelType())
		{
		case ADDITIVE:
			Logger.printUserLog("Additive model is used.");
			for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
			{
				coeffModels[snpIdx] = addModel;
			}
			break;
		case DOMINANCE:
			Logger.printUserLog("Dominance model is used.");
			for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
			{
				coeffModels[snpIdx] = domModel;
			}
			break;
		case RECESSIVE:
			Logger.printUserLog("Recessive model is used.");
			for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
			{
				coeffModels[snpIdx] = recModel;
			}
			break;
		case FILE:
			Logger.printUserLog("Model of each SNP is recorded in file '" + profCmdArgs.getCoeffModelFile() + "'.");
			readCoeffModelFile(profCmdArgs.getCoeffModelFile(), snps, coeffModels);
		}
		
		return coeffModels;
	}
	
	private void readCoeffModelFile(String fileName, SNP[] snps, CoeffModel[] coeffModels)
	{
		CoeffModel addModel = new AdditiveCoeffModel();
		CoeffModel domModel = new DominanceCoeffModel();
		CoeffModel recModel = new RecessiveCoeffModel();
		HashMap<String, CoeffModel> mapSNP2Model = new HashMap<String, CoeffModel>();
		
		BufferedReader reader = BufferedReader.openTextFile(fileName, "model");
		String[] tokens;
		while ((tokens = reader.readTokens(2)) != null)
		{
			String snpName = tokens[0];
			String modelName = tokens[1]; 
			CoeffModel prevModel = null;
			
			if (modelName.equalsIgnoreCase("ADD"))
			{
				prevModel = mapSNP2Model.put(snpName, addModel);
			}
			else if (modelName.equalsIgnoreCase("DOM"))
			{
				prevModel = mapSNP2Model.put(snpName, domModel);
			}
			else if (modelName.equalsIgnoreCase("REC"))
			{
				prevModel = mapSNP2Model.put(snpName, recModel);
			}
			else
			{
				reader.errorPreviousLine("Unrecognized model '" + modelName + "'. Model name must be 'ADD', 'DOM' or 'REC'");
			}
			
			if (prevModel != null)
			{
				reader.errorPreviousLine("SNP '" + snpName + "' is duplicated.");
			}
		}
		reader.close();
		
		for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
		{
			String snpName = snps[snpIdx].getName();
			CoeffModel coeffModel = mapSNP2Model.get(snpName);
			if (coeffModel == null)
			{
				Logger.printUserError("The model file '" + fileName + "' doesn't contain SNP '" + snpName + "'.");
				System.exit(1);
			}
			coeffModels[snpIdx] = coeffModel;
		}
	}

	private void printSNPFilterResult(ScoreFile scoreFile, QRange[] qRanges, FilteredSNPs filteredSNPs)
	{
		Logger.printUserLog("Number of monomorphic loci (removed): " + filteredSNPs.getNumMonoLoci());
		Logger.printUserLog("Number of ambiguous loci (A/T or C/G) " + (profCmdArgs.getIsKeepATGC() ? "detected: " : "removed: ") + filteredSNPs.getNumAmbiguousLoci());

		// Allele Matching Schemes
		Logger.printUserLog("Number of Scheme I predictors: predictor alleles were A1: " + filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE1));
		Logger.printUserLog("Number of Scheme II predictors: predictor alleles were A2: " + filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE2));
		if (profCmdArgs.getIsAutoFlip())
		{
			Logger.printUserLog("Number of Scheme III predictors: predictor alleles were flipped A1: " + filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE1_FLIPPED));
			Logger.printUserLog("Number of Scheme IV predictors: predictor alleles were flipped A2: " + filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE2_FLIPPED));
		}
		Logger.printUserLog("Number of score alleles matching none in the data file: " + filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_NONE));

		if (profCmdArgs.getQScoreFile() != null)
		{
			Logger.printUserLog("Number of loci having no q-scores: " + filteredSNPs.getNumLociNoQScore());
			for (int i = 0; i < filteredSNPs.getNumLocusGroups(); i++)
			{
				Logger.printUserLog("Number of loci within the range: " + qRanges[i].getLowerBound() + ", " + qRanges[i].getUpperBound() + " is " + filteredSNPs.getNumInLocusGroup(i));
			}
		}
		
		for (int traitIdx = 0; traitIdx < scoreFile.getNumberOfTraits(); ++traitIdx)
		{
			Logger.printUserLog("Number of loci having no score for trait " + scoreFile.getTrait(traitIdx) + ": " + filteredSNPs.getNumLociNoScore(traitIdx));
		}
	}

	private HashMap<String, Float> readQScores()
	{
		String qScoreFile = profCmdArgs.getQScoreFile();
		if (qScoreFile == null)
		{
			return null;
		}
		
		HashMap<String, Float> qScores = new HashMap<String, Float>();

		BufferedReader reader = BufferedReader.openTextFile(qScoreFile, "q-score");
		String[] tokens;
		
		while ((tokens = reader.readTokens(2)) != null)
		{
			if (!ConstValues.isNA(tokens[1]))
			{
				try
				{
					if (qScores.containsKey(tokens[0]))
					{
						Logger.printUserLog("Warning: Marker '" + tokens[0] +"' duplicated in '" + profCmdArgs.getQScoreFile() + "'" + ", the first instance used, others skipped.");
						continue;
					}
					qScores.put(/* locusName = */ tokens[0], Float.parseFloat(tokens[1]));
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + tokens[1] + "' is not a valid q-score.");
				}
			}
		}
		reader.close();

		Logger.printUserLog("Number of q-scores: " + qScores.size());
		
		return qScores;
	}

	private QRange[] readQRanges()
	{
		String qRangeFile = profCmdArgs.getQRangeFile();
		if (qRangeFile == null)
		{
			return null;
		}
		
		ArrayList<QRange> qRanges = new ArrayList<QRange>();

		BufferedReader reader = BufferedReader.openTextFile(qRangeFile, "q-score-range");
		String[] tokens;
		
		while ((tokens = reader.readTokens(3)) != null)
		{
			try
			{
				float lowerBound, upperBound;
				lowerBound = ConstValues.isNA(tokens[1]) ? Float.MIN_VALUE : Float.parseFloat(tokens[1]);
				upperBound = ConstValues.isNA(tokens[2]) ? Float.MAX_VALUE : Float.parseFloat(tokens[2]);
				
				for (QRange qRange : qRanges)
				{
					if (qRange.getName().equals(tokens[0]))
					{
						Logger.printUserError("QRange '" + qRange.getName() + "' appears more than once in the qrange file '" + profCmdArgs.getQRangeFile() + "'.");
						System.exit(1);
					}
				}
				
				qRanges.add(new QRange(/* name = */tokens[0], lowerBound, upperBound));
			}
			catch (NumberFormatException e)
			{
				reader.errorPreviousLine("An invalid q-score value is detected.");
			}
		}
		reader.close();

		Logger.printUserLog("Number of q-ranges: " + qRanges.size());

		return qRanges.toArray(new QRange[0]);
	}
	
	private Data initData()
	{
		PLINKParser plinkParser = PLINKParser.parse(profCmdArgs);
		if (plinkParser == null)
		{
			if (profCmdArgs.isScale()) 
			{
				Logger.printUserError("scale option is not only available for Mach format.");
				System.exit(0);
			}
			return MachData.create(profCmdArgs.getMachDosageFile(),
			                       profCmdArgs.getMachInfoFile(),
			                       profCmdArgs.getMachDosageBatch(),
			                       profCmdArgs.getMachInfoBatch());
		}
		else if (profCmdArgs.isScale())
		{

			SampleFilter sf = new SampleFilter(plinkParser.getPedigreeData(), plinkParser.getMapData());
			SumStatQC ssQC = new SumStatQC(plinkParser.getPedigreeData(), plinkParser.getMapData(), sf);
			GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
			freq = PopStat.calAlleleFrequency(gm, gm.getNumMarker());
			
			if (profCmdArgs.getScaleFile() != null)
			{
				HashMap<String, ScaleMaf> scaleMaf = readScaleFile();
				if (scaleMaf.size() > 0)
				{
					ArrayList<SNP> snplist = plinkParser.getMapData().getMarkerList();
					for (int i = 0; i < snplist.size(); i++)
					{
						SNP snp = snplist.get(i);
						if (scaleMaf.containsKey(snp.getName()))
						{
							ScaleMaf sm = scaleMaf.get(snp.getName());
							switch (FilteredSNPs.getMatchScheme(sm.getA1(), snp.getFirstAllele(), snp.getSecAllele(), profCmdArgs.getIsAutoFlip()))
							{
								case MATCH_ALLELE1:
								case MATCH_ALLELE1_FLIPPED:
									freq[i][0] = sm.getMaf();
									freq[i][1] = 1 - sm.getMaf();
									break;
								case MATCH_ALLELE2:
								case MATCH_ALLELE2_FLIPPED:
									freq[i][0] = 1 - sm.getMaf();
									freq[i][1] = sm.getMaf();
									break;
								default:
									continue;
							}
						}
					}
				}
				else
				{
					profCmdArgs.turneoffScale();
				}
			}
		}
		return new PlinkData(plinkParser);
	}

	private HashMap<String, ScaleMaf> readScaleFile()
	{
		BufferedReader reader = BufferedReader.openTextFile(profCmdArgs.getScaleFile(), "maf");
		HashMap<String, ScaleMaf> scaleMaf = NewIt.newHashMap();
		String[] tokens;
		while ((tokens = reader.readTokensAtLeast(3)) != null)
		{
			if (ConstValues.isNA(tokens[2])) continue;
			scaleMaf.put(tokens[0], new ScaleMaf(tokens[0], tokens[1].charAt(0), Double.parseDouble(tokens[2])));
		}
		return scaleMaf;
	}

	private ProfileCommandArguments profCmdArgs;
	private HashSet<String> extractSCsnp = null;
	private double[][] freq = null;
}
