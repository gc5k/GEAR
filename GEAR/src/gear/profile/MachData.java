package gear.profile;

import gear.family.pedigree.file.SNP;
import gear.util.BufferedReader;
import gear.util.Logger;

import java.io.File;
import java.util.ArrayList;

public class MachData extends Data
{
	public static MachData create(String dosageFile, String infoFile, String dosageBatch, String infoBatch)
	{
		String[] dosageFiles, infoFiles;
		
		if (dosageFile != null)
		{
			assert(infoFile != null);
			
			if (!(new File(dosageFile)).exists())
			{
				Logger.printUserError("The dosage file '" + dosageFile + "' does not exist");
				System.exit(1);
			}
			dosageFiles = new String[1];
			dosageFiles[0] = dosageFile;
			
			infoFiles = new String[1];			
			infoFiles[0] = infoFile;
			if (!(new File(infoFiles[0])).exists())
			{
				Logger.printUserError("The .mlinfo file " + infoFiles[0] + "' does not exist");
				System.exit(1);
			}
		}
		else
		{
			assert(dosageBatch != null && infoBatch != null);
			
			String fileName = null;
			ArrayList<String> fileNameList = new ArrayList<String>();
			
			BufferedReader dosageBatchReader = BufferedReader.openTextFile(dosageBatch, "MaCH dosage batch");
			while ((fileName = dosageBatchReader.readNonEmptyLine()) != null)
			{
				if (!(new File(fileName)).exists())
				{
					Logger.printUserError("The dosage file '" + fileName + "' does not exist.");
					System.exit(1);
				}
				fileNameList.add(fileName);
			}
			dosageFiles = (String[]) fileNameList.toArray(new String[0]);
			
			fileNameList.clear();
			
			gear.util.BufferedReader infoBatchReader = BufferedReader.openTextFile(infoBatch, "MaCH information batch");
			while ((fileName = infoBatchReader.readNonEmptyLine()) != null)
			{
				if (!(new File(fileName)).exists())
				{
					Logger.printUserError("The information file '" + fileName + "' does not exist.");
					System.exit(1);
				}
				fileNameList.add(fileName);
			}
			
			infoFiles = (String[]) fileNameList.toArray(new String[0]);
		}
		
		if (dosageFiles.length == 0)
		{
			Logger.printUserError("No dosage file is provided.");
			System.exit(1);
		}

		if (dosageFiles.length != infoFiles.length)
		{
			String msg = "";
			msg += "The number of dosage files (" + dosageFiles.length + ") and ";
			msg += "the number of information files (" + infoFiles.length + ") are inconsistent";
			Logger.printUserError(msg);
			System.exit(1);
		}
		
		return new MachData(dosageFiles, infoFiles);
	}
	
	public MachData(String[] dosageFiles, String[] infoFiles)
	{
		assert dosageFiles != null;
		assert infoFiles != null;
		assert dosageFiles.length > 0;
		assert dosageFiles.length == infoFiles.length;
		
		this.dosageFiles = dosageFiles;
		this.infoFiles = infoFiles;
		initSNPs();
	}
	
	private void initSNPs()
	{
		ArrayList<SNP> snps = new ArrayList<SNP>();
		numLociInFile = new int[infoFiles.length];
		
		for (int fileIdx = 0; fileIdx < infoFiles.length; ++fileIdx)
		{
			BufferedReader reader = BufferedReader.openTextFile(infoFiles[fileIdx], ".mlinfo");
			
			reader.readTokens(7);  // Ignore the title line
			
			String[] tokens;
			while ((tokens = reader.readTokens(7)) != null)
			{
				String locusName = tokens[0];
				
				if (tokens[1].length() != 1)
				{
					String msg = "";
					msg += "'" + tokens[1] + "' is not a valid allele. ";
					msg += "Notice that every allele must be labled by a character.";
					reader.errorPreviousLine(msg);
				}
				char allele1 = tokens[1].charAt(0);
				
				if (tokens[2].length() != 1)
				{
					String msg = "";
					msg += "'" + tokens[2] + "' is not a valid allele. ";
					msg += "Notice that every allele must be labled by a character.";
					reader.errorPreviousLine(msg);
				}
				char allele2 = tokens[2].charAt(0);
				
				snps.add(new SNP(locusName, allele1, allele2));
				
				++numLociInFile[fileIdx];
			}
			reader.close();
		}
		
		this.snps = snps.toArray(new SNP[0]);
	}
	
	public class Iterator extends Data.Iterator
	{
		@Override
		public boolean next()
		{
			if (dosages != null && ++locusIdxInCurFile < dosages.length)
			{
				return true;
			}
			
			String[] tokens;
			
			while (dosageReader == null || (tokens = dosageReader.readTokens(numLociInFile[dosageFileIdx] + 2)) == null)
			{
				if (dosageReader != null)
				{
					dosageReader.close();
					if (indIdx != indIDs.size() - 1)
					{
						String msg = "";
						msg += "'" + dosageFiles[0] + "' contains " + indIDs.size() + " individual(s), ";
						msg += "but '" + dosageFiles[dosageFileIdx] + "' contains " + (indIdx + 1) + " individual(s).";
						Logger.printUserError(msg);
						System.exit(1);
					}
				}
				
				if (++dosageFileIdx >= dosageFiles.length)
				{
					return false;
				}
				
				dosageReader = BufferedReader.openZipFile(dosageFiles[dosageFileIdx], "dosage");
				indIdx = -1;
			}
			
			++indIdx;
			
			if (dosageFileIdx == 0)
			{
				indIDs.add(tokens[0]);
			}
			else
			{
				if (indIdx >= indIDs.size())
				{
					String msg = "";
					msg += "'" + dosageFiles[0] + "' contains " + indIDs.size() + " individual(s), ";
					msg += "but '" + dosageFiles[dosageFileIdx] + "' contains more than that.";
					Logger.printUserError(msg);
					System.exit(1);
				}
				
				if (!tokens[0].equals(indIDs.get(indIdx)))
				{
					String msg = "";
					msg += "Individual " + (indIdx + 1) + " in '" + dosageFiles[0] + "' is '" + indIDs.get(indIdx) + "', ";
					msg += "but individual " + (indIdx + 1) + " in '" + dosageFiles[dosageFileIdx] + "' is '" + tokens[0] + "'.";
					Logger.printUserError(msg);
					System.exit(1);
				}
			}
			
			String[] id = tokens[0].split("->", 2);
			if (id.length != 2)
			{
				dosageReader.errorPreviousLine("The first column must be in the form of 'FamilyID->IndividualID'");
			}
			famID = id[0];
			indID = id[1];
			
			dosages = new float[numLociInFile[dosageFileIdx]];
			for (int col = 2; col < tokens.length; ++col)
			{
				boolean failedToParse = false;
				
				try
				{
					dosages[col-2] = Float.parseFloat(tokens[col]);
				}
				catch (NumberFormatException e)
				{
					failedToParse = true;
				}
				
				if (failedToParse || dosages[col-2] < 0.0f || dosages[col-2] > 2.0f)
				{
					String msg = "";
					msg += "'" + tokens[col-2] + "' is not a valid dosage. ";
					msg += "Dosages must be floating point values ranging from 0.0 to 2.0.";
					dosageReader.errorPreviousLine(msg);
				}
			}
			
			locusIdxInCurFile = 0;
			
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
			return famID;
		}

		@Override
		public String getIndividualID()
		{
			return indID;
		}

		@Override
		public int getLocusIndex()
		{
			int numLociTraversed = 0;
			for (int i = 0; i < dosageFileIdx; ++i)
			{
				numLociTraversed += numLociInFile[i];
			}
			return numLociTraversed + locusIdxInCurFile;
		}

		@Override
		public float getAllele1Fraction()
		{
			return dosages[locusIdxInCurFile];
		}
		
		private BufferedReader dosageReader;
		private int dosageFileIdx = -1;
		private float[] dosages;
		private int locusIdxInCurFile = -1;
		private int indIdx = -1;
		private String famID;
		private String indID;
		private ArrayList<String> indIDs = new ArrayList<String>();
	}  // class Iterator
	
	@Override
	public SNP[] getSNPs()
	{
		return snps;
	}
	
	@Override
	public Iterator iterator()
	{
		return new Iterator();
	}
	
	private String[] dosageFiles;
	private String[] infoFiles;
	private int[] numLociInFile;
	private SNP[] snps; 
}
