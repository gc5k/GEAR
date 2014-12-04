package gear.subcommands.profile;

import gear.ConstValues;
import gear.util.BufferedReader;
import gear.util.Logger;

import java.util.HashMap;
import java.util.HashSet;

class ScoreFile
{
	protected static ScoreFile readTextFile(String fileName, boolean hasHeaders, HashSet<String> SCsnp, boolean isKeepSC)
	{
		BufferedReader reader = BufferedReader.openTextFile(fileName, "score");
		return new ScoreFile(reader, hasHeaders, SCsnp, isKeepSC);
	}

	protected static ScoreFile readTextFileGZ(String fileName, boolean hasHeaders, HashSet<String> SCsnp, boolean isKeepSC)
	{
		BufferedReader reader = BufferedReader.openGZipFile(fileName, "gz score");
		return new ScoreFile(reader, hasHeaders, SCsnp, isKeepSC);
	}
	
	private ScoreFile(BufferedReader reader, boolean hasHeaders, HashSet<String> SCsnp, boolean isKeepSC)
	{	
		String[] tokens = readFirstLine(reader, hasHeaders);
		int cnt = 0;
		while (tokens != null)
		{
			if (tokens[1].length() != 1)
			{
				reader.errorPreviousLine("'" + tokens[1] + "' is not a character, so it is not a valid allele.");
				continue;
			}

			if (SCsnp == null)
			{
				Score score = new Score(tokens[1].charAt(0), tokens.length - 2);

				if (scores.put(/*locusName*/tokens[0], score) != null)
				{
					Logger.printUserError("SNP '" + tokens[0] + "' appears more than once in the score file '" + reader.getFileName() + "'.");
					System.exit(1);
				}
			
				for (int ii = 2; ii < tokens.length; ++ii)
				{
					if (!ConstValues.isNA(tokens[ii]))
					{
						try
						{
							score.setValue(ii - 2, Float.parseFloat(tokens[ii]));
						}
						catch (NumberFormatException e)
						{
							reader.errorPreviousLine("'" + tokens[ii] + "' is not a floating point number, so it it not a valid score.");
						}
					}
				}
				cnt++;	
			}
			else if ( (isKeepSC && SCsnp.contains(tokens[0])) || (!isKeepSC && !SCsnp.contains(tokens[0])) )
			{
				Score score = new Score(tokens[1].charAt(0), tokens.length - 2);

				if (scores.put(/*locusName*/tokens[0], score) != null)
				{
					Logger.printUserError("SNP '" + tokens[0] + "' appears more than once in the score file '" + reader.getFileName() + "'.");
					System.exit(1);
				}
			
				for (int ii = 2; ii < tokens.length; ++ii)
				{
					if (!ConstValues.isNA(tokens[ii]))
					{
						try
						{
							score.setValue(ii - 2, Float.parseFloat(tokens[ii]));
						}
						catch (NumberFormatException e)
						{
							reader.errorPreviousLine("'" + tokens[ii] + "' is not a floating point number, so it it not a valid score.");
						}
					}
				}
				cnt++;
			}

			tokens = reader.readTokens(tokens.length);

		}
		reader.close();
		Logger.printUserLog(cnt + " scores have been included.");
	}

	/**
	 *  
	 * @param reader
	 * @param hasHeaders
	 * @return the tokens of the first non-header line
	 */
	private String[] readFirstLine(BufferedReader reader, boolean hasHeaders)
	{
		String[] tokens = reader.readTokensAtLeast(3);;
		
		traits = new String[tokens.length - 2];
		
		if (hasHeaders)
		{
			System.arraycopy(tokens, 2, traits, 0, traits.length);
			tokens = reader.readTokens(tokens.length);
		}
		else
		{
			for (int ii = 0; ii < traits.length; ++ii)
			{
				traits[ii] = String.valueOf(ii + 1);
			}
		}
		
		return tokens;
	}
	
	protected int getNumberOfTraits()
	{
		return traits.length;
	}
	
	protected String getTrait(int traitIdx)
	{
		return traits[traitIdx];
	}
	
	protected Score getScore(String snp)
	{
		return scores.get(snp);
	}
	
	private String[] traits;
	private HashMap<String, Score> scores = new HashMap<String, Score>();  // name-to-score
}
