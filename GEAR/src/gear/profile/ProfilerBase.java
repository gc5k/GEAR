package gear.profile;

import gear.CmdArgs;
import gear.profile.struct.QScore;
import gear.profile.struct.ScoreUnit;
import gear.util.Logger;
import gear.util.NewIt;

import java.util.ArrayList;
import java.util.HashMap;

public abstract class ProfilerBase
{
	public ProfilerBase()
	{
		// read score file
		gear.util.BufferedReader scoreReader = new gear.util.BufferedReader(CmdArgs.INSTANCE.getProfileArgs().getScoreFile(), "score");
		ScoreUnit scoreUnit;
		while ((scoreUnit = ScoreUnit.getNextScoreUnit(scoreReader)) != null)
		{
			Score.put(scoreUnit.getSNP(), scoreUnit);
		}
		scoreReader.close();

		Logger.printUserLog("Number of predictors: " + Score.size());

		// read q score file and q range file
		String qScoreFile = CmdArgs.INSTANCE.getProfileArgs().getQScoreFile(),
			   qScoreRangeFile = CmdArgs.INSTANCE.getProfileArgs().getQScoreRangeFile();
		if (qScoreFile != null && qScoreRangeFile != null)
		{
			// q score file
			gear.util.BufferedReader qScoreReader = new gear.util.BufferedReader(qScoreFile, "q-score");
			QScore qScore;
			while ((qScore = QScore.getNextQScore(qScoreReader)) != null)
			{
				QS.put(qScore.getSNP(), qScore);
			}
			qScoreReader.close();

			if (QS.size() == 0)
			{
				Logger.printUserError("Nothing is selected in '" + qScoreFile + "'.");
				System.exit(1);
			} else
			{
				Logger.printUserLog("Number of q-scores: " + QS.size());
			}

			// q range file
			gear.util.BufferedReader qRangeReader = new gear.util.BufferedReader(qScoreRangeFile, "q-score-range");
			ArrayList<ArrayList<String>> QR = NewIt.newArrayList();
			while (true)
			{
				String[] tokens = qRangeReader.readTokens(3);
				if (tokens == null)
				{
					break;
				}
				ArrayList<String> qr = NewIt.newArrayList();
				qr.add(tokens[0]); // range label
				qr.add(tokens[1]); // lower bound
				qr.add(tokens[2]); // upper bound
				QR.add(qr);
			}
			qRangeReader.close();

			if (QR.isEmpty())
			{
				Logger.printUserError("Nothing is selected in '" + qScoreRangeFile + "'.");
			} else
			{
				Logger.printUserLog("Number of q-ranges: " + QR.size());
			}

			q_score_range = new double[QR.size()][2];
			QRName = new String[QR.size()];

			for (int i = 0; i < QR.size(); i++)
			{
				ArrayList<String> qr = QR.get(i);
				QRName[i] = qr.get(0);
				q_score_range[i][0] = Double.parseDouble(qr.get(1));
				q_score_range[i][1] = Double.parseDouble(qr.get(2));
			}

			isQ = true;
		}
	}
	
	public void makeProfile()
	{
		if (isQ)
		{
			multipleProfile();
		}
		else
		{
			singleProfile();
		}
	}
	
	protected abstract void multipleProfile();
	
	protected abstract void singleProfile();
	
	protected HashMap<String, ScoreUnit> Score = new HashMap<String, ScoreUnit>();
	
	protected HashMap<String, QScore> QS = new HashMap<String, QScore>();
	
	protected double[][] q_score_range;
	
	protected String[] QRName;
	
	private boolean isQ = false;
}
