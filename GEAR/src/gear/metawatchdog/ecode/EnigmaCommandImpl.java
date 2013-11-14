package gear.metawatchdog.ecode;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.random.RandomDataImpl;

import gear.CommandArguments;
import gear.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;

public class EnigmaCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		EnigmaCommandArguments enigmaArgs = (EnigmaCommandArguments)cmdArgs;
		
		readRefAlleles(enigmaArgs.getMapFile());
		
		RandomDataImpl rnd = new RandomDataImpl();
		rnd.reSeed(enigmaArgs.getSeed());
		double[][] beta = new double[ref.size()][enigmaArgs.getNumberOfColumns()];
		for(int i = 0; i < beta.length; i++)
		{
			for(int j = 0; j < beta[i].length; j++)
			{
				beta[i][j] = rnd.nextGaussian(0, 1);				
			}
		}
		
		String fileName = enigmaArgs.getOutRoot() + ".enigma";
		PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + fileName + "'.");
		}
		
		DecimalFormat df=new DecimalFormat("##.0000");

		for(int snpIdx = 0; snpIdx < ref.size(); snpIdx++)
		{
			writer.print(ref.get(snpIdx));
			for (int colIdx = 0; colIdx < enigmaArgs.getNumberOfColumns(); colIdx++)
			{
				writer.print("\t" + df.format(beta[snpIdx][colIdx]));
			}
			writer.println();
		}
		writer.close();
		Logger.printUserLog(ref.size() + " SNPs were used for generating Enigma scores.");
	}
	
	private void readRefAlleles(String mapFile)
	{
		BufferedReader reader = BufferedReader.openTextFile(mapFile, "map");
		String[] tokens;
		while ((tokens = reader.readTokensAtLeast(2)) != null)
		{
			ref.add(tokens[0] + "\t" + tokens[1]);
		}
	}
	
	private ArrayList<String> ref = new ArrayList<String>();
}
