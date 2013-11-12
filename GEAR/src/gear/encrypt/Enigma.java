package gear.encrypt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.random.RandomDataImpl;

import gear.CmdArgs;
import gear.util.Logger;
import gear.util.NewIt;

public class Enigma
{
	protected static final String DELIMITER = "\\s+";
	private String refFile;
	private long eseed;
	private int ecol;
	private ArrayList<String> ref = NewIt.newArrayList();
	
	private RandomDataImpl rnd;

	public Enigma()
	{
		refFile = CmdArgs.INSTANCE.refAllele;
		eseed = CmdArgs.INSTANCE.Ecode;
		ecol = CmdArgs.INSTANCE.Ecol;
		rnd = new RandomDataImpl();
		rnd.reSeed(eseed);

		Logger.printUserLog("Enigma parameter setting is as below:");
		Logger.printUserLog("Ecode is: " + eseed);
		Logger.printUserLog("ECol number is: " + ecol);
		Logger.printUserLog("Reference allele file is " + refFile);
	}

	public void Revup()
	{
		readRefAllele();
		double[][] beta = new double[ref.size()][ecol];
		for(int i = 0; i < beta.length; i++)
		{
			for(int j = 0; j < beta[i].length; j++)
			{
				beta[i][j] = rnd.nextGaussian(0, 1);				
			}
		}
		
		PrintWriter ec = null;
		try 
		{
			ec = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
		                                    					+ ".enigma")));
		}
		catch (IOException e)
		{
			Logger.handleException(e,
				"An exception occurred when writing files.");
		}
		
		DecimalFormat df=new DecimalFormat("##.0000");

		for(int i = 0; i < ref.size(); i++)
		{
			ec.print(ref.get(i));
			for (int j = 0; j < ecol; j++)
			{
				if(j < (ecol -1) ) 
				{
					ec.print(df.format(beta[i][j]) + "\t");
				}
				else 
				{
					ec.println(df.format(beta[i][j]));
				}
			}
		}
		ec.close();
		Logger.printUserLog(ref.size() + " SNPs were used for generating Enigma scores.");
	}

	private void readRefAllele() 
	{
		File mapfile = new File(refFile);

		BufferedReader reader = null;
		try
		{
			reader = new BufferedReader(new FileReader(mapfile));
		} catch (IOException e)
		{
			Logger.handleException(e, "Cannot open the map file '" + mapfile
					+ "'.");
		}

		String line = null;
		int fc = 0;

		try
		{
			int c = 0;
			while ((line = reader.readLine()) != null)
			{
				c++;
				String[] tokens = line.split(DELIMITER);
				if (tokens.length < 2)
				{
					Logger.printUserLog("line " + c + " was ignored.");
					fc++;
				}
				// Skip genetic distance field at tokens[2].
				StringBuffer sb = new StringBuffer();
				sb.append(tokens[0]+ "\t" + tokens[1] + "\t");
				ref.add(sb.toString());
			}
			reader.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occured when reading the map file '" + refFile
							+ "'.");
		}
		if(fc>0)
		{
			Logger.printUserLog(fc + " SNPs were failed for generating their Enigma codes.");
		}
	}
}
