package gear.subcommands.oath.oathx;

import java.io.PrintStream;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.synthesize.SynthCommandArguments;
import gear.subcommands.oath.synthesize.SynthCommandImpl;
import gear.util.BufferedReader;
import gear.util.Combination;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OATHXCommandImpl extends CommandImpl 
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		obArgs = (OATHXCommandArguments) cmdArgs;

		readCM();
        Integer[] tmp= new Integer[corMat.length - 1];
        for (int i = 0; i < tmp.length; i++)
        {
        	tmp[i] = i+2;
        }

		String Fout = obArgs.getOutRoot() + ".oath-bus";
		PrintStream OBout = FileUtil.CreatePrintStream(Fout);

        int cnt = 1;
        for (int ii = 1; ii <= tmp.length; ii++)
        {
//            ArrayList<Integer[]> rs=cmn(tmp,ii);
            int[][] rs = Combination.combination(tmp, ii);
            for (int i=0; i<rs.length;i++)  
            {
                ArrayList<String> mod = NewIt.newArrayList();
                mod.add("1");
            	for (int j=0; j<rs[i].length;j++)
            	{
            		mod.add(new Integer(tmp[rs[i][j]]).toString());
            	}

            	SynthCommandArguments synArgs = new SynthCommandArguments();
            	synArgs.setNSSBatch(obArgs.getBatchFile());
            	synArgs.setCMFile(obArgs.getCMFile());
            	synArgs.setKeepBatch(mod.toArray(new String[0]));
            	synArgs.setN((new Integer(obArgs.getN())).toString());
            	if(obArgs.isVerbose())
            	{
            		synArgs.setVerbose();
            	}
            	synArgs.setOutRoot(obArgs.getOutRoot() + "-" + (new Integer(cnt).toString()));
            	SynthCommandImpl synImpl = new SynthCommandImpl();

            	if (cnt == 1)
            	{
            		for (int j = 0; j < synArgs.getNSSFile().length; j++)
            		{
                		OBout.println("#" + (j+1) + " " + synArgs.getNSSFile()[j]);            			
            		}
            	}
            	String s = new String();
            	s = "OATH-bus is running model [#1";
        		OBout.print(obArgs.getOutRoot() + "-" + (new Integer(cnt).toString()) + ".oath: #1" );

        		for (int k = 0; k < rs[i].length; k++)
        		{
        			if (k == 0)
        			{
            			OBout.print("=#" + tmp[rs[i][k]]);
            			s += "=#" + tmp[rs[i][k]];
        			}
        			else
        			{
            			OBout.print("+#" + tmp[rs[i][k]]);
            			s += "+#" + tmp[rs[i][k]];
        			}
        		}
        		Logger.printUserLog(s+']');
        		OBout.println();
            	
            	synImpl.execute(synArgs);

            	cnt++;
            }
        }
        OBout.close();
        Logger.printUserLog("Finding the summary information in '" + obArgs.getOutRoot() + ".oath-bus'.");
	}

	private void readCM()
	{
		BufferedReader reader = null;
		reader = BufferedReader.openTextFile(obArgs.getCMFile(), "Correlation matrix file.");

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		if (tokenLen != obArgs.getNSSFile().length)
		{
			Logger.printUserError("The dimension of the matrix does not match the number of nss files.");
			Logger.printUserError("Gear quitted.");
			System.exit(1);
		}
		corMat = new double[tokenLen][tokenLen];

		int cnt=0;
		do
		{
			for(int i = 0; i < tokens.length; i++)
			{
				corMat[cnt][i] = Double.parseDouble(tokens[i]);
			}
			cnt++;
		} while ((tokens = reader.readTokens()) != null);
		Logger.printUserLog("Reading " +tokenLen +"X" +tokenLen +" correlation matrix from '" + obArgs.getCMFile() + "'.");
	}

	private OATHXCommandArguments obArgs = null;
	private double[][] corMat;
}
