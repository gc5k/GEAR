package gear.subcommands.oath.oathbus;

import java.io.PrintStream;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.nss.NSSCommandArguments;
import gear.subcommands.oath.nss.NSSCommandImpl;
import gear.subcommands.oath.synthesize.SynthCommandArguments;
import gear.subcommands.oath.synthesize.SynthCommandImpl;
import gear.util.Combination;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OATHBusCommandImpl extends CommandImpl 
{

	OATHBusCommandArguments obArgs = null;

	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		OATHBusCommandArguments obArgs = (OATHBusCommandArguments) cmdArgs;
		
		NSSCommandArguments nssArgs = new NSSCommandArguments();
		nssArgs.setBFile(obArgs.getBFile());
		nssArgs.setPhenotypeFile(obArgs.getPhenotypeFile());
		nssArgs.setPhenotypeIndex(obArgs.getMpheno()[0]);
		nssArgs.setCovFile(obArgs.getCovFile());
		nssArgs.setCovNumber(obArgs.getCovNumber());
		nssArgs.setMAF((new Double(obArgs.getMAF())).toString());
		if (obArgs.getKeepFile() != null)
		{
			nssArgs.setKeeFile(obArgs.getKeepFile());
		}
		
		if (obArgs.getChr()!= null)
		{
			nssArgs.setChr(obArgs.getChr().toArray(new String[0]));
		}

		nssArgs.setOutRoot(obArgs.getOutRoot());

		NSSCommandImpl nssImpl = new NSSCommandImpl();
		nssImpl.execute(nssArgs);
		
        Integer[] tmp= new Integer[obArgs.getCovNumber().length];
        for(int i = 0; i < obArgs.getCovNumber().length; i++)
        {
        	tmp[i] = i+2;
        }

		String Fout = obArgs.getOutRoot() + ".oath-bus";
		PrintStream OBout = FileUtil.CreatePrintStream(Fout);

        int cnt = 1;
        for(int ii = 1; ii <= tmp.length; ii++)
        {
//            ArrayList<Integer[]> rs=cmn(tmp,ii);
            int[][] rs = Combination.combination(tmp, ii);
            for(int i=0; i<rs.length;i++)  
            {
                ArrayList<String> mod = NewIt.newArrayList();
                mod.add("1");
            	for(int j=0; j<rs[i].length;j++)  
            	{
            		mod.add(new Integer(tmp[rs[i][j]]).toString());
            	}

            	SynthCommandArguments synArgs = new SynthCommandArguments();
            	synArgs.setNSSBatch(obArgs.getOutRoot()+".list.nss");
            	synArgs.setCMFile(obArgs.getOutRoot()+".m.nss");
            	synArgs.setKeepBatch(mod.toArray(new String[0]));
            	synArgs.setN((new Integer(nssImpl.getN())).toString());
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

}
