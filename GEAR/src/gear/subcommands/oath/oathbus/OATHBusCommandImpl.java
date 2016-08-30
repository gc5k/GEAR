package gear.subcommands.oath.oathbus;

import java.io.PrintStream;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.nss.NSSCommandArguments;
import gear.subcommands.oath.nss.NSSCommandImpl;
import gear.subcommands.oath.synthesize.SynthCommandArguments;
import gear.subcommands.oath.synthesize.SynthCommandImpl;
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
		
		if (obArgs.isChrFlagOn())
		{
			nssArgs.setChr((new Integer(obArgs.getChr())).toString());
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
            ArrayList<Integer[]> rs=cmn(tmp,ii);
            for(int i=0; i<rs.size();i++)  
            {
                ArrayList<String> mod = NewIt.newArrayList();
                mod.add("1");
            	for(int j=0; j<rs.get(i).length;j++)  
            	{
            		mod.add(new Integer(rs.get(i)[j]).toString());
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
        		for (int k = 0; k < rs.get(i).length; k++)
        		{
        			if (k == 0)
        			{
            			OBout.print("=#" + rs.get(i)[k]);
            			s += "=#" + rs.get(i)[k];
        			}
        			else
        			{
            			OBout.print("+#" + rs.get(i)[k]);
            			s += "+#" + rs.get(i)[k];        				
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

    static ArrayList<Object[]> RandomC(Object[] source)  
    {  
        ArrayList<Object[]> result=new ArrayList<Object[]>();  
        if(source.length==1)  
        {  
            result.add(source);          
        }  
        else  
        {  
            Object[] psource=new Object[source.length-1];  
            for(int i=0;i<psource.length;i++)  
            {  
                psource[i]=source[i];  
            }  
            result=RandomC(psource);  
            int len=result.size();//fn组合的长度  
            result.add((new Object[]{source[source.length-1]}));  
            for(int i=0;i<len;i++)  
            {  
                Object[] tmp=new Object[result.get(i).length+1];  
                for(int j=0;j<tmp.length-1;j++)  
                {  
                    tmp[j]=result.get(i)[j];  
                }  
                tmp[tmp.length-1]=source[source.length-1];  
                result.add(tmp);  
            }  
  
        }  
        return result;  
    }  
  
    static ArrayList<Integer []> cmn(Integer[] source,int n)  
    {  
        ArrayList<Integer[]> result=new ArrayList<Integer[]>();  
        if(n==1)  
        {  
            for(int i=0;i<source.length;i++)  
            {  
                result.add(new Integer[]{source[i]});  
            }
        }
        else if(source.length==n)
        {
            result.add(source);
        }
        else 
        {
            Integer[] psource=new Integer[source.length-1];  
            for(int i=0;i<psource.length;i++)  
            {  
                psource[i]=source[i];  
            }  
            result=cmn(psource,n);  
            ArrayList<Integer[]> tmp=cmn(psource,n-1);  
            for(int i=0;i<tmp.size();i++)  
            {  
                Integer[] rs=new Integer[n];  
                for(int j=0; j<n-1; j++)
                {
                    rs[j]=tmp.get(i)[j];
                }
                rs[n-1]=source[source.length-1];
                result.add(rs);
            }
        }
        return result;
    }
}
