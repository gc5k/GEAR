package gear.subcommands.hpc;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;

public final class HpcCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		HpcCommandArguments hpcCmdArgs = (HpcCommandArguments)cmdArgs;
		
		String scriptName = hpcCmdArgs.getJobName() + ".sh";

		PrintWriter pw = null;
		try
		{
			pw = new PrintWriter(scriptName);
		}
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "Cannot create the script file '" + scriptName + "'.");
		}

		pw.println("#$ -cwd");
		pw.println("#$ -l vf=" + (hpcCmdArgs.getRam() + 2) + hpcCmdArgs.getRamUnit());
		pw.println("#$ -l h_vmem=" + (hpcCmdArgs.getRam() + 2) + hpcCmdArgs.getRamUnit());
		pw.println("#$ -N " + hpcCmdArgs.getJobName());
		pw.println("#$ -m eas");
		if (hpcCmdArgs.getEmail() != null)
		{
			pw.println("#$ -M " + hpcCmdArgs.getEmail());
		}

		pw.print("java -jar " + getJarName() + " -Xmx" + hpcCmdArgs.getRam() + hpcCmdArgs.getRamUnit());
		
		String[] nonHpcOptsAndArgs = hpcCmdArgs.getNonHpcOptionsAndArguments();
		for (String s : nonHpcOptsAndArgs)
		{
			pw.print(" " + s);
		}
		pw.println();

		pw.close();

		if (hpcCmdArgs.getIsSubmit())
		{
			Runtime rt = Runtime.getRuntime();
			String cmd = "qsub " + scriptName;
			try
			{
				rt.exec(cmd);
			}
			catch (IOException e)
			{
				Logger.handleException(e, "Failed to execute command '" + cmd + "'.");
			}
		}
	}
	
	private String getJarName()
	{
		String myClassFile = this.getClass().getName().replaceAll("\\.", "/") + ".class";
		URL jarURL = ClassLoader.getSystemResource(myClassFile);
		String urlStr = jarURL.toString();   
		int from = "jar:file:".length();   
		int to = urlStr.indexOf("!/");
		return urlStr.substring(from, to);
	}

}
