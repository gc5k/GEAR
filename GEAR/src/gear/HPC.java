package gear;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import gear.util.Logger;

public class HPC
{
	public static void genScript(String[] args)
	{
		StringBuilder sb = new StringBuilder();
		sb.append(CmdArgs.INSTANCE.getHpcArgs().getName());
		sb.append(".sh");

		PrintWriter pw = null;
		try
		{
			pw = new PrintWriter(sb.toString());
		}
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "Cannot create the script file '" + sb.toString() + "'.");
		}

		pw.println("#$ -cwd");

		int G = Integer.parseInt(CmdArgs.INSTANCE.getHpcArgs().getRam().substring(0, CmdArgs.INSTANCE.getHpcArgs().getRam().length()-1)) + 2;

		pw.println("#$ -l vf=" + G + "G");
		pw.println("#$ -l h_vmem=" + G + "G");
		pw.println("#$ -N " + CmdArgs.INSTANCE.getHpcArgs().getName());
		pw.println("#$ -m eas");
		pw.println("#$ -M " + CmdArgs.INSTANCE.getHpcArgs().getEmail());

		pw.print("java -jar -Xmx" + CmdArgs.INSTANCE.getHpcArgs().getRam() + " ");
		pw.print(HPC.class.getProtectionDomain().getCodeSource().getLocation().getPath() + " ");
		for (int i = 0; i < args.length; i++)
		{
			String arg = args[i];
			if (arg.equals("--shell") || arg.equals("--qsub"))
			{
				continue;
			}
			if (arg.equals("--email") || arg.equals("--ram") || arg.equals("--name"))
			{
				++i;
				continue;
			}
			pw.print(arg + " ");
		}
		pw.println();

		pw.close();

		if (CmdArgs.INSTANCE.getHpcArgs().isQsubSet())
		{
			Runtime rt = Runtime.getRuntime();
			String cmd = "qsub " + sb.toString();
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
}
