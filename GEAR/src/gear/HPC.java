package gear;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.logging.Level;

import gear.util.Logger;

public class HPC {

	public static void GenScript(String[] args) {
		StringBuilder sb = new StringBuilder();
		sb.append(Parameter.INSTANCE.name);
		sb.append(".sh");
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(sb.toString());
		} catch (FileNotFoundException e) {
			Logger.printUserError("Cannot create the script file '" + sb.toString() + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}

		StringBuilder shell = new StringBuilder();
		shell.append("#$ -cwd" + "\n");
		shell.append("#$ -l vf=" + Parameter.INSTANCE.ram + "\n");
		shell.append("#$ -N " + Parameter.INSTANCE.name + "\n");
		shell.append("#$ -m eas" + "\n");
		shell.append("#$ -M " + Parameter.INSTANCE.email + "\n");

		shell.append("java -jar -Xmx" + Parameter.INSTANCE.ram + " ");
		shell.append(HPC.class.getProtectionDomain().getCodeSource().getLocation().getPath() + " ");
		for (int i = 0; i < args.length; i++) {
			if (args[i].compareTo("--shell") == 0 || args[i].compareTo("--qsub") == 0)
				continue;
			shell.append(args[i] + " ");
		}
		shell.append("\n");
		pw.append(shell);
		pw.close();
		if (Parameter.INSTANCE.qsubFlag) {
			Runtime rt = Runtime.getRuntime();
			String cmd = "qsub " + sb.toString();
			try {
				rt.exec(cmd);
			} catch (IOException e) {
				Logger.printUserError("Failed to execute command '" + cmd + "'.");
				Logger.printUserError("Exception Message: " + e.getMessage());
				Logger.getDevLogger().log(Level.SEVERE, "Executing qsub command", e);
				System.exit(1);
			}
		}
	}
	
}
