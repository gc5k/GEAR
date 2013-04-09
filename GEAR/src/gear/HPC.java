package gear;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import gear.util.Logger;

public class HPC {

	public static void genScript(String[] args) {
		StringBuilder sb = new StringBuilder();
		sb.append(Parameter.INSTANCE.getHpcParameter().getName());
		sb.append(".sh");
		
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(sb.toString());
		} catch (FileNotFoundException e) {
			Logger.handleException(e, "Cannot create the script file '" + sb.toString() + "'.");
		}

		pw.println("#$ -cwd");
		pw.println("#$ -l vf=" + Parameter.INSTANCE.getHpcParameter().getRam());
		pw.println("#$ -N " + Parameter.INSTANCE.getHpcParameter().getName());
		pw.println("#$ -m eas");
		pw.println("#$ -M " + Parameter.INSTANCE.getHpcParameter().getEmail());

		pw.print("java -jar -Xmx" + Parameter.INSTANCE.getHpcParameter().getRam() + " ");
		pw.print(HPC.class.getProtectionDomain().getCodeSource().getLocation().getPath() + " ");
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			if (arg.equals("--shell") || arg.equals("--qsub")) {
				continue;
			}
			if (arg.equals("--email") || arg.equals("--ram") || arg.equals("--name")) {
				++i;
				continue;
			}
			pw.print(arg + " ");
		}
		pw.println();
		
		pw.close();
		
		if (Parameter.INSTANCE.getHpcParameter().isQsubSet()) {
			Runtime rt = Runtime.getRuntime();
			String cmd = "qsub " + sb.toString();
			try {
				rt.exec(cmd);
			} catch (IOException e) {
				Logger.handleException(e, "Failed to execute command '" + cmd + "'.");
			}
		}
	}
	
}
