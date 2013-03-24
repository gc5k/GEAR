package hpc;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import parameter.Parameter;

public class HPC {
	String[] p;

	public HPC(String[] par) {
		p = par;
	}

	public void GenScript() {
		StringBuilder sb = new StringBuilder();
		sb.append(Parameter.INSTANCE.name);
		sb.append(".sh");
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(sb.toString());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		StringBuilder shell = new StringBuilder();
		shell.append("#$ -cwd" + "\n");
		shell.append("#$ -l vf=" + Parameter.INSTANCE.ram + "\n");
		shell.append("#$ -N " + Parameter.INSTANCE.name + "\n");
		shell.append("#$ -m eas" + "\n");
		shell.append("#$ -M " + Parameter.INSTANCE.email + "\n");

		shell.append("java -jar -Xmx" + Parameter.INSTANCE.ram + " ");
		shell.append(getClass().getProtectionDomain().getCodeSource().getLocation().getPath() + " ");
		for (int i = 0; i < p.length; i++) {
			if (p[i].compareTo("--shell") == 0 || p[i].compareTo("--qsub") == 0)
				continue;
			shell.append(p[i] + " ");
		}
		shell.append("\n");
		pw.append(shell);
		pw.close();
		if (Parameter.INSTANCE.qsubFlag) {
			Runtime rt = Runtime.getRuntime();
			try {
				rt.exec("qsub " + sb.toString());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
}
