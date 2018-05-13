import java.io.*;
import java.util.Calendar;

public class GenAboutInfo
{
	public static void main(String[] args) throws Exception 
	{
		String packageName = args.length > 0 ?  args[0] : "gear";
		String packageDeclaration = "package " + packageName + ";\n";
		String folder = "src/" + packageName.replace('.', '/') + "/";

		Process proc;
		BufferedReader reader;

		// Get branch name.
		proc = Runtime.getRuntime().exec("git rev-parse --abbrev-ref HEAD");
		reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String branch = reader.readLine();

		// Get commit SHA1 revision.
		proc = Runtime.getRuntime().exec("git rev-list head --max-count 1");
		reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String commit = reader.readLine();
		
		// Has uncommitted changes?
		proc = Runtime.getRuntime().exec("git status --porcelain");
		reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		boolean hasUncommittedChanges = reader.readLine() != null;

		// Write AboutInfo.java.
		FileWriter writer = new FileWriter(folder + "AboutInfo.java");
		writer.write(packageDeclaration);
		writer.write("\n");
		writer.write("public class AboutInfo\n");
		writer.write("{\n");
		writer.write("\tpublic static final String WELCOME_MESSAGE =\n");
		writer.write("\t\t\"****************************************************************\\n\" +\n");
		writer.write("\t\t\" GEAR [GEnetic Analysis Repository]\\n\" +\n");
		writer.write("\t\t\" (C) 2013 Guo-Bo Chen, Zhi-Xiang Zhu\\n\" +\n");
		writer.write("\t\t\" Branch " + branch + "\\n\" +\n");
		writer.write("\t\t\" Commit " + commit + "\\n\" +\n");
		writer.write("\t\t\" Built on " + Calendar.getInstance().getTime() +
				(hasUncommittedChanges ? " with uncommitted changes" : "") + "\\n\" +\n");
		writer.write("\t\t\"****************************************************************\\n\";\n");
		writer.write("}");
		writer.close();
	}
}