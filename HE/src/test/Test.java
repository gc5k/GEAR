package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Calendar;


import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;

import util.NewIt;

import parameter.Parameter;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Test {

	public static StringBuffer LOG = new StringBuffer();
	public static boolean fileFlag = false;
	public static boolean bfileFlag = false;
	public static void main(String[] args) throws IOException {

		System.err.print(Parameter.version);
		String[] scmd = script(args);
		printCommandLine(scmd);
		Parameter.INSTANCE.commandListener(scmd);

		PLINKParser pp = null;
		if (Parameter.fileOption) {
			pp = new PLINKParser(Parameter.pedfile, Parameter.mapfile);
		}
		if (Parameter.bfileOption) {
			pp = new PLINKBinaryParser(Parameter.bedfile, Parameter.bimfile, Parameter.famfile);
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp.Parse();
	
		printLog();
	}

	public static void printCommandLine(String[] args) {

		Calendar calendar = Calendar.getInstance();
		System.err.println("The analysis was implemented at: " + calendar.getTime() + "\n");
		LOG.append(Parameter.version);
		LOG.append("The analysis was implemented at: " + calendar.getTime() + "\n\n");
		System.err.println("The command line in effect: ");
		LOG.append("The command line in effect: \n");
		for (int i = 0; i < args.length; i++) {
			System.err.print(args[i] + " ");
			LOG.append(args[i] + " ");
		}
		System.err.println("\n");
		LOG.append("\n\n");

	}

	public static void printLog() {
		StringBuilder sb = new StringBuilder(Parameter.out);
		sb.append(".log");
		PrintStream pw = null;
		try {
			pw = new PrintStream(sb.toString());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Calendar calendar = Calendar.getInstance();
		System.err.println("\nThe analysis was finished at: " + calendar.getTime() + "\n");
		System.err.println("These above messages were printed into " + sb.toString() + ".\n");
		LOG.append("\nThe analysis was finished at: " + calendar.getTime() + "\n");
		LOG.append("\n");

		pw.append(LOG.toString());
		pw.close();
	}

	public static String[] script(String[] args) {

		String sf = null;
		boolean scriptFlag = false;
		String[] scmd;
		int c = 0;
		for (int i = 0; i < args.length; i++) {
			if (args[i].compareTo("--script") == 0) {
				c = i;
				sf = args[i + 1];
				scriptFlag = true;
				continue;
			} else if (i == c+1) {
				continue;
			}
			if (args[i].compareTo("--bfile")==0 || args[i].compareTo("--bed") == 0 || args[i].compareTo("--bim") == 0 || args[i].compareTo("--fam") == 0) {
				bfileFlag = true;
			}
			if (args[i].compareTo("--file") == 0 || args[i].compareTo("--ped") == 0 || args[i].compareTo("--map") == 0) {
				fileFlag = true;
			}
		}

		if (bfileFlag && fileFlag) {
			System.err.println("specified both text and binary format files.");
			Test.LOG.append("specified both text and binary format files.\n");
			Test.printLog();
			System.exit(0);
		}

		if (scriptFlag) {
			File f = new File(sf);
			if (!f.exists()) {
				System.err.print("could not find " + sf + ".");
				Test.LOG.append("could not find " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(f));
			} catch (IOException E) {
				System.err.print("could not open " + sf + ".");
				Test.LOG.append("could not open " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}

			ArrayList<String> cmd = NewIt.newArrayList();
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split("\\s+");
					for (int i = 0; i < s.length; i++) {
						cmd.add(s[i]);
					}
				}
				reader.close();
			} catch (IOException E) {
				System.err.println("bad lines in " + sf + ".");
				Test.LOG.append("bad lines in " + sf + ".\n");
				Test.printLog();
				System.exit(0);
			}

			scmd = (String[]) cmd.toArray(new String[0]);
		} else {
			scmd = args;
		}
		return scmd;
	}

	public static void savecmd(String[] args) {
		StringBuffer sb = new StringBuffer();

		sb.append(Parameter.out);
		PrintStream ps = null;
		try {
			ps = new PrintStream(sb.toString());
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		for (String s:args) {
			ps.print(s + " ");
		}
		ps.close();
	}
}
