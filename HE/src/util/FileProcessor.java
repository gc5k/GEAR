package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.zip.GZIPInputStream;

import parameter.Parameter;

import test.Test;

public class FileProcessor {
	public static BufferedReader FileOpen(String file) {
		File f = new File(file);
		if (!f.exists()) {
			System.err.println("could not open " + file + ".");
			Test.LOG.append("could not open " + file + ".\n");
			Test.printLog();
			System.exit(0);
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
		} catch (IOException E) {
			System.err.println("could not read " + file + ".");
			System.exit(0);
		}
		return reader;
	}
	
	public static PrintStream CreatePrintStream(String file) {

		PrintStream ps = null;
		try {
			ps = new PrintStream(file);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return ps;

	}
	
	public static BufferedReader ZipFileOpen(String file) {
		FileInputStream fin = null;
		try {
			fin = new FileInputStream(file);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		GZIPInputStream gzis = null;
		try {
			gzis = new GZIPInputStream(fin);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		InputStreamReader xover = new InputStreamReader(gzis);
		BufferedReader is = new BufferedReader(xover);
		return is;
	}
}
