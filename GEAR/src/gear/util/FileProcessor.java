package gear.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.zip.GZIPInputStream;

public class FileProcessor {
	public static BufferedReader FileOpen(String file) {
		File f = new File(file);
		if (!f.exists()) {
			Logger.printUserError("File '" + file + "' does not exist.");
			System.exit(1);
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(new File(file)));
		} catch (IOException e) {
			Logger.printUserError("Cannot open '" + file + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}
		return reader;
	}
	
	public static PrintStream CreatePrintStream(String file) {
		PrintStream ps = null;
		try {
			ps = new PrintStream(file);
		} catch (FileNotFoundException e) {
			Logger.printUserError("File '" + file + "' does not exist.");
			System.exit(1);
		}
		return ps;
	}
	
	public static BufferedReader ZipFileOpen(String file) {
		FileInputStream fin = null;
		try {
			fin = new FileInputStream(file);
		} catch (FileNotFoundException e) {
			Logger.printUserError("File '" + file + "' does not exist.");
			System.exit(1);
		}
		GZIPInputStream gzis = null;
		try {
			gzis = new GZIPInputStream(fin);
		} catch (IOException e) {
			Logger.printUserError("Cannot open the archive '" + file + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}
		InputStreamReader xover = new InputStreamReader(gzis);
		BufferedReader is = new BufferedReader(xover);
		return is;
	}
}
