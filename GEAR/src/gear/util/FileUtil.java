package gear.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FileUtil
{
	public static BufferedReader FileOpen(String file)
	{
		File f = new File(file);
		if (!f.exists())
		{
			Logger.printUserError("File '" + file + "' does not exist.");
			System.exit(1);
		}
		BufferedReader reader = null;
		try
		{
			reader = new BufferedReader(new FileReader(new File(file)));
		} catch (IOException e)
		{
			Logger.handleException(e, "Cannot open '" + file + "'.");
		}
		return reader;
	}

	public static PrintStream CreatePrintStream(String file)
	{
		PrintStream ps = null;
		try
		{
			return new PrintStream(file);
		}
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "Unable to create file '" + file + "'.");
		}
		return ps;
	}

	public static BufferedReader ZipFileOpen(String file)
	{
		FileInputStream fin = null;
		try
		{
			fin = new FileInputStream(file);
		} catch (FileNotFoundException e)
		{
			Logger.handleException(e, "File '" + file + "' does not exist.");
		}
		GZIPInputStream gzis = null;
		try
		{
			gzis = new GZIPInputStream(fin);
		} catch (IOException e)
		{
			Logger.handleException(e, "Cannot open the archive '" + file + "'.");

		}
		InputStreamReader xover = new InputStreamReader(gzis);
		BufferedReader is = new BufferedReader(xover);
		return is;
	}

	public static BufferedWriter ZipFielWriter(String file)
	{
		BufferedWriter writer = null;
		try
		{
			GZIPOutputStream zip = new GZIPOutputStream(new FileOutputStream(
					new File(file)));
			writer = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
		} catch (IOException e)
		{
			Logger.handleException(e, "Cannot create the archive '" + file
					+ "'.");
		}
		return writer;
	}

	public static void exists(String file)
	{
		java.io.File f = new java.io.File(file);
		if (!f.exists())
		{
			Logger.printUserError("File '" + file + "' does not exist.");
			System.exit(1);
		}
	}

}
