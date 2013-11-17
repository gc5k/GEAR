package gear.util;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.logging.Level;
import java.util.zip.GZIPInputStream;

public class BufferedReader
{
	public static BufferedReader openTextFile(String fileName, String fileType)
	{
		java.io.BufferedReader reader = null;
		try
		{
			reader = new java.io.BufferedReader(new java.io.FileReader(fileName));
		}
		catch (Exception e)
		{
			Logger.handleException(e, "Cannot open the " + fileType + " file '" + fileName + "'.");
		}
		return new BufferedReader(reader, fileName, fileType);
	}
	
	public static BufferedReader openGZipFile(String fileName, String fileType)
	{
		java.io.FileInputStream fileInStrm = null;
		try
		{
			fileInStrm = new java.io.FileInputStream(fileName);
		}
		catch (FileNotFoundException e)
		{
			Logger.handleException(e, "File '" + fileName + "' does not exist.");
		}
		
		GZIPInputStream gzip = null;
		try
		{
			gzip = new GZIPInputStream(fileInStrm);
		}
		catch (IOException e)
		{
			Logger.handleException(e, "Cannot open the archive '" + fileName + "'.");
		}
		
		java.io.InputStreamReader inStrmReader = new java.io.InputStreamReader(gzip);
		java.io.BufferedReader bufferedReader = new java.io.BufferedReader(inStrmReader);
		return new BufferedReader(bufferedReader, fileName, fileType);
	}
	
	private BufferedReader(java.io.BufferedReader reader, String fileName, String fileType)
	{
		innerReader = reader; 
		this.fileName = fileName;
		this.fileType = fileType;
		curLineNum = 1;
	}

	public void close()
	{
		try
		{
			innerReader.close();
		}
		catch (java.io.IOException e)
		{
			String msg = "An I/O exception occurred when closing the " + fileType + " file '" + fileName + "'.";
			Logger.printUserError(msg);
			Logger.getDevLogger().log(Level.WARNING, msg, e);
		}
	}

	public String readLine()
	{
		String line = null;
		
		try
		{
			line = innerReader.readLine();
		}
		catch (java.io.IOException e)
		{
			String msg = "";
			msg += "An I/O exception occurred when reading to line " + curLineNum;
			msg += " of the " + fileType + " file '" + fileName + "'.";
			Logger.handleException(e, msg);
		}
		
		if (line != null)
		{
			++curLineNum;
		}
		
		return line;
	}

	public String readNonEmptyLine()
	{
		String line = null;
		while ((line = readLine()) != null)
		{
			line = line.trim();
			if (!line.isEmpty())
			{
				break;
			}
		}
		return line;
	}

	/**
	 * @return split tokens on the next non-empty line
	 */
	public String[] readTokens()
	{
		String[] tokens = null;

		while (true)
		{
			String line = readNonEmptyLine();

			if (line == null)
			{
				tokens = null;
				break;
			}

			tokens = line.trim().split("\\s+");

			if (tokens.length > 0)
			{
				break;
			}
		}

		return tokens;
	}

	public String[] readTokens(int expectedNumCols)
	{
		String[] tokens = readTokens();
		if (tokens != null && tokens.length != expectedNumCols)
		{
			String msg = "";
			msg += "The format of the " + fileType + " file '" + fileName + "' is incorrect: ";
			msg += "The file should consists of " + expectedNumCols + " column(s) at line " + (curLineNum - 1);
			msg += ", but this line actually contains " + tokens.length + " column(s).";
			Logger.printUserError(msg);
			System.exit(1);
		}
		return tokens;
	}
	
	public String[] readTokensAtLeast(int minNumCols)
	{
		String[] tokens = readTokens();
		if (tokens != null && tokens.length < minNumCols)
		{
			String msg = "";
			msg += "The format of the " + fileType + " file '" + fileName + "' is incorrect: ";
			msg += "The file should consists of at least " + minNumCols + " column(s) at line " + (curLineNum - 1);
			msg += ", but this line actually contains only " + tokens.length + " column(s).";
			Logger.printUserError(msg);
			System.exit(1);
		}
		return tokens;
	}

	public String getFileName()
	{
		return fileName;
	}

	public int getCurLineNum()
	{
		return curLineNum;
	}
	
	public void error(int lineNo, String msg)
	{
		Logger.printUserError(fileType + " file '" + getFileName() + "', line " + lineNo + ": " + msg);
		System.exit(1);
	}
	
	public void error(String msg)
	{
		error(getCurLineNum(), msg);
	}
	
	public void errorPreviousLine(String msg)
	{
		error(getCurLineNum() - 1, msg);
	}
	
	public void warning(int lineNo, String msg)
	{
		Logger.printUserWarning(fileType + " file '" + getFileName() + "', line " + lineNo + ": " + msg);
	}
	
	public void warning(String msg)
	{
		warning(getCurLineNum(), msg);
	}
	
	public void warningPreviousLine(String msg)
	{
		warning(getCurLineNum() - 1, msg);
	}

	private java.io.BufferedReader innerReader;

	private int curLineNum;

	private String fileName;

	private String fileType;
}
