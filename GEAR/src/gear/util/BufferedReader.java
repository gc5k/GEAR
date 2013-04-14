package gear.util;

public class BufferedReader {
	
	public BufferedReader(String fileName, String fileType) {
		try {
			innerReader = new java.io.BufferedReader(new java.io.FileReader(fileName));
		}
		catch (Exception e)
		{
			Logger.handleException(e, "Cannot open the " + fileType + " file '" + fileName + "'.");
		}
		this.fileName = fileName;
		this.fileType = fileType;
		curLineNum = 1;
	}
	
	public String readLine() {
		String line = null;
		try {
			line = innerReader.readLine();
		}
		catch (java.io.IOException e) {
			Logger.handleException(e, "An I/O exception occurred when reading to line " + curLineNum + " of the " + fileType + " file '" + fileName + "'.");
		}
		if (line != null) {
			++curLineNum;
		}
		return line;
	}
	
	/**
	 * @return split tokens on the next non-empty line
	 */
	public String[] readTokens() {
		String[] tokens = null;
		
		while (true) {
			String line = readLine();
			
			if (line == null) {
				tokens = null;
				break;
			}
			
			tokens = line.trim().split("\\s+");
			
			if (tokens.length > 0) {
				break;
			}
		}
		
		return tokens;
	}
	
	public String[] readTokens(int expectedNumCols) {
		String[] tokens = readTokens();
		if (tokens != null && tokens.length != expectedNumCols) {
			Logger.printUserError(
				"The format of the " + fileType + " file '" + fileName + "' is incorrect: " +
				"A " + fileType + " file should consists of " + expectedNumCols + " columns, but line " + curLineNum + " contains " + tokens.length + " column(s)."
			);
			System.exit(1);
		}
		return tokens;
	}
	
	public String getFileName() {
		return fileName;
	}
	
	public int getCurLineNum() {
		return curLineNum;
	}
	
	private java.io.BufferedReader innerReader;
	
	private int curLineNum;
	
	private String fileName;
	
	private String fileType;

}
