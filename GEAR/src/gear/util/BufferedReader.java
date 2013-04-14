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
