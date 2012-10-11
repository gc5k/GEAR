package gui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

public class DataFile {

	//

	private String format;
	private File file;
	private boolean hasHeader;

	//

	public String getFormat() {
		return format;
	}

	public void setFormat(String format) {
		this.format = format;
	}

	public File getFile() {
		return file;
	}

	public void setFile(File file) {
		this.file = file;
	}

	public boolean hasFile() {
		if (file == null) {
			return false;
		}
		if (!file.exists()) {
			return false;
		}
		if (!file.canRead()) {
			return false;
		}
		return true;
	}

	public boolean hasHeader() {
		return hasHeader;
	}

	public void setHasHeader(boolean hasHeader) {
		this.hasHeader = hasHeader;
	}

	public Vector[] getData(long rowLimit) throws IOException {
		if (!hasFile()) {
			return null;
		}
		Vector dataVector = new Vector();
		Vector headVector = new Vector();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		if (hasHeader) {
			line = br.readLine();
			String[] cells = line.split("\\s");
			headVector.addAll(Arrays.asList(cells));
		}
		long rowIndex = 0;
		int colMax = 0;
		while ((line = br.readLine()) != null) {
			if (rowIndex++ == rowLimit) {
				break;
			}
			String[] cells = line.split("\\s");
			colMax = Math.max(colMax, cells.length);
			Vector rowVector = new Vector(Arrays.asList(cells));
			dataVector.addElement(rowVector);
		}
		br.close();
		if (!hasHeader) {
			String head = "A";
			for (int i = 0; i < colMax; i++) {
				headVector.addElement(head);
				head = Radix26.next(head);
			}
		}
		return new Vector[] { dataVector, headVector };
	}

	public static void main(String[] args) {
	}

}
