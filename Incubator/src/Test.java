import org.apache.commons.lang3.ArrayUtils;

import java.io.FileWriter;
import java.io.BufferedWriter;

public class Test {

	public static void main(String[] args) {

		try {
			// Create file
			FileWriter fstream = new FileWriter("out.txt", true);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write("Hello Java");
			// Close the output stream
			out.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}
}
