package test;

import java.io.IOException;

public class Test2 {

	public static void main(String[] args) throws IOException {
		test2();
	}

	static void test1() throws IOException {
		Test.main(new String[] { "--bfile", "C:/Users/Lei/Documents/GMDR/poly" });
	}

	static void test2() throws IOException {
		Test.main(new String[] { "--bfile", "C:/Users/Lei/Documents/GMDR/example2011" });
	}
}
