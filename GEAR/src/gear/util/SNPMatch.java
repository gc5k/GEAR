package gear.util;

import gear.Parameter;

public class SNPMatch {

	public static boolean IsBiallelic(char a1, char a2) {
		boolean f = true;
		if (a1 == Parameter.INSTANCE.missing_allele.charAt(0) || a1 == Parameter.INSTANCE.missing_allele.charAt(0)) {
			f = false;
		}
		return f;
	}

	public static String Flip(String a) {
		String f = "A";
		if(a.compareTo("A") == 0) {
			f = "T";
		} else if(a.compareTo("C") == 0) {
			f = "G";
		} else if(a.compareTo("G") == 0) {
			f = "C";
		} else if (a.compareTo("T") == 0) {
			f = "A";
		} else {
			f = "N";
		}
		return f;		
	}
	
	public static char Flip(char a) {
		char f = 'A';
		if(a == 'A') {
			f = 'T';
		} else if(a == 'C') {
			f = 'G';
		} else if(a == 'G') {
			f = 'C';
		} else if (a == 'T') {
			f = 'A';
		} else {
			f = 'N';
		}
		
		return f;
	}
	
	public static boolean Confusion(String a, String b) {
		boolean f = false;
		char a1 = a.charAt(0);
		char a2 = b.charAt(0);
		if(a1>a2) {
			char t = a1;
			a1 = a2;
			a2 = t;
		}

		if(a1=='A' && a2=='T') {
			f = true;
		}
		if(a1=='C' && a2 =='G') {
			f = true;
		}
		
		return f;
	}

	
	public static boolean Confusion(char a1, char a2) {
		boolean f = false;
		if(a1>a2) {
			char t = a1;
			a1 = a2;
			a2 = t;
		}
		
		if(a1=='A' && a2=='T') {
			f = true;
		}
		if(a1=='C' && a2 =='G') {
			f = true;
		}
		
		return f;
	}
	
	public static boolean IsBiallelic(char a1, char a2, char b1, char b2) {
		char La1 = a1;
		char La2 = a2;
		char Lb1 = b1;
		char Lb2 = b2;

		boolean f = true;
		if((La1 - La2) > 0) {
			char t = La1;
			La1 = La2;
			La2 = t;
		}
		StringBuilder a = new StringBuilder(La1);
		a.append(La2);
		String g1 = a.toString();
		
		if((Lb1-Lb2) > 0) {
			char t = Lb1;
			Lb1 = Lb2;
			Lb2 = t;
		}
		StringBuilder b = new StringBuilder(Lb1);
		b.append(Lb2);
		String g2 = b.toString();

		if(La1!=La2 && Lb1!=Lb2) { //both bialelic
			if (g1.compareTo("AC") == 0) {
				if (g2.compareTo("AC") == 0 || g2.compareTo("GT") == 0) {
					return true;
				}
			} else if (g1.compareTo("AG") == 0) {
				if (g2.compareTo("AG") == 0 || g2.compareTo("CT") == 0) {
					return true;
				}
			} else if (g1.compareTo("AT") == 0) {
				if (g2.compareTo("AT") == 0) {
					return true;
				}
			} else if (g1.compareTo("CG") == 0) {
				if (g2.compareTo("CG") == 0) {
					return true;
				}
			} else if (g1.compareTo("CT") == 0) {
				if (g2.compareTo("CT") == 0 || g2.compareTo("AG") == 0) {
					return true;
				}
			} else {
				if (g2.compareTo("GT") == 0 || g2.compareTo("AC") == 0) {
					return true;
				}
			}
		} else {
			return false;
		}
		return f;
	}
}
