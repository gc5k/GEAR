package gui;

import java.util.Collections;
import java.util.Vector;

public class Radix26 {

	public static String next(String s) {
		Vector<Character> v = new Vector<Character>();
		for (char c : s.toCharArray()) {
			v.addElement(c);
		}
		Collections.reverse(v);
		for (int i = 0; i <= v.size(); i++) {
			if (i == v.size()) {
				v.addElement('A');
				break;
			}
			char c = v.get(i);
			if (c == 'Z') {
				v.setElementAt('A', i);
				continue;
			} else {
				v.setElementAt(++c, i);
				break;
			}
		}
		Collections.reverse(v);
		StringBuilder sb = new StringBuilder();
		for (Character c : v) {
			sb.append(c);
		}
		return sb.toString();
	}

	public static void main(String[] args) {
		String s = "A";
		for (int i = 0; i < 1000; i++) {
			System.out.println(i + "\t" + s);
			s = next(s);
		}
	}

}
