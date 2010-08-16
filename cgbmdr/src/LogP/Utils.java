package LogP;

public class Utils {

	public static String makeKey(int[] element) {
		String key = new String();
		for (int i = 0; i < element.length; i++) {
			key = key + Integer.toString(element[i]);
			if (i < element.length-1) {
				key = key + "-";
			}
		}
		return key;
	}
}
