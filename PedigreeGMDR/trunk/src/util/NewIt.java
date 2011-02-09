package util;

import java.util.HashMap;
import java.util.ArrayList;
public class NewIt {
	public static <K, V> HashMap<K, V> newHashMap() {
		return new HashMap<K, V>();
	}
	public static <E> ArrayList<E> newArrayList() {
		return new ArrayList<E>();
	}
}
