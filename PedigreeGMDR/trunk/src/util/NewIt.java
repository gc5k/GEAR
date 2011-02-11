package util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.ArrayList;

public class NewIt {
	public static <K, V> HashMap<K, V> newHashMap() {
		return new HashMap<K, V>();
	}
	public static <E> HashSet<E> newHashSet() {
		return new HashSet<E>();
	}
	public static <K, V> Hashtable<K, V> newHashtable() {
		return new Hashtable<K, V>();
	}
	public static <K, V> TreeMap<K, V> newTreeMap() {
		return new TreeMap<K, V>();
	}
	public static <E> ArrayList<E> newArrayList() {
		return new ArrayList<E>();
	}
}
