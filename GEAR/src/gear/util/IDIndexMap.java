package gear.util;

import java.util.ArrayList;
import java.util.HashMap;

public class IDIndexMap<E>
{
	public boolean add(E id)
	{
		if (id2idx.containsKey(id))
		{
			return false;
		}
		id2idx.put(id, idArray.size());
		idArray.add(id);
		return true;
	}
	
	public int getIndex(E id)
	{
		Integer idx = id2idx.get(id);
		return idx == null ? -1 : idx;
	}
	
	public E getID(int idx)
	{
		return idArray.get(idx);
	}
	
	public int getNumberOfEntries()
	{
		if (id2idx.size() != idArray.size())
		{
			Logger.fatalInternalBug("Inconsistent state in IDIndexMap.");
		}
		return idArray.size();
	}
	
	public void swapEntries(int idx1, int idx2)
	{
		E id1 = idArray.get(idx1), id2 = idArray.get(idx2);
		id2idx.put(id1, idx2);
		id2idx.put(id2, idx1);
		idArray.set(idx2, id1);
		idArray.set(idx1, id2);
	}
	
	private HashMap<E, Integer> id2idx = new HashMap<E, Integer>();
	private ArrayList<E> idArray = new ArrayList<E>();
}
