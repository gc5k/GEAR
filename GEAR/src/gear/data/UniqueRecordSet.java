package gear.data;

import gear.util.Logger;

import java.util.ArrayList;
import java.util.Hashtable;

public class UniqueRecordSet<R extends UniqueRecord>
{
	public void put(R record)
	{
		if (recordMap.containsKey(record.getID()))
		{
			Logger.fatalInternalBug("'" + record.getID() + "' is added twice to UniqueRecordSet.");
		}
		recordMap.put(record.getID(), record);
		recordList.add(record);
	}
	
	public R get(String id)
	{
		return recordMap.get(id);
	}
	
	public R get(int index)
	{
		return recordList.get(index);
	}
	
	public boolean has(String id)
	{
		return get(id) != null;
	}
	
	public int size()
	{
		return recordList.size();
	}
	
	/*
	 * UniqueRecordSet is very similar with gear.util.IDIdxMap. But I can't
	 * reuse IDIdxMap, because it needs the map to be of type <R, Integer>.
	 */
	private Hashtable<String, R> recordMap = new Hashtable<String, R>();
	private ArrayList<R> recordList = new ArrayList<R>();
}
