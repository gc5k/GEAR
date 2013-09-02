package gear;

import java.util.NavigableSet;
import java.util.TreeSet;

public abstract class Command
{
	public abstract String getName();
	
	public boolean hasAlias(String alias)
	{
		return aliases.contains(alias);
	}
	
	protected void addAlias(String alias)
	{
		aliases.add(alias);
	}
	
	@Override
	public boolean equals(Object obj)
	{
		String s = null;
		
		if (obj instanceof String)
		{
			s = (String)obj;
		}
		else if (obj instanceof Command)
		{
			s = ((Command)obj).getName();
		}
		
		return s != null && (this.getName().equals(s) || this.hasAlias(s));
	}
	
	public abstract String getDescription();
	public abstract String getLongDescription();
	public abstract void printHelp();
	public abstract void execute(String[] args);
	
	private NavigableSet<String> aliases = new TreeSet<String>();
}
