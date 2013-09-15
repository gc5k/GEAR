package gear.hpc;

import gear.CommandArguments;

public class HpcCommandArguments extends CommandArguments
{
	public String getJobName()
	{
		return jobName;
	}
	
	public void setJobName(String jobName)
	{
		this.jobName = jobName;
	}
	
	public int getRam()
	{
		return ram;
	}
	
	public void setRam(int ram)
	{
		this.ram = ram;
	}
	
	public String getRamUnit()
	{
		return "G";
	}
	
	public String getEmail()
	{
		return email;
	}
	
	public void setEmail(String email)
	{
		this.email = email;
	}
	
	public boolean getIsSubmit()
	{
		return isSubmit;
	}
	
	public void setIsSubmit(boolean isSubmit)
	{
		this.isSubmit = isSubmit;
	}
	
	public void setNonHpcOptionsAndArguments(String[] nonHpcOptsAndArgs)
	{
		this.nonHpcOptsAndArgs = new String[nonHpcOptsAndArgs.length];
		System.arraycopy(nonHpcOptsAndArgs, 0, this.nonHpcOptsAndArgs, 0, nonHpcOptsAndArgs.length);
	}
	
	public String[] getNonHpcOptionsAndArguments()
	{
		String[] result = new String[nonHpcOptsAndArgs.length];
		System.arraycopy(nonHpcOptsAndArgs, 0, result, 0, nonHpcOptsAndArgs.length);
		return result;
	}
	
	private String jobName;
	private int ram;
	private String email;
	private boolean isSubmit;
	private String[] nonHpcOptsAndArgs;
}
