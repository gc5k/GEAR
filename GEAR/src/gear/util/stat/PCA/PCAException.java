package gear.util.stat.PCA;

/*
 * Exception thrown by the PCA class
 */
public class PCAException extends Exception 
{
	// constructor signatures all match constructors of the Exception class
	
	public PCAException()
	{
		super();
	}
	
	public PCAException(String message)
	{
		super(message);
	}
	
	public PCAException(String message, Throwable cause)
	{
		super(message,cause);
	}
	
	public PCAException(Throwable cause)
	{
		super(cause);
	}
}