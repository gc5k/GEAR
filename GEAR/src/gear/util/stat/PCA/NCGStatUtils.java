package gear.util.stat.PCA;

import org.apache.commons.math.stat.StatUtils;

/*
 * Various simple statistical methods for the ncg.statistics package
 * 
 * Author : Peter Foley
 *
 */

public class NCGStatUtils
{
	
	//*************************************************************************
	// Name    : standardize
	// 
	// Purpose : returns z scores of data array
	// 
	// Notes   : returns a zero length array if the the input is null
	//           This method does NOT modify the input array
	// 
	//*************************************************************************
	public static double[] standardize(final double[] data)
	{
		
		double[] dataTransform = null;
		
		if ( data != null )
		{
		
			double dataMean = StatUtils.mean(data);
			double dataStdDev = Math.sqrt(StatUtils.variance(data));
			
			dataTransform = new double[data.length];
			
			for (int i = 0; i < data.length; i++ )
			{
				dataTransform[i] = (data[i] - dataMean) / dataStdDev;
			}
		} 
		else
		{
			dataTransform = new double[0];
		}
		return dataTransform;
	}
	
	//*************************************************************************
	// Name    : standardize
	// 
	// Purpose : returns z scores of columns of data array
	// 
	// Notes   : rowOrder is set to true if the first dimension of the input
	//           array contains the rows (observations). If rowOrder is set to
	//           false the first dimension of the input array refers to the 
	//           columns (Variables). 
	//           returns a zero length 2d array of doubles if the the input is null
	//           This method does NOT modify the input array
	// 
	//*************************************************************************
	public static double[][] standardize(double[][] data, final boolean rowOrder)
	{
		
		double[][] dataZScores = null;
		if (data != null)
		{
			
			// if rowOrder is true then the first dimension refers to the 
			// rows - need to transpose this matrix so that the first dimension
			// refers to the columns (variables)
			if (rowOrder == true)
			{ 
					data = transpose(data);
			} 
			
			int numVars = data.length;
			
			dataZScores = new double[numVars][0];
			
			// for each column in the data set, standardize it
			for ( int j=0; j<data.length; j++ )
			{
				dataZScores[j] = standardize(data[j]);
			}
			
			// if rowOrder is true then transpose the dataZScores array
			// to make the first dimension refer to the rows again - 
			// in other words it should match the input data array
			if (rowOrder == true)
			{
				dataZScores = transpose(dataZScores);
			}
		}
		else
		{
			dataZScores = new double[0][0];
		}
		return dataZScores;
	}
	
	//*************************************************************************
	// Name    : transpose
	// 
	// Purpose : transpose a 2d array of doubles 'data'
	// 
	// Notes   : returns a zero length array if the the input is null
	//           This method does NOT modify the input array
	// 
	//*************************************************************************
	public static double[][] transpose(double[][] data)
	{
				
		double[][] dataTransposed = new double[0][0];

		if ( data != null )
		{
			int numRows = data.length;			
			if ( numRows > 0 )
			{	
				int numCols = data[0].length;
				dataTransposed = new double[numCols][numRows];
				for ( int i = 0; i < numRows; i++ )
				{
					for ( int j = 0; j < numCols; j++ )
					{
						dataTransposed[j][i] = data[i][j];
					}
				}
			}
		}	
		return dataTransposed;
	}
}
