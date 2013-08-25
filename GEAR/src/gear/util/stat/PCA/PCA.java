package gear.util.stat.PCA;

/*
 * Implement Principal Components Analysis (PCA)
 * 
 * Author : Peter Foley, 07.10.2010
 * 
 * Method taken from 'Multivariate Statistical Methods : A Primer'
 * by Bryan F.J. Manly (3rd edition), Chapman & Hall/CRC 2005.
 * See chapter 6 for details
 */

import java.util.logging.Logger;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.EigenDecomposition;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.InvalidMatrixException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.Covariance;

public class PCA {
	
	// attributes
	private transient RealMatrix observations = null;
	private transient RealMatrix principalComponents = null;
	private transient RealMatrix eigenVectors = null;
	private transient RealVector eigenValues = null;
	
	//logger object
	protected final static Logger logger = 
		Logger.getLogger(PCA.class.getName());
	
	// empty constructor
	public PCA() {};
	
	// reset all member variables
	public void reset() {
		
		observations = null;
		principalComponents = null;
		eigenVectors = null;
		eigenValues = null;
		
		// run the java garbage collector
		Runtime.getRuntime().gc();
	}
	
	//*************************************************************************
	// Name    : validateObservations
	// 
	// Purpose : check to see if the observations attribute has been set
	// 
	// Notes   : throws a new PCAException if observations is not set
	// 
	//*************************************************************************
	public void validateObservations() throws PCAException {	
		if ( observations == null ) {
			throw new PCAException("input observations not set");
		}
	}
	
	//*************************************************************************
	// Name    : validatePrincipalComponents
	// 
	// Purpose : check to see if the principalComponents attribute has been set
	// 
	// Notes   : throws a new PCAException if principalComponents is not set
	// 
	//*************************************************************************
	public void validatePrincipalComponents() throws PCAException {	
		if ( principalComponents == null ) {
			throw new PCAException("output principalcomponents not set");
		}
	}
	
	//*************************************************************************
	// Name    : validateEigenValues
	// 
	// Purpose : check to see if the eigenValues attribute has been set
	// 
	// Notes   : throws a new PCAException if eigenValues is not set
	// 
	//*************************************************************************
	public void validateEigenValues() throws PCAException {	
		if ( eigenValues == null ) {
			throw new PCAException("output eigen values not set");
		}
	}
	
	//*************************************************************************
	// Name    : validateEigenVectors
	// 
	// Purpose : check to see if the eigenVectors attribute has been set
	// 
	// Notes   : throws a new PCAException if eigenVectors is not set
	// 
	//*************************************************************************
	public void validateEigenVectors() throws PCAException {	
		if ( eigenVectors == null ) {
			throw new PCAException("eigen vectors not set");
		}
	}
	

	//*************************************************************************
	// Name    : setObservations
	// 
	// Purpose : set the observations matrix
	// 
	// Notes   : rowOrder is set to true if the first dimension of the input
	//           array contains the rows (observations). If rowOrder is set to
	//           false the first dimension of the input array refers to the 
	//           columns. If standardize is set to true, z scores are computed
	//           for each variable. Each row of the observations matrix 
	//           refers to an observation and each column, a variable.
	// 
	//*************************************************************************
	public void setObservations(double[][] observations,
					boolean rowOrder, boolean standardize) {
		try {
			
			// standardize the observations array if required
			if ( standardize == true ) {
				observations = NCGStatUtils.standardize(observations,rowOrder);
			}
			
			if (rowOrder == true) {
				this.observations =  
					MatrixUtils.createRealMatrix(observations);
			} else {
				this.observations =  
					MatrixUtils.createRealMatrix(observations).transpose();
			}
			
		} catch (IllegalArgumentException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
			this.observations = null;
		} catch (NullPointerException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
			this.observations = null;
		}
				
	}
	
	//*************************************************************************
	// Name    : getObservations
	// 
	// Purpose : return a copy of the observations matrix as a 2d array of doubles
	// 
	// Notes   : throws a new PCAException if observations is not set
	// 
	//*************************************************************************
	public double[][] getObservations() throws PCAException {		
		validateObservations();
		return observations.getData();		
	}
	
	//*************************************************************************
	// Name    : getPrincipalComponents
	// 
	// Purpose : returns a copy of the principalComponents matrix as a
	//           2d array of doubles
	// 
	// Notes   : throws a new PCAException if principalComponents is not set
	// 
	//*************************************************************************
	public double[][] getPrincipalComponents() throws PCAException {		
		validatePrincipalComponents();
		return principalComponents.getData();		
	}
	
	//*************************************************************************
	// Name    : getEigenValues
	// 
	// Purpose : returns a copy of the (real) eigenValues vector
	// 
	// Notes   : throws a new PCAException if eigenValues is not set
	// 
	//*************************************************************************
	public double[] getEigenValues() throws PCAException {		
		validateEigenValues();
		return eigenValues.getData();		
	}
	
	
	//*************************************************************************
	// Name    : getEigenVectors
	// 
	// Purpose : returns a copy of the EigenVectors matrix as a 2d array of doubles
	// 
	// Notes   : throws a new PCAException if eigenVectors is not set
	// 
	//*************************************************************************
	public double[][] getEigenVectors() throws PCAException {		
		validateEigenVectors();
		return eigenVectors.getData();		
	}
		
	//*************************************************************************
	// Name    : getCovarianceMatrix
	// 
	// Purpose : returns the covariance matrix for the observations matrix
	// 
	// Notes   : throws a new PCAException if observations is not set
	//           returns null if an error occurs
	// 
	//*************************************************************************
	private RealMatrix getCovarianceMatrix() throws PCAException {
		
		// validate the observations matrix
		validateObservations();
		
		RealMatrix covMatrix = null;
			
		try {
			
			if ( observations.getColumnDimension() > 1) {
				
				// compute covariance matrix if we have more than 1 attribute
				Covariance c = new Covariance(observations);
				covMatrix = c.getCovarianceMatrix();
			
			} else {
				
				// if we only have one attribute calculate the variance instead
				covMatrix = MatrixUtils.createRealMatrix(1,1);
				covMatrix.setEntry(0, 0, StatUtils.variance(observations.getColumn(0)));	
			
			}
		} catch (MatrixIndexException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
			covMatrix = null;
		} catch (IllegalArgumentException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
			covMatrix = null;
		} catch (NullPointerException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
			covMatrix = null;
		}
		
		return covMatrix;
	}
	
	//*************************************************************************
	// Name    : transform
	// 
	// Purpose : compute the eigen values, eigen vectors for the covariance matrix
	//           of the observations matrix and transform observations 
	// 
	// Notes   : throws a new PCAException if observations is not set
	// 
	//*************************************************************************
	// do the transformation
	public void transform() throws PCAException {
		
		// check to see if the observations matrix has been set
		validateObservations();
		
		// compute the covariance matrix
		RealMatrix covMatrix = getCovarianceMatrix();
		
		try {
			
			// get the eigenvalues and eigenvectors of the covariance matrixE
			EigenDecomposition eDecomp = new EigenDecompositionImpl(covMatrix,0.0);
			
			// set the eigenVectors matrix
			// the columns of the eigenVectors matrix are the eigenVectors of
			// the covariance matrix
			eigenVectors = eDecomp.getV();
			
			// set the eigenValues vector
			eigenValues = new ArrayRealVector(eDecomp.getRealEigenvalues());
			
			//transform the data
			principalComponents = observations.multiply(eigenVectors);
			
		} catch (InvalidMatrixException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
		} catch (NullPointerException e) {
			logger.severe(e.toString() + " : " + e.getMessage());
			e.printStackTrace();
		}	
	}
}