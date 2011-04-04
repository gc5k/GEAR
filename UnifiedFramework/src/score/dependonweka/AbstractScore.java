package score.dependonweka;

/**
 * A public interface, is inherited by all other generalized model methods.
 * @author Guo-Bo Chen, Zhejiang University, chenguobo@zju.edu.cn
 *
 */

public interface AbstractScore {
	/**
	 * Set the ridge of the model, a virtual method.
	 * @param r 
	 */
	void SetRidge (double r);
	/**
	 * Calculate score after the implementation of generalized linear regression with given parameters. 
	 */
	void CalculateScore();
	double[] GetScore();
	double[] getCoefficients();
}
