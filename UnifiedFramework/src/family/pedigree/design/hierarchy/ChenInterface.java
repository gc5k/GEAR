package family.pedigree.design.hierarchy;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public interface ChenInterface {
	public static int Linear = 0;
	public static int Logistc = 1;

	public void setSeed(long s);
	public String[] getMarkerName();
	public String[] getScoreName();
	public byte[][] getGenotype();
	public byte[] getStatus();
	public double[] getScore();
	public double[][] getScore2();
	public double[] getPermutedScore(boolean nested);
	public void print2MDRFormat(String f);
}
