package family.pedigree.design.hierarchy;

public interface ChenInterface {
	public static int Linear = 0;
	public static int Logistc = 1;
	
	public void generateScore(int pheIdx, int[] covIdx, int method);
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
