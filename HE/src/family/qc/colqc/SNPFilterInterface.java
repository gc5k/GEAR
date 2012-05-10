package family.qc.colqc;

public interface SNPFilterInterface {

	public void Select();
	public int[] getWorkingSNP();
	public int[] getBgSeq();
	public int[] getWSeq();
	public int[][] getWSeq2();

}
