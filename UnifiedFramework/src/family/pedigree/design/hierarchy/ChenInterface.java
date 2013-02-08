package family.pedigree.design.hierarchy;

import java.util.ArrayList;

import family.mdr.data.PersonIndex;
import family.pedigree.file.MapFile;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public interface ChenInterface {
	public static int Linear = 0;
	public static int Logistc = 1;

	public void setSeed(long s);
	public void recoverSeed();
	public double[] getStatus();
	public void getPermutedScore(boolean nested);
	public void RecoverScore();
	public MapFile getMapFile();
	public int getNumberMarker();
	public int SampleSize();
	public ArrayList<PersonIndex> getSample();
}
