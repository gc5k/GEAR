package pscontrol.hierarchy;

import java.util.ArrayList;

import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public interface ChenInterface {

	public MapFile getMapFile();
	public int getNumberMarker();
	public int SampleSize();
	public ArrayList<PersonIndex> getSample();
}
