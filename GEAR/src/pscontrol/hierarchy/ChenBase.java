package pscontrol.hierarchy;

import java.util.ArrayList;
import java.util.Random;


import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import gear.Parameter;
import gear.util.NewIt;


/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public abstract class ChenBase implements ChenInterface {

	protected long seed = 2011;
	public static Random rnd = new Random();
	protected MapFile MapData;
	protected PedigreeFile PedData;

	protected int qualified_Unrelated;
	protected int qualified_Sib;
	protected int[] numSib;
	// protected double[] permuted_score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records

	protected int[] subsetMarker;

	public ChenBase(PedigreeFile ped, MapFile map) {
		rnd.setSeed(Parameter.INSTANCE.nontransSeed);
		PedData = ped;
		MapData = map;

		PersonTable = NewIt.newArrayList();
	}

	@Override
	public int getNumberMarker() {
		return MapData.getMarkerNumber();
	}

	@Override
	public MapFile getMapFile() {
		return MapData;
	}

	public int SampleSize() {
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample() {
		return PersonTable;
	}
}
