package simulation.qc.rowqc;

import java.util.ArrayList;
import java.util.Hashtable;

import family.pedigree.Hukou;
import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import family.qc.rowqc.SampleFilter;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class RealDataSimulationQC {

	protected MapFile MapData;
	protected PedigreeFile PedData;
	private SampleFilter sf;
	// protected double[] permuted_score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
	private ArrayList<Hukou> HukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	public RealDataSimulationQC(PedigreeFile ped, MapFile map, SampleFilter sf) {
		PedData = ped;
		MapData = map;
		this.sf = sf;
		PersonTable = sf.getSample();
		HukouBook = sf.getHukouBook();
		QC();
	}

	public void QC() {
		qualification();
	}

	private void qualification() {

	}

	public int getNumberMarker() {
		return MapData.getMarkerNumber();
	}

	public MapFile getMapFile() {
		return MapData;
	}

	public int SampleSize() {
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample() {
		return PersonTable;
	}

	public ArrayList<Hukou> getHukouBook() {
		return HukouBook;
	}
}
