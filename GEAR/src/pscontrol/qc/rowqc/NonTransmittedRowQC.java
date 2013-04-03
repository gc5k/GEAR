package pscontrol.qc.rowqc;

import java.util.ArrayList;

import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.qc.rowqc.SampleFilter;

public class NonTransmittedRowQC {
	protected MapFile MapData;
	protected PedigreeFile PedData;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records

	protected int[][] num_qualified;//
	protected boolean[][] keep;
	
	public NonTransmittedRowQC(PedigreeFile ped, MapFile map, SampleFilter sf) {
		MapData = map;
		PedData = ped;
	}
	
}
