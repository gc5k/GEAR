package gear.pscontrol.qc.rowqc;

import java.util.ArrayList;

import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.qc.rowqc.SampleFilter;

public class NonTransmittedRowQC
{
	protected MapFile MapData;
	protected PedigreeFile PedData;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records

	protected int[][] num_qualified;//
	protected boolean[][] keep;

	public NonTransmittedRowQC(PedigreeFile ped, MapFile map, SampleFilter sf)
	{
		MapData = map;
		PedData = ped;
	}

}
