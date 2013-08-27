package gear.sumstat.qc.rowqc;

import java.util.ArrayList;

import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.qc.rowqc.SampleFilter;

public class SumStatQC
{

	private SampleFilter samFilter;
	protected MapFile MapData;
	protected PedigreeFile PedData;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records

	protected int[][] num_qualified;//
	protected boolean[][] keep;

	public SumStatQC(PedigreeFile ped, MapFile map, SampleFilter sf)
	{

		samFilter = sf;
		PedData = ped;
		MapData = map;
		PersonTable = sf.getSample();

		QC();

	}

	public void QC()
	{

	}

	public ArrayList<PersonIndex> getSample()
	{
		return PersonTable;
	}

}
