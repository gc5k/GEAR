package sumstat.qc.rowqc;

import java.util.ArrayList;

import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.qc.rowqc.SampleFilter;

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
