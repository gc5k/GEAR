package gear.simulation.qc.rowqc;

import java.util.ArrayList;

import gear.family.pedigree.Hukou;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.qc.rowqc.SampleFilter;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class RealDataSimulationQC
{

	protected MapFile MapData;
	protected PedigreeFile PedData;
	// protected double[] permuted_score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
	private ArrayList<Hukou> HukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	public RealDataSimulationQC(PedigreeFile ped, MapFile map, SampleFilter sf)
	{
		PedData = ped;
		MapData = map;
		PersonTable = sf.getSample();
		HukouBook = sf.getHukouBook();
		QC();
	}

	public void QC()
	{
		qualification();
	}

	private void qualification()
	{

	}

	public int getNumberMarker()
	{
		return MapData.getMarkerNumber();
	}

	public MapFile getMapFile()
	{
		return MapData;
	}

	public int SampleSize()
	{
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample()
	{
		return PersonTable;
	}

	public ArrayList<Hukou> getHukouBook()
	{
		return HukouBook;
	}
}
