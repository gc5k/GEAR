package admixture.chromosome;

public class FamilySingleChromosome {
	protected int[][][] p_g; // parental chromosomes;
							 // p_g[0] for dad, p_g[1] for mom
	protected int[][][] o_g; // offspring chromosomes;

	protected int modCount = 0;
	
	public FamilySingleChromosome(int[][][] p, int[][][] o) {
		p_g = p;
		o_g = o;
	}
	
	public FamilySingleChromosome(int id, int num_kid) {
		p_g = new int[2][][];
		o_g = new int[num_kid][][];
	}

	public void AddFatherChr(int[][] g) {
		p_g[0] = g;
	}

	public void AddMotherChr(int[][] g) {
		p_g[1] = g;
	}
	
	public void AddOffspringChr(int[][] g) {
		try {
			if(modCount == o_g.length) {
				throw new Exception("Family " + " no space for more than " + modCount + "kids");
			}
			o_g[modCount++] = g;
		} catch(Exception E) {
			E.printStackTrace(System.err);
		}
	}

}
