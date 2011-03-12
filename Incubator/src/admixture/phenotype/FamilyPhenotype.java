package admixture.phenotype;

public class FamilyPhenotype {
	private double[][] p_phe;
	private int[] p_status;
	private double[][] o_phe;
	private int[] o_status;

	private int num_kid;
	private int FamID;

	public FamilyPhenotype(int FI, double[][] p_p, int[] p_s, double[][] o_p, int[] o_s ) {
		FamID = FI;
		p_phe = p_p;
		p_status = p_s;
		o_phe = o_p;
		o_status = o_s;
	}

	public int[] getOffspringStatus() {
		return o_status;
	}

	public int getNumberOffspring() {
		return o_status.length;
	}

	public int getNumberAffectedOffspring() {
		int s = 0; 
		for (int i = 0; i < o_status.length; i++) {
			s += o_status[i];
		}
		return s;
	}

	public String getStringParentPhenotype(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		sb.append(p_status[idx]+" ");
		for(int i = 0; i < p_phe[idx].length; i++) {
			sb.append(p_phe[idx][i]+" ");
		}
		return sb.toString();
	}

	public String getStringOffspringPhenotype(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		sb.append(o_status[idx]+" ");
		for(int i = 0; i < o_phe[idx].length; i++) {
			sb.append(o_phe[idx][i]+" ");
		}
		return sb.toString();
	}

	public void print() {
		System.out.println("FamID " + FamID);
		for(int i = 0; i < p_phe.length; i++) {
			System.out.println("Parent " + i + " status: " + p_status[i]);
			for(int j = 0; j < p_phe[i].length; j++) {
				System.out.print(p_phe[i][j] + " ");
			}
			System.out.println();
		}

		for(int i = 0; i < o_phe.length; i++) {
			System.out.println("Kid " + i + " status: " + o_status[i]);
			for(int j = 0; j < o_phe[i].length; j++) {
				System.out.print(o_phe[i][j] + ", ");
			}
			System.out.println();
		}
	}
}
