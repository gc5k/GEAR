package statistics.FisherExactTest;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
//--------------------
//--------------------
//----| a | b | n1
//-----------------
//----| c | d | n2
//-----------------
//----| n3| n4| N
//-------------------
public class FisherExactTestExtended {
	private double n1;//cases
	private double n2;//controls
	private double n3;//affected
	private double n4;//unaffacted
	private double N;//sample size
	private double[] f;
	private boolean complete;
	double[][] cum_P;
	double[][] P;
	int[] n3_offset_for_a;

	private boolean debug = false;

	public FisherExactTestExtended(int cases, int controls) {
		n1 = cases;
		n2 = controls;
		N = n1 + n2;
		f = new double[(int) N + 1];
		f[0] = 0.0;
		for (int i = 1; i <= N; i++) {
			f[i] = f[i - 1] + Math.log(i);
		}
	}

	public void calculateExtensiveTable() {
		complete = false;
		cum_P = new double[(int) N - 1][];
		P = new double[(int) N - 1][];
		n3_offset_for_a = new int[(int) N - 1];
		for (int i = 0; i < N - 1; i++) {
			n3 = i + 1;
			n4 = N - n3;
			int a = (int) Math.max(n3 - n2, Math.ceil((n1 * n3 / N)));
			int b = (int) n1 - a;
			int c = (int) n3 - a;
			int d = (int) n2 - c;
			int bound = Math.min((int) n1, (int) n3);
			n3_offset_for_a[i] = a;
			cum_P[i] = new double[bound - a + 1];
			P[i] = new double[bound - a + 1];
			if (debug) {
				System.out.println(" a=" + a + "bound=" + bound + " n1*n3/N=" + n1 * n3 / N);
				System.out.println("n1=" + n1 + ", n2=" + n2 + ", n3=" + n3 + ", n4=" + n4);
			}

			for (int j = n3_offset_for_a[i]; j <= bound; j++) {
				P[i][j - n3_offset_for_a[i]] = getP(a++, b--, c--, d++);
			}
			cum_P[i][cum_P[i].length - 1] = P[i][P[i].length - 1];
			for (int j = P[i].length - 2; j >= 0; j--) {
				cum_P[i][j] = cum_P[i][j + 1] + P[i][j];
			}
		}
	}

	public void calculateCompleteExtensiveTable() {
		complete = true;
		cum_P = new double[(int) N + 1][];
		P = new double[(int) N + 1][];
		n3_offset_for_a = new int[(int) N + 1];
		for (int i = 0; i < N + 1; i++) {
			n3 = i;
			n4 = N - n3;
			int a = Math.max(0, (int) (n3 - n2));
			int b = (int) n1 - a;
			int c = (int) n3 - a;
			int d = (int) n2 - c;
			int bound = Math.min((int) n1, (int) n3);
			n3_offset_for_a[i] = a;
			cum_P[i] = new double[bound - a + 1];
			P[i] = new double[bound - a + 1];
			if (debug) {
				System.out.println(" a=" + a + "bound=" + bound + " n1*n3/N=" + n1 * n3 / N);
				System.out.println("n1=" + n1 + ", n2=" + n2 + ", n3=" + n3 + ", n4=" + n4);
			}

			for (int j = n3_offset_for_a[i]; j <= bound; j++) {
				P[i][j - n3_offset_for_a[i]] = getP(a++, b--, c--, d++);
			}
			cum_P[i][cum_P[i].length - 1] = P[i][P[i].length - 1];
			for (int j = P[i].length - 2; j >= 0; j--) {
				cum_P[i][j] = cum_P[i][j + 1] + P[i][j];
			}
		}
	}

	public final void calculatePTable(boolean C) {
		complete = C;
		int dim = (int) Math.ceil((N+1)/2);
		cum_P = new double[dim][];
		P = new double[dim][];
		n3_offset_for_a = new int[dim];
		for (int i = 0; i < dim; i++) {
			n3 = i;
			n4 = N - n3;
			int a = Math.max(0, (int) (n3 - n2));
			int b = (int) n1 - a;
			int c = (int) n3 - a;
			int d = (int) n2 - c;
			int bound = Math.min((int) n1, (int) n3);
			n3_offset_for_a[i] = complete ? Math.max(0, (int) (n3 - n2)) : (int) Math.max(n3 - n2, Math.ceil((n1 * n3 / N)));
			cum_P[i] = new double[bound - a + 1];
			P[i] = new double[bound - a + 1];
			if (debug) {
				System.out.println(" a=" + a + "bound=" + bound + " n1*n3/N=" + n1 * n3 / N);
				System.out.println("n1=" + n1 + ", n2=" + n2 + ", n3=" + n3 + ", n4=" + n4);
			}

			for (int j = Math.max(0, (int) (n3 - n2)); j <= bound; j++) {
				P[i][j - Math.max(0, (int) (n3 - n2))] = getP(a++, b--, c--, d++);
			}
			cum_P[i][cum_P[i].length - 1] = P[i][P[i].length - 1];
			for (int j = P[i].length - 2; j >= 0; j--) {
				cum_P[i][j] = cum_P[i][j + 1] + P[i][j];
			}
		}
	}

	public final double getP(int a, int b, int c, int d) {
		int n = a + b + c + d;
		if (debug) {
			System.out.println(a + "," + b + "," + c + "," + d + "--");
			System.out.println((a + b) + "," + (c + d) + "," + (a + c) + "," + (b + d) + "-");
		}
		if (n > N) {
			return Double.NaN;
		}
		double p;
		p = (f[a + b] + f[c + d] + f[a + c] + f[b + d]) - (f[a] + f[b] + f[c] + f[d] + f[n]);
		return Math.exp(p);
	}

	public double[][] getCumulativePmatrix() {
		return cum_P;
	}

	public double[][] getPmatrix() {
		return P;
	}

	public double getPvalue(int n_3, int a) {
		int mod = cum_P.length;
		int x1;
		if (n_3 < mod) {
			x1 = n_3;
		} else {
			int res = n_3 % mod;
			int o = (int) (N + 1) % 2;
			x1 = mod - res -1 - o;
		}
		int x2 = a;
		if( a - x2 < 0 || a - x2 > cum_P[x1].length || x1 < 0 || x1 > cum_P.length) {
			System.err.println("n3 ="+n_3 +", a=" + a +" is abnormal" );
			System.exit(0);
		}
		return cum_P[x1][x2]/cum_P[x1][n3_offset_for_a[x1]];
	}

	public int[] getOffset() {
		return n3_offset_for_a;
	}

	static public void main(String[] args) throws IOException {
		FisherExactTestExtended FETE = new FisherExactTestExtended(3, 4);
		FETE.calculatePTable(false);
		System.out.println(FETE.getPvalue(5, 2));
		double[][] p = FETE.getCumulativePmatrix();
		int[] offset = FETE.getOffset();
		PrintWriter pedout = new PrintWriter(new File("p.txt"));
		for (int i = 0; i < p.length; i++) {
			pedout.print(offset[i] + " ");
			for (int j = 0; j < p[i].length; j++) {
				pedout.print(p[i][j] + " ");
			}
			pedout.println();
		}
		pedout.close();
	}
}
