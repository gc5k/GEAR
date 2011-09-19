package statistics.FisherExactTest;

public class FastFisherExactTest {

	private double Pobs;
	private double Pexcess;
	private double Pdificit;
	private double Phwe;
	private double Pw;
	
	private int N;
	private int NAB;
	private int NA;
	public FastFisherExactTest(int N, int NAB, int NA) {
		this.N = N;
		this.NAB = NAB;
		this.NA = NA;

		Pobs = base();

		Pexcess();
		Pdificit();
		Phwe += Pobs;
	}

	private double base() {
		double f = 0;
		f += Math.log(2) * NAB;
		for (int i = 1; i <= N; i++) {
			f += Math.log(i);
		}
		for (int i = 1; i <= NA; i++) {
			f += Math.log(i);
		}

		int NB = 2 * N - NA;
		for (int i = 1; i <= NB; i++) {
			f += Math.log(i);
		}

		int NAA = (NA - NAB) / 2;
		for (int i = 1; i <= NAA; i++) {
			f -= Math.log(i);
		}

		for (int i = 1; i <= NAB; i++) {
			f -= Math.log(i);
		}

		int NBB = (2 * N - NA - NAB) / 2;
		for (int i = 1; i <= NBB; i++) {
			f -= Math.log(i);
		}

		for (int i = 1; i <= 2 * N; i++) {
			f -= Math.log(i);
		}

		return Math.exp(f);
	}

	public void Pexcess() {
		Pexcess = Pobs;
		double f = Pobs;
		for (int i = NAB; i <= NA - 2 ; i += 2) {
			
			int NAA = (NA - i) / 2;
			int NBB = (2 * N - NA - i) / 2;
			f = f * 4 * NAA * NBB / ((i + 2) * (i + 1));
			Pexcess += f;
			if ( Pobs >= f) {
				Phwe += f;
			}
			Pw += f;
		}
	}
	
	public void Pdificit() {
		Pdificit = Pobs;
		double f = Pobs;
		for (int i = NAB; i >= 0; i -= 2) {
			int NAA = (NA - i) / 2;
			int NBB = (2 * N - NA - i) / 2;
			f = f * i * (i - 1) / (4 * (NAA + 1) * (NBB + 1));
			Pdificit += f;
			if ( Pobs >= f) {
				Phwe += f;
			}
			Pw += f;
		}
	}
	
	public double getPexcess() {
		return Pexcess;
	}

	public double getPdificit() {
		return Pdificit;
	}

	public double HDP() {
		return Phwe;
	}

	public double getPW() {
		return Pw + Pobs;
	}
	
	public double getPobs() {
		return Pobs;
	}

	public static void main(String[] args) {
		FastFisherExactTest f0 = new FastFisherExactTest(20, 0, 8);
		System.out.println(f0.getPobs() + " " + f0.getPW() + " " + f0.HDP());
		FastFisherExactTest f2 = new FastFisherExactTest(20, 2, 8);
		System.out.println(f2.getPobs() + " " + f2.getPW() + " " + f2.HDP());
		FastFisherExactTest f4 = new FastFisherExactTest(20, 4, 8);
		System.out.println(f4.getPobs() + " " + f4.getPW() + " " + f4.HDP());
		FastFisherExactTest f6 = new FastFisherExactTest(20, 6, 8);
		System.out.println(f6.getPobs() + " " + f6.getPW() + " " + f6.HDP());
		FastFisherExactTest f8 = new FastFisherExactTest(20, 8, 8);
		System.out.println(f8.getPobs() + " " + f8.getPW() + " " + f8.HDP());
	}
}
