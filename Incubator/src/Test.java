import jsc.distributions.DiscreteUniform;

public class Test {

	public static void main(String[] args) {
		DiscreteUniform U = new DiscreteUniform(1, 1);
		for (int i = 0; i < 10; i++) {
			System.out.println(U.random());
		}
	}
}
