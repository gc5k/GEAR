package gear.util.stat;

public class CommonMath {
	public static int calculateTriangleSize(int rows) {
		return (rows * rows + rows) >> 1;
	}
}
