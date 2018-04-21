package gear.subcommands;

import java.util.ArrayList;
import java.util.Arrays;
import gear.util.Logger;
import gear.util.NewIt;

public class CommandSplitter {
	private final String[] cmd;
	private String opt = "GEAR";

	private int intMin;
	private boolean isIntMin;

	private int intMax;
	private boolean isIntMax;

	private double doubleMin;
	private boolean isDoubleMin;

	private double doubleMax;
	private boolean isDoubleMax;

	private CommandSplitter(Builder builder) {
		this.cmd = builder.cmd;
		this.opt = builder.opt;
		this.intMin = builder.intMin;
		this.intMax = builder.intMax;
		this.isIntMin = builder.isIntMin;
		this.isIntMax = builder.isIntMax;	
		this.doubleMax = builder.doubleMax;
		this.doubleMin = builder.doubleMin;
		this.isDoubleMax = builder.isDoubleMax;
		this.isDoubleMin = builder.isDoubleMin;
	}

	public int[] ParseToInt() {
		ArrayList<String> intData = ParseToStr();
		int[] dat = new int[intData.size()];
		for(int i = 0; i < intData.size(); i++) {
			dat[i] = Integer.parseInt(intData.get(i)) - 1;
		}
		Arrays.sort(dat);
		return dat;
	}

	public ArrayList<String> ParseToStr() {
		ArrayList<String> S = NewIt.newArrayList();
		for(int i = 0; i < cmd.length; i++) {
			if (cmd[i].contains("-")) {
				String[] s = cmd[i].split("-");
				try {
					int s0 = Math.min(Integer.parseInt(s[0]), Integer.parseInt(s[1]));
					int s1 = Math.max(Integer.parseInt(s[0]), Integer.parseInt(s[1]));
					if (isIntMin) {
						tooSmallInt(s0, cmd[i]);
					}
					if (isIntMax) {
						tooBigInt(s1, cmd[i]);
					}
					for(int j = s0; j <= s1; j++) {
						S.add( (new Integer(j)).toString());
					}
				} catch (NumberFormatException e) {
					Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
					System.exit(1);
				}
			} else {
				try {
					int s = Integer.parseInt(cmd[i]);
					if (isIntMin) tooSmallInt(s, cmd[i]);
					if (isIntMax) tooBigInt(s, cmd[i]);
					S.add(cmd[i]);
				} catch (NumberFormatException e) {
					Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
					System.exit(1);
				}
			}
		}
		return S;
	}

	public double[][] ParseToDouble() {
		double[][] rg = new double[cmd.length][2];
		for(int i = 0; i < cmd.length; i++) {
			if (cmd[i].contains("-")) {
				String[] s = cmd[i].split("-");
				try {
					double s0 = 0, s1 = 0;
					if (isDoubleMin) {
						s0 = Double.parseDouble(s[0]); 
						tooSmallDouble(s0, cmd[i]);
					}
					if (isDoubleMax) {
						s1 = Double.parseDouble(s[1]);
						tooBigDouble(s1, cmd[i]);
					}
					rg[i][0] = s0; rg[i][1] = s1;
				} catch (NumberFormatException e) {
					Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
					System.exit(1);
				}
			} else {
				Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
				System.exit(1);				
			}
		}
		return rg;
	}
	
	private void tooSmallDouble(double s, String c) {
		if (s < doubleMin) {
			Logger.printUserLog(s + " is smaller than " + doubleMin  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}
	}

	private void tooBigDouble(double b, String c) {
		if (b > doubleMax) {
			Logger.printUserLog(b + " is greater than " + doubleMax  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}
	}

	private void tooSmallInt(int s, String c) {
		if (s < intMin) {
			Logger.printUserLog(s + " is smaller than " + intMin  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}
	}

	private void tooBigInt(int b, String c) {
		if (b > intMax) {
			Logger.printUserLog(b + " is greater than " + intMax  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}
	}

	public static class Builder {
		private final String[] cmd;

		private String opt = null;
		private int intMin = -1;
		private boolean isIntMin = false;

		private int intMax = 1000;
		private boolean isIntMax = false;

		private double doubleMin;
		private boolean isDoubleMin;

		private double doubleMax;
		private boolean isDoubleMax;

		public Builder(String[] cmd) {
			this.cmd = cmd;
		}

		public Builder OPT(String opt) {
			this.opt = opt;
			return this;
		}
		
		public Builder IntMin(int min) {
			this.intMin = min;
			isIntMin = true;
			return this;
		}

		public Builder IntMax(int max) {
			this.intMax = max;
			isIntMax = true;
			return this;
		}

		public Builder DoubleMin(double min) {
			this.doubleMin = min;
			isDoubleMin = true;
			return this;
		}

		public Builder DoubleMax(double max) {
			this.doubleMax = max;
			isDoubleMax = true;
			return this;
		}

		public CommandSplitter create() {
			return new CommandSplitter(this);
		}
	}
}
