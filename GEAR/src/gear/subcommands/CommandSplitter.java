package gear.subcommands;

import java.util.ArrayList;
import java.util.Arrays;
import gear.util.Logger;
import gear.util.NewIt;

public class CommandSplitter {
	private final String[] cmd;
	private String opt = "GEAR";

	private int min;
	private boolean isMin;

	private int max;
	private boolean isMax;

	private CommandSplitter(Builder builder) {
		this.cmd = builder.cmd;
		this.opt = builder.opt;
		this.min = builder.min;
		this.max = builder.max;
		this.isMin = builder.isMin;
		this.isMax = builder.isMax;		
	}

	public int[] ParseToInt() {
		ArrayList<String> intData = Parse();
		int[] dat = new int[intData.size()];
		for(int i = 0; i < intData.size(); i++) {
			dat[i] = Integer.parseInt(intData.get(i)) - 1;
		}
		Arrays.sort(dat);
		return dat;
	}

	public ArrayList<String> Parse() {
		ArrayList<String> S = NewIt.newArrayList();
		for(int i = 0; i < cmd.length; i++) {
			if (cmd[i].contains("-")) {
				String[] s = cmd[i].split("-");
				try {
					int s0 = Integer.parseInt(s[0]); 
					int s1 = Integer.parseInt(s[1]);
					if (isMin) tooSmall(s0, cmd[i]);
					if (isMax) tooBig(s1, cmd[i]);
					for(int j = Math.min(s0, s1); j <= Math.max(s0, s1); j++) {
						S.add( (new Integer(j)).toString());
					}
				} catch (NumberFormatException e) {
					Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
					System.exit(1);
				}
			} else {
				try {
					int s = Integer.parseInt(cmd[i]);
					if (isMin) tooSmall(s, cmd[i]);
					if (isMax) tooBig(s, cmd[i]);
					S.add(cmd[i]);
				} catch (NumberFormatException e) {
					Logger.printUserLog("The value of " + cmd[i] + " is invalid for --" + opt + ". GEAR quit.");
					System.exit(1);
				}
			}
		}
		return S;
	}

	private void tooSmall(int s, String c) {
		if (s < min) {
			Logger.printUserLog(s + " is smaller than " + min  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}

	}

	private void tooBig(int b, String c) {
		if (b > max) {
			Logger.printUserLog(b + " is greater than " + max  + " in --" + c + ". GEAR quit.");
			System.exit(1);
		}
	}
	
	
	public static class Builder {
		private final String[] cmd;

		private String opt = null;
		private int min = -1;
		private boolean isMin = false;
		
		private int max = 1000;
		private boolean isMax = false;

		public Builder(String[] cmd) {
			this.cmd = cmd;
		}

		public Builder opt(String opt) {
			this.opt = opt;
			return this;
		}
		
		public Builder min(int min) {
			this.min = min;
			isMin = true;
			return this;
		}

		public Builder max(int max) {
			this.max = max;
			isMax = true;
			return this;
		}

		public CommandSplitter create() {
			return new CommandSplitter(this);
		}
	}
}
