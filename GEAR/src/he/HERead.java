package he;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import gear.util.FileProcessor;
import gear.util.Logger;

import parameter.Parameter;

public class HERead {
	private final String delim = "\\s+";
	protected boolean[] flag;
	HashMap<String, Integer> ID;
	protected String grmFile;
	protected String grmID;
	protected String keepFile;
	protected String phenoFile;
	protected String output;
	protected int[] mpheno;
	protected int perm;
	protected boolean permFlag = false;

	protected boolean reverse;
	protected boolean k_button;
	protected double k;
	protected parameter.HEType heType;

	protected HashMap<String, Integer> ID2Idx;

	protected double yyProd;

	protected double[][] XtX;
	protected double[] XtY;
	protected double[][] y;
	protected int dim;
	protected double P;
	protected boolean isCC = false;

	protected Lambda lambda;
	StringBuffer sb = new StringBuffer();

	public HERead() {
		grmFile = Parameter.INSTANCE.getHEParameter().getGrm();
		grmID = Parameter.INSTANCE.getHEParameter().getGrmId();
		keepFile = Parameter.INSTANCE.keepFile;
		phenoFile = Parameter.INSTANCE.getHEParameter().getPheno();
		mpheno = Parameter.INSTANCE.getHEParameter().getMPheno();
		reverse = Parameter.INSTANCE.reverse;
		k_button = Parameter.INSTANCE.k_button;
		k = Parameter.INSTANCE.k;
		output = Parameter.INSTANCE.out;
		heType = Parameter.INSTANCE.getHEParameter().getType();
		permFlag = Parameter.INSTANCE.permFlag;
		perm = Parameter.INSTANCE.perm;

		XtX = new double[mpheno.length + 1][mpheno.length + 1];
		XtY = new double[mpheno.length + 1];

		String line;

		// *********************************** read grm id
		BufferedReader reader = FileProcessor.FileOpen(grmID);

		int i2 = 0;
		ID2Idx = new HashMap<String, Integer>();
		try {
			while ((line = reader.readLine()) != null) {
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				ID2Idx.put(sb.toString(), i2++);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		flag = new boolean[ID2Idx.size()];
		Arrays.fill(flag, false);

		// *********************************** read pheno file
		reader = FileProcessor.FileOpen(phenoFile);

		y = new double[flag.length][mpheno.length + 1];

		try {
			while ((line = reader.readLine()) != null) {
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				if (ID2Idx.containsKey(sb.toString())) {
					int ii = ID2Idx.get(sb.toString());
					boolean f = true;
					y[ii][0] = 1;
					for (int j = 0; j < mpheno.length; j++) {
						if (Parameter.INSTANCE.isNA(s[1 + mpheno[j]])) {
							f = false;
							break;
						} else {
							y[ii][j + 1] = Double.parseDouble(s[1 + mpheno[j]]);
						}
					}
					flag[ii] = f;
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// ************************keep
		if (keepFile != null) {
			reader = FileProcessor.FileOpen(keepFile);
			boolean[] ff = new boolean[flag.length];
			Arrays.fill(ff, false);
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder(s[0] + "." + s[1]);
					if (ID2Idx.containsKey(sb.toString())) {
						int ii = ID2Idx.get(sb.toString());
						ff[ii] = true;
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			for (int i = 0; i < ff.length; i++) {
				flag[i] &= ff[i];
			}
		}

		int Len = 0;
		for (int i = 0; i < flag.length; i++)
			if (flag[i]) Len++;
		dim = Len*(Len-1)/2;
		
		// ************************************standardising
		double[] ss = new double[y[0].length-1];
		double[] ssx = new double[y[0].length-1];
		for (int i = 0; i < flag.length; i++) {
			if(!flag[i]) continue; 
			for(int j = 0; j < ss.length; j++) {
				ss[j] += y[i][j+1];
				ssx[j] += y[i][j+1] * y[i][j+1];
			}
		}
		double[] sd = new double[ssx.length];
		for (int i = 0; i < sd.length; i++) {
			ss[i] /= Len;
			sd[i] = Math.sqrt((ssx[i] - Len * ss[i] * ss[i])/(Len-1));
		}

		if (Parameter.INSTANCE.scale) {
			Logger.printUserLog("Standardising phentoype.");
			for (int i = 0; i < flag.length; i++) {
				if(!flag[i]) continue;
				for (int j = 1; j < y[i].length; j++) {
					y[i][j] = (y[i][j] - ss[j-1])/sd[j-1];
				}
			}
		}
	}
	
	public static void main(String[] args) {
		Parameter.INSTANCE.commandListener(args);
	}

}
