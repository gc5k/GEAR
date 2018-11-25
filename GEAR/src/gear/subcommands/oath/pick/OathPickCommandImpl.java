package gear.subcommands.oath.pick;

import java.io.PrintStream;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.BufferedReader;

public class OathPickCommandImpl extends CommandImpl {

	OathPickCommandArguments pkArgs = null;
	ArrayList<String> pList = null;

	@Override
	public void execute(CommandArguments cmdArgs) {
		pkArgs = (OathPickCommandArguments) cmdArgs;
		pList = pkArgs.getPickList();
		readPList();
	}

	private void readPList() {

		String[] title = null;
		double[][] beta = null;
		double[][] se = null;

		ArrayList<ArrayList<String>> Info = NewIt.newArrayList();

		for (int i = 0; i < pList.size(); i++) {
			Logger.printUserLog("Reading '" + pList.get(i) + "'.");
			BufferedReader br = BufferedReader.openGZipFile(pList.get(i), "oath nss");

			String[] val;
			int tokenSize = 7;

			if (i == 0) {
				tokenSize = 9;
				title = br.readTokensAtLeast(tokenSize);

				ArrayList<String> bArray = NewIt.newArrayList();
				ArrayList<String> sArray = NewIt.newArrayList();

				int cnt = 0;
				while ((val = br.readTokensAtLeast(tokenSize)) != null) {
					bArray.add(val[7]);
					sArray.add(val[8]);

					ArrayList<String> info = NewIt.newArrayList();
					if (cnt == 0) {
						info.add(val[0]);
						info.add(val[1]);
						info.add(val[2]);
						info.add(val[3]);
						info.add(val[4]);
						info.add(val[5]);
						info.add(val[6]);
					}
					Info.add(info);
				}

				beta = new double[bArray.size()][pList.size()];
				se = new double[bArray.size()][pList.size()];

				for (int j = 0; j < bArray.size(); j++) {
					beta[j][i] = Double.parseDouble(bArray.get(j));
					se[j][i] = Double.parseDouble(sArray.get(j));
				}
			} else {
				val = br.readTokensAtLeast(tokenSize); // read title

				int cnt = 0;
				while ((val = br.readTokensAtLeast(tokenSize)) != null) {
					beta[cnt][i] = Double.parseDouble(val[5]);
					se[cnt][i] = Double.parseDouble(val[6]);
					cnt++;
				}
			}
		}

		PrintStream infoOut = FileUtil.CreatePrintStream(pkArgs.getOutRoot() + ".oath.info");
		PrintStream betaOut = FileUtil.CreatePrintStream(pkArgs.getOutRoot() + ".oath.beta");
		PrintStream seOut = FileUtil.CreatePrintStream(pkArgs.getOutRoot() + ".oath.se");

		for (int i = 0; i < 7; i++) {
			if (i != 6) {
				infoOut.print(title[i] + " ");
			} else {
				infoOut.println(title[i]);
			}
		}

		for (int i = 0; i < pList.size(); i++) {
			if (i < (pList.size() - 1)) {
				betaOut.print("Beta" + i + " ");
				seOut.print("SE" + i + " ");
			} else {
				betaOut.println("Beta" + i);
				seOut.println("SE" + i);
			}
		}

		for (int i = 0; i < beta.length; i++) {
			ArrayList<String> info = Info.get(i);

			for (int j = 0; j < info.size(); j++) {
				if (j < (info.size()-1)) {
					infoOut.print(info.get(j) + " ");					
				} else {
					infoOut.println(info.get(j));
				}
			}

			for (int j = 0; j < beta[i].length; j++) {
				if (j < (beta[i].length - 1)) {
					betaOut.print(beta[i][j] + " ");
					seOut.print(se[i][j] + " ");
				} else {
					betaOut.println(beta[i][j]);
					seOut.println(se[i][j]);
				}
			}
		}
		infoOut.close();
		betaOut.close();
		seOut.close();
		Logger.printUserLog("\n");
		Logger.printUserLog("Saving loci information to '" + pkArgs.getOutRoot() + ".info'.");
		Logger.printUserLog("Saving beta to '" + pkArgs.getOutRoot() + ".beta'.");
		Logger.printUserLog("Saving se to '" + pkArgs.getOutRoot() + ".se'.");
	}
}
