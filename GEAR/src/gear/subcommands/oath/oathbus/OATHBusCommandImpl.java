package gear.subcommands.oath.oathbus;

import java.io.PrintStream;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.nss.NSSCommandArguments;
import gear.subcommands.oath.nss.NSSCommandImpl;
import gear.subcommands.oath.pick.OathPickCommandArguments;
import gear.subcommands.oath.pick.OathPickCommandImpl;
import gear.subcommands.oath.synthesize.SynthCommandArguments;
import gear.subcommands.oath.synthesize.SynthCommandImpl;
import gear.util.Combination;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OATHBusCommandImpl extends CommandImpl {

	private OATHBusCommandArguments obArgs = null;
	private NSSCommandArguments nssArgs = null;
	private OathPickCommandArguments pkArgs = null;
	private int N;

	@Override
	public void execute(CommandArguments cmdArgs) {
		obArgs = (OATHBusCommandArguments) cmdArgs;
		step1();
		step2();
		step3();

		writeMod();

		if (!obArgs.isKeepInter()) {
			Logger.printUserLog("Clean up intermediate files.");
			ArrayList<String> fList = pkArgs.getPickList();
			for(int i = 1; i < fList.size(); i++) {
				java.io.File myFilePath = new java.io.File(fList.get(i));
		        myFilePath.delete();				
			}

			java.io.File bFile = new java.io.File(obArgs.getOutRoot() + ".oath-bus");
			bFile.delete();

			java.io.File lbFile = new java.io.File(obArgs.getOutRoot() + ".oath-list");
			lbFile.delete();

			java.io.File lFile = new java.io.File(obArgs.getOutRoot() + ".list.nss");
			lFile.delete();

			java.io.File mFile = new java.io.File(obArgs.getOutRoot() + ".m.nss");
			mFile.delete();

			java.io.File pheFile = new java.io.File(obArgs.getOutRoot() + ".p." + (obArgs.getSelectedPhenotype()[0]+1)+".nss.gz");
			pheFile.delete();

			int[] covIdx = obArgs.getCovNumber();
			for(int i = 0; i < covIdx.length; i++) {
				java.io.File covFile = new java.io.File(obArgs.getOutRoot() + ".c." + (covIdx[i]+1)+".nss.gz");
				covFile.delete();				
			}

		}
	}
	
	private void writeMod() {
		Integer[] tmp = new Integer[obArgs.getCovNumber().length];
		ArrayList<ArrayList<Integer>> ModCode = NewIt.newArrayList();
		ArrayList<Integer> M1 = NewIt.newArrayList();
		for (int i = 0; i < tmp.length; i++) {
			M1.add(0);
		}
		ModCode.add(M1);

		for (int i = 1; i <= tmp.length; i++) {
			int[][] rs = Combination.combination(tmp, i);
			for (int j = 0; j < rs.length; j++) {
				int[] MD = new int[tmp.length];
				for (int k = 0; k < rs[j].length; k++) {
					MD[rs[j][k]] = 1;
				}

				ArrayList<Integer> M2 = NewIt.newArrayList();
				for (int k = 0; k < MD.length; k++) {
					M2.add(MD[k]);
				}
				ModCode.add(M2);
			}
		}

		PrintStream Mout = FileUtil.CreatePrintStream(obArgs.getOutRoot() + ".oath.mod");
		for (int i = 0; i < ModCode.size(); i++) {
			ArrayList<Integer> mm = ModCode.get(i);
			for (int j = 0; j < mm.size(); j++) {
				if (j < (mm.size() - 1)) {
					Mout.print(mm.get(j)+" ");
				} else {
					Mout.println(mm.get(j)+" ");
				}
			}
		}
		Mout.close();
		Logger.printUserLog("Saving model information in '"+ obArgs.getOutRoot() + ".oath.mod'.");
	}
	
	private void step1() {
		//step 1
		nssArgs = new NSSCommandArguments();
		nssArgs.setBFile(obArgs.getBFile());
		nssArgs.setPhenotypeFile(obArgs.getPhenotypeFile());
		nssArgs.setPhenotypeIndex(obArgs.getSelectedPhenotype(0)+1); //need to add 1 here, because the phenotype reading has "-1" option.
		nssArgs.setCovFile(obArgs.getCovFile());
		nssArgs.setCovNumber(obArgs.getCovNumber());
		nssArgs.setMAF((new Double(obArgs.getMAF())).toString());
		if (obArgs.isKeepFile())
			nssArgs.setKeepFile(obArgs.getKeepFile());

		if (obArgs.isChr())
			nssArgs.setChr(obArgs.getChr().toArray(new String[0]));

		nssArgs.setOutRoot(obArgs.getOutRoot());

		NSSCommandImpl nssImpl = new NSSCommandImpl();
		nssImpl.execute(nssArgs);
		N = nssImpl.getN();
	}

	private void step2() {
		//step 2
		Integer[] tmp = new Integer[obArgs.getCovNumber().length];
		for (int i = 0; i < obArgs.getCovNumber().length; i++) {
			tmp[i] = i + 2;
		}

		String Fout = obArgs.getOutRoot() + ".oath-bus";
		PrintStream OBout = FileUtil.CreatePrintStream(Fout);

		String OathList = obArgs.getOutRoot() +".oath-list";
		PrintStream OathListOut = FileUtil.CreatePrintStream(OathList);

		OathListOut.println(obArgs.getOutRoot() + ".p." + (obArgs.getSelectedPhenotype(0)+1) +".nss.gz");

		int cnt = 1;

		for (int ii = 1; ii <= tmp.length; ii++) {
			// ArrayList<Integer[]> rs=cmn(tmp,ii);
			int[][] rs = Combination.combination(tmp, ii);
			for (int i = 0; i < rs.length; i++) {
				ArrayList<String> mod = NewIt.newArrayList();
				mod.add("1");
				for (int j = 0; j < rs[i].length; j++) {
					mod.add(new Integer(tmp[rs[i][j]]).toString());
				}

				SynthCommandArguments synArgs = new SynthCommandArguments();
				synArgs.setGZ(true);
				synArgs.setNSSBatch(obArgs.getOutRoot() + ".list.nss");
				synArgs.setCMFile(obArgs.getOutRoot() + ".m.nss");
				synArgs.setKeepBatch(mod.toArray(new String[0]));
				synArgs.setN((new Integer(N)).toString());
				synArgs.setOutRoot(obArgs.getOutRoot() + "-" + (new Integer(cnt).toString()));

				SynthCommandImpl synImpl = new SynthCommandImpl();

				if (cnt == 1) {
					for (int j = 0; j < synArgs.getNSSFile().length; j++) {
						OBout.println("#" + (j + 1) + " " + synArgs.getNSSFile()[j]);
					}
				}

				String s = new String();
				s = "OATH-bus is running model [#1";
				OBout.print(obArgs.getOutRoot() + "-" + (new Integer(cnt).toString()) + ".oath.gz: #1");

				for (int k = 0; k < rs[i].length; k++) {
					if (k == 0) {
						OBout.print("=#" + tmp[rs[i][k]]);
						s += "=#" + tmp[rs[i][k]];
					} else {
						OBout.print("+#" + tmp[rs[i][k]]);
						s += "+#" + tmp[rs[i][k]];
					}
				}
				Logger.printUserLog(s + ']');
				OBout.println();

				synImpl.execute(synArgs);

				OathListOut.println(obArgs.getOutRoot() + "-" + (new Integer(cnt).toString())+".oath.gz");
				cnt++;
			}
		}
		OBout.close();
		
		OathListOut.close();

		Logger.printUserLog("Finding the summary information in '" + obArgs.getOutRoot() + ".oath-bus'.");

	}

	private void step3() {
		pkArgs = new OathPickCommandArguments();
		pkArgs.setPickList(obArgs.getOutRoot()+".oath-list");
		pkArgs.setOutRoot(obArgs.getOutRoot());
		OathPickCommandImpl pkImpl = new OathPickCommandImpl();
		pkImpl.execute(pkArgs);

	}
}