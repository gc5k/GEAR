package gear.subcommands.projectedpc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import gear.ConstValues;

//import java.io.PrintStream;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwas.EigenGWASCommandArguments;
import gear.subcommands.eigengwas.EigenGWASCommandImpl;
import gear.subcommands.exsnp2.ExSNP2CommandArguments;
import gear.subcommands.exsnp2.ExSNP2CommandImpl;
import gear.subcommands.grm.GRMCommandArguments;
import gear.subcommands.grm.GRMCommandImpl;
import gear.subcommands.profile.CoeffModelType;
import gear.subcommands.profile.ProfileCommandArguments;
import gear.subcommands.profile.ProfileCommandImpl;
import gear.subcommands.qpca.QPCACommandArguments;
import gear.subcommands.qpca.QPCACommandImpl;
import gear.util.FileUtil;
//import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class ProjectedPCCommandImpl extends CommandImpl {

	private ProjectedPCCommandArguments proArgs;

	@Override
	public void execute(CommandArguments cmdArgs) {
		proArgs = (ProjectedPCCommandArguments) cmdArgs;

		ExSNP2CommandArguments exArgs = new ExSNP2CommandArguments();
		exArgs.setBFile(proArgs.getBFile());
		exArgs.setBatch(proArgs.getBatch());
		exArgs.setOutRoot(proArgs.getOutRoot());

		if (proArgs.isExtractFile()) {
			exArgs.setExtractFile(proArgs.getExtractFile());
		} else if (proArgs.isExcludeFile()) {
			exArgs.setExcludeFile(proArgs.getExcludeFile());
		} else if (proArgs.isChr()) {
			exArgs.setChr(proArgs.getChr().toArray(new String[0]));
		} else if (proArgs.isNotChr()) {
			exArgs.setNotChr(proArgs.getNotChr().toArray(new String[0]));
		}

		if (proArgs.isMAF()) {
			exArgs.setMAF((new Double(proArgs.getMAF())).toString());
		}
		if (proArgs.isMaxMAF()) {
			exArgs.setMaxMAF((new Double(proArgs.getMaxMAF())).toString());
		}
		if (proArgs.isGENO()) {
			exArgs.setGENO((new Double(proArgs.getGENO())).toString());
		}
		if (proArgs.isZeroVar()) {
			exArgs.setZeroVar();
		}
		if (proArgs.isMAFRange()) {
			exArgs.setMAFRange2(proArgs.getMAFRange());
		}

		if (proArgs.isKeepFile()) {
			exArgs.setKeepFile(proArgs.getKeepFile());
		} else if (proArgs.isRemoveFile()) {
			exArgs.setRemoveFile(proArgs.getRemoveFile());
		} else if (proArgs.isKeepFamFile()) {
			exArgs.setKeepFamFile(proArgs.getKeepFamFile());
		} else if (proArgs.isRemoveFamFile()) {
			exArgs.setRemoveFamFile(proArgs.getRemoveFamFile());
		}

		Logger.printUserLog("---------Extracting SNPs---------");
		ExSNP2CommandImpl exsnp2Impl = new ExSNP2CommandImpl();
		exsnp2Impl.execute(exArgs);

		
		GRMCommandArguments gArgs = new GRMCommandArguments();
		gArgs.setBFile(proArgs.getBFile());
		gArgs.setGZ();
		gArgs.setOutRoot(proArgs.getOutRoot());
		if (proArgs.isInbred()) {
			gArgs.setInbred();
		}
		// if (proArgs.isGUI()) {
		// gArgs.setGUI();
		// }

		gArgs.setExtractFile(proArgs.getOutRoot()+".snp");

		if (proArgs.isKeepFile()) {
			gArgs.setKeepFile(proArgs.getKeepFile());
		} else if (proArgs.isRemoveFile()) {
			gArgs.setRemoveFile(proArgs.getRemoveFile());
		} else if (proArgs.isKeepFamFile()) {
			gArgs.setKeepFamFile(proArgs.getKeepFamFile());
		} else if (proArgs.isRemoveFamFile()) {
			gArgs.setRemoveFamFile(proArgs.getRemoveFamFile());
		}

		Logger.printUserLog("");
		Logger.printUserLog("---------Generating GRM---------");
		GRMCommandImpl gImpl = new GRMCommandImpl();
		gImpl.execute(gArgs);

		// PrintStream gui_file = null;
		// if (proArgs.isGUI()) {
		// gui_file = FileUtil.CreatePrintStream(proArgs.getOutRoot()+".gui");
		// gui_file.println("s1");
		// }

		Logger.printUserLog("");
		Logger.printUserLog("---------Generating eigenvectors---------");

		QPCACommandArguments qpcaArgs = new QPCACommandArguments();
		qpcaArgs.setGrmGZ(proArgs.getOutRoot() + ".grm.gz");
		qpcaArgs.setGrmID(proArgs.getOutRoot() + ".grm.id");
		qpcaArgs.setEV(Integer.toString(proArgs.getEV()));
		qpcaArgs.setOutRoot(proArgs.getOutRoot());
		if (proArgs.isKeepFile()) {
			qpcaArgs.setKeepFile(proArgs.getKeepFile());
		} else if (proArgs.isRemoveFile()) {
			qpcaArgs.setRemoveFile(proArgs.getRemoveFile());
		}

		QPCACommandImpl qpcaImpl = new QPCACommandImpl();
		qpcaImpl.execute(qpcaArgs);
//
//		// if (eArgs.isGUI()) {
//		// gui_file.println("s2");
//		// }
//
		for (int i = 1; i <= proArgs.getEV(); i++) {
			EigenGWASCommandArguments eigenArgs = new EigenGWASCommandArguments();
			eigenArgs.setBFile(proArgs.getBFile());
			eigenArgs.setPhenotypeFile(proArgs.getOutRoot() + ".eigenvec");
			eigenArgs.setPhenotypeIndex(i);
			eigenArgs.setOutRoot(proArgs.getOutRoot() + "." + i);
			// if (proArgs.isGUI()) {
			// eigenArgs.setGUI();
			// }
			eigenArgs.setTAB();

			eigenArgs.setExtractFile(proArgs.getOutRoot()+".snp");

			if (proArgs.isKeepFile()) {
				eigenArgs.setKeepFile(proArgs.getKeepFile());
			} else if (proArgs.isRemoveFile()) {
				eigenArgs.setRemoveFile(proArgs.getRemoveFile());
			} else if (proArgs.isKeepFamFile()) {
				eigenArgs.setKeepFamFile(proArgs.getKeepFamFile());
			} else if (proArgs.isRemoveFamFile()) {
				eigenArgs.setRemoveFamFile(proArgs.getRemoveFamFile());
			}

			Logger.printUserLog("");
			Logger.printUserLog("---------Running EigenGWAS (Additive model) for the " + i + "th eigenvector.");
			EigenGWASCommandImpl eigenImpl = new EigenGWASCommandImpl();
			eigenImpl.execute(eigenArgs);
			// if (eArgs.isGUI()) {
			// gui_file.println("s3-" + i);
			// }
			Logger.printUserLog("Saved EigenGWAS results in '" + eigenArgs.getOutRoot() + ".egwas'.");
		}

		////make score file
		ArrayList< ArrayList<String> >  snpTab = NewIt.newArrayList();
		for (int i = 1; i <= proArgs.getEV(); i++) {
			BufferedReader eFile = FileUtil.FileOpen(proArgs.getOutRoot() + "." + i + ".egwas");
			String line;
			int cnt = 0;
			try {
				while ((line = eFile.readLine()) != null) {
					String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);

					cnt++;
					if(cnt==1) continue;

					if (i == 1) {
						ArrayList<String> sreader = NewIt.newArrayList();
						for (int j = 0; j < s.length; j++) {
							sreader.add(s[j]);
						}
						snpTab.add(sreader);
					} else {
						ArrayList<String> sreader = snpTab.get((cnt-2));
						sreader.add(s[2]);
						snpTab.set(cnt-2, sreader);
					}
				}
				eFile.close();
			} catch (IOException e) {
				Logger.handleException(e, "An exception occurred when reading '"
						+ proArgs.getOutRoot() + "." + i + ".egwas" + "'.");
			}
		}

		PrintStream proS = FileUtil.CreatePrintStream(proArgs.getOutRoot() + ".score");
		proS.print("SNP\tRefAllele");
		for(int i = 1; i <= proArgs.getEV(); i++) {
			proS.print("\tBeta"+i);
		}
		proS.println();

		for(int i = 0; i < snpTab.size(); i++) {
			ArrayList<String> s = snpTab.get(i);
			proS.print(s.get(0)+"\t" + s.get(1));	
			for(int j = 2; j < s.size(); j++) {
				proS.print("\t" + s.get(j));
			}
			proS.println();
		}
		proS.close();

		///generating score

		// if (eArgs.isGUI()) {
		// gui_file.close();
		// }

		for (int i = 0; i < proArgs.getAllBedFiles().size(); i++) {
			Logger.printUserLog("");
			Logger.printUserLog("------Generating profile score for " + proArgs.getAllBedFiles().get(i) + ".*-------");
			ProfileCommandArguments profArgs = new ProfileCommandArguments();

			profArgs.setCoeffModelType(CoeffModelType.ADDITIVE);
			profArgs.setIsSameAsPlink(true);
			profArgs.setBFile(proArgs.getAllBedFiles().get(i));
			profArgs.setScoreFile(proArgs.getOutRoot()+".score");
			profArgs.setHasScoreHeader(true);
			profArgs.setIsWeighted(false);
			profArgs.setResultFile(proArgs.getAllBedFiles().get(i));

			ProfileCommandImpl profImpl = new ProfileCommandImpl();
			profImpl.execute(profArgs);
			Logger.printUserLog("");
		}
	}
}
