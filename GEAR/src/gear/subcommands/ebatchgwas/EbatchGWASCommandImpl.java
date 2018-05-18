package gear.subcommands.ebatchgwas;

import java.io.PrintStream;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwas.EigenGWASCommandArguments;
import gear.subcommands.eigengwas.EigenGWASCommandImpl;
//import gear.subcommands.eigengwasdom.EigenGWASDomCommandArguments;
//import gear.subcommands.eigengwasdom.EigenGWASDomCommandImpl;
//import gear.subcommands.eigengwasepi.EigenGWASEpiCommandArguments;
import gear.subcommands.grm.GRMCommandArguments;
import gear.subcommands.grm.GRMCommandImpl;
import gear.subcommands.qpca.QPCACommandArguments;
import gear.subcommands.qpca.QPCACommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class EbatchGWASCommandImpl extends CommandImpl {
	@Override
	public void execute(CommandArguments cmdArgs) {
		EbatchGWASCommandArguments eArgs = (EbatchGWASCommandArguments) cmdArgs;

		GRMCommandArguments gArgs = new GRMCommandArguments();
		gArgs.setBFile(eArgs.getBFile());
		gArgs.setGZ();
		gArgs.setOutRoot(eArgs.getOutRoot());
		if (eArgs.isInbred()) {
			gArgs.setInbred();
		}
		if (eArgs.isGUI()) {
			gArgs.setGUI();
		}

		if (eArgs.isExtractFile()) {
			gArgs.setExtractFile(eArgs.getExtractFile());
		} else if (eArgs.isExcludeFile()) {
			gArgs.setExcludeFile(eArgs.getExcludeFile());
		} else if (eArgs.isChr()) {
			gArgs.setChr(eArgs.getChr().toArray(new String[0]));
		} else if (eArgs.isNotChr()) {
			gArgs.setNotChr(eArgs.getNotChr().toArray(new String[0]));
		}
		
		if (eArgs.isMAF()) {
			gArgs.setMAF((new Double(eArgs.getMAF())).toString());
		}
		if (eArgs.isMaxMAF()) {
			gArgs.setMaxMAF((new Double(eArgs.getMaxMAF())).toString());
		} 
		if (eArgs.isGENO()) {
			gArgs.setGENO((new Double(eArgs.getGENO())).toString());
		}
		if (eArgs.isZeroVar()) {
			gArgs.setZeroVar();
		}
		if (eArgs.isMAFRange()) {
			gArgs.setMAFRange2(eArgs.getMAFRange());
		}

		if (eArgs.isKeepFile()) {
			gArgs.setKeepFile(eArgs.getKeepFile());
		} else if (eArgs.isRemoveFile()) {
			gArgs.setRemoveFile(eArgs.getRemoveFile());
		} else if (eArgs.isKeepFamFile()) {
			gArgs.setKeepFamFile(eArgs.getKeepFamFile());
		} else if (eArgs.isRemoveFamFile()) {
			gArgs.setRemoveFamFile(eArgs.getRemoveFamFile());
		}

		GRMCommandImpl gImpl = new GRMCommandImpl();
		gImpl.execute(gArgs);

		PrintStream gui_file = null;

		if (eArgs.isGUI()) {
			gui_file = FileUtil.CreatePrintStream(eArgs.getOutRoot()+".gui");
			gui_file.println("s1");
		}

		Logger.printUserLog("\n---------Generating eigenvectors---------");

		QPCACommandArguments qpcaArgs = new QPCACommandArguments();
		qpcaArgs.setGrmGZ(eArgs.getOutRoot() + ".grm.gz");
		qpcaArgs.setGrmID(eArgs.getOutRoot() + ".grm.id");
		qpcaArgs.setEV(Integer.toString(eArgs.getEV()));
		qpcaArgs.setOutRoot(eArgs.getOutRoot());
		if (eArgs.isKeepFile()) {
			qpcaArgs.setKeepFile(eArgs.getKeepFile());
		} else if (eArgs.isRemoveFile()) {
			qpcaArgs.setRemoveFile(eArgs.getRemoveFile());
		}

		QPCACommandImpl qpcaImpl = new QPCACommandImpl();
		qpcaImpl.execute(qpcaArgs);

		if (eArgs.isGUI()) {
			gui_file.println("s2");
		}
		
		if (!eArgs.isDom() && !eArgs.isEpi()) {
			for (int i = 1; i <= eArgs.getEV(); i++) {
				Logger.printUserLog("\n---------Running EigenGWAS (Additive model) for the " + i + "th eigenvector.");
				EigenGWASCommandArguments eigenArgs = new EigenGWASCommandArguments();
				eigenArgs.setBFile(eArgs.getBFile());
				eigenArgs.setPhenotypeFile(eArgs.getOutRoot() + ".eigenvec");
				eigenArgs.setPhenotypeIndex(i);
				eigenArgs.setOutRoot(eArgs.getOutRoot() + "." + i);
				if (eArgs.isGUI()) {
					eigenArgs.setGUI();
				}

				if (eArgs.isExtractFile()) {
					eigenArgs.setExtractFile(eArgs.getExtractFile());
				} else if (eArgs.isExcludeFile()) {
					eigenArgs.setExcludeFile(eArgs.getExcludeFile());
				} else if (eArgs.isChr()) {
					eigenArgs.setChr(eArgs.getChr().toArray(new String[0]));
				} else if (eArgs.isNotChr()) {
					eigenArgs.setNotChr(eArgs.getNotChr().toArray(new String[0]));
				}

				if (eArgs.isMAF()) {
					eigenArgs.setMAF((new Double(eArgs.getMAF())).toString());
				}
				if (eArgs.isMaxMAF()) {
					eigenArgs.setMaxMAF((new Double(eArgs.getMaxMAF())).toString());
				}
				if (eArgs.isGENO()) {
					eigenArgs.setGENO((new Double(eArgs.getGENO())).toString());
				}
				if (eArgs.isZeroVar()) {
					eigenArgs.setZeroVar();
				}
				if (eArgs.isMAFRange()) {
					eigenArgs.setMAFRange2(eArgs.getMAFRange());
				}

				if (eArgs.isKeepFile()) {
					eigenArgs.setKeepFile(eArgs.getKeepFile());
				} else if (eArgs.isRemoveFile()) {
					eigenArgs.setRemoveFile(eArgs.getRemoveFile());
				}else if (eArgs.isKeepFamFile()) {
					eigenArgs.setKeepFamFile(eArgs.getKeepFamFile());
				} else if (eArgs.isRemoveFamFile()) {
					eigenArgs.setRemoveFamFile(eArgs.getRemoveFamFile());
				}

				EigenGWASCommandImpl eigenImpl = new EigenGWASCommandImpl();
				eigenImpl.execute(eigenArgs);
				if (eArgs.isGUI()) {
					gui_file.println("s3-" + i);
				}
				Logger.printUserLog("Saved EigenGWAS results in '" + eigenArgs.getOutRoot() + ".egwas'.");
			}
		}
		
		if (eArgs.isGUI()) {
			gui_file.close();
		}
		// else if (eArgs.isEpi())
		// {
		// for(int i = 1; i <= eArgs.getEV(); i++)
		// {
		// Logger.printUserLog("\n---------Running EigenGWAS (Add+Dom+AA model) for the
		// "+i+"th eigenvector.");
		// EigenGWASEpiCommandArguments eigenEpiArgs = new
		// EigenGWASEpiCommandArguments();
		// eigenEpiArgs.setBFile(eArgs.getBFile());
		// eigenEpiArgs.setPhenotypeFile(eArgs.getOutRoot()+".eigenvec");
		// eigenEpiArgs.setPhentypeIndex(i);
		// eigenEpiArgs.setOutRoot(eArgs.getOutRoot() + "." + i);
		// if (eArgs.isInbred())
		// {
		// eigenEpiArgs.setInbred();
		// }
		//
		// EigenGWASDomImpl eigenDomImpl = new EigenGWASDomImpl();
		// eigenDomImpl.execute(eigenEpiArgs);
		// Logger.printUserLog("Saved EigenGWAS results in '"+eigenEpiArgs.getOutRoot()
		// + ".egwasepi'.");
		// }
		// }
		// else if (eArgs.isDom())
		// {
		// for(int i = 1; i <= eArgs.getEV(); i++)
		// {
		// Logger.printUserLog("\n---------Running EigenGWAS (Additive+Dominance model)
		// for the "+i+"th eigenvector.");
		// EigenGWASDomCommandArguments eigenDomArgs = new
		// EigenGWASDomCommandArguments();
		// eigenDomArgs.setBFile(eArgs.getBFile());
		// eigenDomArgs.setPhenotypeFile(eArgs.getOutRoot()+".eigenvec");
		// eigenDomArgs.setPhentypeIndex(i);
		// eigenDomArgs.setOutRoot(eArgs.getOutRoot() + "." + i);
		//
		// EigenGWASDomImpl eigenDomImpl = new EigenGWASDomImpl();
		// eigenDomImpl.execute(eigenDomArgs);
		// Logger.printUserLog("Saved EigenGWAS results in '"+eigenDomArgs.getOutRoot()
		// + ".egwasd'.");
		// }
		// }
	}
}
