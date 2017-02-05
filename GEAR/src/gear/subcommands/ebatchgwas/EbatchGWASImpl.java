package gear.subcommands.ebatchgwas;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwas.EigenGWASArguments;
import gear.subcommands.eigengwas.EigenGWASImpl;
import gear.subcommands.grm.GRMArguments;
import gear.subcommands.grm.GRMImpl;
import gear.subcommands.qpca.QPCACommandArguments;
import gear.subcommands.qpca.QPCACommandImpl;
import gear.util.Logger;

public class EbatchGWASImpl extends CommandImpl 
{
	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		EbatchGWASArguments eArgs = (EbatchGWASArguments) cmdArgs;

		GRMArguments gArgs = new GRMArguments();
		gArgs.setBFile(eArgs.getBFile());
		gArgs.setGZ();
		GRMImpl gImpl = new GRMImpl();
		gImpl.execute(gArgs);

		Logger.printUserLog("\n---------Generating eigenvectors---------");

		QPCACommandArguments qpcaArgs = new QPCACommandArguments();
		qpcaArgs.setGrmGZ(eArgs.getOutRoot() + ".grm.gz");
		qpcaArgs.setGrmID(eArgs.getOutRoot() + ".grm.id");
		qpcaArgs.setEV(Integer.toString(eArgs.getEV()));
		qpcaArgs.setOutRoot(eArgs.getOutRoot());
		QPCACommandImpl qpcaImpl = new QPCACommandImpl();
		qpcaImpl.execute(qpcaArgs);
		Logger.printUserLog("Saving the top "+ eArgs.getEV() + " eigenvectors in '" + eArgs.getOutRoot() + ".eigenvec'.");

		for(int i = 1; i <= eArgs.getEV(); i++)
		{
			Logger.printUserLog("\n---------Running EigenGWAS for the "+i+"th eigenvector.");
			EigenGWASArguments eigenArgs = new EigenGWASArguments();
			eigenArgs.setBFile(eArgs.getBFile());
			eigenArgs.setPhenotypeFile(eArgs.getOutRoot()+".eigenvec");
			eigenArgs.setPhentypeIndex(i);
			eigenArgs.setOutRoot(eArgs.getOutRoot() + "." + i);
			
			EigenGWASImpl eigenImpl = new EigenGWASImpl();
			eigenImpl.execute(eigenArgs);
			Logger.printUserLog("Saved EigenGWAS results in '"+eigenArgs.getOutRoot() + ".egwas'.");
		}
		
	}
}
