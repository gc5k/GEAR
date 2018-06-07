package gear.subcommands.fst;

import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;

public class FSTCommandImpl extends CommandImpl {

	private GenotypeMatrix pGM;
	private FSTCommandArguments fstArgs;
	@Override
	public void execute(CommandArguments cmdArgs) {
		fstArgs = (FSTCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(fstArgs);
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);
		
		Logger.printUserLog("Estimating Fst for " + fstArgs.getBFile());
	}

}
