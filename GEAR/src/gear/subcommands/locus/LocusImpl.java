package gear.subcommands.locus;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import gear.family.pedigree.file.PedigreeFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.pop.PopStat;

public class LocusImpl extends CommandImpl
{
	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private double[] allelevar;
	private PedigreeFile pf;
	private ArrayList<SNP> snpList;
	private LocusArguments locusArgs;
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		this.locusArgs = (LocusArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(this.locusArgs);
		pf = pp.getPedigreeData();
		snpList = pp.getMapData().getMarkerList();

		SampleFilter sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());
		SumStatQC ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(),
				sf);
		GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
		G = gm;
		numMarker = G.getNumMarker();
		allelefreq = new double[numMarker][3];
		
		allelefreq = PopStat.calAlleleFrequency(G, numMarker);
		allelevar = PopStat.calGenoVariance(G, numMarker);
		printResult();

	}

	private void printResult()
	{
		DecimalFormat fmt = new DecimalFormat("0.0000");
		DecimalFormat fmtp = new DecimalFormat("0.00E000");
		PrintStream LocusPrint = FileUtil.CreatePrintStream(this.locusArgs.getOutRoot() + ".locus");
		LocusPrint.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tFreq\tVar");

		for(int i = 0; i < snpList.size(); i++)
		{
			SNP snp = (SNP) snpList.get(i);
			LocusPrint.print(snp.getName() + "\t"+snp.getChromosome()+"\t"+snp.getPosition()+"\t"+snp.getFirstAllele()+ "\t"+snp.getSecAllele()+"\t");
			if (Math.abs(allelefreq[i][0]) > 0.0001)
			{
				LocusPrint.print(fmt.format(allelefreq[i][0])+"\t");
			}
			else
			{
				LocusPrint.print(fmtp.format(allelefreq[i][0])+"\t");				
			}
			if (Math.abs(allelevar[i]) > 0.0001)
			{
				LocusPrint.println(fmt.format(allelevar[i]));
			}
			else
			{
				LocusPrint.println(fmtp.format(allelevar[i]));				
			}
		}
		LocusPrint.close();
	}
}
