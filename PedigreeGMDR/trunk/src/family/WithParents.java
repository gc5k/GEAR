package family;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Pattern;

import edu.mit.wi.pedfile.PedFileException;

import score.CalEngineException;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.GMDRData;
import family.pedigree.GMDRPhenoFileException;
import family.pedigree.MDRPedFileException;

public class WithParents {

	public static void main(String[] args) throws IOException {

		PedigreeParameter pr = new PedigreeParameter();
		if(args.length > 0) {
			pr.read(args[0]);
		} else {
			pr.setIsLouAlgorithm(false);
			pr.setUsingFounderGenotype(true);
			pr.setPedigreeFile("0.ped");
			pr.setPhenotypeFile("0.phe");
			pr.setConvertedPedigreeFile("Converted_0.ped");
			pr.setConvertedPhenotypeFile("Converted_0.phe");
			pr.setFamilyIDFile("Family_ID_0.txt");
			pr.setCovariateIndex("2");
			pr.setPhenotypeIndex("0");
			pr.setScoreBuildMethod(2);
			pr.setAdjustScore(true);
			pr.setScoreBuildWithFounder(true);
			pr.setReplication(10);
			pr.setSeed(10);
		}
		GMDRData GD = new GMDRData(pr.usingFounderGenotype(), pr.isLouAlgorithm());
		AbstractGenoDistribution.rnd = new Random(pr.seed);
		boolean isRabinowitzProc = false;
		GD.RevvingUp(pr.getPedigreeFile(), pr.getPhenotypeFile());
		GD.RabinowitzApproach();
		GD.CreateTable(isRabinowitzProc);
		if (pr.getScoreBuildMethod() >= 0) {
			GD.buildScore(pr.getPhenotypeIndex(), pr.getCovarianteIndex(), pr.isAdjustScore(), pr.getScoreBuildMethod(), pr.getScoreBuildWithFounder());
		} else {
			GD.fetchScore(pr.getPhenotypeIndex(), pr.getScoreBuildWithFounder());
		}
		GD.PrintGMDR(pr.getConvertedPedigreeFile(), pr.getConvertedPhenotypeFile(), false);
		isRabinowitzProc = true;
		for (int i = 0; i < pr.getReplication(); i++) {
			String opfN = "Rabin_" + Integer.toString(i) + ".txt";
			GD.RabinowitzApproach();
			GD.CreateTable(isRabinowitzProc);
			GD.PrintGMDR(opfN, null, true);
		}
	}
}
