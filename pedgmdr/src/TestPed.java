
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import score.CalEngineException;
import edu.mit.wi.pedfile.PedFileException;
import family.GMDRData;
import family.GMDRPhenoFileException;
import family.MDRPedFileException;

/**
 * It's a beta version of the PedGMDR, which produces nontransmitted genotypes and scores. TestPed is mainly developed
 * to realize the first step of PedGMDR, and the files produced by TestPed can be exported to GMDR.jar.
 * 
 * @author Guo-Bo Chen, Zhejiang University, gc5k@zju.edu.cn
 */
public class TestPed {

    /**
     * It reads parameters from the configuration file, which is given in the command line. When parameters are provided
     * correctly, the intermediate files will be produced by invoking those relevant classes. An IOException is thrown
     * when the number of parameters is not correct.
     * 
     * @param args
     *            a configuration file is required to provide parameters.
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        File config = new File(args[0]);
        ArrayList<String> PedFileStrings = new ArrayList(8);
        BufferedReader Reader = new BufferedReader(new FileReader(config));
        String Line;
        while ((Line = Reader.readLine()) != null) {
            if (Line.length() == 0)// skip blank Lines
            {
                continue;
            }
            if (Line.startsWith("#"))// skip comments
            {
                continue;
            }
            PedFileStrings.add(Line);
        }
        /*
         * There are two circumstances, 1) When there are 8 parameters, build score with picked covariates and a
         * phenotype; 2) when there are 6 parameters, a score column has been listed in the phenotype file already.
         */
        int ParamNum = PedFileStrings.size();
        if (ParamNum != 8) {
            throw new IOException("Problems in reading paremeters.");
        }
        boolean IsPedigree = true;// false for case-control design; true for
        // family based design

        int WhichDataSet = 0; // 0 for all; 1 for affected only ; need add a
        // control at the GMDR software to accommodate
        // this feature

        boolean Adjustment = true;
        int[] SelectedPheIndex = null;
        int Reg = 0; // 0 for linear regression; 1 for logistic

        GMDRData GD = new GMDRData(IsPedigree);
        File PedFile = new File((String) PedFileStrings.get(0));// load the pedigree file

        File PhenoFile = new File((String) PedFileStrings.get(1));// load the phenotype file

        String OutPedFile = new String((String) PedFileStrings.get(2));// it is used to store the new content
        // after transformation of the pedigree
        // file

        String OutScoreFile = new String((String) PedFileStrings.get(3));// it is used to store the
        // score files

        String OutIDFile = new String((String) PedFileStrings.get(4));// it is used to store the IDs of selected
        // individuals.

        // indexes of the covariats selected.
        String[] Cov = PedFileStrings.get(5).split("[,\\s]++");
        int SelectedCovIndex[] = stringTOinteger(Cov, 1);
        if (SelectedCovIndex != null && SelectedCovIndex.length > 1) {
            for (int j = 0; j < SelectedCovIndex.length; j++) {
                if (SelectedCovIndex[j] < 0) {
                    throw new IOException(
                            "The index of the selected predictor cannot be " + Cov[j]);
                }
            }
        }
        if (SelectedCovIndex.length == 1 && SelectedCovIndex[0] < 0) {
            Adjustment = false;
        }

        String[] Phe = PedFileStrings.get(6).split("[,\\s]++");
        SelectedPheIndex = stringTOinteger(Phe, 1);
        for (int j = 0; j < SelectedPheIndex.length; j++) {
            if (SelectedPheIndex[j] < -1) {
                throw new IOException(
                        "The index of the selected response cannot be " + Phe[j]);
            }
            for (int jj = 0; jj < SelectedCovIndex.length; jj++) {
                if (SelectedPheIndex[j] == SelectedCovIndex[jj]) {
                    throw new IOException(
                            "The " + Phe[j] + "th phenotype cannot work as predictor and response simultaneously.");
                }
            }
        }

        String[] RegParam = PedFileStrings.get(7).split("[,\\s]++");
        int[] regression = stringTOinteger(RegParam, 0); // 0 for statistics
        // built outside the
        // package; 1 for
        // linear
        // regression; 2 for
        // logistic

        Reg = regression[0];
        if (Reg < 0 || Reg > 2) {
            throw new IOException(
                    "incorrect parameter for generalized linear model.");
        }
        try {
            GD.InitialPedFile(PedFile);// initial Pedfile

        } catch (MDRPedFileException e) {
            System.err.println("Pedigree File Exception.");
            e.printStackTrace(System.err);
        } catch (PedFileException e) {
            System.err.println("Pedigree File checking Excetpion.");
            e.printStackTrace(System.err);
        }
        GD.Allele2Genotype();
        GD.GenotypeImputation();
        GD.NonTransmittedGenoType();

        try {
            GD.InitialPhenoFile(PhenoFile);
        } catch (GMDRPhenoFileException e) {
            System.err.println("Phenotype File Exception.");
            e.printStackTrace(System.err);
        }
        GD.Match();
        GD.realCreateTable();

        try {
            if (Reg > 0) {
                GD.buildScore(SelectedPheIndex, SelectedCovIndex, Adjustment,
                        Reg);
            } else {
                GD.fetchScore(SelectedPheIndex[0]);
            }
            GD.realPrintGMDR(OutPedFile, OutScoreFile, OutIDFile, WhichDataSet);
        } catch (CalEngineException e) {
            e.printStackTrace(System.err);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }

        ArrayList PopulationStatistics = GD.getPopulationStatistics();// get
        // population
        // statistics

        ArrayList MendError = GD.getMendError();// get Mendelian Error information

        Hashtable FamInformative = GD.getFamilyInformative();// get family
    // informativeness
        /*
     * Vector unmatched = GD.getUnMatched(); Vector matched = GD.getMatched(); System.out.println("unMatched: " +
     * unmatched.size()); System.out.println("Matched: " + matched.size());
     * System.out.println("Families: "+GD.PedData.getNumFamilies()); Enumeration fkeys = FamInformative.keys();
     * while( fkeys.hasMoreElements() ) { String fid = (String) fkeys.nextElement(); Boolean ft = (Boolean)
     * FamInformative.get( fid ); if( !ft.booleanValue() ) { System.out.println(fid); } } Iterator it =
     * PopulationStatistics.iterator(); for(; it.hasNext(); ) { Vector item = (Vector) it.next(); Iterator iit =
     * item.iterator(); for(; iit.hasNext(); ) { System.out.print(iit.next()+"\t"); } System.out.println(); }
     * System.out.println( "No Mendelian Error = " + MendError.size() ); Iterator It = MendError.iterator(); for(;
     * It.hasNext(); ) { MendErrorTrace MET = (MendErrorTrace) It.next();
     * System.out.println(MET.FID+" "+MET.PID+" "+MET.Marker); }
     */

    }

    /**
     * Covert a string array to an integer array with an offset.
     * 
     * @param s         The String array to be converted
     * @param offset    The offset
     * @return Converted integer array
     */
    private static int[] stringTOinteger(String[] s, int offset) {
        if (s == null) {
            return null;
        }
        int[] d = new int[s.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = Integer.parseInt(s[i]) - offset;
        }
        return d;
    }
}