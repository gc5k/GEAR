
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


import edu.mit.wi.pedfile.PedFileException;

import score.CalEngineException;

import family.GMDRData;
import family.GMDRPhenoFileException;
import family.MDRPedFileException;

public class TestRabinowitz {

    public static void main(String[] args) throws IOException{
        
        File config = new File(args[0]);
        ArrayList<String> PedFileStrings = new ArrayList(9);
        BufferedReader Reader = new BufferedReader(new FileReader(config));
        String Line;
        while ((Line = Reader.readLine()) != null) {
            if (Line.length() == 0) {// skip blank Lines
                continue;
            }
            if (Line.startsWith("#")) {// skip comments
                continue;
            }
            PedFileStrings.add(Line);
        }
        /*
         * There are two circumstances, 1) When there are 9 parameters, build score with picked covariates and a
         * phenotype; 2) when there are 6 parameters, a score column has been listed in the phenotype file already.
         */
        int ParamNum = PedFileStrings.size();
        if (ParamNum != 9) {
            throw new IOException("Problems in reading paremeters." + ParamNum + "\n");
        }
        boolean IsPedigree = true;// false for case-control design; true for
        // family based design

        boolean Adjustment = true;
        int Reg = 0; // 0 for linear regression; 1 for logistic

        GMDRData GD = new GMDRData(IsPedigree);
        File PedFile = new File((String) PedFileStrings.get(0));// load the pedigree file

        File PhenoFile = new File((String) PedFileStrings.get(1));// load the phenotype file

        String OutPedFile = new String((String) PedFileStrings.get(2));// it is used to store the new content
        // after transformation of the pedigree
        // file

        String OutScoreFile = new String((String) PedFileStrings.get(3));// it is used to store the
        // score files

        // indexes of the covariats selected.
        String[] Cov = PedFileStrings.get(4).split("[,\\s]++");
        int CovIdx[] = stringTOinteger(Cov, 1);
        if (CovIdx != null && CovIdx.length > 1) {
            for (int j = 0; j < CovIdx.length; j++) {
                if (CovIdx[j] < 0) {
                    throw new IOException(
                            "The index of the selected predictor cannot be " + Cov[j]);
                }
            }
        }
        if (CovIdx.length == 1 && CovIdx[0] < 0) {
            Adjustment = false;
        }

        String[] Phe = PedFileStrings.get(5).split("[,\\s]++");
        int[] PheIdx = stringTOinteger(Phe, 1);
        for (int j = 0; j < PheIdx.length; j++) {
            if (PheIdx[j] < 0) {
                throw new IOException(
                        "The index of the selected response cannot be " + Phe[j]);
            }
            for (int jj = 0; jj < CovIdx.length; jj++) {
                if (PheIdx[j] == CovIdx[jj]) {
                    throw new IOException(
                            "The " + Phe[j] + "th phenotype cannot work as predictor and response simultaneously.");
                }
            }
        }
        Adjustment = (((String) PedFileStrings.get(6)).compareTo("true") == 0)? true:false;
        String[] RegParam = PedFileStrings.get(7).split("[,\\s]++");
        int[] regression = stringTOinteger(RegParam, -1); 

        Reg = regression[0];
        // 0 for statistics
        // built outside the
        // package; 1 for
        // linear
        // regression; 2 for
        // logistic        
        int replication = Integer.parseInt((String) PedFileStrings.get(8));

        boolean includeFounder = false;

        try {
            GD.InitialPedFile(PedFile);//initial Pedfile
        } catch (MDRPedFileException E) {
            E.printStackTrace(System.err);
        } catch (PedFileException E) {
            E.printStackTrace(System.err);
        }
        GD.Allele2Genotype();
        GD.GenotypeImputation();
        try {
            GD.InitialPhenoFile(PhenoFile);
        } catch (GMDRPhenoFileException e) {
            System.err.println("Phenotype File Exception.");
        }
        GD.Match();
        GD.RabinowitzApproach();
        GD.realCreateTable();
        try {
            if (Reg >= 0) {
                GD.buildScore2(PheIdx, CovIdx, Adjustment, Reg, includeFounder);
            } else {
                GD.fetchScore(PheIdx[0]);
            }
            GD.RabinowitzPrintGMDR(OutPedFile, OutScoreFile, false);
        } catch (CalEngineException e) {
            e.printStackTrace(System.err);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
        for (int i = 0; i < replication; i++) {
            String opfN = "Rabin_" + Integer.toString(i)+".txt";
            GD.RabinowitzApproach();
            GD.RabinowitzCreateTable();
            try {
                GD.RabinowitzPrintGMDR(opfN, null, true);
            } catch (CalEngineException e) {
                e.printStackTrace(System.err);
            } catch (Exception e) {
                e.printStackTrace(System.err);
            }
            System.out.println("Rabinowitz approach is simulating " + i);
        }
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
