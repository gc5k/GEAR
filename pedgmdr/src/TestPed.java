
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Pattern;
import score.CalEngineException;
import edu.mit.wi.pedfile.PedFileException;
import family.GMDRData;
import family.GMDRPhenoFileException;
import family.MDRPedFileException;

/**
 * It's a beta version of the PedGMDR, which produces nontransmitted genotypes and scores. TestPed is mainly developed
 * to realize the first step of PedGMDR, and the files produced by TestPed can be exported to GMDR.jar.
 * 
 * @author Guo-Bo Chen, Zhejiang University, chenguobo@gmail.com
 */
public class TestPed {

	public static class Parameter {
		protected ArrayList<String> lines;
		protected BufferedReader buffer;
		protected String filename;
		protected boolean isPedigree; // 1
		protected String ped_file;  //2
		protected String phe_file;  //3
		protected String converted_ped_file;  //4
		protected String converted_phe_file;  //5
		protected String id_file;  //6
		protected int[] cov_idx;  //7 starts from 1
		protected int[] phe_idx;  //8 starts from 1
		protected int method; //9; 1 for linear, 2 for logistic
		protected boolean adjustment; //10
		protected boolean includeFounder; //11
		protected int replication; //12

		public Parameter() {
			lines = new ArrayList();
		}

		public void read(String file) throws IOException {
			filename = file;
			buffer = new BufferedReader(new FileReader(new File(file)));
			sweepComments();
			parseValue();
		}

		public void sweepComments() throws IOException {
			boolean flag = true;
			String line;
			while ((line = buffer.readLine()) != null) {
				if (Pattern.matches("^//.*", line)) {// empty line
					continue;
				} else {
					lines.add(line);
				}
			}
		}

		public void parseValue() {
			isPedigree = Boolean.parseBoolean(lines.get(0));
			ped_file = lines.get(1);
			phe_file = lines.get(2);
			converted_ped_file = lines.get(3);
			converted_phe_file = lines.get(4);
			id_file = lines.get(5);
			String[] c = lines.get(6).split(",");
			cov_idx = new int[c.length];
			for( int i = 0; i < c.length; i++) {
				cov_idx[i] = Integer.parseInt(c[i]) - 1;
			}
			String[] p = lines.get(7).split(",");
			phe_idx = new int[p.length];
			for( int i = 0; i < p.length; i++) {
				phe_idx[i] = Integer.parseInt(p[i]) - 1;
			}
			method = Integer.parseInt(lines.get(8));
			adjustment = Boolean.parseBoolean(lines.get(9));
			includeFounder = Boolean.parseBoolean(lines.get(10));
			replication = Integer.parseInt(lines.get(11));
		}
	}
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
    	
		Parameter pr = new Parameter();
		if(args.length > 0) {
			pr.read(args[0]);
		} else {
			pr.isPedigree = true;
			pr.ped_file = "0.ped";
			pr.phe_file = "0.phe";
			pr.converted_ped_file = "Converted_0.ped";
			pr.converted_phe_file = "Converted_0.phe";
			pr.id_file = "Family_ID.txt";
			pr.cov_idx = new int[1];
			pr.cov_idx[0] = 1;
			pr.phe_idx = new int[1];
			pr.phe_idx[0] = 0;
			pr.method = 1;
			pr.adjustment = true;
			pr.includeFounder = false;
			pr.replication = 10;
		}

		GMDRData GD = new GMDRData(pr.isPedigree);
        int WhichDataSet = 0; // 0 for all; 1 for affected only ; need add a

        try {
            GD.InitialPedFile(pr.ped_file);
        } catch (MDRPedFileException e) {
            System.err.println("Pedigree File Exception.");
            e.printStackTrace(System.err);
        } catch (PedFileException e) {
            System.err.println("Pedigree File checking Excetpion.");
            e.printStackTrace(System.err);
        }
        GD.NonTransmittedGenoType();
        try {
            GD.InitialPhenoFile(pr.phe_file);
        } catch (GMDRPhenoFileException e) {
            System.err.println("Phenotype File Exception.");
            e.printStackTrace(System.err);
        }
        GD.Match();
        GD.realCreateTable();

        try {
            if (pr.method > 0) {
                GD.buildScore2(pr.phe_idx, pr.cov_idx, pr.adjustment, pr.method, pr.includeFounder);
            } else {
                GD.fetchScore(pr.phe_idx[0]);
            }
            GD.realPrintGMDR(pr.converted_ped_file, pr.converted_phe_file, pr.id_file, WhichDataSet);
        } catch (CalEngineException e) {
            e.printStackTrace(System.err);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }

		for (int i = 0; i < pr.replication; i++) {
			String opfN = "Lou_" + Integer.toString(i) + ".txt";
			try {
				GD.PrintNullGMDR(opfN, WhichDataSet, (new Long (i)).longValue());
			} catch(Exception E) {
				E.printStackTrace(System.err);
			}
		}
    }
}