package family;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;
import edu.mit.wi.pedfile.PedFileException;

import score.CalEngineException;

import family.pedigree.GMDRData;
import family.pedigree.GMDRPhenoFileException;
import family.pedigree.MDRPedFileException;

public class InnerTest {

	public static class Parameter {
		protected ArrayList<String> lines;
		protected BufferedReader buffer;
		protected String filename;
		protected boolean isPedigree; // 1
		protected String ped_file;  // 2
		protected String phe_file;  // 3
		protected String converted_ped_file;  // 4
		protected String converted_phe_file;  // 5
		protected String id_file;  // 6
		protected int[] cov_idx;  // 7 starts from 1
		protected int[] phe_idx;  // 8 starts from 1
		protected int method; // 9; 1 for linear, 2 for logistic
		protected boolean adjustment; // 10
		protected boolean includeFounder; // 11
		protected int replication; // 12

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
				if (Pattern.matches("^\\s*//*.*", line)) {// empty line
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

    public static void main(String[] args) throws IOException{
		String[] arg = new String[1];
		arg[0] = "config.txt";
		Parameter pr = new Parameter();
		if(arg.length > 0) {
			pr.read(arg[0]);
			pr.sweepComments();
			pr.parseValue();
		} else {
			pr.isPedigree = true;
			pr.ped_file = "0.ped";
			pr.phe_file = "0.phe";
			pr.converted_ped_file = "Converted_0.ped";
			pr.converted_phe_file = "Converted_0.phe";
			pr.id_file = "Family_ID.txt";
			pr.cov_idx = new int[1];
			pr.cov_idx[0] = 2;
			pr.phe_idx = new int[1];
			pr.phe_idx[0] = 1;
			pr.method = 1;
			pr.adjustment = true;
			pr.includeFounder = false;
			pr.replication = 1000;
		}

		GMDRData GD = new GMDRData(pr.isPedigree);

        boolean includeFounder = false;

        try {
            GD.InitialPedFile(pr.ped_file);
        } catch (MDRPedFileException E) {
            E.printStackTrace(System.err);
        }
        GD.Allele2Genotype();
        try {
            GD.InitialPhenoFile(pr.phe_file);
        } catch (GMDRPhenoFileException e) {
            System.err.println("Phenotype File Exception.");
        }
        GD.Match();
        GD.RabinowitzApproach();
        GD.realCreateTable();
        try {
            if (pr.method >= 0) {
                GD.buildScore(pr.phe_idx, pr.cov_idx, pr.adjustment, pr.method, pr.includeFounder);
            } else {
                GD.fetchScore(pr.phe_idx[0]);
            }
            GD.RabinowitzPrintGMDR(pr.converted_ped_file, pr.converted_phe_file, false);
        } catch (CalEngineException e) {
            e.printStackTrace(System.err);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
        for (int i = 0; i < pr.replication; i++) {
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
}
