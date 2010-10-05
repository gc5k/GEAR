
import java.io.*;
import java.util.ArrayList;
import java.util.regex.Pattern;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.moore.LinearMergeSearch;
import mdr.moore.Permutation;
import mdr.moore.PermutationResult;
import mdr.SavedModels;
import mdr.moore.AbstractMergeSearch;

/**
 *
 * @author Guo-Bo Chen
 */
public class main {

	public static class Parameter {
		protected ArrayList<String> lines;
		protected BufferedReader buffer;		
		protected String filename;
		protected String ped_file;
		protected String phe_file;
		protected int interaction_from;
		protected int interaction_end;
		protected int[] scr_idx;
		protected int interval;
		protected long seed;
		protected int partition_method;
		protected boolean isMooreMDR;
		protected int replication_permutation;
		protected int search_method;

		public Parameter() {
			
		}

		public void read(String file) throws IOException {
			filename = file;
			lines = new ArrayList();
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
			ped_file = lines.get(0);
			phe_file = lines.get(1);

			interaction_from = Integer.parseInt(lines.get(2));
			interaction_end = Integer.parseInt(lines.get(3));
			String[] c = lines.get(4).split(",");
			scr_idx = new int[c.length];
			for( int i = 0; i < c.length; i++) {
				scr_idx[i] = Integer.parseInt(c[i]) - 1;
			}
			interval = Integer.parseInt(lines.get(5));
			seed = Long.parseLong(lines.get(6));
			partition_method = Integer.parseInt(lines.get(7));
			isMooreMDR = Boolean.parseBoolean(lines.get(8));
			replication_permutation = Integer.parseInt(lines.get(9));
			search_method = Integer.parseInt(lines.get(10));
		}		
	}

    public static void main(String[] args) throws IOException {
    	Parameter pr = new Parameter();
    	if (args.length > 0) {
    		pr.read(args[0]);
    	} else {
    		pr.ped_file = "ped0.txt"; 
        	pr.phe_file = "phe0.txt";
        	pr.interaction_from = 1;
        	pr.interaction_end = 2;
        	pr.scr_idx = new int[1];
        	pr.scr_idx[0] = 1 - 1;
        	pr.interval = 10;
        	pr.seed = 7;
        	pr.partition_method = 0;
        	pr.isMooreMDR = false;
        	pr.replication_permutation = 100;
        	pr.search_method = 0;
    	}

//        pFile = null;
        DataFile dr = new DataFile(pr.ped_file, pr.phe_file, pr.scr_idx);

        Subdivision sd = new Subdivision(pr.interval, pr.seed, dr);
        sd.RandomPartition();

        CombinationGenerator cg = new CombinationGenerator(pr.interaction_from, pr.interaction_end, dr.getMarkerNum());
        cg.generateCombination();
        AbstractMergeSearch as = new LinearMergeSearch(dr, sd, cg, pr.scr_idx.length, dr.getOffset(), pr.isMooreMDR);
        for (int i = pr.interaction_from; i <= pr.interaction_end; i++) {
            as.search(i);
            as.summarise();
        }

        for (int i = pr.interaction_from; i <= pr.interaction_end; i++) {
            SavedModels sm = as.getBestSavedModelAtOrder(new Integer(i));
            Permutation pm = new Permutation(dr, sm, dr.getOffset(), sd, pr.replication_permutation);
            pm.evaluate();
            PermutationResult permu_res = pm.getResult(0);
            permu_res.summarise();
            System.out.println(permu_res);
        }
    }
}
