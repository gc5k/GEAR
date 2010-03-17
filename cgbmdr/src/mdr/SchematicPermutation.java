package mdr;

import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Random;
import im.population.IMPopulation;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class SchematicPermutation {

    long seed;
    IMPopulation imp;
    ArrayList permutationIDs = new ArrayList();

    public SchematicPermutation(IMPopulation Imp) {
        imp = Imp;
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            permutationIDs.add(new Integer(i));
        }
    }

    public ArrayList getShuffledIDs(long sd) {
        Random rnd = new Random(sd);

        int ogi = imp.MarkerNumber();
        int env = imp.getEnvironment();
        int replication = imp.getReplication();

        Pattern cross = Pattern.compile("if2", Pattern.CASE_INSENSITIVE);
        Matcher match = cross.matcher(imp.getCrossParameter());
        if (match.matches()) {
            ArrayList PI = new ArrayList();
            for (int i = 0; i < env; i++) {
                for (int j = 0; j < ogi; j++) {
                    PI.add(new Integer(ogi * i * replication + j));
                }
            }

            Collections.shuffle(PI, rnd);
            for (int i = 0; i < env; i++) {
                for (int j = 0; j < ogi; j++) {
                    Integer I = (Integer) PI.get(i * ogi + j);
                    for (int k = 0; k < replication; k++) {
                        int idx = i * replication * ogi + k * ogi + j;
                        int v = I.intValue() + k * ogi;
                        permutationIDs.set(idx, new Integer(v));
                    }
                }
            }
        } else {
            if (permutationIDs != null) {
                Collections.shuffle(permutationIDs, rnd);
            }
        }
        return permutationIDs;
    }
}
