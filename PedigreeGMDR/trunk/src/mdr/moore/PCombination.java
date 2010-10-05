
package mdr.moore;

import java.util.HashMap;

/**
 *
 * @author Guo-Bo Chen
 */
public class PCombination extends HashMap {

    protected int numTraits;

    public PCombination(int sI) {
        numTraits = sI;
    }

    public void addPSuite(String key, PSuite s) {
        put(key, s);
    }
}
