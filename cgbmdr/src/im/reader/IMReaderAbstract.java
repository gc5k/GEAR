package im.reader;

import java.io.BufferedReader;
import java.io.IOException;

import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen
 */
public abstract class IMReaderAbstract {

    protected String file;
    protected String fileID;
    protected boolean isIndividualMarkers;
    protected ChromosomeParameter chrp;
    protected CrossParameter crp;
    protected BufferedReader buffer;
    protected String[][] markerName;
    protected double[][] distance;
    protected double[][] recombination;
    protected String[][][] marker;
    protected String[] phenoTitle;
    protected String[][] pheno;
    protected String[][] ancillaryPheno;

    public static class ChromosomeParameter {

        private String type;
        private int function;
        private String Unites;
        private int chromosome;
        private int maximum;
        private String named;

        public ChromosomeParameter() {
        }

        public void setParameter(String[] parameter) {
            type = parameter[IMConstant.chr_type];
            function = Integer.parseInt(parameter[IMConstant.chr_function]);
            Unites = parameter[IMConstant.chr_Units];
            chromosome = Integer.parseInt(parameter[IMConstant.chr_chromosomes]);
            maximum = Integer.parseInt(parameter[IMConstant.chr_maximum]);
            named = parameter[IMConstant.chr_named];
        }

        public String getType() {
            return type;
        }

        public int getFunction() {
            return function;
        }

        public String getUnites() {
            return Unites;
        }

        public int getChromosome() {
            return chromosome;
        }

        public int getMaximum() {
            return maximum;
        }

        public String getNamed() {
            return named;
        }
    }

    public static class CrossParameter {

        private int sampleSize;
        private String cross;
        private int trait;
        private int otrait;
        private String missingTrait;
        private String Case;
        private String[][] markercoding;

        public CrossParameter() {
        }

        public void setParameter(String[] parameter, String[] mc) {
            sampleSize = Integer.parseInt(parameter[IMConstant.crs_samplesize]);
            cross = parameter[IMConstant.crs_cross];
            trait = Integer.parseInt(parameter[IMConstant.crs_traits]);
            otrait = Integer.parseInt(parameter[IMConstant.crs_otraits]);
            missingTrait = parameter[IMConstant.crs_missingtrait];
            Case = parameter[IMConstant.crs_case];
            markercoding = new String[mc.length][3];
            for (int i = 0; i < markercoding.length; i++) {
                String[] m = mc[i].split("\\s+");
                System.arraycopy(m, 0, markercoding[i], 0, m.length);
            }
        }

        public int getSampleSize() {
            return sampleSize;
        }

        public String getCross() {
            return cross;
        }

        public int getTrait() {
            return trait;
        }

        public int getOTrait() {
            return otrait;
        }

        public String getMissingTrait() {
            return missingTrait;
        }

        public String getCase() {
            return Case;
        }

        public String[][] getMarkerCoding() {
            return markercoding;
        }
    }

    public IMReaderAbstract(String f) {
        file = f;
    }

    public void setMarker(ArrayList m) {
        if(!isIndividualMarkers) {
            marker = new String[m.size()][][];
            String[][] mc = crp.getMarkerCoding();
            for (int i = 0; i < m.size(); i++) {
                String mk = (String) m.get(i);
                String[] mkr = mk.split("\\s+");
                int idx = 0;
                marker[i] = new String[distance.length][];
                for (int j = 0; j < distance.length; j++) {
                    marker[i][j] = new String[distance[j].length];
                    for (int k = 0; k < marker[i][j].length; k++) {//read markers; if a marker's coding is out of definition, set it missed.
                        int v;
                        for (v = 0; v < mc.length; v++) {
                            if (mkr[k + idx].compareTo(mc[v][2]) == 0) {
                                break;
                            }
                        }
                        if (v < mc.length) {
                            marker[i][j][k] = mkr[k + idx];
                        } else {
                            marker[i][j][k] = mc[IMConstant.Missing][2];
                        } 
                    }
                    idx += distance[j].length;
                }
            }
        }
        else {
            // add something here when reading in marker by marker format
        }
    }

    public void setPhenotypeInformation(String pt, ArrayList pheValue) {
        String[] t = pt.split("\\s+");
        phenoTitle = new String[t.length];
        System.arraycopy(t, 0, phenoTitle, 0, t.length);

        ancillaryPheno = new String[pheValue.size()][IMConstant.ancillary_End];
        pheno = new String[pheValue.size()][phenoTitle.length - IMConstant.ancillary_End];
        for (int i = 0; i < pheValue.size(); i++) {
            String s = (String) pheValue.get(i);
            t = s.split("\\s+");
            for (int j = 0; j < IMConstant.ancillary_End; j++) {
                ancillaryPheno[i][j] = t[j];
            }
            for (int j = IMConstant.ancillary_End; j < t.length; j++) {
                pheno[i][j - IMConstant.ancillary_End] = t[j];
            }
        }
    }

    public String[][][] getMarkers() {
        return marker;
    }

    public double[][] getDistance() {
        return distance;
    }

    public double[][] getRecombination() {
        return recombination;
    }

    public String[][] getMarkerName() {
        return markerName;
    }

    public int MarkerNumber() {
        int s = 0;
        for(int i = 0; i < markerName.length; i++) {
            s += markerName[i].length;
        }
        return s;
    }

    public int IndividualNumber() {
        return pheno.length;
    }
    public String[][] getAncillaryPheno() {
        return ancillaryPheno;
    }

    public String[][] getPhenotype() {
        return pheno;
    }

    public ChromosomeParameter getChromosomeParameter() {
        return chrp;
    }

    public CrossParameter getCrossParameter() {
        return crp;
    }

    public String getFile() {
        return file;
    }

    public abstract void readData() throws IOException;
}
