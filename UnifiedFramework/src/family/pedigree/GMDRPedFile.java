package family.pedigree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.Vector;

import edu.mit.wi.pedfile.MarkerResult;
import edu.mit.wi.pedfile.PedFile;
import edu.mit.wi.pedfile.PedFileException;

public class GMDRPedFile extends PedFile {

    /**
     * GMDRPedFile derived from pedfile
     */
    private int[] GMDRmarkerRatings;
    private int[] GMDRdups;
    private ArrayList markerInfor;
    private Vector pedgreeStrings;
    private ArrayList pedFileStrings;
    private String titleLine;
    private File pedfile;

    /**
     * Constructor of GMDRPedFile
     */
    public GMDRPedFile() {
        super();

    }

    /**
     * Initial of GMDRPedFile object (File infile) pedfile
     */
    public void Initial(File infile) throws GMDRPedFileException, IOException {
        pedFileStrings = new ArrayList();
        BufferedReader reader = new BufferedReader(new FileReader(infile));
        pedfile = infile;
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.length() == 0) {
                // skip blank lines
                continue;
            }
            if (line.startsWith("#")) {
                // skip comments
                continue;
            }
            pedFileStrings.add(line);
        }
        pedgreeStrings = (Vector) pedFileStrings.clone();
        titleLine = (String) pedgreeStrings.get(0);
        pedgreeStrings.remove(0);
        int numLines = pedFileStrings.size();

        if (numLines < 2) {
            throw new GMDRPedFileException(
                    "Pedgree data format error: empty file");
        }
        StringTokenizer tokenizer = new StringTokenizer(titleLine, "\n\t\" \"");
        int numTokens = tokenizer.countTokens();

        if (numTokens < 7) {
            throw new GMDRPedFileException(
                    "Pedgree format error: the title line is incorrect");
        }

        // reading the title line:get the marker number
        int c = 0;
        ArrayList temp = new ArrayList();
        while (tokenizer.hasMoreTokens()) {
            if (c++ < 6) {
                tokenizer.nextToken();
            } else {
                String marker = (String) tokenizer.nextToken();
                temp.add(marker);
            }
        }
        markerInfor = (ArrayList) temp.clone();
    }

    public Vector getPedgreeStrings() {
        return pedgreeStrings;
    }

    public ArrayList getMarkerInfor() {
        return markerInfor;
    }

    public String getTitleLine() {
        return titleLine;
    }

    /**
     * get pedigress file
     * 
     * @return pedfile
     */
    public File getFile() {
        return pedfile;
    }

    public void GMDRparseLinkage() throws GMDRPedFileException,
            PedFileException {
        try {
            parseLinkage(pedgreeStrings);
        } catch (PedFileException e) {
            System.out.println("PedFileException!!");
        }
        String firstLine = (String) pedgreeStrings.get(0);
        StringTokenizer tokenizer = new StringTokenizer(firstLine, "\n\t\" \"");
        if ((tokenizer.countTokens() - 6) / 2 != markerInfor.size()) {
            throw new GMDRPedFileException("The header line is incorrect!");
        }
    }

    public static void main(String[] args[]) {
        File gf = new File("Title.ped");
        GMDRPedFile gped = new GMDRPedFile();
    }

    public ArrayList getGMDRTable() {
        ArrayList tableData = new ArrayList();
        int numResults = getResults().size();
        GMDRmarkerRatings = new int[numResults];
        GMDRdups = new int[numResults];
        ArrayList header = new ArrayList();
        header.add(new String("Position"));
        header.add(new String("Marker"));
        header.add(new String("ObsHet"));
        header.add(new String("PreHet"));
        header.add(new String("HWpvalue"));
        header.add(new String("GenoPercent"));
        header.add(new String("FamTrioNum"));
        header.add(new String("MendErrNum"));
        header.add(new String("MAF"));
        header.add(new String("MinAllele"));
        tableData.add(header);
        for (int i = 0; i < numResults; i++) {
            ArrayList tempVect = new ArrayList();
            MarkerResult currentResult = (MarkerResult) (getResults().get(i));
            tempVect.add((new Integer(i + 1)).toString());

            tempVect.add((String) markerInfor.get(i));
            tempVect.add((new Double(currentResult.getObsHet())).toString());
            tempVect.add((new Double(currentResult.getPredHet())).toString());
            tempVect.add((new Double(currentResult.getHWpvalue())).toString());
            tempVect.add((new Double(currentResult.getGenoPercent())).toString());
            tempVect.add((new Integer(currentResult.getFamTrioNum())).toString());
            tempVect.add((new Integer(currentResult.getMendErrNum())).toString());
            tempVect.add((new Double(currentResult.getMAF())).toString());
            tempVect.add(currentResult.getMinorAllele());

            // these values are never displayed, just kept for bookkeeping
            tableData.add(tempVect.clone());
        }
        return tableData;
    }
}