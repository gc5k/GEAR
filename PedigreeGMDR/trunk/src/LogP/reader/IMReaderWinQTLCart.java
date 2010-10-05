package LogP.reader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;

import im.IMToolKit;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */ 
public class IMReaderWinQTLCart extends IMReaderAbstract {

    private String line;
    private ArrayList<String> lines;

    public IMReaderWinQTLCart(String file) {
        super(file);
        lines = new ArrayList();
    }

    public void readData() throws IOException {
        buffer = new BufferedReader(new FileReader(new File(file)));
        sweepComments();
        int index = 0;
        index = readFileID(index);
        index = readChromosomeParameter(index);
        index = readMap(index);
        index = readCrossParameter(index);
        index = readMarker(index);
        index = readPhenotype(index);
    }

    //step 2
    private int readFileID(int idx) throws IOException {
        Pattern pattern = Pattern.compile("(\\s*#fileid\\s*)(.+)\\s*", Pattern.CASE_INSENSITIVE);
        while (idx < lines.size()) {
            line = lines.get(idx++);
            Matcher matcher = pattern.matcher(line);
            if (matcher.matches()) {
                fileID = matcher.group(2);
                break;
            }
        }
        return idx;
    }

    //step 3
    private int readChromosomeParameter(int idx) throws IOException {
        Pattern[] p = {
            Pattern.compile("\\s*-type\\s+(\\w+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-function\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-Units\\s+(\\w+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-chromosome\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-maximum\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-named\\s+(\\w+)\\s*", Pattern.CASE_INSENSITIVE)
        };
        Pattern pstop = Pattern.compile("\\s*-start.*", Pattern.CASE_INSENSITIVE);
        Matcher matcher;
        String[] ChrP = new String[6];
        while (idx < lines.size()) {
            line = lines.get(idx++);
            for (int i = 0; i < p.length; i++) {
                matcher = p[i].matcher(line);
                if (matcher.matches()) {
                    ChrP[i] = matcher.group(1);
                    break;
                }
            }
            matcher = pstop.matcher(line);
            if (matcher.matches()) {
                break;
            }
        }
        chrp = new IMReaderAbstract.ChromosomeParameter();
        chrp.setParameter(ChrP);
        return idx;
    }

    //stpe 4
    private int readMap(int idx) throws IOException {
        Pattern p_named = Pattern.compile("y|yes", Pattern.CASE_INSENSITIVE);
        Matcher m_named = p_named.matcher(chrp.getNamed());
        Pattern p_chr = Pattern.compile("-Chromosome.*", Pattern.CASE_INSENSITIVE);
        Matcher m_chr;
        if (m_named.matches()) {
            markerName = new String[chrp.getChromosome()][];
        }
        distance = new double[chrp.getChromosome()][];

        Pattern p_stop = Pattern.compile("\\s*-stop.*");
        Matcher m_stop;
        ArrayList mk = null;
        int chrIdx = -1;
        while (idx < lines.size()) {
            line = lines.get(idx++);
            m_stop = p_stop.matcher(line);
            m_chr = p_chr.matcher(line);
            if (m_chr.matches() || m_stop.matches()) {//start a chromosome
                chrIdx++;
                if (chrIdx > 0) {
                    if (m_named.matches()) {//read marker name and distance
                        markerName[chrIdx - 1] = new String[mk.size()];
                        distance[chrIdx - 1] = new double[mk.size()];
                        for (int i = 0; i < mk.size(); i++) {
                            String[] name_dist = ((String) (mk.get(i))).split("\\s+");
                            markerName[chrIdx - 1][i] = name_dist[0];
                            if (chrp.getType().startsWith("d") || chrp.getType().startsWith("D")) {
                                distance[chrIdx - 1][i] = Double.parseDouble(name_dist[1]);
                            } else {
                                if(i == 0) {
                                    distance[chrIdx - 1][i] = Double.parseDouble(name_dist[1]);
                                } else {
                                    distance[chrIdx - 1][i] = distance[chrIdx - 1][i-1] + Double.parseDouble(name_dist[1]);
                                }
                            }
                        }
                    } else {//read distance only
                        distance[chrIdx - 1] = new double[mk.size()];
                        for (int i = 0; i < mk.size(); i++) {
                            if (chrp.getType().startsWith("p") || chrp.getType().startsWith("P")) {
                                distance[chrIdx - 1][i] = Double.parseDouble((String) mk.get(0));
                            } else {
                                if(i == 0) {
                                    distance[chrIdx - 1][i] = Double.parseDouble((String) mk.get(0));
                                } else {
                                    distance[chrIdx - 1][i] = distance[chrIdx - 1][i-1] + Double.parseDouble((String) mk.get(0));
                                }
                            }
                        }
                    }
                }
                if (m_stop.matches()) {
                    break;
                }
                mk = new ArrayList();
                continue;
            }
            mk.add(line);
        }
        double unit = (chrp.getUnites().compareTo("m") == 0 || chrp.getUnites().compareTo("M") == 0) ? 1 : 100;
        recombination = new double[distance.length][];
        for( int i = 0; i < distance.length; i++) {
            recombination[i] = new double[distance[i].length];
            for( int j = 0; j < distance[i].length; j++) {
                distance[i][j] /=unit;
            }
            recombination[i][0] = 0.5;
            for( int j = 1; j < distance[i].length; j++) {
                recombination[i][j] = IMToolKit.Felsenstein(distance[i][j] - distance[i][j-1], chrp.getFunction());
            }
        }
        return idx;
    }

    //step 5
    private int readCrossParameter(int idx) throws IOException {
        Pattern[] p = {
            Pattern.compile("\\s*-samplesize\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-cross\\s+([\\d|\\w]+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-traits\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-otraits\\s+(\\d+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-missingtrait\\s+(\\S+)\\s*", Pattern.CASE_INSENSITIVE),
            Pattern.compile("\\s*-case\\s+(\\w+)\\s*", Pattern.CASE_INSENSITIVE)
        };
        Matcher m_p;

        Pattern[] PT = {
            Pattern.compile("AA\\s+2\\s+(.+)"),
            Pattern.compile("Aa\\s+1\\s+(.+)"),
            Pattern.compile("aa\\s+0\\s+(.+)"),
            Pattern.compile("A-\\s+12\\s+(.+)"),
            Pattern.compile("a-\\s+10\\s+(.+)"),
            Pattern.compile("--\\s+-1\\s+(.+)")
        };
        Matcher m_pt;

        Pattern pTranslationTable = Pattern.compile("-TranslationTable", Pattern.CASE_INSENSITIVE);
        Matcher m_TT;
        Pattern pstop = Pattern.compile("\\s*-start\\s+markers*", Pattern.CASE_INSENSITIVE);
        Matcher m_stop;

        String[] CrP = new String[6];
        String[] TranslationTable = new String[6];
        int idxTT = 0;
        boolean isTT = false;
        while (idx < lines.size()) {
            line = lines.get(idx++);
            m_stop = pstop.matcher(line);
            if (m_stop.matches()) {//finish read CrossParameters
                break;
            }
            m_TT = pTranslationTable.matcher(line);
            if (m_TT.matches()) {//turn to translation table
                isTT = true;
                continue;
            }
            if (!isTT) {
                for (int i = 0; i < p.length; i++) {
                    m_p = p[i].matcher(line);
                    if (m_p.matches()) {
                        CrP[i] = m_p.group(1);
                        break;
                    }
                }
            } else {
                for (int i = 0; i < PT.length; i++) {
                    m_pt = PT[i].matcher(line);
                    if (m_pt.matches()) {
                        TranslationTable[i] = m_pt.group(0);
                        break;
                    }
                }
            }
        }
        crp = new IMReaderAbstract.CrossParameter();
        crp.setParameter(CrP, TranslationTable);
        return idx;
    }

    //step 6;
    private int readMarker(int idx) throws IOException {
        Pattern p_stop1 = Pattern.compile("-stop\\s+markers\\s*", Pattern.CASE_INSENSITIVE);
        Matcher m_stop1;

        Pattern p_stop2 = Pattern.compile("-stop\\s+individuals\\s+markers\\s*", Pattern.CASE_INSENSITIVE);
        Matcher m_stop2;

        Pattern p_data = Pattern.compile("\\d+\\s+(.*)");
        Matcher m_data;
        ArrayList Marker = new ArrayList();
        while (idx < lines.size()) {
            line = lines.get(idx++);
            m_stop1 = p_stop1.matcher(line);
            if (m_stop1.matches()) {
                isIndividualMarkers = false;
                break;
            }
            m_stop2 = p_stop2.matcher(line);
            if (m_stop2.matches()) {
                isIndividualMarkers = true;
                break;
            }
            m_data = p_data.matcher(line);
            if (m_data.matches()) {
                Marker.add(m_data.group(1));
            }
        }
        setMarker(Marker);
        return idx;
    }

    //step 7
    private int readPhenotype(int idx) throws IOException {
        Pattern p_start = Pattern.compile("\\s*-start\\s+traits\\s*", Pattern.CASE_INSENSITIVE);
        Matcher m_start;
        Pattern p_stop = Pattern.compile("\\s*-stop\\s+traits\\s*", Pattern.CASE_INSENSITIVE);
        Matcher m_stop;
        ArrayList phenotype = new ArrayList();
        String titleline = null;
        boolean flag = false;
        int pheIdx = 0;
        while (idx < lines.size()) {
            line = lines.get(idx++);
            m_stop = p_stop.matcher(line);
            if (m_stop.matches()) {
                break;
            }
            if (!flag) {
                m_start = p_start.matcher(line);
                if (m_start.matches()) {
                    flag = true;
                }
            } else {
                if (pheIdx > 0) {
                    phenotype.add(line);
                } else {
                    titleline = line;
                }
                pheIdx++;
            }
        }
        setPhenotypeInformation(titleline, phenotype);
        return idx;
    }

    //step 1
    private void sweepComments() throws IOException {
        boolean flag = true;
        Pattern pslc = Pattern.compile("(.*)(//.*)"); //compile single line comment
        while ((line = buffer.readLine()) != null) {
            if (Pattern.matches("\\s*", line)) {//empty line
                continue;
            }
            if (Pattern.matches("^/\\*.*", line)) {//multiple line comments
                flag = false;
            } else if (Pattern.matches("^\\s*\\*/.*", line)) {
                flag = true;
                continue;
            }
            if (flag) {
                Matcher matcher = pslc.matcher(line);
                String s1 = line;
                if (matcher.matches()) {//single line comment
                    s1 = matcher.group(1);
                }
                Pattern prs = Pattern.compile("(\\s)+$");
                Matcher mrs = prs.matcher(s1);
                String s2 = mrs.replaceFirst("");//remove the meaningless whitespaces in the end of matched region.
                Pattern pfs = Pattern.compile("^(\\s)+");
                Matcher mfs = pfs.matcher(s2);
                String s3 = mfs.replaceFirst("");
                lines.add(s3);
            }
        }
    }
}
