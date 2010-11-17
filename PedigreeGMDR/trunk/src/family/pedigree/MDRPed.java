package family.pedigree;

/*
 * $Id: PedFile.java,v 3.21 2006/05/12 18:01:28 jmaller Exp $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Vector;

import org._3pq.jgrapht.graph.SimpleGraph;

import publicAccess.PublicData;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.imputation.*;
import edu.mit.wi.haploview.Chromosome;
import edu.mit.wi.haploview.SNP;
import edu.mit.wi.pedfile.Family;
import edu.mit.wi.pedfile.Individual;
import edu.mit.wi.pedfile.MarkerResult;
import edu.mit.wi.pedfile.PedFileException;
import edu.mit.wi.pedparser.PedParser;
import edu.mit.wi.pedparser.PedigreeException;

/**
 * Handles input and storage of Pedigree files this class is not thread safe (untested). modified from original Pedfile
 * and checkdata classes by Hui Gong
 * 
 * @author Julian Maller
 * @author Guo-Bo Chen
 */
public class MDRPed {

	private Hashtable famNonInformative;
    private Hashtable famInformative;
    private Hashtable families;
    private Hashtable familystructure;
    private ArrayList axedPeople = new ArrayList();

    // stores the individuals found by parse() in allIndividuals. this is useful for outputting Pedigree information to
    // a file of another type.
    private ArrayList allIndividuals;

    // stores the individuals chosen by pedparser
    private ArrayList unrelatedIndividuals;
    private ArrayList results = null;
    private String[][] hminfo;
    // bogusParents is true if someone in the file referenced a parent not in the file
    private boolean bogusParents = false;
    private boolean haploidHets = false;
    private static Hashtable hapMapTranslate;
    private int[] markerRatings;
    private int[] dups;
    private HashSet whitelist;
    private ArrayList markerInfor;
    private ArrayList pedigrees;
    private ArrayList pedFileStrings;
    private ArrayList menderrortrace;
    private String titleLine;
    private File pedfile;
    private double MissingThreshold = 1;
    static int MissingAllele = 0;

    // public static String MissingGenotype="00";
    public MDRPed() {

        // hardcoded hapmap info
        this.famInformative = new Hashtable();
        this.famNonInformative = new Hashtable();
        this.families = new Hashtable();
        this.familystructure = new Hashtable();
        this.menderrortrace = new ArrayList();
        hapMapTranslate = new Hashtable(90, 1);
        hapMapTranslate.put("NA10846", "1334 NA10846 NA12144 NA12145 1 0");
        hapMapTranslate.put("NA12144", "1334 NA12144 0 0 1 0");
        hapMapTranslate.put("NA12145", "1334 NA12145 0 0 2 0");
        hapMapTranslate.put("NA10847", "1334a NA10847 NA12146 NA12239 2 0");
        hapMapTranslate.put("NA12146", "1334a NA12146 0 0 1 0");
        hapMapTranslate.put("NA12239", "1334a NA12239 0 0 2 0");
        hapMapTranslate.put("NA07029", "1340 NA07029 NA06994 NA07000 1 0");
        hapMapTranslate.put("NA06994", "1340 NA06994 0 0 1 0");
        hapMapTranslate.put("NA07000", "1340 NA07000 0 0 2 0");
        hapMapTranslate.put("NA07019", "1340a NA07019 NA07022 NA07056 2 0");
        hapMapTranslate.put("NA07022", "1340a NA07022 0 0 1 0");
        hapMapTranslate.put("NA07056", "1340a NA07056 0 0 2 0");
        hapMapTranslate.put("NA07048", "1341 NA07048 NA07034 NA07055 1 0");
        hapMapTranslate.put("NA07034", "1341 NA07034 0 0 1 0");
        hapMapTranslate.put("NA07055", "1341 NA07055 0 0 2 0");
        hapMapTranslate.put("NA06991", "1341a NA06991 NA06993 NA06985 2 0");
        hapMapTranslate.put("NA06993", "1341a NA06993 0 0 1 0");
        // hapMapTranslate.put("NA06993.dup", "dup NA06993.dup 0 0 1 0");
        hapMapTranslate.put("NA06985", "1341a NA06985 0 0 2 0");
        hapMapTranslate.put("NA10851", "1344 NA10851 NA12056 NA12057 1 0");
        hapMapTranslate.put("NA12056", "1344 NA12056 0 0 1 0");
        hapMapTranslate.put("NA12057", "1344 NA12057 0 0 2 0");
        hapMapTranslate.put("NA07348", "1345 NA07348 NA07357 NA07345 2 0");
        hapMapTranslate.put("NA07357", "1345 NA07357 0 0 1 0");
        hapMapTranslate.put("NA07345", "1345 NA07345 0 0 2 0");
        hapMapTranslate.put("NA10857", "1346 NA10857 NA12043 NA12044 1 0");
        hapMapTranslate.put("NA12043", "1346 NA12043 0 0 1 0");
        hapMapTranslate.put("NA12044", "1346 NA12044 0 0 2 0");
        hapMapTranslate.put("NA10859", "1347 NA10859 NA11881 NA11882 2 0");
        hapMapTranslate.put("NA11881", "1347 NA11881 0 0 1 0");
        hapMapTranslate.put("NA11882", "1347 NA11882 0 0 2 0");
        hapMapTranslate.put("NA10854", "1349 NA10854 NA11839 NA11840 2 0");
        hapMapTranslate.put("NA11839", "1349 NA11839 0 0 1 0");
        hapMapTranslate.put("NA11840", "1349 NA11840 0 0 2 0");
        hapMapTranslate.put("NA10856", "1350 NA10856 NA11829 NA11830 1 0");
        hapMapTranslate.put("NA11829", "1350 NA11829 0 0 1 0");
        hapMapTranslate.put("NA11830", "1350 NA11830 0 0 2 0");
        hapMapTranslate.put("NA10855", "1350a NA10855 NA11831 NA11832 2 0");
        hapMapTranslate.put("NA11831", "1350a NA11831 0 0 1 0");
        hapMapTranslate.put("NA11832", "1350a NA11832 0 0 2 0");
        hapMapTranslate.put("NA12707", "1358 NA12707 NA12716 NA12717 1 0");
        hapMapTranslate.put("NA12716", "1358 NA12716 0 0 1 0");
        hapMapTranslate.put("NA12717", "1358 NA12717 0 0 2 0");
        hapMapTranslate.put("NA10860", "1362 NA10860 NA11992 NA11993 1 0");
        hapMapTranslate.put("NA11992", "1362 NA11992 0 0 1 0");
        hapMapTranslate.put("NA11993", "1362 NA11993 0 0 2 0");
        // hapMapTranslate.put("NA11993.dup", "dup NA11993.dup 0 0 2 0");
        hapMapTranslate.put("NA10861", "1362a NA10861 NA11994 NA11995 2 0");
        hapMapTranslate.put("NA11994", "1362a NA11994 0 0 1 0");
        hapMapTranslate.put("NA11995", "1362a NA11995 0 0 2 0");
        hapMapTranslate.put("NA10863", "1375 NA10863 NA12264 NA12234 2 0");
        hapMapTranslate.put("NA12264", "1375 NA12264 0 0 1 0");
        hapMapTranslate.put("NA12234", "1375 NA12234 0 0 2 0");
        hapMapTranslate.put("NA10830", "1408 NA10830 NA12154 NA12236 1 0");
        hapMapTranslate.put("NA12154", "1408 NA12154 0 0 1 0");
        hapMapTranslate.put("NA12236", "1408 NA12236 0 0 2 0");
        hapMapTranslate.put("NA10831", "1408a NA10831 NA12155 NA12156 2 0");
        hapMapTranslate.put("NA12155", "1408a NA12155 0 0 1 0");
        hapMapTranslate.put("NA12156", "1408a NA12156 0 0 2 0");
        // hapMapTranslate.put("NA12156.dup", "dup NA12156.dup 0 0 2 0");
        hapMapTranslate.put("NA10835", "1416 NA10835 NA12248 NA12249 1 0");
        hapMapTranslate.put("NA12248", "1416 NA12248 0 0 1 0");
        // hapMapTranslate.put("NA12248.dup", "dup NA1248.dup 0 0 1 0");
        hapMapTranslate.put("NA12249", "1416 NA12249 0 0 2 0");
        hapMapTranslate.put("NA10838", "1420 NA10838 NA12003 NA12004 1 0");
        hapMapTranslate.put("NA12003", "1420 NA12003 0 0 1 0");
        // hapMapTranslate.put("NA12003.dup", "dup NA12003.dup 0 0 1 0");
        hapMapTranslate.put("NA12004", "1420 NA12004 0 0 2 0");
        hapMapTranslate.put("NA10839", "1420a NA10839 NA12005 NA12006 2 0");
        hapMapTranslate.put("NA12005", "1420a NA12005 0 0 1 0");
        hapMapTranslate.put("NA12006", "1420a NA12006 0 0 2 0");
        hapMapTranslate.put("NA12740", "1444 NA12740 NA12750 NA12751 2 0");
        hapMapTranslate.put("NA12750", "1444 NA12750 0 0 1 0");
        hapMapTranslate.put("NA12751", "1444 NA12751 0 0 2 0");
        hapMapTranslate.put("NA12752", "1447 NA12752 NA12760 NA12761 1 0");
        hapMapTranslate.put("NA12760", "1447 NA12760 0 0 1 0");
        hapMapTranslate.put("NA12761", "1447 NA12761 0 0 2 0");
        hapMapTranslate.put("NA12753", "1447a NA12753 NA12762 NA12763 2 0");
        hapMapTranslate.put("NA12762", "1447a NA12762 0 0 1 0");
        hapMapTranslate.put("NA12763", "1447a NA12763 0 0 2 0");
        hapMapTranslate.put("NA12801", "1454 NA12801 NA12812 NA12813 1 0");
        hapMapTranslate.put("NA12812", "1454 NA12812 0 0 1 0");
        hapMapTranslate.put("NA12813", "1454 NA12813 0 0 2 0");
        hapMapTranslate.put("NA12802", "1454a NA12802 NA12814 NA12815 2 0");
        hapMapTranslate.put("NA12814", "1454a NA12814 0 0 1 0");
        hapMapTranslate.put("NA12815", "1454a NA12815 0 0 2 0");
        hapMapTranslate.put("NA12864", "1459 NA12864 NA12872 NA12873 1 0");
        hapMapTranslate.put("NA12872", "1459 NA12872 0 0 1 0");
        hapMapTranslate.put("NA12873", "1459 NA12873 0 0 2 0");
        hapMapTranslate.put("NA12865", "1459a NA12865 NA12874 NA12875 2 0");
        hapMapTranslate.put("NA12874", "1459a NA12874 0 0 1 0");
        hapMapTranslate.put("NA12875", "1459a NA12875 0 0 2 0");
        hapMapTranslate.put("NA12878", "1463 NA12878 NA12891 NA12892 2 0");
        hapMapTranslate.put("NA12891", "1463 NA12891 0 0 1 0");
        hapMapTranslate.put("NA12892", "1463 NA12892 0 0 2 0");
        hapMapTranslate.put("NA18526", "chi1 NA18526 0 0 2 0");
        hapMapTranslate.put("NA18524", "chi2 NA18524 0 0 1 0");
        hapMapTranslate.put("NA18529", "chi3 NA18529 0 0 2 0");
        hapMapTranslate.put("NA18558", "chi4 NA18558 0 0 1 0");
        hapMapTranslate.put("NA18532", "chi5 NA18532 0 0 2 0");
        hapMapTranslate.put("NA18561", "chi6 NA18561 0 0 1 0");
        hapMapTranslate.put("NA18942", "jap1 NA18942 0 0 2 0");
        hapMapTranslate.put("NA18940", "jap2 NA18940 0 0 1 0");
        hapMapTranslate.put("NA18951", "jap3 NA18951 0 0 2 0");
        hapMapTranslate.put("NA18943", "jap4 NA18943 0 0 1 0");
        hapMapTranslate.put("NA18947", "jap5 NA18947 0 0 2 0");
        hapMapTranslate.put("NA18944", "jap6 NA18944 0 0 1 0");
        hapMapTranslate.put("NA18562", "chi7 NA18562 0 0 1 0");
        hapMapTranslate.put("NA18537", "chi8 NA18537 0 0 2 0");
        hapMapTranslate.put("NA18603", "chi9 NA18603 0 0 1 0");
        hapMapTranslate.put("NA18540", "chi10 NA18540 0 0 2 0");
        hapMapTranslate.put("NA18605", "chi11 NA18605 0 0 1 0");
        hapMapTranslate.put("NA18542", "chi12 NA18542 0 0 2 0");
        hapMapTranslate.put("NA18945", "jap7 NA18945 0 0 1 0");
        hapMapTranslate.put("NA18949", "jap8 NA18949 0 0 2 0");
        hapMapTranslate.put("NA18948", "jap9 NA18948 0 0 1 0");
        hapMapTranslate.put("NA18952", "jap10 NA18952 0 0 1 0");
        hapMapTranslate.put("NA18956", "jap11 NA18956 0 0 2 0");
        hapMapTranslate.put("NA18545", "chi13 NA18545 0 0 2 0");
        hapMapTranslate.put("NA18572", "chi46 NA18572 0 0 1 0");
        hapMapTranslate.put("NA18547", "chi15 NA18547 0 0 2 0");
        hapMapTranslate.put("NA18609", "chi16 NA18609 0 0 1 0");
        hapMapTranslate.put("NA18550", "chi17 NA18550 0 0 2 0");
        hapMapTranslate.put("NA18608", "chi18 NA18608 0 0 1 0");
        hapMapTranslate.put("NA18964", "jap12 NA18964 0 0 2 0");
        hapMapTranslate.put("NA18953", "jap13 NA18953 0 0 1 0");
        hapMapTranslate.put("NA18968", "jap14 NA18968 0 0 2 0");
        hapMapTranslate.put("NA18959", "jap15 NA18959 0 0 1 0");
        hapMapTranslate.put("NA18969", "jap16 NA18969 0 0 2 0");
        hapMapTranslate.put("NA18960", "jap17 NA18960 0 0 1 0");
        hapMapTranslate.put("NA18552", "chi19 NA18552 0 0 2 0");
        hapMapTranslate.put("NA18611", "chi20 NA18611 0 0 1 0");
        hapMapTranslate.put("NA18555", "chi21 NA18555 0 0 2 0");
        hapMapTranslate.put("NA18564", "chi22 NA18564 0 0 2 0");
        hapMapTranslate.put("NA18961", "jap18 NA18961 0 0 1 0");
        hapMapTranslate.put("NA18972", "jap19 NA18972 0 0 2 0");
        hapMapTranslate.put("NA18965", "jap20 NA18965 0 0 1 0");
        hapMapTranslate.put("NA18973", "jap21 NA18973 0 0 2 0");
        hapMapTranslate.put("NA18966", "jap22 NA18966 0 0 1 0");
        hapMapTranslate.put("NA18975", "jap23 NA18975 0 0 2 0");
        hapMapTranslate.put("NA18566", "chi23 NA18566 0 0 2 0");
        hapMapTranslate.put("NA18563", "chi24 NA18563 0 0 1 0");
        hapMapTranslate.put("NA18570", "chi25 NA18570 0 0 2 0");
        hapMapTranslate.put("NA18612", "chi26 NA18612 0 0 1 0");
        hapMapTranslate.put("NA18571", "chi27 NA18571 0 0 2 0");
        hapMapTranslate.put("NA18620", "chi28 NA18620 0 0 1 0");
        hapMapTranslate.put("NA18976", "jap24 NA18976 0 0 2 0");
        hapMapTranslate.put("NA18967", "jap25 NA18967 0 0 1 0");
        hapMapTranslate.put("NA18978", "jap26 NA18978 0 0 2 0");
        hapMapTranslate.put("NA18970", "jap27 NA18970 0 0 1 0");
        hapMapTranslate.put("NA18980", "jap28 NA18980 0 0 2 0");
        hapMapTranslate.put("NA18995", "jap29 NA18995 0 0 1 0");
        hapMapTranslate.put("NA18621", "chi29 NA18621 0 0 1 0");
        hapMapTranslate.put("NA18594", "chi30 NA18594 0 0 2 0");
        // hapMapTranslate.put("NA18594.dup", "dup 0 0 0 0 0");
        // hapMapTranslate.put("NA18603.dup", "dup 0 0 0 0 0");
        // hapMapTranslate.put("NA18609.dup", "dup 0 0 0 0 0");
        // hapMapTranslate.put("NA18951.dup", "dup 0 0 0 0 0");
        // hapMapTranslate.put("NA18995.dup", "dup 0 0 0 0 0");
        hapMapTranslate.put("NA18622", "chi31 NA18622 0 0 1 0");
        hapMapTranslate.put("NA18573", "chi32 NA18573 0 0 2 0");
        hapMapTranslate.put("NA18623", "chi33 NA18623 0 0 1 0");
        hapMapTranslate.put("NA18576", "chi34 NA18576 0 0 2 0");
        hapMapTranslate.put("NA18971", "jap30 NA18971 0 0 1 0");
        hapMapTranslate.put("NA18981", "jap31 NA18981 0 0 2 0");
        hapMapTranslate.put("NA18974", "jap32 NA18974 0 0 1 0");
        hapMapTranslate.put("NA18987", "jap33 NA18987 0 0 2 0");
        hapMapTranslate.put("NA18990", "jap34 NA18990 0 0 1 0");
        hapMapTranslate.put("NA18991", "jap35 NA18991 0 0 2 0");
        hapMapTranslate.put("NA18577", "chi35 NA18577 0 0 2 0");
        hapMapTranslate.put("NA18624", "chi36 NA18624 0 0 1 0");
        hapMapTranslate.put("NA18579", "chi37 NA18579 0 0 2 0");
        hapMapTranslate.put("NA18632", "chi38 NA18632 0 0 1 0");
        hapMapTranslate.put("NA18582", "chi39 NA18582 0 0 2 0");
        hapMapTranslate.put("NA18633", "chi40 NA18633 0 0 1 0");
        hapMapTranslate.put("NA18994", "jap36 NA18994 0 0 1 0");
        hapMapTranslate.put("NA18992", "jap37 NA18992 0 0 2 0");
        hapMapTranslate.put("NA18997", "jap38 NA18997 0 0 2 0");
        hapMapTranslate.put("NA18996", "jap39 NA18996 0 0 1 0");
        hapMapTranslate.put("NA18635", "chi41 NA18635 0 0 1 0");
        hapMapTranslate.put("NA18592", "chi42 NA18592 0 0 2 0");
        hapMapTranslate.put("NA18636", "chi43 NA18636 0 0 1 0");
        hapMapTranslate.put("NA18593", "chi44 NA18593 0 0 2 0");
        hapMapTranslate.put("NA18637", "chi45 NA18637 0 0 1 0");
        hapMapTranslate.put("NA19000", "jap40 NA19000 0 0 1 0");
        hapMapTranslate.put("NA18998", "jap41 NA18998 0 0 2 0");
        hapMapTranslate.put("NA19005", "jap42 NA19005 0 0 1 0");
        hapMapTranslate.put("NA18999", "jap43 NA18999 0 0 2 0");
        hapMapTranslate.put("NA19007", "jap44 NA19007 0 0 1 0");
        hapMapTranslate.put("NA19003", "jap45 NA19003 0 0 2 0");
        hapMapTranslate.put("NA19012", "jap46 NA19012 0 0 1 0");
        hapMapTranslate.put("NA18500", "Yoruba004 NA18500 NA18501 NA18502 1 0");
        hapMapTranslate.put("NA18501", "Yoruba004 NA18501 0 0 1 0");
        hapMapTranslate.put("NA18502", "Yoruba004 NA18502 0 0 2 0");
        hapMapTranslate.put("NA18503", "Yoruba005 NA18503 NA18504 NA18505 1 0");
        hapMapTranslate.put("NA18504", "Yoruba005 NA18504 0 0 1 0");
        hapMapTranslate.put("NA18505", "Yoruba005 NA18505 0 0 2 0");
        hapMapTranslate.put("NA18506", "Yoruba009 NA18506 NA18507 NA18508 1 0");
        hapMapTranslate.put("NA18507", "Yoruba009 NA18507 0 0 1 0");
        hapMapTranslate.put("NA18508", "Yoruba009 NA18508 0 0 2 0");
        hapMapTranslate.put("NA18860", "Yoruba012 NA18860 NA18859 NA18858 1 0");
        hapMapTranslate.put("NA18859", "Yoruba012 NA18859 0 0 1 0");
        hapMapTranslate.put("NA18858", "Yoruba012 NA18858 0 0 2 0");
        hapMapTranslate.put("NA18515", "Yoruba013 NA18515 NA18516 NA18517 1 0");
        hapMapTranslate.put("NA18516", "Yoruba013 NA18516 0 0 1 0");
        hapMapTranslate.put("NA18517", "Yoruba013 NA18517 0 0 2 0");
        hapMapTranslate.put("NA18521", "Yoruba016 NA18521 NA18522 NA18523 1 0");
        hapMapTranslate.put("NA18522", "Yoruba016 NA18522 0 0 1 0");
        hapMapTranslate.put("NA18523", "Yoruba016 NA18523 0 0 2 0");
        hapMapTranslate.put("NA18872", "Yoruba017 NA18872 NA18871 NA18870 1 0");
        hapMapTranslate.put("NA18871", "Yoruba017 NA18871 0 0 1 0");
        hapMapTranslate.put("NA18870", "Yoruba017 NA18870 0 0 2 0");
        hapMapTranslate.put("NA18854", "Yoruba018 NA18854 NA18853 NA18852 1 0");
        hapMapTranslate.put("NA18853", "Yoruba018 NA18853 0 0 1 0");
        hapMapTranslate.put("NA18852", "Yoruba018 NA18852 0 0 2 0");
        hapMapTranslate.put("NA18857", "Yoruba023 NA18857 NA18856 NA18855 1 0");
        hapMapTranslate.put("NA18856", "Yoruba023 NA18856 0 0 1 0");
        hapMapTranslate.put("NA18855", "Yoruba023 NA18855 0 0 2 0");
        hapMapTranslate.put("NA18863", "Yoruba024 NA18863 NA18862 NA18861 1 0");
        hapMapTranslate.put("NA18862", "Yoruba024 NA18862 0 0 1 0");
        hapMapTranslate.put("NA18861", "Yoruba024 NA18861 0 0 2 0");
        hapMapTranslate.put("NA18914", "Yoruba028 NA18914 NA18913 NA18912 1 0");
        hapMapTranslate.put("NA18913", "Yoruba028 NA18913 0 0 1 0");
        hapMapTranslate.put("NA18912", "Yoruba028 NA18912 0 0 2 0");
        hapMapTranslate.put("NA19094", "Yoruba040 NA19094 NA19092 NA19093 2 0");
        hapMapTranslate.put("NA19092", "Yoruba040 NA19092 0 0 1 0");
        hapMapTranslate.put("NA19093", "Yoruba040 NA19093 0 0 2 0");
        hapMapTranslate.put("NA19103", "Yoruba042 NA19103 NA19101 NA19102 1 0");
        hapMapTranslate.put("NA19101", "Yoruba042 NA19101 0 0 1 0");
        hapMapTranslate.put("NA19102", "Yoruba042 NA19102 0 0 2 0");
        hapMapTranslate.put("NA19139", "Yoruba043 NA19139 NA19138 NA19137 1 0");
        hapMapTranslate.put("NA19138", "Yoruba043 NA19138 0 0 1 0");
        hapMapTranslate.put("NA19137", "Yoruba043 NA19137 0 0 2 0");
        hapMapTranslate.put("NA19202", "Yoruba045 NA19202 NA19200 NA19201 2 0");
        hapMapTranslate.put("NA19200", "Yoruba045 NA19200 0 0 1 0");
        hapMapTranslate.put("NA19201", "Yoruba045 NA19201 0 0 2 0");
        hapMapTranslate.put("NA19173", "Yoruba047 NA19173 NA19171 NA19172 1 0");
        hapMapTranslate.put("NA19171", "Yoruba047 NA19171 0 0 1 0");
        hapMapTranslate.put("NA19172", "Yoruba047 NA19172 0 0 2 0");
        hapMapTranslate.put("NA19205", "Yoruba048 NA19205 NA19203 NA19204 1 0");
        hapMapTranslate.put("NA19203", "Yoruba048 NA19203 0 0 1 0");
        hapMapTranslate.put("NA19204", "Yoruba048 NA19204 0 0 2 0");
        hapMapTranslate.put("NA19211", "Yoruba050 NA19211 NA19210 NA19209 1 0");
        hapMapTranslate.put("NA19210", "Yoruba050 NA19210 0 0 1 0");
        hapMapTranslate.put("NA19209", "Yoruba050 NA19209 0 0 2 0");
        hapMapTranslate.put("NA19208", "Yoruba051 NA19208 NA19207 NA19206 1 0");
        hapMapTranslate.put("NA19207", "Yoruba051 NA19207 0 0 1 0");
        hapMapTranslate.put("NA19206", "Yoruba051 NA19206 0 0 2 0");
        hapMapTranslate.put("NA19161", "Yoruba056 NA19161 NA19160 NA19159 1 0");
        hapMapTranslate.put("NA19160", "Yoruba056 NA19160 0 0 1 0");
        hapMapTranslate.put("NA19159", "Yoruba056 NA19159 0 0 2 0");
        hapMapTranslate.put("NA19221", "Yoruba058 NA19221 NA19223 NA19222 2 0");
        hapMapTranslate.put("NA19223", "Yoruba058 NA19223 0 0 1 0");
        hapMapTranslate.put("NA19222", "Yoruba058 NA19222 0 0 2 0");
        hapMapTranslate.put("NA19120", "Yoruba060 NA19120 NA19119 NA19116 1 0");
        hapMapTranslate.put("NA19119", "Yoruba060 NA19119 0 0 1 0");
        hapMapTranslate.put("NA19116", "Yoruba060 NA19116 0 0 2 0");
        hapMapTranslate.put("NA19142", "Yoruba071 NA19142 NA19141 NA19140 1 0");
        hapMapTranslate.put("NA19141", "Yoruba071 NA19141 0 0 1 0");
        hapMapTranslate.put("NA19140", "Yoruba071 NA19140 0 0 2 0");
        hapMapTranslate.put("NA19154", "Yoruba072 NA19154 NA19153 NA19152 1 0");
        hapMapTranslate.put("NA19153", "Yoruba072 NA19153 0 0 1 0");
        hapMapTranslate.put("NA19152", "Yoruba072 NA19152 0 0 2 0");
        hapMapTranslate.put("NA19145", "Yoruba074 NA19145 NA19144 NA19143 1 0");
        hapMapTranslate.put("NA19144", "Yoruba074 NA19144 0 0 1 0");
        hapMapTranslate.put("NA19143", "Yoruba074 NA19143 0 0 2 0");
        hapMapTranslate.put("NA19129", "Yoruba077 NA19129 NA19128 NA19127 2 0");
        hapMapTranslate.put("NA19128", "Yoruba077 NA19128 0 0 1 0");
        hapMapTranslate.put("NA19127", "Yoruba077 NA19127 0 0 2 0");
        hapMapTranslate.put("NA19132", "Yoruba101 NA19132 NA19130 NA19131 2 0");
        hapMapTranslate.put("NA19130", "Yoruba101 NA19130 0 0 1 0");
        hapMapTranslate.put("NA19131", "Yoruba101 NA19131 0 0 2 0");
        hapMapTranslate.put("NA19100", "Yoruba105 NA19100 NA19098 NA19099 2 0");
        hapMapTranslate.put("NA19098", "Yoruba105 NA19098 0 0 1 0");
        hapMapTranslate.put("NA19099", "Yoruba105 NA19099 0 0 2 0");
        hapMapTranslate.put("NA19194", "Yoruba112 NA19194 NA19192 NA19193 1 0");
        hapMapTranslate.put("NA19192", "Yoruba112 NA19192 0 0 1 0");
        hapMapTranslate.put("NA19193", "Yoruba112 NA19193 0 0 2 0");
        hapMapTranslate.put("NA19240", "Yoruba117 NA19240 NA19239 NA19238 2 0");
        hapMapTranslate.put("NA19239", "Yoruba117 NA19239 0 0 1 0");
        hapMapTranslate.put("NA19238", "Yoruba117 NA19238 0 0 2 0");

        // data entered for the additional perlegen samples
        hapMapTranslate.put("NA10844", "NA10844 NA10844 0 0 2 0");
        hapMapTranslate.put("NA17134", "NA17134 NA17134 0 0 2 0");
        hapMapTranslate.put("NA17115", "NA17115 NA17115 0 0 1 0");
        hapMapTranslate.put("NA12560", "NA12560 NA12560 0 0 1 0");
        hapMapTranslate.put("NA17137", "NA17137 NA17137 0 0 2 0");
        hapMapTranslate.put("NA17747", "NA17747 NA17747 0 0 2 0");
        hapMapTranslate.put("NA12547", "NA12547 NA12547 0 0 1 0");
        hapMapTranslate.put("NA17138", "NA17138 NA17138 0 0 2 0");
        hapMapTranslate.put("NA17116", "NA17116 NA17116 0 0 2 0");
        hapMapTranslate.put("NA17753", "NA17753 NA17753 0 0 1 0");
        hapMapTranslate.put("NA17102", "NA17102 NA17102 0 0 1 0");
        hapMapTranslate.put("NA10858", "NA10858 NA10858 0 0 1 0");
        hapMapTranslate.put("NA17106", "NA17106 NA17106 0 0 1 0");
        hapMapTranslate.put("NA17114", "NA17114 NA17114 0 0 1 0");
        hapMapTranslate.put("NA12548", "NA12548 NA12548 0 0 2 0");
        hapMapTranslate.put("NA17109", "NA17109 NA17109 0 0 1 0");
        hapMapTranslate.put("NA17740", "NA17740 NA17740 0 0 2 0");
        hapMapTranslate.put("NA17139", "NA17139 NA17139 0 0 2 0");
        hapMapTranslate.put("NA07349", "NA07349 NA07349 0 0 1 0");
        hapMapTranslate.put("NA17761", "NA17761 NA17761 0 0 2 0");
        hapMapTranslate.put("NA10853", "NA10853 NA10853 0 0 2 0");
        hapMapTranslate.put("NA17737", "NA17737 NA17737 0 0 1 0");
        hapMapTranslate.put("NA17744", "NA17744 NA17744 0 0 2 0");
        hapMapTranslate.put("NA10845", "NA10845 NA10845 0 0 1 0");
        hapMapTranslate.put("NA17140", "NA17140 NA17140 0 0 2 0");
        hapMapTranslate.put("NA17105", "NA17105 NA17105 0 0 1 0");
        hapMapTranslate.put("NA17752", "NA17752 NA17752 0 0 2 0");
        hapMapTranslate.put("NA17746", "NA17746 NA17746 0 0 2 0");
        hapMapTranslate.put("NA17135", "NA17135 NA17135 0 0 2 0");
        hapMapTranslate.put("NA17742", "NA17742 NA17742 0 0 1 0");
        hapMapTranslate.put("NA17104", "NA17104 NA17104 0 0 1 0");
        hapMapTranslate.put("NA17741", "NA17741 NA17741 0 0 2 0");
        hapMapTranslate.put("NA17738", "NA17738 NA17738 0 0 2 0");
        hapMapTranslate.put("NA17201", "NA17201 NA17201 0 0 1 0");
        hapMapTranslate.put("NA17745", "NA17745 NA17745 0 0 2 0");
        hapMapTranslate.put("NA17749", "NA17749 NA17749 0 0 1 0");
        hapMapTranslate.put("NA17133", "NA17133 NA17133 0 0 2 0");
        hapMapTranslate.put("NA10842", "NA10842 NA10842 0 0 1 0");
        hapMapTranslate.put("NA17736", "NA17736 NA17736 0 0 1 0");
        hapMapTranslate.put("NA17107", "NA17107 NA17107 0 0 1 0");
        hapMapTranslate.put("NA10852", "NA10852 NA10852 0 0 2 0");
        hapMapTranslate.put("NA17756", "NA17756 NA17756 0 0 2 0");
        hapMapTranslate.put("NA17735", "NA17735 NA17735 0 0 2 0");
        hapMapTranslate.put("NA10848", "NA10848 NA10848 0 0 1 0");
        hapMapTranslate.put("NA10850", "NA10850 NA10850 0 0 2 0");
        hapMapTranslate.put("NA17110", "NA17110 NA17110 0 0 2 0");
        hapMapTranslate.put("NA17111", "NA17111 NA17111 0 0 1 0");
        hapMapTranslate.put("NA17136", "NA17136 NA17136 0 0 2 0");
        hapMapTranslate.put("NA17755", "NA17755 NA17755 0 0 1 0");
        hapMapTranslate.put("NA06990", "NA06990 NA06990 0 0 2 0");
        hapMapTranslate.put("NA17733", "NA17733 NA17733 0 0 2 0");
        hapMapTranslate.put("NA17103", "NA17103 NA17103 0 0 1 0");
        hapMapTranslate.put("NA17739", "NA17739 NA17739 0 0 2 0");
        hapMapTranslate.put("NA17108", "NA17108 NA17108 0 0 1 0");
        hapMapTranslate.put("NA17759", "NA17759 NA17759 0 0 1 0");
        hapMapTranslate.put("NA17112", "NA17112 NA17112 0 0 2 0");
        hapMapTranslate.put("NA17113", "NA17113 NA17113 0 0 2 0");
        hapMapTranslate.put("NA10843", "NA10843 NA10843 0 0 2 0");
        hapMapTranslate.put("NA17743", "NA17743 NA17743 0 0 1 0");
        hapMapTranslate.put("NA17757", "NA17757 NA17757 0 0 2 0");
        hapMapTranslate.put("NA17754", "NA17754 NA17754 0 0 2 0");
        hapMapTranslate.put("NA17734", "NA17734 NA17734 0 0 2 0");
    }

    public void Allele2Genotype() throws MDRPedFileException {
        Enumeration famList = this.families.keys();
        MendErrorTrace met;
        while (famList.hasMoreElements()) {
            String fid = (String) famList.nextElement();

            Family Fam = getFamily(fid);
            Enumeration sibList = Fam.getMemberList();

            FamilyStruct FamStr = getFamilyStruct(fid);

            while (sibList.hasMoreElements()) {
                String iid = (String) sibList.nextElement();
                try {
                    Individual ind = Fam.getMember(iid);
                    Person per = FamStr.getPerson(iid);
                    ArrayList VM = new ArrayList();
                    for (int i = 0; i < ind.getNumMarkers(); i++) {
                        if (!FamStr.getPersons().containsKey(per.getPersonID())) {
                            throw new MDRPedFileException("Person " + per.getPersonID() + " in family " + ind.getFamilyID() + " can't be matched.");
                        }
                        String marker;
                        if (ind.getZeroed(i)) {
                            marker = new String(PublicData.MissingGenotype);
                            met = new MendErrorTrace(fid, iid,
                                    (String) markerInfor.get(i), i);
                            menderrortrace.add(met);
                        } else {
                            byte m1 = ind.getAllele(i, 0);
                            byte m2 = ind.getAllele(i, 1);
                            if (m2 < m1) {
                                byte m3 = m1;
                                m1 = m2;
                                m2 = m3;
                            }
                            marker = new String((new Byte(m1).toString()) + (new Byte(m2).toString()));
                        }
                        VM.add(marker);
                    }
                    per.setGenotype(VM);
                } catch (MDRPedFileException e) {
                    System.err.println("no such person with ID " + iid + " in family " + fid + "in ped files.");
                } catch (PedFileException e) {
                    System.err.println("no such individual with ID " + iid + " in family " + fid + "in ped file.");
                }
            }
        }
        GenotypeSummary(true);
    }

    public void GenotypeImputation() throws MDRPedFileException {
        Enumeration fsList = this.familystructure.keys();
        ArrayList error_report = new ArrayList();
        while (fsList.hasMoreElements()) {
            FamilyStruct fs = (FamilyStruct) getFamilyStruct((String) fsList.nextElement());
            for (int i = 0; i < getNumMarkers(); i++) {
                Imputation(fs, i);
            }
        }
        GenotypeSummary(false);
        for (Iterator it = error_report.iterator(); it.hasNext();) {
            System.err.println((String) it.next());
        }
    }

    protected void Imputation(FamilyStruct fs, int genoIdx) throws MDRPedFileException {
        GenoSet gSet = fs.getObservedGenoSet(genoIdx);
        AbstractImputation ai;
        if ((gSet.getNumUntypedChildren() > 0) && (gSet.getNumUntypedChildren() < gSet.getNumChildren())) {
            if (gSet.getNumTypedParents() == 2) {
                ai = new GenotypedParents(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
            } else if (gSet.getNumParents() == 1) {
                String pg = (String) ((TreeMap) gSet.getfullchildrenGenoMap()).firstKey();
                if (!AbstractGenoDistribution.isHeterozygous(pg)) {
                    ai = new OneHomozygousParent(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
                } else {
                    ai = new OneHeterozygousParent(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
                }
            } else {
                ai = new UngenotypedParents(gSet.getchildrenGenoMap());
            }

            Enumeration perList = fs.getPersonList();
            while (perList.hasMoreElements()) {
                Person per;
                per = fs.getPerson((String) perList.nextElement());
                String genotype = per.getGenotype(genoIdx);
                if (fs.hasAncestor(per.getPersonID())) {
                    if (genotype.compareTo(PublicData.MissingGenotype) == 0) {
                        String ImputedGenotype = ai.RandomAssign();
                        per.setGenotype(genoIdx, ImputedGenotype);
                        System.out.println(fs.getFamilyStructName() + " " + per.getPersonID() + " " + genoIdx);
                    }
                }
            }
        }
    }

    protected void GenotypeSummary(boolean IsObservedGenoSet) throws MDRPedFileException {
        Enumeration fsList = this.familystructure.keys();
        ArrayList error_report = new ArrayList();
        while (fsList.hasMoreElements()) {
            FamilyStruct fs = (FamilyStruct) getFamilyStruct((String) fsList.nextElement());
            Boolean b = new Boolean(true);
            famInformative.put(fs.getFamilyStructName(), b);
            TreeMap Ps;
            TreeMap Ks;
            GenoSet gSet;
            for (int i = 0; i < getNumMarkers(); i++) {
                Ps = new TreeMap();
                Ks = new TreeMap();
                Enumeration perList = fs.getPersonList();
                while (perList.hasMoreElements()) {
                    Person per = fs.getPerson((String) perList.nextElement());
                    String genotype = per.getGenotype(i);
                    if (fs.hasAncestor(per.getPersonID())) {
                        if (Ks.containsKey(genotype)) {
                            Integer c = ((Integer) Ks.get(genotype));
                            int v = (c.intValue());
                            v++;
                            c = new Integer(v);
                            Ks.put(genotype, c);
                        } else {
                            Integer c = new Integer(1);
                            Ks.put(new String(genotype), c);
                        }
                    } else {
                        if (Ps.containsKey(genotype)) {
                            Integer c = ((Integer) Ps.get(genotype));
                            int v = (c.intValue());
                            v++;
                            c = new Integer(v);
                            Ps.put(genotype, c);
                        } else {
                            Integer c = new Integer(1);
                            Ps.put(new String(genotype), c);
                        }
                    }
                }
                gSet = new GenoSet(Ps, Ks, i);
                if (IsObservedGenoSet) {
                    fs.addObservedGenoSet(gSet);
                } else {
                    if (gSet.getchildrenGenoMap().size() == 0) {
                        error_report.add(new String("Family " + fs.getFamilyStructName() + " Locus " + (String) markerInfor.get(i)));
                        b = new Boolean(false);
                        famInformative.put(fs.getFamilyStructName(), b);
                    }
                    fs.addImputedGenoSet(gSet);
                }
            }
        }
    }

    /**
     * gets the allIndividuals ArrayList
     */
    public ArrayList getAllIndividuals() {
        return allIndividuals;
    }

    public ArrayList getUnrelatedIndividuals() {
        return unrelatedIndividuals;
    }

    /**
     * @return enumeration containing a list of familyID's in the families hashtable
     */
    public Enumeration getFamList() {
        return this.families.keys();
    }

    public Enumeration getFamStrList() {
        return this.familystructure.keys();
    }

    public String[] getFamListSorted() {
    	Enumeration famstrList = this.familystructure.keys();
		String[] FID = new String[getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			FID[ind++] = (String) famstrList.nextElement();
		}
		Arrays.sort(FID);
		return FID;
    }

    public String[] getInformativeFamListSorted() {
    	Enumeration famstrList = this.familystructure.keys();
		String[] FID = new String[getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			String F = (String) famstrList.nextElement();
			if ( ((Boolean) (famInformative.get(F))).booleanValue()) {
				FID[ind++] = F;
			}
		}
		Arrays.sort(FID);
		return FID;
    }

    public String[] getUninformativeFamListSorted() {
    	Enumeration famstrList = this.familystructure.keys();
		String[] FID = new String[getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			String F = (String) famstrList.nextElement();
			if ( !((Boolean) (famInformative.get(F))).booleanValue()) {
				FID[ind++] = F;
			}
		}
		Arrays.sort(FID);
		return FID;
    }

    public Hashtable getFamInformative() {
        return this.famInformative;
    }

    /**
     * @param familyID
     *            id of desired family
     * @return Family identified by familyID in families hashtable
     */
    public Family getFamily(String familyID) {
        return (Family) this.families.get(familyID);
    }

    public FamilyStruct getFamilyStruct(String familystrID) {
        return (FamilyStruct) this.familystructure.get(familystrID);
    }

    /**
     * @return the number of Family objects in the families hashtable
     */
    public int getNumFamilies() {
        return this.families.size();
    }

    /**
     * Setting the threshould of missinggenopercentage.
     * 
     * @param threshould
     *            (0~1);
     */
    public void setMissingThreshold(double threshold) {
        if (threshold >= 0 || threshold <= 1) {
            MissingThreshold = threshold;
        }
    }

    /**
     * get the MissingThreshould value;
     * 
     * @return
     */
    public double getMissingThreshold() {
        return MissingThreshold;
    }

    /**
     * this method iterates through each family in Hashtable families and adds up the number of individuals in total
     * across all families
     * 
     * @return the total number of individuals in all the family objects in the families hashtable
     */
    public int getNumIndividuals() {
        Enumeration famEnum = this.families.elements();
        int total = 0;
        while (famEnum.hasMoreElements()) {
            Family fam = (Family) famEnum.nextElement();
            total += fam.getNumMembers();
        }
        return total;
    }

    /**
     * finds the first individual in the first family and returns the number of markers for that individual
     * 
     * @return the number of markers
     */
    public int getNumMarkers() {
        Enumeration famList = this.families.elements();
        int numMarkers = 0;
        while (famList.hasMoreElements()) {
            Family fam = (Family) famList.nextElement();
            Enumeration indList = fam.getMemberList();
            Individual ind = null;
            while (indList.hasMoreElements()) {
                try {
                    ind = fam.getMember((String) indList.nextElement());
                } catch (PedFileException pfe) {
                }
                numMarkers = ind.getNumMarkers();
                if (numMarkers > 0) {
                    return numMarkers;
                }
            }
        }
        return 0;
    }

    public File GetPedFile() {
        return pedfile;
    }

    public void Initial(File infile) throws MDRPedFileException, IOException {
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
        pedigrees = (ArrayList) pedFileStrings.clone();
        titleLine = (String) pedigrees.get(0);
        pedigrees.remove(0);
        int numLines = pedFileStrings.size();

        if (numLines < 2) {
            throw new MDRPedFileException(
                    "Pedgree data format error: empty file");
        }
        StringTokenizer tokenizer = new StringTokenizer(titleLine, "\n\t\" \"");
        int numTokens = tokenizer.countTokens();

        if (numTokens < 7) {
            throw new MDRPedFileException(
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
        markerInfor = new ArrayList(temp);
    // System.out.println(markerInfor);
    }

    /**
     * Taking in a pedigree file in the form of a vector of strings and parses it.
     * The data parsed is stored in families in the member hashtable families.
     * Note that the "Linkage" here is the relationship between relatives in a
     * pedigree, but is not that term of genetics.
     */
    public void parseLinkage() throws MDRPedFileException {
        int colNum = markerInfor.size() * 2 + 6;

        int numMarkers = 0;
        boolean genoError = false;
        int numLines = pedigrees.size();
        if (numLines == 0) {
            throw new MDRPedFileException("Data format error: empty file");
        }
        Individual ind;
        Person per;
        PseudoPerson pseudoper;
        this.allIndividuals = new ArrayList();

        for (int k = 0; k < numLines; k++) {
            StringTokenizer tokenizer = new StringTokenizer((String) pedigrees.get(k), "\n\t\" \"");
            int numTokens = tokenizer.countTokens();

            if (numTokens % 2 == 0) {
                numMarkers = (numTokens - 6) / 2;
                if (numMarkers != markerInfor.size()) {
                    new MDRPedFileException(
                            "Mismatched Colunm in pedfile. line " + (k + 2));
                }
            } else {
                new MDRPedFileException("Mismatched Column in pedfile. line " + (k + 2));
            }

            if (colNum != numTokens) {
                // this line has a different number of columns
                // should send some sort of error message
                throw new MDRPedFileException(
                        "Column number mismatch in pedfile. line " + (k + 2));
            }

            ind = new Individual(numMarkers);
            per = new Person(numMarkers);
            pseudoper = new PseudoPerson(numMarkers);
            if (numTokens < 6) {
                throw new MDRPedFileException(
                        "Incorrect number of fields on line " + (k + 2));
            }

            if (tokenizer.hasMoreTokens()) {
                String FID = (String) tokenizer.nextToken().trim();
                String IID = (String) tokenizer.nextToken().trim();
                String DID = (String) tokenizer.nextToken().trim();
                String MID = (String) tokenizer.nextToken().trim();

                ind.setFamilyID(new String(FID));
                ind.setIndividualID(new String(IID));
                ind.setDadID(new String(DID));
                ind.setMomID(new String(MID));

                per.setFamilyID(new String(FID));
                per.setPersonID(new String(IID));
                per.setDadID(new String(DID));
                per.setMomID(new String(MID));

                pseudoper.setFamilyID(new String(FID));
                pseudoper.setPseudoPersonID(new String(IID));
                pseudoper.setDadID(new String(DID));
                pseudoper.setMomID(new String(MID));

                try {
                    int Gender = Integer.parseInt(tokenizer.nextToken().trim());
                    int Status = Integer.parseInt(tokenizer.nextToken().trim());
                    ind.setGender(Gender);
                    ind.setAffectedStatus(Status);

                    per.setGender(Gender);
                    per.setAffectedStatus(Status);

                    pseudoper.setGender(Gender);
                    pseudoper.setAffectedStatus(Status);
                } catch (NumberFormatException nfe) {
                    throw new MDRPedFileException(
                            "Pedfile error: invalid gender or affected status on line " + (k + 2));
                }

                byte genotype1;
                byte genotype2;
                while (tokenizer.hasMoreTokens()) {
                    try {
                        String alleleA = tokenizer.nextToken();
                        String alleleB = tokenizer.nextToken();
                        int[] checker1, checker2;
                        checker1 = checkGenotype(alleleA);
                        checker2 = checkGenotype(alleleB);
                        if (checker1[1] != checker2[1]) {
                            genoError = !genoError;
                        }

                        if (genoError) {
                            throw new MDRPedFileException(
                                    "File input error on line " + (k + 2) + ", marker " + (ind.getNumMarkers() + 2) + ".\nFor any marker, an individual's genotype must be only letters or only numbers.");
                        }
                        if (checker1[0] == MissingAllele || checker2[0] == MissingAllele) {
                            checker1[0] = checker2[0] = MissingAllele;
                        }
                        genotype1 = (byte) checker1[0];
                        genotype2 = (byte) checker2[0];
                        ind.addMarker(genotype1, genotype2);
                        per.addMarker(genotype1, genotype2);
                    } catch (NumberFormatException nfe) {
                        throw new MDRPedFileException(
                                "Pedigree file input error: invalid genotype on line " + (k + 2));
                    }
                }

                // check if the family exists already in the Hashtable
                Family fam = (Family) this.families.get(ind.getFamilyID());
                if (fam == null) {
                    // it doesnt exist, so create a new Family object
                    fam = new Family(ind.getFamilyID());
                }

                if (fam.getMembers().containsKey(ind.getIndividualID())) {
                    throw new MDRPedFileException("Individual " + ind.getIndividualID() + " in family " + ind.getFamilyID() + " appears more than once.");
                }

                fam.addMember(ind);
                this.families.put(ind.getFamilyID(), fam);
                this.allIndividuals.add(ind);

                // check if the family exists already in the Hashtable
                FamilyStruct famstr = (FamilyStruct) this.familystructure.get(per.getFamilyID());
                if (famstr == null) {
                    // it doesn't exist, so create a new FamilyStruct object
                    famstr = new FamilyStruct(per.getFamilyID());
                    famstr.setNumMarkers(numMarkers);
                }

                if (famstr.getPersons().containsKey(per.getPersonID()) || famstr.getPseudoPersons().containsKey(
                        pseudoper.getPseudoPersonID())) {
                    throw new MDRPedFileException("Person" + per.getPersonID() + " in family " + ind.getFamilyID() + " appears more than once.");
                }

                famstr.addPerson(per);
                famstr.addPseudoPerson(pseudoper);
                this.familystructure.put(per.getFamilyID(), famstr);
            }
        }

        // now we check if anyone has a reference to a parent who isn't in the file, and if so, we remove the reference
        for (int i = 0; i < allIndividuals.size(); i++) {
            Individual currentInd = (Individual) allIndividuals.get(i);
            Hashtable curFam = ((Family) (families.get(currentInd.getFamilyID()))).getMembers();
            if (!currentInd.getDadID().equals("0") && !(curFam.containsKey(currentInd.getDadID()))) {
                currentInd.setDadID("0");
                bogusParents = true;
            }
            if (!currentInd.getMomID().equals("0") && !(curFam.containsKey(currentInd.getMomID()))) {
                currentInd.setMomID("0");
                bogusParents = true;
            }
        }
    }

    public void parseHapMap(ArrayList lines, ArrayList hapsData)
            throws MDRPedFileException {
        int colNum = -1;
        int numLines = lines.size();
        if (numLines < 2) {
            throw new MDRPedFileException(
                    "Hapmap data format error: empty file");
        }
        if (hapsData != null) {
            String indName;
            for (int i = 0; i < hapsData.size(); i++) {
                StringTokenizer hd = new StringTokenizer((String) hapsData.get(i));
                if (hd.countTokens() < 6) {
                    throw new MDRPedFileException(
                            "Hapmap data format error: pedigree data on line " + (i + 1) + ".");
                }
                if (hd.countTokens() > 7) {
                    throw new MDRPedFileException(
                            "Hapmap data format error: pedigree data on line " + (i + 1) + ".");
                }
                hd.nextToken();
                indName = hd.nextToken();
                hapMapTranslate.put(indName, (String) hapsData.get(i));
            }
        }
        Individual ind;

        this.allIndividuals = new ArrayList();

        // enumerate indivs
        StringTokenizer st = new StringTokenizer((String) lines.get(0),
                "\n\t\" \"");
        int numMetaColumns = 0;
        boolean doneMeta = false;
        boolean genoErrorB = false;
        while (!doneMeta && st.hasMoreTokens()) {
            String thisfield = st.nextToken();
            numMetaColumns++;
            // first indiv ID will be a string beginning with "NA"
            if (thisfield.startsWith("NA")) {
                doneMeta = true;
            }
        }
        numMetaColumns--;

        st = new StringTokenizer((String) lines.get(0), "\n\t\" \"");
        for (int i = 0; i < numMetaColumns; i++) {
            st.nextToken();
        }
        ArrayList namesIncludingDups = new ArrayList();
        StringTokenizer dt;
        while (st.hasMoreTokens()) {
            // todo: sort out how this used to work. now it's counting the header line so we subtract 1
            ind = new Individual(numLines - 1);

            String name = st.nextToken();
            namesIncludingDups.add(name);
            if (name.endsWith("dup")) {
                // skip dups (i.e. don't add 'em to ind array)
                continue;
            }
            String details = (String) hapMapTranslate.get(name);
            if (details == null) {
                throw new MDRPedFileException("Hapmap data format error: " + name);
            }
            dt = new StringTokenizer(details, "\n\t\" \"");
            ind.setFamilyID(dt.nextToken().trim());
            ind.setIndividualID(dt.nextToken().trim());
            ind.setDadID(dt.nextToken().trim());
            ind.setMomID(dt.nextToken().trim());
            try {
                ind.setGender(Integer.parseInt(dt.nextToken().trim()));
                ind.setAffectedStatus(Integer.parseInt(dt.nextToken().trim()));
            } catch (NumberFormatException nfe) {
                throw new MDRPedFileException(
                        "File error: invalid gender or affected status for indiv " + name);
            }

            // check if the family exists already in the Hashtable
            Family fam = (Family) this.families.get(ind.getFamilyID());
            if (fam == null) {
                // it doesnt exist, so create a new Family object
                fam = new Family(ind.getFamilyID());
            }
            fam.addMember(ind);
            this.families.put(ind.getFamilyID(), fam);
            this.allIndividuals.add(ind);
        }

        // start at k=1 to skip header which we just processed above.
        hminfo = new String[numLines - 1][];
        for (int k = 1; k < numLines; k++) {
            StringTokenizer tokenizer = new StringTokenizer((String) lines.get(k));
            // reading the first line
            if (colNum < 0) {
                // only check column number count for the first line
                colNum = tokenizer.countTokens();
            }
            if (colNum != tokenizer.countTokens()) {
                // this line has a different number of columns
                // should send some sort of error message
                // TODO: add something which stores number of markers for all lines and checks that they're consistent
                throw new MDRPedFileException(
                        "Line number mismatch in input file. line " + (k + 1));
            }

            if (tokenizer.hasMoreTokens()) {
                hminfo[k - 1] = new String[2];
                for (int skip = 0; skip < numMetaColumns; skip++) {
                    // meta-data crap
                    String s = tokenizer.nextToken().trim();

                    // get marker name, chrom and pos
                    if (skip == 0) {
                        hminfo[k - 1][0] = s;
                    }
                    if (skip == 2) {
                        String dc = Chromosome.getDataChrom();
                        if (dc != null && !dc.equals("none")) {
                            if (!dc.equalsIgnoreCase(s)) {
                                throw new MDRPedFileException(
                                        "Hapmap file format error on line " + (k + 1) + ":\n The file appears to contain multiple chromosomes:" + "\n" + dc + ", " + s);
                            }
                        } else {
                            Chromosome.setDataChrom(s);
                        }
                    }
                    if (skip == 3) {
                        hminfo[k - 1][1] = s;
                    }
                    if (skip == 5) {
                        Chromosome.setDataBuild(s);
                    }
                }
                int index = 0;
                int indexIncludingDups = -1;
                while (tokenizer.hasMoreTokens()) {
                    String alleles = tokenizer.nextToken();

                    indexIncludingDups++;
                    // we've skipped the dups in the ind array, so we skip their genotypes
                    if (((String) namesIncludingDups.get(indexIncludingDups)).endsWith("dup")) {
                        continue;
                    }

                    ind = (Individual) allIndividuals.get(index);
                    int[] checker1, checker2;
                    try {
                        checker1 = checkGenotype(alleles.substring(0, 1));
                        checker2 = checkGenotype(alleles.substring(1, 2));
                    } catch (NumberFormatException nfe) {
                        throw new MDRPedFileException(
                                "Invalid genotype on individual " + ind.getIndividualID() + ".");
                    }
                    if (checker1[1] != checker2[1]) {
                        genoErrorB = !genoErrorB;
                    }
                    byte allele1 = (byte) checker1[0];
                    byte allele2 = (byte) checker2[0];
                    ind.addMarker(allele1, allele2);
                    if (genoErrorB) {
                        throw new MDRPedFileException(
                                "File input error: individual " + ind.getIndividualID() + ", marker " + this.hminfo[ind.getNumMarkers() - 1][0] + ".\nFor any marker, an individual's genotype must be only letters or only numbers.");
                    }
                    index++;
                }
            }
        }
    }

    public void parseHapsFile(ArrayList individs) throws MDRPedFileException {
        // This method is used to parse haps files which now go through similar processing to ped files.
        String currentLine;
        byte[] genos = new byte[0];
        String ped, indiv;
        int numLines = individs.size();
        if (numLines == 0) {
            throw new MDRPedFileException("Data format error: empty file");
        }

        Individual ind = null;
        this.allIndividuals = new ArrayList();
        int lineCount = 0;
        int numTokens = 0;
        ArrayList chromA = new ArrayList();
        ArrayList chromB = new ArrayList();
        boolean hapsEven = false;
        boolean hapsError = false;

        for (int i = 0; i < numLines; i++) {
            lineCount++;
            currentLine = (individs.get(i)).toString();
            if (currentLine.length() == 0) {
                continue;
            }
            StringTokenizer st = new StringTokenizer(currentLine);
            // first two tokens are expected to be ped, indiv
            if (st.countTokens() > 2) {
                ped = st.nextToken();
                indiv = st.nextToken();
            } else {
                throw new MDRPedFileException("Genotype file error:\nLine " + lineCount + " appears to have fewer than 3 columns.");
            }
            if (hapsEven) {
                ind = new Individual(st.countTokens());
                ind.setFamilyID(ped);
                ind.setIndividualID(indiv);
                ind.setDadID("");
                ind.setMomID("");
                ind.setGender(0);
                ind.setAffectedStatus(0);
            }

            // all other tokens are loaded into a vector (they should all be genotypes)
            genos = new byte[st.countTokens()];

            int q = 0;

            if (numTokens == 0) {
                numTokens = st.countTokens();
            }
            if (numTokens != st.countTokens()) {
                throw new MDRPedFileException("Genotype file error:\nLine " + lineCount + " appears to have an incorrect number of entries");
            }
            // Allowed for A/C/G/T input in Haps files.
            while (st.hasMoreTokens()) {
                String thisGenotype = (String) st.nextElement();
                if (!hapsEven) {
                    chromA.add(thisGenotype);
                } else {
                    chromB.add(thisGenotype);
                }
                if (thisGenotype.equals("h")) {
                    genos[q] = 9;
                } else if (thisGenotype.equalsIgnoreCase("A")) {
                    genos[q] = 1;
                } else if (thisGenotype.equalsIgnoreCase("C")) {
                    genos[q] = 2;
                } else if (thisGenotype.equalsIgnoreCase("G")) {
                    genos[q] = 3;
                } else if (thisGenotype.equalsIgnoreCase("T")) {
                    genos[q] = 4;
                } else {
                    try {
                        genos[q] = Byte.parseByte(thisGenotype);
                    } catch (NumberFormatException nfe) {
                        throw new MDRPedFileException(
                                "Genotype file input error:\ngenotype value \"" + thisGenotype + "\" on line " + lineCount + " not allowed.");
                    }
                }
                // Allele values other then 0-4 or 9 generate exceptions.
                if ((genos[q] < 0 || genos[q] > 4) && (genos[q] != 9)) {
                    throw new MDRPedFileException(
                            "Genotype file input error:\ngenotype value \"" + genos[q] + "\" on line " + lineCount + " not allowed.");
                }
                q++;
            }

            if (hapsEven) {
                for (int m = 0; m < chromA.size(); m++) {
                    if (((String) chromA.get(m)).equals("h")) {
                        chromA.set(m, "9");
                    } else if (((String) chromA.get(m)).equalsIgnoreCase("A")) {
                        chromA.set(m, "1");
                        hapsError = !hapsError;
                    } else if (((String) chromA.get(m)).equalsIgnoreCase("C")) {
                        chromA.set(m, "2");
                        hapsError = !hapsError;
                    } else if (((String) chromA.get(m)).equalsIgnoreCase("G")) {
                        chromA.set(m, "3");
                        hapsError = !hapsError;
                    } else if (((String) chromA.get(m)).equalsIgnoreCase("T")) {
                        chromA.set(m, "4");
                        hapsError = !hapsError;
                    }
                    if (((String) chromB.get(m)).equals("h")) {
                        chromB.set(m, "9");
                    } else if (((String) chromB.get(m)).equalsIgnoreCase("A")) {
                        chromB.set(m, "1");
                        hapsError = !hapsError;
                    } else if (((String) chromB.get(m)).equalsIgnoreCase("C")) {
                        chromB.set(m, "2");
                        hapsError = !hapsError;
                    } else if (((String) chromB.get(m)).equalsIgnoreCase("G")) {
                        chromB.set(m, "3");
                        hapsError = !hapsError;
                    } else if (((String) chromB.get(m)).equalsIgnoreCase("T")) {
                        chromB.set(m, "4");
                        hapsError = !hapsError;
                    }
                    if (hapsError) {
                        throw new MDRPedFileException(
                                "File input error: Individual " + ind.getFamilyID() + " strand " + ind.getIndividualID() + ", marker " + (m + 1) + ".\nFor any marker, an individual's genotype must be only letters or only numbers.");
                    }
                    byte allele1 = Byte.parseByte(chromA.get(m).toString());
                    byte allele2 = Byte.parseByte(chromB.get(m).toString());
                    ind.addMarker(allele1, allele2);
                }
                // check if the family exists already in the Hashtable
                Family fam = (Family) this.families.get(ind.getFamilyID());
                if (fam == null) {
                    // it doesnt exist, so create a new Family object
                    fam = new Family(ind.getFamilyID());
                }

                if (fam.getMembers().containsKey(ind.getIndividualID())) {
                    throw new MDRPedFileException("Individual " + ind.getIndividualID() + " in family " + ind.getFamilyID() + " appears more than once.");
                }

                fam.addMember(ind);
                this.families.put(ind.getFamilyID(), fam);
                this.allIndividuals.add(ind);
                chromA = new ArrayList();
                chromB = new ArrayList();
            }
            hapsEven = !hapsEven;
        }
        if (hapsEven) {
            // we're missing a line here
            throw new MDRPedFileException(
                    "Genotype file appears to have an odd number of lines.\n" + "Each individual is required to have two chromosomes");
        }

    }

    public int[] checkGenotype(String allele) throws MDRPedFileException {
        // This method cleans up the genotype checking process for hap map and ped files & allows for both numerical and
        // alphabetical input.
        int[] genotype = new int[2];

        if (allele.equalsIgnoreCase("N")) {
            genotype[0] = 0;
        } else if (allele.equalsIgnoreCase("A")) {
            genotype[0] = 1;
        } else if (allele.equalsIgnoreCase("C")) {
            genotype[0] = 2;
        } else if (allele.equalsIgnoreCase("G")) {
            genotype[0] = 3;
        } else if (allele.equalsIgnoreCase("T")) {
            genotype[0] = 4;
        } else {
            genotype[0] = Integer.parseInt(allele.trim());
            genotype[1] = 1;
        }

        return genotype;
    }

    public ArrayList check() throws MDRPedFileException, PedFileException {

        // Before we perform the check we want to prune out individuals with too
        // much missing data or trios which contain individuals with too much
        // missing data

        Iterator fitr = families.values().iterator();
        ArrayList useable = new ArrayList();
        while (fitr.hasNext()) {
            Family curFam = (Family) fitr.next();
            Enumeration indIDEnum = curFam.getMemberList();
            Vector victor = new Vector();
            while (indIDEnum.hasMoreElements()) {
                victor.add(curFam.getMember((String) indIDEnum.nextElement()));
            }

            PedParser pp = new PedParser();
            try {
                SimpleGraph sg = pp.buildGraph(victor, getMissingThreshold());
                ArrayList indStrings = new ArrayList(pp.parsePed(sg));
                if (indStrings != null) {
                    Iterator sitr = indStrings.iterator();
                    while (sitr.hasNext()) {
                        useable.add(curFam.getMember((String) sitr.next()));
                    }
                }
            } catch (PedigreeException pe) {
                String pem = pe.getMessage();
                if (pem.indexOf("one parent") != -1) {
                    indIDEnum = curFam.getMemberList();
                    while (indIDEnum.hasMoreElements()) {
                        curFam.getMember((String) indIDEnum.nextElement()).setReasonImAxed(pem);
                    }
                } else {
                    throw new MDRPedFileException(pem + "\nin family " + curFam.getFamilyName());
                }
            }
        }

        unrelatedIndividuals = useable;

        ArrayList indList = (ArrayList) allIndividuals.clone();
        Individual currentInd;
        Family currentFamily;

        // deal with individuals who are missing too much data
        for (int x = 0; x < indList.size(); x++) {
            currentInd = (Individual) indList.get(x);
            currentFamily = getFamily(currentInd.getFamilyID());

            if (currentInd.getGenoPC() < 1 - getMissingThreshold()) {
                allIndividuals.remove(currentInd);
                axedPeople.add(currentInd);
                currentInd.setReasonImAxed("% Genotypes: " + new Double(currentInd.getGenoPC() * 100).intValue());
                currentFamily.removeMember(currentInd.getIndividualID());
                if (currentFamily.getNumMembers() == 0) {
                    // if everyone in a family is gone, we remove it from the list
                    families.remove(currentInd.getFamilyID());
                }
            } else if (!useable.contains(currentInd)) {
                axedPeople.add(currentInd);
                if (currentInd.getReasonImAxed() == null) {
                    currentInd.setReasonImAxed("Not a member of maximum unrelated subset.");
                }
            }
        }
        if (useable.size() == 0) {
            // todo: this should be more specific about the problems.
            throw new MDRPedFileException(
                    "File contains zero valid individuals.");
        }

        CheckData cd = new CheckData(this);
        ArrayList results = cd.check();
        this.results = results;
        return results;
    }

    public String[][] getHMInfo() {
        return hminfo;
    }

    public ArrayList getResults() {
        return results;
    }

    public void setResults(ArrayList res) {
        results = res;
    }

    public ArrayList getAxedPeople() {
        return axedPeople;
    }

    public boolean isBogusParents() {
        return bogusParents;
    }

    public ArrayList getTableData() {
        ArrayList tableData = new ArrayList();
        int numResults = results.size();
        markerRatings = new int[numResults];
        dups = new int[numResults];
        ArrayList title = new ArrayList(10);
        title.add(new String("Position"));
        title.add(new String("Marker"));
        title.add(new String("ObsHet"));
        title.add(new String("PreHet"));
        title.add(new String("HWpvalue"));
        title.add(new String("GenoPercent"));
        title.add(new String("FamTrio"));
        title.add(new String("MendError"));
        title.add(new String("MAF"));
        title.add(new String("MinorAllele"));
        tableData.add(title);
        for (int i = 0; i < numResults; i++) {
            ArrayList tempVect = new ArrayList(10);
            MarkerResult currentResult = (MarkerResult) results.get(i);
            tempVect.add(new Integer(i + 1));

            tempVect.add(markerInfor.get(i));
            tempVect.add(new Double(currentResult.getObsHet()));
            tempVect.add(new Double(currentResult.getPredHet()));
            tempVect.add(new Double(currentResult.getHWpvalue()));
            tempVect.add(new Double(currentResult.getGenoPercent()));
            tempVect.add(new Integer(currentResult.getFamTrioNum()));
            tempVect.add(new Integer(currentResult.getMendErrNum()));
            tempVect.add(new Double(currentResult.getMAF()));
            tempVect.add(currentResult.getMinorAllele());

            markerRatings[i] = currentResult.getRating();
            tableData.add(tempVect.clone());
        }

        return tableData;
    }

    public ArrayList getMarkerInformation() {
        return markerInfor;
    }

    public ArrayList getMarkerInformation(int[] subsetMarker) {
    	if(subsetMarker.length == markerInfor.size() || subsetMarker == null) {
    		return markerInfor;
    	} else {
    		ArrayList mk = new ArrayList();
    		for(int i = 0; i < subsetMarker.length; i++) {
    			mk.add((String) markerInfor.get(subsetMarker[i]));
    		}
    		return mk;
    	}
    }
    public int[] getMarkerRatings() {
        return markerRatings;
    }

    public int[] getDups() {
        return dups;
    }

    public ArrayList getColumnNames() {
        ArrayList c = new ArrayList();
        c = new ArrayList();
        c.add("#");
        if (Chromosome.getUnfilteredMarker(0).getName() != null) {
            c.add("Name");
            c.add("Position");
        }
        c.add("ObsHET");
        c.add("PredHET");
        c.add("HWpval");
        c.add("%Geno");
        c.add("FamTrio");
        c.add("MendErr");
        c.add("MAF");
        c.add("M.A.");
        c.add("Rating");
        return c;
    }

    public void saveCheckDataToText(File outfile) throws IOException {
        FileWriter checkWriter = null;
        if (outfile != null) {
            checkWriter = new FileWriter(outfile);
        } else {
            throw new IOException("Error saving checkdata to file.");
        }

        ArrayList names = getColumnNames();
        int numCols = names.size();
        StringBuffer header = new StringBuffer();
        for (int i = 0; i < numCols; i++) {
            header.append(names.get(i)).append("\t");
        }
        header.append("\n");
        checkWriter.write(header.toString());

        ArrayList tableData = getTableData();
        for (int i = 0; i < tableData.size(); i++) {
            StringBuffer sb = new StringBuffer();
            ArrayList row = (ArrayList) tableData.get(i);
            // don't print the true/false vals in last column
            for (int j = 0; j < numCols - 1; j++) {
                sb.append(row.get(j)).append("\t");
            }
            // print BAD if last column is false
            if (((Boolean) row.get(numCols - 1)).booleanValue()) {
                sb.append("\n");
            } else {
                sb.append("BAD\n");
            }
            checkWriter.write(sb.toString());
        }
        checkWriter.close();
    }

    public void setWhiteList(HashSet whiteListedCustomMarkers) {
        whitelist = whiteListedCustomMarkers;
    }

    public boolean isWhiteListed(SNP snp) {
        return whitelist.contains(snp);
    }

    public boolean isHaploidHets() {
        return haploidHets;
    }

    public void setHaploidHets(boolean haploidHets) {
        this.haploidHets = haploidHets;
    }

    public void m() throws MDRPedFileException {
        Enumeration famList = this.families.keys();
        MendErrorTrace met;
        while (famList.hasMoreElements()) {
            String fid = (String) famList.nextElement();

            Family Fam = getFamily(fid);
            Enumeration sibList = Fam.getMemberList();

            FamilyStruct FamStr = getFamilyStruct(fid);

            while (sibList.hasMoreElements()) {
                String iid = (String) sibList.nextElement();
                try {
                    Individual ind = Fam.getMember(iid);
                    Person per = FamStr.getPerson(iid);
                    ArrayList VM = new ArrayList();
                    for (int i = 0; i < ind.getNumMarkers(); i++) {
                        if (!FamStr.getPersons().containsKey(per.getPersonID())) {
                            throw new MDRPedFileException("Person " + per.getPersonID() + " in family " + ind.getFamilyID() + " can't be matched.");
                        }
                        String marker;
                        if (ind.getZeroed(i)) {
                            marker = new String(PublicData.MissingGenotype);
                            met = new MendErrorTrace(fid, iid,
                                    (String) markerInfor.get(i), i);
                            menderrortrace.add(met);
                        } else {
                            byte m1 = ind.getAllele(i, 0);
                            byte m2 = ind.getAllele(i, 1);
                            if (m2 < m1) {
                                byte m3 = m1;
                                m1 = m2;
                                m2 = m3;
                            }
                            marker = new String((new Byte(m1).toString()) + (new Byte(m2).toString()));
                        }
                        VM.add(marker);
                    }
                    per.setGenotype(VM);
                } catch (MDRPedFileException e) {
                    System.err.println("no such person with ID " + iid + " in family " + fid + "in ped files.");
                } catch (PedFileException e) {
                    System.err.println("no such individual with ID " + iid + " in family " + fid + "in ped file.");
                }
            }
        }

        Enumeration fsList = this.familystructure.keys();
        ArrayList error_report = new ArrayList();
        while (fsList.hasMoreElements()) {
            FamilyStruct fs = (FamilyStruct) getFamilyStruct((String) fsList.nextElement());
            TreeMap Ps;
            TreeMap Ks;
            GenoSet gSet;
            for (int i = 0; i < getNumMarkers(); i++) {
                Ps = new TreeMap();
                Ks = new TreeMap();
                Enumeration perList = fs.getPersonList();
                while (perList.hasMoreElements()) {
                    Person per = fs.getPerson((String) perList.nextElement());
                    String genotype = per.getGenotype(i);
                    if (fs.hasAncestor(per.getPersonID())) {
                        if (Ks.containsKey(genotype)) {
                            Integer c = ((Integer) Ks.get(genotype));
                            int v = (c.intValue());
                            v++;
                            c = new Integer(v);
                            Ks.put(genotype, c);
                        } else {
                            Integer c = new Integer(1);
                            Ks.put(new String(genotype), c);
                        }
                    } else {
                        if (Ps.containsKey(genotype)) {
                            Integer c = ((Integer) Ps.get(genotype));
                            int v = (c.intValue());
                            v++;
                            c = new Integer(v);
                            Ps.put(genotype, c);
                        } else {
                            Integer c = new Integer(1);
                            Ps.put(new String(genotype), c);
                        }
                    }
                }

                gSet = new GenoSet(Ps, Ks, i);
                if (gSet.getchildrenGenoMap().size() == 0) {
                    error_report.add(new String("Family " + fs.getFamilyStructName() + " Locus " + (String) markerInfor.get(i)));
                }
                fs.addObservedGenoSet(gSet);
            }
        }
    /*
     * if( error_report.size() > 0 ) { Iterator it = error_report.iterator(); for( ; it.hasNext(); ) {
     * System.err.println( (String) it.next() ); } throw new MDRPedFileException("empty set of children genosets");
     * }
     */
    }

    public void RabinowitzApproach(boolean forNontransmitted, int[] subsetMarker) {
        Enumeration fsList = this.familystructure.keys();
        boolean informative;
        String fid;
        while (fsList.hasMoreElements()) {
            fid = (String) fsList.nextElement();
            FamilyStruct fs = (FamilyStruct) getFamilyStruct(fid);
            if (!((Boolean) famInformative.get(fid)).booleanValue()) {
                System.err.println("Omitted" + fs.getFamilyStructName());
                continue;
            }
            try {
            	if (forNontransmitted) {
            		fs.NontransmittedProc(markerInfor, subsetMarker);
            	} else {
            		fs.RabinowitzProc(markerInfor, subsetMarker);
            	}
            } catch (FamilyStructException E) {
                System.err.println("Exception in family " + fs.getFamilyStructName() + " when in RabinowitzApproach.");
            }
        }
    }

    /**
     * return the record of menderrortrace(); Allele2Genotype should be invoked precedingly;
     * 
     * @return
     */
    public ArrayList getMendErrorTrace() {
        return menderrortrace;
    }

    public void printFamilyStruct() {
        String[] FID = new String[famInformative.size()];
        int idx = 0;
        Set keys = famInformative.keySet();
        for (Iterator i = keys.iterator(); i.hasNext();) {
            FID[idx++] = (String) i.next();
        }
        Arrays.sort(FID);

        for (int i = 0; i < FID.length; i++) {
            String fid = FID[i];
            FamilyStruct fs = getFamilyStruct(fid);
            Hashtable htPerson = fs.getPseudoPersons();
            Enumeration pid = htPerson.keys();
            String[] PID = new String[htPerson.size()];
            idx = 0;
            while (pid.hasMoreElements()) {
                PID[idx++] = (String) pid.nextElement();
            }
            Arrays.sort(PID);
            for (int j = 0; j < PID.length; j++) {
                PseudoPerson ps = (PseudoPerson) htPerson.get(PID[j]);
                System.out.print(ps.getPseudoPersonID() + "\t");
                ArrayList mk = (ArrayList) ps.getGenotype();
                for (int k = 0; k < mk.size(); k++) {
                    System.out.print(mk.get(k) + "\t");
                }
                System.out.println();
            }
        }
//        Enumeration fkeys = famInformative.keys();
//        String fid;
//        while (fkeys.hasMoreElements()) {
//            fid = (String) fkeys.nextElement();
//            FamilyStruct fs = getFamilyStruct(fid);
//            Hashtable psperson = fs.getPseudoPersons();
//            Enumeration pid = psperson.keys();
//            while (pid.hasMoreElements()) {
//                PseudoPerson ps = (PseudoPerson) psperson.get((String) pid.nextElement());
//                System.out.print(ps.getPseudoPersonID() + "\t");
//                ArrayList mk = (ArrayList) ps.getGenotype();
//                for (int i = 0; i < mk.size(); i++) {
//                    System.out.print(mk.get(i) + "\t");
//                }
//                System.out.println();
//            }
//        }
    }
}