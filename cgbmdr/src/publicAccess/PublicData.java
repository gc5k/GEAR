package publicAccess;

public class PublicData {

    public static boolean UniversalThreshould = true;
    public static int MissingAllele = 0;
    public static int tieValue = 1;
    public static double epsilon = 0.001;
    public static double tooBigRatio = 5;
    public static String MissingPhenotypeValue = ".";
    public static String delim = "\t";
    public static String seperator = ",";
    public static String MissingGenotype = "00";
    public static String IgnorFile = "ND";
    
    public static int LinearSearch = 0;
    public static int TreeSearch = 1;
    public static int DynamicSearch = 2;
    public static int GenotypeSearch = 3;

    public static int NumOfStatistics = 2;
    public static int TestingAccuIdx = 0;
    public static int TrainingAccuIdx = 1;

    public static int RandomPartition = 0;
    public static int UnpairedPartition = 1;
    public static int PairedPartition = 2;

}
