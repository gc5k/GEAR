package family.pedigree;

public class MendErrorTrace {

    public String FID;
    public String PID;
    public String Marker;
    public int index;

    MendErrorTrace(String FID, String PID, String Marker, int index) {
        this.FID = new String(FID);
        this.PID = new String(PID);
        this.Marker = new String(Marker);
        this.index = index;
    }
}
