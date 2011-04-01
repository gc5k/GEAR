package score.dependonweka;

import weka.core.*;
import weka.core.matrix.Matrix;
import weka.classifiers.functions.LinearRegression;

/**
 * It implements AbstractScore.  It is utilized when the linear regression should be used to calculate score 
 * @author Guo-Bo Chen, Zhejiang University, gc5k@zju.edu.cn
 */
public class LinearScore implements AbstractScore {

    private double[] Residual;//residual score
    private Matrix Response;//Y
    private Matrix Observe;//X
    private Attribute[] AttriX;//Attribute X
    private Attribute AttriY;// Attribute Y
    private Instances LS;//dataset
    private double ridge;//ridge, default 0;
    boolean intercept;//intercept, default true;
    private double[] coefficients;
    /**
    0, "M5 method",
    (step through the attributes removing the one
    with the smallest standardised coefficient until no improvement is observed
    in the estimate of the error given by the Akaike
    information criterion)
    1, "No attribute selection",
    2, "Greedy method"
    a greedy selection using the Akaike information metric.
     */
    private int Selection;//

    /**
     * Construct the object, and initialize the predictors and response.
     * @param X Predictors
     * @param Y Response
     */
    public LinearScore(double[][] X, double[][] Y) {
        ridge = 0;
        ridge = 0;
        Selection = 1;
        intercept = true;
        Observe = new Matrix(X);
        Response = new Matrix(Y);
        Residual = new double[Response.getRowDimension()];
        CreateInstance();
    }

    /**
     * Set model selection method.
     * @param i should be no less than 0 and no greater than 2.
     */
    public void SetModelSelectionMethod(int i) {
        if (i >= 0 && i <= 2) {
            Selection = i;
        }
    }

    /**
     * Constructing model with the given predictors and response.
     */
    private void CreateInstance() {
        AttriX = new Attribute[Observe.getColumnDimension()];
        AttriY = new Attribute("Response");
        FastVector attributes = new FastVector(Observe.getColumnDimension() + Response.getColumnDimension());
        for (int i = 0; i < AttriX.length; i++) {
            AttriX[i] = new Attribute("X" + (String) (new Integer(i).toString()));
            attributes.addElement(AttriX[i]);
        }
        attributes.addElement(AttriY);
        LS = new Instances("LS", attributes, Observe.getRowDimension());
        LS.setClassIndex(AttriY.index());

        for (int i = 0; i < Observe.getRowDimension(); i++) {
            Instance ins = new Instance(Observe.getColumnDimension() + Response.getColumnDimension());
            for (int j = 0; j < Observe.getColumnDimension(); j++) {
                ins.setValue(AttriX[j], Observe.get(i, j));
            }
            ins.setValue(AttriY, Response.get(i, 0));
            ins.setDataset(LS);
            LS.add(ins);
        }
    }

    public void SetRidge(double r) {
        ridge = r;
    }

    public void SetIntercept(boolean i) {
        intercept = i;
    }

    public double[] GetScore() {
        return Residual;
    }

    public double[] getCoefficients() {
        LinearRegression Reg = new LinearRegression();
        SelectedTag TS = Reg.getAttributeSelectionMethod();
        SelectedTag T = new SelectedTag(Selection, TS.getTags());
        //System.out.println(((Tag)T.getSelectedTag()).getReadable());
        Reg.setAttributeSelectionMethod(T);

        try {
            Reg.buildClassifier(LS);
        } catch (Exception err) {
            err.printStackTrace();
        }
        System.out.println(Reg);
        double[] coe = Reg.coefficients();	//something weird for linear regression in weka is the order of coefficients are:
        //for Y=b0+b1*X1+b2*X2+e
        //the return coefficients are b1, b2, 0(for e), b0.
        double[] co = new double[coe.length - 1];
        System.arraycopy(coe, 0, co, 1, coe.length - 2);
        co[0] = coe[coe.length - 1];
        return co;
    }

    public void CalculateScore() {
        LinearRegression Reg = new LinearRegression();
        SelectedTag TS = Reg.getAttributeSelectionMethod();
        SelectedTag T = new SelectedTag(Selection, TS.getTags());
        //System.out.println(((Tag)T.getSelectedTag()).getReadable());
        Reg.setAttributeSelectionMethod(T);

        try {
            Reg.buildClassifier(LS);
        } catch (Exception err) {
            err.printStackTrace();
        }
        System.out.println(Reg);
        double[] co = Reg.coefficients();//something weird for linear regression in weka is the order of coefficients are:
        //for Y=b0+b1*X1+b2*X2+e
        //the return coefficients are b1, b2, 0(for e), b0.
        for (int j = 0; j < Observe.getRowDimension(); j++) {
            double ex = 0;
            for (int k = 0; k < Observe.getColumnDimension(); k++) {
                ex += Observe.get(j, k) * co[k];
            }
            ex += co[co.length - 1];
            Residual[j] = Response.get(j, 0) - ex;
        }
    }

    public static double[] stringTOdouble(String[] s) {
        if (s == null) {
            return null;
        }
        double[] d = new double[s.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = Double.parseDouble(s[i]);
        }
        return d;
    }

    public static void main(String[] args) {
        double[][] vals = {{5.6},
            {12},
            {11}
        };

        Matrix A = new Matrix(vals);
        double[][] response = {{5.6}, {12.2}, {11.3}};
        Matrix Y = new Matrix(response);
        LinearScore QS = new LinearScore(vals, response);
        QS.SetModelSelectionMethod(1);
        QS.SetIntercept(false);
        QS.CalculateScore();

        double[] residual = QS.GetScore();
        for (int i = 0; i < residual.length; i++) {
            System.out.println(residual[i]);
        }
    }
}
