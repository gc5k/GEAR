package score.dependonweka;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class CalEngine {

    /**
     * @param args
     */
    static private int Linear = 1;
    static private int Logistic = 2;
    private double[][] Covariable;
    private double[][] Traits;
    int NumofTrait;
    private double[][] res;
    private boolean adjustment;

    public CalEngine(final double[][] cov, final double[][] trait, boolean adjust) throws CalEngineException {
        adjustment = adjust;

        if (adjustment == true) {
            if (cov == null || trait == null || cov.length != trait.length) {
                throw new CalEngineException("unequal column sizes of covariables and traits.");
            }

            Covariable = new double[cov.length][cov[0].length];
            for (int i = 0; i < Covariable.length; i++) {
                System.arraycopy(cov[i], 0, Covariable[i], 0, cov[i].length);
            }
        } else {
            if (trait == null) {
                throw new CalEngineException("empty trait.");
            }
        }

        Traits = new double[trait.length][trait[0].length];
        for (int i = 0; i < trait.length; i++) {
            System.arraycopy(trait[i], 0, Traits[i], 0, trait[i].length);
        }

        NumofTrait = Traits[0].length;
        res = new double[Traits.length][Traits[0].length];
    }

    public boolean isAdjusted() {
        return adjustment;
    }

    public double[] getCoefficents(int[] covInd, int pheInd, int model) throws CalEngineException {
        if (pheInd > (Traits.length - 1)) {
            throw new CalEngineException("Selected phenotype does not exist.");
        }
        AbstractScore GS;
        double[] co;
        if (adjustment) {
            double[][] Y = new double[Traits.length][1];
            double[][] X = new double[Covariable.length][covInd.length];
            for (int k = 0; k < Traits.length; k++) {
                Y[k][0] = Traits[k][pheInd];
            }
            for (int k = 0; k < Covariable.length; k++) {
                for (int kk = 0; kk < covInd.length; kk++) {
                    X[k][kk] = Covariable[k][covInd[kk]];
                }
            }

            if (model == Linear) {
                GS = new LinearScore(X, Y);
            } else {
                GS = new LogisticScore(X, Y);
            }
            co = GS.getCoefficients();
        } else {
            double[][] Y = new double[Traits.length][1];
            for (int k = 0; k < Traits.length; k++) {
                Y[k][0] = Traits[k][pheInd];
            }
            if (model == Linear) {
                GS = new LinearMeanScore(Y);
            } else {
                GS = new LogisticMeanScore(Y);
            }
            co = GS.getCoefficients();
        }
        return co;
    }

    public double[][] GeneralScore(int index)//index=0 if method is linear regression;
    //index=1 if method is logistic regression;
    {
        AbstractScore[] GS = new AbstractScore[NumofTrait];
        if (adjustment) {
            for (int i = 0; i < NumofTrait; i++) {
                double[][] Y = new double[Traits.length][1];
                for (int k = 0; k < Traits.length; k++) {
                    Y[k][0] = Traits[k][i];
                }
                if (index == Linear) {
                    GS[i] = new LinearScore(Covariable, Y);
                } else {
                    GS[i] = new LogisticScore(Covariable, Y);
                }
                GS[i].CalculateScore();
                double[] s = GS[i].GetScore();
                for (int k = 0; k < s.length; k++) {
                    res[k][i] = s[k];
                }
            }
        } else {
            for (int i = 0; i < NumofTrait; i++) {
                double[][] Y = new double[Traits.length][1];

                for (int k = 0; k < Traits.length; k++) {
                    Y[k][0] = Traits[k][i];
                }

                if (index == Linear) {
                    GS[i] = new LinearMeanScore(Y);
                } else {
                    GS[i] = new LogisticMeanScore(Y);
                }
                GS[i].CalculateScore();
                double[] s = GS[i].GetScore();
                for (int k = 0; k < s.length; k++) {
                    res[k][i] = s[k];
                }
            }
        }
        return res;
    }

    public static void main(String[] args) {
        System.out.println("CalEngine!");

        double[][] Covs = {{1, 3}, {1, 5}, {1, 4}, {1, 1}, {1, 3}, {1, 3}, {1, 1}, {1, 3}, {1, 1}, {1, 3}};
//		double[][] Covs = null;
        double[][] trait = {{1}, {1}, {0}, {0}, {0}, {0}, {0}, {1}, {1}, {0}};
//		double[][] Covs = 
//		{{3.614492, 4.745349},{4.299145, 7.057928},{7.485585, 2.483148},{2.850424, 3.561575},
//		{6.75514, 4.667191},{5.02019, 6.043154},{2.76813, 4.214377},{3.852767, 3.469689},
//		{7.476863, 4.471005},{2.562243, 1.061865}};
		/*
        double[][] trait = 
        { {1, 20.7128223, 14.856477}, {1, 23.5494665, 25.602254}, {1, 23.7693979, 23.162992},
        {0, 13.6457537, 15.451552}, {0, 21.1541691, 23.939506}, {0, 17.9310275, 18.33978},
        {0, 12.5772315, 9.993274}, {0, 12.3063936, 12.33613},{0, 27.4465537, 26.473369},
        {1, 0.9049123,	5.728055}};
         */
        boolean adjust = false;
        CalEngine CE;
        try {
            CE = new CalEngine(Covs, trait, adjust);//Covs are the data from Covariable data table from "score page";
            //trait are the data from Trait data table "score page";
            int[] covInd = {0};
            int pheInd = 0;
            int model = Logistic;
            double co[] = CE.getCoefficents(covInd, pheInd, model);
            for (int i = 0; i < co.length; i++) {
                System.out.println("coefficients  " + " " + co[i]);
            }
            double[][] res = CE.GeneralScore(model); //res are the double[][] return to you, and you have 
        //to fill the data into the score table in "filter" page	
/*
        for( int i=0; i < res.length; i++)
        {
        for( int j=0; j < res[i].length; j++)
        {
        System.out.print(res[i][j]+"\t");
        }
        System.out.println();
        }
         */
        } catch (CalEngineException e) {
            System.err.println("PedFile initial exception.");
            e.printStackTrace(System.err);
        }
    }
}
