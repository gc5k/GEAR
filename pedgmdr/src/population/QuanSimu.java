package population;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.math.random.RandomDataImpl;

import score.CalEngine;
import score.CalEngineException;
import family.GMDRData;
import family.MDRPedFileException;

public class QuanSimu {
	int replication;
	int samplesize;
	int numLocus;
	int dimension;
	int numFC;
	String cell[];
	int numKids=2;
	int numAffected=1;

	String Allele[][]; 
	double frequency[][];
	double recombination[];
	ArrayList Parents;
	ArrayList Children;
	ArrayList PCovariate;
	ArrayList PObservation;
	ArrayList CCovariate;
	ArrayList CObservation;
	ArrayList Traits;
	ArrayList PTraits;
	ArrayList PQuan;
	ArrayList CQuan;

	double error;
	double variation;
	double intercept;

	double covariable;
	double geneAffect;
	double threshold;

	String model;
	RandomDataImpl randomData;	
	double missingrate;
	
	QuanSimu (int replication, int samplesize, int numLocus, String m)
	{
		this.replication = replication;
		this.samplesize = samplesize;
		this.numLocus = numLocus;
		this.model = new String(m);

		Parents = new ArrayList();
		Children = new ArrayList();
		CCovariate = new ArrayList();
		CObservation = new ArrayList();
		PCovariate = new ArrayList();
		PObservation = new ArrayList();
		Traits = new ArrayList();
		PTraits = new ArrayList();
		PQuan = new ArrayList();
		CQuan = new ArrayList();
		
		error=0;
		variation=0;
		intercept=0;
		geneAffect=0;
		threshold=0;
		randomData = new RandomDataImpl();
		randomData.reSeed(20316017);
		missingrate = 0;
	}
	
	public int affection(double obs)
	{
		if( model.compareTo("B") == 0 )//binary
		{
			double prob = Math.exp(obs)/(1+Math.exp(obs));
			if( Math.random() < prob )
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			if( obs > threshold )
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}
/*
	public void create()
	{
		String pair[][];
		ArrayList ptemp;
		ArrayList ptrait;
		ArrayList ctemp;
		ArrayList ctrait;
		ArrayList pcovariate;
		ArrayList ccovariate;
		ArrayList pobs;
		ArrayList cobs;
		String P[][];
		String child[][];
		int i=0;
		int TP=0;
		while( i < samplesize )
		{
			ptemp = new ArrayList();
			ptrait = new ArrayList();
			pcovariate = new ArrayList();
			pobs = new ArrayList();
			for( int j=0; j<2; j++ )
			{
				pair = new String[numLocus][2];
				for( int l=0; l<numLocus; l++)
				{
					for( int k=0; k<2; k++)
					{
						double rd = Math.random();
						int index = 0;
						double fq = frequency[l][index];
						while( rd > fq )
						{
							fq+=frequency[l][++index];
						}
						pair[l][k]=Allele[l][index];
					}
				}
				ptemp.add(pair);
				Integer genestatus = getStatus(pair);
				double covariate;
				double obs;
				if( model.compareTo("B") == 0 )
				{
					covariate = random(0, variation);
					obs = genestatus.doubleValue() * geneAffect + covariate*covariable+ intercept;
				}
				else
				{
					covariate = random(0, variation);
					obs = genestatus.doubleValue() * geneAffect + covariate*covariable + intercept + random(0,error);					
				}

				Integer status = new Integer(affection(obs));

				pcovariate.add( new Double(covariate) );
				ptrait.add(status);
				pobs.add(new Double(obs));
			}

			ctemp = new ArrayList();
			ctrait = new ArrayList();
			ccovariate = new ArrayList();
			cobs = new ArrayList();
			int affect=0;
			for( int ii = 0; ii<numKids; ii++)
			{
				child = new String[numLocus][2];
				for( int j=0; j<2; j++)
				{
					P = (String[][]) ptemp.get(j );
					int chr=0;
					for( int k=0; k<numLocus; k++)
					{
						double rd = Math.random();
						if( rd>recombination[k] )
						{
							chr = 1-chr;
						}
						child[k][j] = P[k][chr];
					}
				}

				Integer genestatus = getStatus(child);
				double covariate;
				double obs;

				if( model.compareTo("B") == 0 )
				{
					covariate = random(0, variation);
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept;
				}
				else
				{
					covariate = random(0, variation) > 0 ? 1:0;
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept + random(0, error);
				}

				Integer status = new Integer(affection(obs));

				affect+=status.intValue();
				ccovariate.add(new Double(covariate));
				ctrait.add(status);
				ctemp.add(child);
				cobs.add(new Double(obs));
			}
			Children.add(ctemp);
			Parents.add(ptemp);
			PCovariate.add(pcovariate);
			CCovariate.add(ccovariate);
			PTraits.add(ptrait);
			Traits.add(ctrait);
			PQuan.add(pobs);
			CQuan.add(cobs);
			i++;
		}
	}
*/
	
	public void create2(int[] kid, int[] affkid, int CovDis)//0 for bernoulli distribution;
	{
		String pair[][];
		ArrayList ptemp;
		ArrayList ptrait;
		ArrayList ctemp;
		ArrayList ctrait;
		ArrayList pcovariate;
		ArrayList pobservation;
		ArrayList ccovariate;
		ArrayList cobservation;
		String P[][];
		String child[][];
		int SIB[] = new int[kid.length];
		int AFF[] = new int[affkid.length];
		System.arraycopy(kid, 0, SIB, 0, kid.length);
		System.arraycopy(affkid, 0, AFF, 0, affkid.length);
		int i=0;
		int TP=0;
		int cgb=0;
		while( i < samplesize )
		{
			cgb++;
			ptemp = new ArrayList();
			ptrait = new ArrayList();
			pcovariate = new ArrayList();
			pobservation = new ArrayList();
			for( int j=0; j<2; j++ )
			{
				pair = new String[numLocus][2];
				for( int l=0; l<numLocus; l++)
				{
					for( int k=0; k<2; k++)
					{
						double rd = Math.random();
						int index = 0;
						double fq = frequency[l][index];
						while( rd > fq )
						{
							fq+=frequency[l][++index];
						}
						pair[l][k]=Allele[l][index];
					}
				}
				ptemp.add(pair);
				Integer genestatus = getStatus(pair);
				double obs;
				double covariate;
				if ( CovDis == 0 )
				{
					covariate = ( Math.random() < 0.5 ) ? 0 : 1 ;
				}
				else
				{
					covariate = random(0, variation);
				}
				if( model.compareTo("B") == 0 )
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept;
				}
				else
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept + random(0, error);					
				}
				Integer status = new Integer(affection(obs));
				pcovariate.add( new Double(covariate) );
				pobservation.add( new Double(obs) );
				ptrait.add(status);
			}

			ctemp = new ArrayList();
			ctrait = new ArrayList();
			ccovariate = new ArrayList();
			cobservation = new ArrayList();
			int affect=0;

			int sib;
			if( i < 200 )//family type1
			{
				numAffected=AFF[0];
				sib=SIB[0];
			}
			else if ( i < 400 )
			{
				numAffected=AFF[1];
				sib=SIB[1];
			}
			else
			{
				numAffected=AFF[2];
				sib=SIB[2];
			}

			for( int ii = 0; ii<sib; ii++)
			{
				child = new String[numLocus][2];
				for( int j=0; j<2; j++)
				{
					P = (String[][]) ptemp.get(j);
					int chr=0;
					for( int k=0; k<numLocus; k++)
					{
						double rd = Math.random();
						if( rd>recombination[k] )
						{
							chr = 1-chr;
						}
						child[k][j] = P[k][chr];
					}
				}

				Integer genestatus = getStatus(child);
				double obs;
				double covariate;
				if ( CovDis == 0 )
				{
					covariate = ( Math.random() < 0.5 ) ? 0 : 1 ;
				}
				else
				{
					covariate = random(0, variation);
				}
				if( model.compareTo("B") == 0 )
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept;
				}
				else
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept + random(0, error);
				}

				Integer status = new Integer(affection(obs));

				affect+=status.intValue();
				if( status.intValue() == 1 && genestatus.intValue() == 1 )
				{
					TP++;
				}
				ccovariate.add(new Double(covariate));
				cobservation.add(new Double(obs));
				ctrait.add(status);
				ctemp.add(child);
			}

			if( affect >= numAffected )
			{
				Children.add(ctemp);
				Parents.add(ptemp);
				PCovariate.add(pcovariate);
				PObservation.add(pobservation);
				CCovariate.add(ccovariate);
				CObservation.add(cobservation);
				PTraits.add(ptrait);
				Traits.add(ctrait);
				i++;				
			}
		}
		System.out.println("CGB="+cgb);
		System.out.println("TP="+TP);
	}

	public void create3(int[] kid, int[] affkid, int CovDis)//0 for bernoulli distribution
	{
		String pair[][];
		ArrayList ptemp;
		ArrayList ptrait;
		ArrayList ctemp;
		ArrayList ctrait;
		ArrayList pcovariate;
		ArrayList ccovariate;
		String P[][];
		String child[][];
		int SIB[] = new int[kid.length];
		int AFF[] = new int[affkid.length];
		System.arraycopy(kid, 0, SIB, 0, kid.length);
		System.arraycopy(affkid, 0, AFF, 0, affkid.length);
		int i=0;
		int TP=0;
		Parents.clear();
		Children.clear();
		while( i < samplesize )
		{
			ptemp = new ArrayList();
			ptrait = new ArrayList();
			pcovariate = new ArrayList();
			for( int j=0; j<2; j++ )
			{
				pair = new String[numLocus][2];
				for( int l=0; l<numLocus; l++)
				{
					for( int k=0; k<2; k++)
					{
						double rd = Math.random();
						int index = 0;
						double fq = frequency[l][index];
						while( rd > fq )
						{
							fq+=frequency[l][++index];
						}
						pair[l][k]=Allele[l][index];
					}
				}
				ptemp.add(pair);
				Integer genestatus = getStatus(pair);
				double obs;
				double covariate;
				if ( CovDis == 0 )
				{
					covariate = ( Math.random() < 0.5 ) ? 0 : 1 ;
				}
				else
				{
					covariate = random(0, variation);
				}
				
				if( model.compareTo("B") == 0 )
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept;
				}
				else
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept + random(0,error);					
				}
				Integer status = new Integer(affection(obs));
				pcovariate.add( new Double(covariate) );
				ptrait.add(status);
			}

			ctemp = new ArrayList();
			ctrait = new ArrayList();
			ccovariate = new ArrayList();
			int affect=0;

			int sib;
			if( i < 200 )//family type1
			{
				numAffected=AFF[0];
				sib=SIB[0];
			}
			else if ( i < 400 )
			{
				numAffected=AFF[1];
				sib=SIB[1];
			}
			else
			{
				numAffected=AFF[2];
				sib=SIB[2];
			}

			for( int ii = 0; ii<sib; ii++)
			{
				child = new String[numLocus][2];
				for( int j=0; j<2; j++)
				{
					P = (String[][]) ptemp.get(j);
					int chr=0;
					for( int k=0; k<numLocus; k++)
					{
						double rd = Math.random();
						if( rd>recombination[k] )
						{
							chr = 1-chr;
						}
						child[k][j] = P[k][chr];
					}
				}

				Integer genestatus = getStatus(child);
				double obs;
				double covariate;
				if ( CovDis == 0 )
				{
					covariate = ( Math.random() < 0.5 ) ? 0 : 1 ;
				}
				else
				{
					covariate = random(0, variation);
				}
				if( model.compareTo("B") == 0 )
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept;
				}
				else
				{
					obs = genestatus.doubleValue()*geneAffect + covariate*covariable + intercept + random(0, error);
				}

				Integer status = new Integer(affection(obs));
				affect+=status.intValue();
				if( status.intValue() == 1 && genestatus.intValue() == 1 )
				{
					TP++;
				}
				ccovariate.add(new Double(covariate));
				ctrait.add(status);
				ctemp.add(child);
			}
			Children.add(ctemp);
			Parents.add(ptemp);
			i++;
		}
		System.out.println("TP="+TP);
	}
	


	public ArrayList getParentCovariate()
	{
		return PCovariate;
	}

	public ArrayList getChildrenCovariate()
	{
		return CCovariate;
	}

	public ArrayList getTrait()
	{
		return Traits;
	}

	public ArrayList getPTrait()
	{
		return PTraits;
	}

	public ArrayList getChildren()
	{
		return Children;
	}
	
	public Integer getStatus(String[][] geno)
	{
		String genostr = new String();
		int dim = cell[0].length()/2;
		for( int i=0; i<dim; i++)
		{
			if( geno[i][0].compareTo(geno[i][1]) > 0 )
			{
				genostr+=geno[i][1]+geno[i][0];
			}
			else
			{
				genostr+=geno[i][0]+geno[i][1];
			}
		}

		for( int i=0; i<cell.length; i++)
		{
			if( genostr.compareTo(cell[i]) == 0 )
			{
				if( i==1 && cell.length == 3 )
				{
					return new Integer(1);
				}
				else
				{
					return new Integer(2);
				}
			}
		}
		return new Integer(0);
	}
	public ArrayList getParents()
	{
		return Parents;
	}
	

	public void print(String ped, String phe, int[][] miss) throws IOException
	{
	    PrintWriter pedout = new PrintWriter(new File(ped));
	    PrintWriter pheout = new PrintWriter(new File(phe));
	    
		pedout.print("FID\tID\tFA\tM0\tSex\tAffection\t");
		for( int i=0; i<Allele.length; i++)
		{
			pedout.print("M"+(i+1)+"\t");
		}
		pedout.println();

		pheout.println("FIP\tID\tpheno\tcov");
		int Missing[][] = new int[miss.length][miss[0].length];
		for( int i=0; i<miss.length; i++)
		{
			System.arraycopy(miss[i], 0, Missing[i], 0, miss[i].length);
		}
		for(int i=0; i<samplesize; i++)
		{
			int []MR;
			if( i<200 )
			{
				MR=Missing[0];
			}
			else if( i<400 )
			{
				MR=Missing[1];
			}
			else
			{
				MR=Missing[2];
			}
			int FID=2000+i;
			int ID=FID*100;
			int FA=ID;
			int MO=ID+1;
			ArrayList temp = (ArrayList) Parents.get(i);
			ArrayList pt_temp = (ArrayList) PTraits.get(i);
			ArrayList pc_temp = (ArrayList) PCovariate.get(i);
			ArrayList pob_temp = (ArrayList) PObservation.get(i);
			String pair1[][] = (String[][]) temp.get(0);
			String pair2[][] = (String[][]) temp.get(1);

			pedout.print(FID+"\t"+ID+ "\t0\t0\t0\t"+ pt_temp.get(0) +"\t");
			
			for( int j=0; j<Allele.length; j++)
			{
				double rd = Math.random();
				if( rd > MR[0] )
				{
					pedout.print(pair1[j][0]+" "+pair1[j][1]+"\t");
				}
				else
				{
					pedout.print("0 0\t");
				}
			}
			pedout.println();
			pheout.println(FID+"\t"+ID+"\t"+pob_temp.get(0)+"\t"+pc_temp.get(0));
			ID++;

			pedout.print(FID+"\t"+ID+ "\t0\t0\t1\t"+ pt_temp.get(1) +"\t");
			for( int j=0; j<Allele.length; j++)
			{
				double rd = Math.random();
				if( rd > MR[1] )
				{
					pedout.print(pair2[j][0]+" "+pair2[j][1]+"\t");
				}
				else
				{
					pedout.print("0 0\t");
				}
			}
			pedout.println();
			pheout.println(FID+"\t"+ID+"\t"+pob_temp.get(1)+"\t"+pc_temp.get(1));
			ID++;

		    ArrayList children = (ArrayList) Children.get(i);
		    ArrayList t_temp = (ArrayList) Traits.get(i);
		    ArrayList cp_temp = (ArrayList) CCovariate.get(i);
		    ArrayList cob_temp = (ArrayList) CObservation.get(i);
		    for( int j=0; j<children.size(); j++)
		    {
		    	String child[][] = (String[][]) children.get(j);
		    	pedout.print(FID + "\t" + ID + "\t" + FA + "\t" + MO + "\t1\t" + t_temp.get(j) + "\t");

		    	for( int k=0; k<child.length; k++)
		    	{
		    		pedout.print(child[k][0]+" "+child[k][1]+"\t");
		    	}
				pedout.println();
		    	pheout.println(FID+"\t"+ID+"\t"+cob_temp.get(j)+"\t"+cp_temp.get(j));
		    	ID++;
		    }
		}
		pedout.close();
		pheout.close();
	}

	public void printGMDR(String TDTped, String TDTpheM, String TDTphe, double[][] Mres, double[][] score, double[][] trait, ArrayList Marker, int datatype, int pick) throws IOException
	{
	    PrintWriter TDTpedout = new PrintWriter(new File(TDTped));
	    PrintWriter TDTpheout = new PrintWriter(new File(TDTphe));
	    PrintWriter TDTpheoutM = new PrintWriter(new File(TDTpheM));
	    
	    ArrayList marker;
	    
	    for( int i=0; i<((ArrayList) Marker.get(0)).size(); i++)
	    {
	    	TDTpedout.print("M"+(1+i)+"\t");
	    }
	    TDTpedout.println("class");

	    for( int i=0; i<Marker.size(); i+=2)
	    {
	    	marker = (ArrayList) Marker.get(i);
	    	for( int j=0; j<marker.size(); j++)
	    	{
	    		TDTpedout.print(marker.get(j)+"\t");
	    	}
	    	TDTpedout.println((new Double(trait[i/2][0])).intValue());
	    	marker = (ArrayList) Marker.get(i+1);
	    	for( int j=0; j<marker.size(); j++)
	    	{
	    		TDTpedout.print(marker.get(j)+"\t");
	    	}
	    	TDTpedout.println((1-(new Double(trait[i/2][0])).intValue()));
	    }
	    TDTpedout.close();

	    for( int i=0; i<score[0].length; i++)
	    {
	    	TDTpheout.print("score"+(i+1)+"\t");
	    }
	    TDTpheout.println();

	    for( int i=0; i<score.length; i++)
	    {
	    	for( int j=0; j<score[i].length; j++)
	    	{
	    		TDTpheout.print(score[i][j]+"\t");
	    	}
	    	TDTpheout.println();
	    	for( int j=0; j<score[i].length; j++)
	    	{
	    		TDTpheout.print((-1)*score[i][j]+"\t");
	    	}
	    	TDTpheout.println();
	    }
	    TDTpheout.close();
	    
	    for( int i=0; i<Mres[0].length; i++)
	    {
	    	TDTpheoutM.print("score"+(i+1)+"\t");
	    }
	    TDTpheoutM.println();
	    
	    for( int i=0; i<Mres.length; i++)
	    {
	    	for( int j=0; j<Mres[i].length; j++)
	    	{
	    		TDTpheoutM.print(Mres[i][j]+"\t");
	    	}
	    	TDTpheoutM.println();
	    	for( int j=0; j<Mres[i].length; j++)
	    	{
	    		TDTpheoutM.print((-1)*Mres[i][j]+"\t");
	    	}
	    	TDTpheoutM.println();
	    }
	    TDTpheoutM.close();
	}
	
	public double random(double intc, double var)
	{
		return randomData.nextGaussian(intc, var);
	}
	
	public void setNumKids(int i)
	{
		numKids=i;
	}

	public void setMissingrate(double rate)
	{
		missingrate = rate;
	}
	
	public void setNumAffected(int na)
	{
		numAffected = na;
	}
	
	public void setSeed(int seed)
	{
		randomData.reSeed(seed);
	}

	public void setAffect(double gene)
	{
		geneAffect = gene;
	}

	public void setParameters(double intc, double cov, double var, double err)
	{
		intercept = intc;
		covariable = cov;
		variation = var;
		error = err;
	}
	
	public void setThreshold(double t)
	{
		threshold = t;
	}

	public void setAllele ( String[][] allele, double[][] freq, double[] recombination )
	{
		Allele = new String[allele.length][2];
		frequency = new double[freq.length][2];
		this.recombination = new double[recombination.length];

		System.arraycopy(recombination, 0, this.recombination, 0, recombination.length);

		for( int i=0; i<Allele.length; i++)
		{
			System.arraycopy(allele[i], 0, Allele[i], 0, allele[i].length);
			System.arraycopy(freq[i], 0, this.frequency[i], 0, freq[i].length);
		}
	}
	
	public void setCell(String[] Cell)
	{
		cell = new String[Cell.length];
		System.arraycopy(Cell, 0, cell, 0, Cell.length);
	}

	public static void main(String[] args)
	{
		double intc = 0;
		double var = 1; //3.16
		double err = 1; //residual
		double gene = 0.25; //gene affect;
		double cov = 1;
		double threshold = 2.05;
		String model = "C";
		int[] Kid={2, 4, 3};
		int[] AffKid={2, 2, 2};
		int[][] MissingRate={{0,0}, {1,1},{0,1}};
		
		int kids = 4;
		int AffectNum = 1;
		int SIMU = 200;
		int Linear = 0;
		int samplesize = 600;
		int CovDis = 1; //0 for bernoulli distribution; 1 for normal distribution
		double rate =1;

		int datatype = 0;//0 for all; 1 for either affected or unaffected;
		int pick = 0;//0 for affected and 1 for unaffected, when datatype = 1;	
		String allele[][] = {
							{"1", "2"},{"1", "2"},{"1", "2"},{"1", "2"},{"1", "2"},
							{"1", "2"},{"1", "2"},{"1", "2"},{"1", "2"},{"1", "2"}
			  			  };
		double freq[][] = {
							{0.5, 0.5},{0.5, 0.5},{0.5, 0.5},{0.5, 0.5},{0.5, 0.5},
							{0.5, 0.5},{0.5, 0.5},{0.5, 0.5},{0.5, 0.5},{0.5, 0.5}
						};
		double rec[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

		String CELL[] = //{"12"};
			 {"1111", "1212", "2222"};
//			{"1112", "1211", "1222", "2212"};
			//{"221211", "221112", "122211", "112212", "121122", "111222", "121212"};
			/*
			{
							"22221111", "22112211", "22111122",
							"11222211", "11221122", "11112222",
							"22121211", "22121112", "22111212",
							"12221211", "12221112", "11221212",
							"12122211", "12112212", "11122212",
							"12121122", "12111222", "11121222",
					       	"12121212"
			};*/

		for( int i=0; i<SIMU; i++)
		{
			System.out.println("Simulation:"+i);
			String Ped= (new Integer(i)).toString()+".ped";
			String Phe= (new Integer(i)).toString()+".phe";
			QuanSimu qs=new QuanSimu(100, samplesize, rec.length, model);
			qs.setSeed(i);
			qs.setAffect(gene);
			qs.setNumAffected(AffectNum);//Affect Num of Sibs
			
			qs.setAllele(allele, freq, rec);
			qs.setCell(CELL);
			qs.setParameters(intc, cov, var, err);
			qs.setThreshold(threshold);
			qs.setNumKids(kids);
			qs.setMissingrate(rate);
			qs.create2(Kid, AffKid, CovDis);
//			qs.create3(Kid, AffKid, CovDis);

			try
			{
				qs.print(Ped, Phe, MissingRate);
			}
			catch (IOException e)
			{
				System.err.println("io");
			}
			
//////////////////////////////////////////////////////////////
			boolean filetype = true;
			GMDRData GD = new GMDRData(filetype);
			File pf = new File(Ped);
			File gf = new File(Phe);
			try
			{
				GD.Initial(pf, gf);
				GD.Match();
				GD.PedData.Allele2Genotype();
				GD.PedData.NonTransmittedGenoType();
				GD.CreateTable();
			}
			catch (MDRPedFileException e)
			{
				System.err.println("test");
			}
			catch (IOException e)
			{
				System.err.println("IOException.");
			}
			double trait[][] = GD.getStatus();
			double Covs[][] = GD.getCovariates();
			double covs[][] = new double[Covs.length][1];
			for( int iii =0 ; iii<Covs.length; iii++)
			{
				covs[iii][0] = Covs[iii][1];
			}

		    String TDTPed= new String("Ped"+(new Integer(i)).toString()+".txt");
		    String TDTPhe= new String("Phe"+(new Integer(i)).toString()+".txt");
		    String TDTPheM = new String("PheM"+(new Integer(i)).toString()+".txt");
		    ArrayList Table = GD.getTDTGenoTable();
		    System.out.println( Table.size() );

		    int CovIndex[] = {1};
		    int PhenoIndex = 0;
		    boolean adjust = true;

		    try
		    {
				adjust = true;
				CalEngine CE = new CalEngine(covs, trait, adjust);
				double[][] res = CE.GeneralScore(Linear);
				
		    	GD.CalculateScore(CovIndex, PhenoIndex, adjust, Linear);

		    	ArrayList PITable = GD.getPersonIndexTable();
		    	ArrayList FullPI = GD.getFullPersonIndexTable();
		    	System.out.println("FullPI="+FullPI.size());    
//		    	double[][] res = GD.pickupScore(PITable);

		    	CovIndex[0]=1;
		    	PhenoIndex=0;
		    	adjust = false;
		    	GD.CalculateScore(CovIndex, PhenoIndex, adjust, Linear);
//		    	double[][] Mres = GD.pickupScore(PITable);
		    	double[][] Mres=GD.getCovariates();

		    	qs.printGMDR(TDTPed, TDTPheM, TDTPhe, Mres, res, trait, Table, datatype, pick);
		    }
		    catch (CalEngineException e)
		    {
		    	System.err.println("Calculation Engine Exception");
		    	e.printStackTrace(System.err);
		    }
		    catch (IOException e)
		    {
		    	System.err.print("IO Excepiton in printGMDR" );
		    }
			System.out.println("finished the first round: "+i);
		}
	}
}