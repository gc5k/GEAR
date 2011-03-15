package admixture.phenotype;

import admixture.AdmixtureConstant;

public class QualityControl {
	private int N_AffectedKids;

	private int scheme; //0 for family 
	private CC cc;
	static class CC {
		int cases;
		int controls;
		int count_cases;
		int count_controls;
		public CC(int cs, int cl) {
			cases = cs;
			controls = cl;
		}
		
		public boolean AcceptSubject (int c) {
			boolean flag = true;
			if(c>0) {
				if(count_cases < cases) {
					count_cases++;
				} else {
					flag = false;
				}
			} else {
				if(count_controls < controls) {
					count_controls++;
				} else {
					flag = false;
				}
			}
			return flag;
		}
	}

	public QualityControl(int ak, int sc) {
		N_AffectedKids = ak;
		scheme = sc;
	}

	public QualityControl(int cs, int cl, int sc) {
		cc = new CC(cs, cl);
		scheme = sc;
	}

	public boolean Accept(FamilyPhenotype fp) {
		boolean flag = false;
		switch (scheme) {
			case AdmixtureConstant.FamilyExactAffected: flag = Exact(fp); break;
			case AdmixtureConstant.FamilyMoreThanAffected: flag = MoreThan(fp); break;
			case AdmixtureConstant.CaseControl: flag = AcceptSubject(fp); break;
			case AdmixtureConstant.No_selection: flag = true; break; 
		}
		return flag;
	}

	private boolean AcceptSubject(FamilyPhenotype fp) {
		int s = fp.getNumberAffectedOffspring();
		return cc.AcceptSubject(s);
	}

	private boolean MoreThan(FamilyPhenotype fp) {
		int s = fp.getNumberAffectedOffspring();
		return s >= N_AffectedKids ? true:false;
	}
	
	private boolean Exact(FamilyPhenotype fp) {
		int s = fp.getNumberAffectedOffspring();
		return s == N_AffectedKids ? true:false;
	}

}
