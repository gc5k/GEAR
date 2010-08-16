package im;

import java.util.ArrayList;

import im.population.IMPopulation;
import im.IntervalPriorProbability;
import publicAccess.ToolKit;

/**
 * 
 * @author Guo-Bo Chen
 */
public class IMBMatrix {

	IMPopulation imp;
	GenomeScan gs;
	int[] SNPIdx;
	int[][] ChrInt;
	ArrayList id;
	int order;

	public IMBMatrix(GenomeScan g, IMPopulation im) {
		gs = g;
		imp = im;
	}

	public double[][] getICIMMatrixAtPoint(String com, double[][] coeff,
			int interval) {
		SNPIdx = ToolKit.StringToIntArray(com);
		double[][] matrix = new double[imp.IndividualNumber()][coeff.length+1];
		ChrInt();
		IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
		for (int i = 0; i < SNPIdx.length; i++) {
			iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
		}
		for (int i = 0; i < imp.IndividualNumber(); i++) {
			matrix[i][0] = 1;
			int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[0][0],
					ChrInt[0][1]);
			for (int j = 0; j < coeff.length; j++) {
				for (int k = 0; k < iip[0].NumQTLtypes(); k++) {
					matrix[i][1+j] += iip[0].PriorProbabilityAt(IIPIdx, interval, k) * coeff[j][k];
				}
			}
		}
		return matrix;
	}

	public double[][] getCIMMatrixAtPoint(ArrayList selectedMarker, String com,
			double[][] coeff, int interval) {
		SNPIdx = ToolKit.StringToIntArray(com);
		double[][] matrix;
		if (selectedMarker == null) {
			matrix = new double[imp.IndividualNumber()][imp.MarkerNumber()
					- (SNPIdx.length * coeff.length + 1) + 2];
			ChrInt();
			IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
			for (int i = 0; i < SNPIdx.length; i++) {
				iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
			}
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				int markerindex = 0;
				int c = coeff.length * SNPIdx.length + 1;
				matrix[i][0] = 1;
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int jj = 0; jj < SNPIdx.length; jj++) {
						int IIPIdx = gs.getIPPRowIndexForIndividual(i,
								ChrInt[jj][0], ChrInt[jj][1]);
						for (int jjj = 0; jjj < coeff.length; jjj++) {
							for (int k = 0; k < iip[jj].NumQTLtypes(); k++) {
								matrix[i][1 + jj * coeff.length + jjj] += iip[jj]
										.PriorProbabilityAt(IIPIdx, interval, k)
										* coeff[jj][k];
							}
						}
					}
					for (int k = 0; k < imp.MarkerNumber(j); k++) {
						boolean flag = true;
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (markerindex == SNPIdx[kk]
									|| markerindex == (SNPIdx[kk] + 1)) {
								flag = false;
								break;
							}
						}
						if (flag) {
							matrix[i][c++] = imp.MarkerAt(i, j, k);
						}
						markerindex++;
					}
				}
			}
		} else {
			int m = 0;
			for (int i = 0; i < selectedMarker.size(); i++) {
				int mi = ((Integer) selectedMarker.get(i)).intValue();
				if ((SNPIdx[0] + 1) == mi || (SNPIdx[0]) == mi) {
					m++;
				}
			}
			matrix = new double[imp.IndividualNumber()][selectedMarker.size()
					- m + 2];
			ChrInt();
			IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
			for (int i = 0; i < SNPIdx.length; i++) {
				iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
			}
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				int c = coeff.length * SNPIdx.length + 1;
				matrix[i][0] = 1;
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int jj = 0; jj < SNPIdx.length; jj++) {
						int IIPIdx = gs.getIPPRowIndexForIndividual(i,
								ChrInt[jj][0], ChrInt[jj][1]);
						for (int jjj = 0; jjj < coeff.length; jjj++) {
							for (int k = 0; k < iip[jj].NumQTLtypes(); k++) {
								matrix[i][1 + jj * coeff.length + jjj] += iip[jj]
										.PriorProbabilityAt(IIPIdx, interval, k)
										* coeff[jj][k];
							}
						}
					}
					for (int k = 0; k < selectedMarker.size(); k++) {
						int mi = ((Integer) selectedMarker.get(k)).intValue();
						boolean flag = true;
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (mi == SNPIdx[kk]
									|| mi == (SNPIdx[kk] + 1)) {
								flag = false;
								break;
							}
						}
						if (flag) {
							matrix[i][c++] = imp.MarkerAt(i, j, mi);
						}
					}
				}
			}
		}
		return matrix;
	}

	public double[][] getNullCIMMatrix(ArrayList selectedMarker, String com) {
		SNPIdx = ToolKit.StringToIntArray(com);
		double[][] matrix;
		if (selectedMarker == null) {
			matrix = new double[imp.IndividualNumber()][imp.MarkerNumber()
					- (SNPIdx.length + 1) + 1];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 1;
				int markerindex = 0;
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int k = 0; k < imp.MarkerNumber(j); k++) {
						boolean flag = true;
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (markerindex == SNPIdx[kk]
									|| markerindex == (SNPIdx[kk] + 1)) {
								flag = false;
							}
						}
						if (flag) {
							matrix[i][c++] = imp.MarkerAt(i, j, k);
						}
						markerindex++;
					}
				}
			}
		} else {
			int m = 0;
			for (int i = 0; i < selectedMarker.size(); i++) {
				if ((SNPIdx[0] + 1) == ((Integer) selectedMarker.get(i))
						.intValue()
						|| (SNPIdx[0]) == ((Integer) selectedMarker.get(i))
								.intValue()) {
					m++;
				}
			}
			matrix = new double[imp.IndividualNumber()][selectedMarker.size()
					- (m) + 1];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 1;
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int k = 0; k < selectedMarker.size(); k++) {
						int mi = ((Integer) selectedMarker.get(k)).intValue();
						boolean flag = true;
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (mi == SNPIdx[kk] || mi == (SNPIdx[kk] + 1)) {
								flag = false;
							}
						}
						if (flag) {
							matrix[i][c++] = imp.MarkerAt(i, j, mi);
						}
					}
				}
			}
		}
		return matrix;
	}

	public double[][] getFullMatrix(ArrayList selectedMarker) {
		double[][] matrix;
		if (selectedMarker == null) {
			matrix = new double[imp.IndividualNumber()][imp.MarkerNumber() + 1];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 0;
				for (int j = 0; j < imp.ChromosomeNumberOriginal(); j++) {
					for (int k = 0; k < imp.MarkerNumber(j); k++) {
						matrix[i][++c] = imp.MarkerAt(i, j, k);
					}
				}
			}
		} else {
			matrix = new double[imp.IndividualNumber()][selectedMarker.size() + 1];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 0;
				for (int j = 0; j < imp.ChromosomeNumberOriginal(); j++) {
					for (int k = 0; k < selectedMarker.size(); k++) {
						matrix[i][++c] = imp.MarkerAt(i, j, (Integer) selectedMarker.get(k));
					}
				}
			}
		}
		return matrix;
	}

	public double[][] getPPMatrix(String com, double[][] coeff) {
		SNPIdx = ToolKit.StringToIntArray(com);
		double[][] matrix = new double[imp.IndividualNumber()][SNPIdx.length
				* (coeff.length) + 1];
		ChrInt();
		IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
		for (int i = 0; i < SNPIdx.length; i++) {
			iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
		}
		for (int i = 0; i < imp.IndividualNumber(); i++) {
			matrix[i][0] = 1;
			for (int j = 0; j < SNPIdx.length; j++) {
				int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[j][0],
						ChrInt[j][1]);
				for (int jj = 0; jj < coeff.length; jj++) {
					for (int k = 0; k < iip[j].NumQTLtypes(); k++) {
						matrix[i][1 + j * coeff.length + jj] += iip[j]
								.MidPriorProbability(IIPIdx, k)
								* coeff[j][k];
					}
				}
			}
		}
		return matrix;
	}

	public double[][] getPPMatrixAtPoint(String com, int w, double[][] coeff) {
		SNPIdx = ToolKit.StringToIntArray(com);
		double[][] matrix = new double[imp.IndividualNumber()][SNPIdx.length
				* (coeff.length) + 1];
		ChrInt();
		IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
		for (int i = 0; i < SNPIdx.length; i++) {
			iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
		}
		for (int i = 0; i < imp.IndividualNumber(); i++) {
			matrix[i][0] = 1;
			for (int j = 0; j < SNPIdx.length; j++) {
				int IIPIdx = gs.getIPPRowIndexForIndividual(i, ChrInt[j][0],
						ChrInt[j][1]);
				for (int jj = 0; jj < coeff.length; jj++) {
					for (int k = 0; k < iip[j].NumQTLtypes(); k++) {
						matrix[i][1 + j * coeff.length + jj] += iip[j]
								.PriorProbabilityAt(IIPIdx, w, k)
								* coeff[jj][k];
					}
				}
			}
		}
		return matrix;
	}

	public void setOrder(int o) {
		order = o;
		ChrInt = new int[order][2];
	}

	public String getString() {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < ChrInt.length; i++) {
			sb.append(ChrInt[i][0] + "," + ChrInt[i][1] + "\t");
		}
		return sb.toString();
	}

	private void ChrInt() {
		int c = 0;
		int idx = 0;
		for (int i = 0; i < imp.ChromosomeNumber(); i++) {
			for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
				if (c == SNPIdx[idx]) {
					ChrInt[idx][0] = i;
					ChrInt[idx][1] = j;
					idx++;
					if (idx == SNPIdx.length) {
						break;
					}
				}
				c++;
			}
			if (idx == SNPIdx.length) {
				break;
			}
		}
	}
}
