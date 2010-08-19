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
		ChrInt();
		double[][] matrix;
		if (selectedMarker == null) {
			matrix = new double[imp.IndividualNumber()][1 + imp.MarkerNumber()
					- (SNPIdx.length * 2) + (SNPIdx.length * coeff.length)];

			IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
			for (int i = 0; i < SNPIdx.length; i++) {

				iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
			}
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				int c = coeff.length * SNPIdx.length + 1;
				matrix[i][0] = 1;
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
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int k = 0; k < imp.MarkerNumber(j); k++) {
						boolean flag = true;
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (j == ChrInt[kk][0] && (k == ChrInt[kk][1] || k == (ChrInt[kk][1] + 1))) {
								continue;
							}
							matrix[i][c++] = imp.MarkerAt(i, j, k);
						}
					}
				}
			}
		} else {
			int m = 0;
			ChrInt();
			for (int i = 0; i < SNPIdx.length; i++) {
				for (int j = 0; j < selectedMarker.size(); j++) {
					ArrayList mp = (ArrayList) selectedMarker.get(j);
					int chr = ((Integer) mp.get(0)).intValue();
					int mk = ((Integer) mp.get(1)).intValue();
					if(chr == ChrInt[i][0] && (ChrInt[i][1] == mk || mk == (ChrInt[i][1] + 1))) {
						m++;
					}
				}
			}
			matrix = new double[imp.IndividualNumber()][1 + selectedMarker.size()
					- m + SNPIdx.length * coeff.length];
			IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
			for (int i = 0; i < SNPIdx.length; i++) {
				iip[i] = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
			}
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				int c = coeff.length * SNPIdx.length + 1;
				matrix[i][0] = 1;
				int index = 0;
				for (int jj = 0; jj < SNPIdx.length; jj++) {
					int IIPIdx = gs.getIPPRowIndexForIndividual(i,
							ChrInt[jj][0], ChrInt[jj][1]);
					for (int jjj = 0; jjj < coeff.length; jjj++) {
						for (int k = 0; k < iip[jj].NumQTLtypes(); k++) {
							matrix[i][1 + jj * coeff.length + jjj] += iip[jj].PriorProbabilityAt(IIPIdx, interval, k) * coeff[jjj][k];
						}
					}
				}
				for (int j = 0; j < selectedMarker.size(); j++) {
					ArrayList mp = (ArrayList) selectedMarker.get(j);
					int chr = ((Integer) mp.get(0)).intValue();
					int mk = ((Integer) mp.get(1)).intValue();
					boolean flag = true;
					for (int k = 0; k < SNPIdx.length; k++) {
						if(chr == ChrInt[k][0] && (mk == ChrInt[k][1] || mk == (ChrInt[k][1] + 1))) {
							flag = false;
							break;
						}
					}
					if(flag) {
						matrix[i][c++] = imp.MarkerAt(i, chr, mk);
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
					- (SNPIdx.length * 2) + 1];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 1;
				for (int j = 0; j < imp.ChromosomeNumber(); j++) {
					for (int k = 0; k < imp.MarkerNumber(j); k++) {
						boolean flag = true;
						ChrInt();
						for (int kk = 0; kk < SNPIdx.length; kk++) {
							if (j == ChrInt[kk][0] && (k == ChrInt[kk][1] ||  k == (ChrInt[kk][1] + 1))) {
								continue;
							}
							matrix[i][c++] = imp.MarkerAt(i, j, k);
						}
					}
				}
			}
		} else {
			int m = 0;
			ChrInt();
			for (int i = 0; i < SNPIdx.length; i++) {
				for (int j = 0; j < selectedMarker.size(); j++) {
					ArrayList mp = (ArrayList) selectedMarker.get(j);
					int chr = ((Integer) mp.get(0)).intValue();
					int mk = ((Integer) mp.get(1)).intValue();
					if(chr == ChrInt[i][0] && (ChrInt[i][1] == mk || mk == (ChrInt[i][1] + 1))) {
						m++;
					}
				}
			}
			matrix = new double[imp.IndividualNumber()][1 + selectedMarker.size() - m];
			for (int i = 0; i < imp.IndividualNumber(); i++) {
				matrix[i][0] = 1;
				int c = 1;
				for (int j = 0; j < selectedMarker.size(); j++) {
					ArrayList mp = (ArrayList) selectedMarker.get(j);
					int chr = ((Integer) mp.get(0)).intValue();
					int mk = ((Integer) mp.get(1)).intValue();
					boolean flag = true;
					for (int k = 0; k < SNPIdx.length; k++) {
						if(chr == ChrInt[k][0] && (mk == ChrInt[k][1] || mk == (ChrInt[k][1] + 1))) {
							flag = false;
							break;
						}
					}
					if(flag) {
						matrix[i][c++] = imp.MarkerAt(i, chr, mk);
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
				for (int j = 0; j < selectedMarker.size(); j++) {
					ArrayList mp = (ArrayList) selectedMarker.get(j);
					int chr = ((Integer) mp.get(0)).intValue();
					int mk = ((Integer) mp.get(1)).intValue();
					matrix[i][j+1] = imp.MarkerAt(i, chr, mk);
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
