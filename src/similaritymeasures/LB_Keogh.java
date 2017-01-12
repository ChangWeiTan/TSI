/*******************************************************************************
 * Copyright (C) 2017 Chang Wei Tan
 * 
 * This file is part of TSI.
 * 
 * TSI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * 
 * TSI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with TSI.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
package similaritymeasures;

import java.util.ArrayList;

import tools.Tools;

public class LB_Keogh {
	/* LB_Keogh
	 * This is a class for Lower-Bound Keogh
	 * Refer to http://www.cs.ucr.edu/~eamonn/LB_Keogh.htm for more information
	 * 
	 * Last modified: 12/01/2017
	 */
	private final Tools tools = new Tools();	// tools for the functions
	private double[][] W;						// upper and lower envelope 
	private long distComputation = 0;			// distance computations
	
	public final double[][] envelope(final double[] Q, final double[] C, final int w) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Candidate time series
		 * 	w	: Warping window
		 * Output:
		 * 	W	: Upper and lower envelope 
		 */
		final int nq = Q.length;				// length of query
		final int nc = C.length;				// length of compare
		final int minLength = Math.min(nq,nc);	// minimum length of the two
		W = new double[2][minLength];			// initialize envelope
		int jStart, jStop, i, j;				// counters
		
		for (i = 0; i < minLength; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > minLength) ? minLength : i+w+1;
			for (j = jStart; j < jStop; j++) {
				W[0][i] = Math.max(W[0][i], Q[j]);
				W[1][i] = Math.min(W[1][i], Q[j]);
			}
		}
		distComputation = minLength;
		
		return W;
	}
	
	public final double[][] envelope(final double[] Q, final ArrayList<double[]> C, final int w) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Candidate time series
		 * 	w	: Warping window
		 * Output:
		 * 	W	: Upper and lower envelope 
		 */
		final int nbCandidate = C.size();
		final int nq = Q.length;				// length of query
		final int nc = C.get(0).length;			// length of candidate
		final int minLength = Math.min(nq,nc);	// minimum length of the two
		W = new double[2][minLength];			// initialize envelope
		int jStart, jStop, i, j, k;				// counters
		
		for (i = 0; i < minLength; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > minLength) ? minLength : i+w+1;
			for (j = jStart; j < jStop; j++) {
				for (k = 0; k < nbCandidate; k++) {
					W[0][i] = Math.max(W[0][i], C.get(k)[j]);
					W[1][i] = Math.min(W[1][i], C.get(k)[j]);
				}
			}
		}
		distComputation = minLength;
		
		return W;
	}
	
	public final double[][] envelope(final double[] Q, final double[] C, final int w, final double[][] Wedge) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 		: Query time series
		 * 	C		: Candidate time series
		 * 	w		: Warping window
		 * 	Wedge	: Preset upper and lower envelope 
		 * Output:
		 * 	W		: Upper and lower envelope
		 */
		final int nq = Q.length;				// length of query
		final int nc = C.length;				// length of compare
		final int minLength = Math.min(nq,nc);	// minimum length of the two
		int jStart, jStop, i, j;				// counters 
		
		for (i = 0; i < minLength; i++) {
			Wedge[0][i] = Double.NEGATIVE_INFINITY;
			Wedge[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > minLength) ? minLength : i+w+1;
			for (j = jStart; j < jStop; j++) {
				Wedge[0][i] = Math.max(Wedge[0][i], Q[j]);
				Wedge[1][i] = Math.min(Wedge[1][i], Q[j]);
			}
		}
		distComputation = minLength;
		W = Wedge;
		
		return W;
	}
	
	public final double[][] envelope(final double[] Q, final ArrayList<double[]> C, final int w, final double[][] W) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Candidate time series
		 * 	w	: Warping window
		 *  W	: Preset upper and lower envelope 
		 * Output:
		 * 	W	: Upper and lower envelope 
		 */
		final int nbCandidate = C.size();
		final int nq = Q.length;				// length of query
		final int nc = C.get(0).length;			// length of candidate
		final int minLength = Math.min(nq,nc);	// minimum length of the two
		int jStart, jStop, i, j, k;				// counters
		
		for (i = 0; i < minLength; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > minLength) ? minLength : i+w+1;
			for (j = jStart; j < jStop; j++) {
				for (k = 0; k < nbCandidate; k++) {
					W[0][i] = Math.max(W[0][i], C.get(k)[j]);
					W[1][i] = Math.min(W[1][i], C.get(k)[j]);
				}
			}
		}
		distComputation = minLength;
		
		return W;
	}
	
	public final double[][] envelope(final double[] Q, final int w) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	w	: Warping window
		 * Output:
		 * 	W	: Upper and lower envelope
		 */
		final int nq = Q.length;	// length of query
		W = new double[2][nq];		// initialize envelope
		int jStart, jStop, i, j;	// counters
		
		for (i = 0; i < nq; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > nq) ? nq : i+w+1;
			for (j = jStart; j < jStop; j++) {
				W[0][i] = Math.max(W[0][i], Q[j]);
				W[1][i] = Math.min(W[1][i], Q[j]);
			}
			
		}
		distComputation = nq;
		
		return W;
	}
	
	public final double[][] envelope(final ArrayList<double[]> C, final int w) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	w	: Warping window
		 * Output:
		 * 	W	: Upper and lower envelope
		 */
		final int nbCandidate = C.size();
		final int nc = C.get(0).length;	// length of query
		W = new double[2][nc];		// initialize envelope
		int jStart, jStop, i, j, k;	// counters
		
		for (i = 0; i < nc; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > nc) ? nc : i+w+1;
			for (j = jStart; j < jStop; j++) {
				for (k = 0; k < nbCandidate; k++) {
					W[0][i] = Math.max(W[0][i], C.get(k)[j]);
					W[1][i] = Math.min(W[1][i], C.get(k)[j]);
				}
			}
			
		}
		distComputation = nc;
		
		return W;
	}
	
	public final double[][] envelope(final double[] Q, final int w, final double[][] Wedge) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 		: Query time series
		 * 	w		: Warping window
		 * 	Wedge	: Preset upper and lower envelope
		 * Output:
		 * 	W		: Upper and lower envelope
		 */
		final int nq = Q.length;	// length of query
		int jStart, jStop, i, j;	// counters
		
		for (i = 0; i < nq; i++) {
			Wedge[0][i] = Double.NEGATIVE_INFINITY;
			Wedge[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > nq) ? nq : i+w+1;
			for (j = jStart; j < jStop; j++) {
				Wedge[0][i] = Math.max(Wedge[0][i], Q[j]);
				Wedge[1][i] = Math.min(Wedge[1][i], Q[j]);
			}
			
		}
		distComputation = nq;
		W = Wedge;
		
		return W;
	}
	
	public final double[][] envelope(final ArrayList<double[]> C, final int w, final double[][] W) {
		/* Computes the upper and lower envelope 
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	w	: Warping window
		 * 	W	: Preset upper and lower envelope
		 * Output:
		 * 	W	: Upper and lower envelope
		 */
		final int nbCandidate = C.size();
		final int nc = C.get(0).length;	// length of query
		int jStart, jStop, i, j, k;		// counters
		
		for (i = 0; i < nc; i++) {
			W[0][i] = Double.NEGATIVE_INFINITY;
			W[1][i] = Double.POSITIVE_INFINITY;
			jStart = (i - w < 0) ? 0 : i-w;
			jStop = (i + w + 1 > nc) ? nc : i+w+1;
			for (j = jStart; j < jStop; j++) {
				for (k = 0; k < nbCandidate; k++) {
					W[0][i] = Math.max(W[0][i], C.get(k)[j]);
					W[1][i] = Math.min(W[1][i], C.get(k)[j]);
				}
			}
			
		}
		distComputation = nc;
		
		return W;
	}
	
	public final double compute(final double[][] W, final double[] C) {
		/* Compute
		 * This sub-routine computes Lower Bound Keogh
		 * Inputs:
		 * 	W	: Upper and lower envelope
		 * 	C	: Compare time series
		 * Output:
		 * 	lb	: Lower-bound distance 
		 */
		final int minLength = Math.min(W[1].length, C.length);
		double res = 0;
		int i;
		
		for (i = 0; i < minLength; i++) {
			if (C[i] < W[1][i]) {
				res += tools.distanceTo(W[1][i], C[i]);
			} else if (W[0][i] < C[i]) {
				res += tools.distanceTo(W[0][i], C[i]);
			}
		}
		distComputation = minLength;
		
		return res;
	}
	
	public final double compute(final double[] C) {
		/* Compute
		 * This sub-routine computes Lower Bound Keogh
		 * Inputs:
		 * 	C	: Candidate time series
		 * Output:
		 * 	lb	: Lower-bound distance 
		 */
		final int minLength = Math.min(W[1].length, C.length);
		double res = 0;
		int i;
		
		for (i = 0; i < minLength; i++) {
			if (C[i] < W[1][i]) {
				res += tools.distanceTo(W[1][i], C[i]);
			} else if (W[0][i] < C[i]) {
				res += tools.distanceTo(W[0][i], C[i]);
			}
		}
		distComputation = minLength;
		
		return res;
	}
	
	public final double computeEA(final double[][] W, final double[] Q, final double r) {
		/* Computes Lower Bound Keogh with Early Abandon
		 * Inputs:
		 * 	W	: Upper and lower envelope
		 * 	C	: Candidate time series
		 * Output:
		 * 	lb	: Lower-bound distance 
		 */
		final int minLength = Math.min(W[1].length, Q.length);
		double res = 0;
		int i;
		
		for (i = 0; i < minLength; i++) {
			if (Q[i] < W[1][i]) {
				res += tools.distanceTo(W[1][i], Q[i]);
			} else if (W[0][i] < Q[i]) {
				res += tools.distanceTo(W[0][i], Q[i]);
			}
			if (res > r) {
				distComputation = i;
				return Double.POSITIVE_INFINITY;
			}
		}
		distComputation = minLength;
		
		return res;
	}
	
	public final long getDistComputation() {
		return distComputation;
	}
}
