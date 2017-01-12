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

import tools.Tools;

public class DTW {
	/* DTW
	 * This is a class for Dynamic Time Warping (DTW)
	 * 
	 * Last modified: 12/01/2017
	 */
	private static final Tools tools = new Tools();					// tools for the functions
	private static final int maxSeqLen = 3000;						// maximum sequence length
	private static double[][] D = new double[maxSeqLen][maxSeqLen];	// cost matrix
	private static long distComputation = 0;						// distance computations
	
	public final double compute(final double[] Q, final double[] C) {
		/* Computes standard DTW distance
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Compare time series
		 * Output:
		 * 	dtw	: DTW distance 
		 */
		final int nq = Q.length;	// length of query
		final int nc = C.length;	// length of compare
		int i, j;					// for loops counter
		double cost, temp;			// temporary variables
				
		D[0][0] = tools.distanceTo(Q[0], C[0]);
				
		for (i = 1; i < nq; i++) {
			D[i][0] = D[i-1][0] + tools.distanceTo(Q[i], C[0]);
		}
		for (j = 1; j < nc; j++) {
			D[0][j] = D[0][j-1] + tools.distanceTo(Q[0], C[j]);
		}
		
		for (i = 1; i < nq; i++) {
			for (j = 1; j < nc; j++) {
				cost = tools.distanceTo(Q[i], C[j]);
				temp = tools.Min3(D[i-1][j-1],D[i-1][j],D[i][j-1]);
				D[i][j] = cost + temp;
			}
		}
		distComputation = nq * nc;
		
		return D[nq-1][nc-1];
	}
	
	public final double compute(final double[] Q, final double[] C, final double[][] D) {
		/* Computes standard DTW distance
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Compare time series
		 * 	D	: Cost matrix
		 * Output:
		 * 	dtw	: DTW distance 
		 */
		final int nq = Q.length;	// length of query
		final int nc = C.length;	// length of compare
		int i, j;					// for loops counter
		double cost, temp;			// temporary variables
		
		D[0][0] = tools.distanceTo(Q[0], C[0]);
				
		for (i = 1; i < nq; i++) {
			D[i][0] = D[i-1][0] + tools.distanceTo(Q[i], C[0]);
		}
		for (j = 1; j < nc; j++) {
			D[0][j] = D[0][j-1] + tools.distanceTo(Q[0], C[j]);
		}
		
		for (i = 1; i < nq; i++) {
			for (j = 1; j < nc; j++) {
				cost = tools.distanceTo(Q[i], C[j]);
				temp = tools.Min3(D[i-1][j-1],D[i-1][j],D[i][j-1]);
				D[i][j] = cost + temp;
			}
		}
		distComputation = nq * nc;
		
		return D[nq-1][nc-1];
	}
	
	public final double compute(final double[] Q, final double[] C, final int W) {
		/* Computes standard DTW distance
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Compare time series
		 * 	W	: Warping window
		 * Output:
		 * 	dtw	: DTW distance 
		 */
		final int nq = Q.length;					// length of query
		final int nc = C.length;					// length of compare
		final int w = Math.max(W, Math.abs(nq-nc));	// warping window
		final int mq = Math.min(nq, 1+w);			// boundary of the cost matrix
		final int mc = Math.min(nc, 1+w);				
		int i, j, jStart, jStop;					// counters
		double cost, temp;							// temporary variables
				
		D[0][0] = tools.distanceTo(Q[0], C[0]);
		distComputation = 1;
		
		for (i = 1; i < mq; i++) {
			D[i][0] = D[i-1][0] + tools.distanceTo(Q[i], C[0]);
		}
		if (i < nq)
			D[i][0] = Double.POSITIVE_INFINITY;
		
		for (j = 1; j < mc; j++) {
			D[0][j] = D[0][j-1] + tools.distanceTo(Q[0], C[j]);
		}
		if (j < nc)
			D[0][j] = Double.POSITIVE_INFINITY;
		
		distComputation += (mq + mc); 
		
		for (i = 1; i < nq; i++) {
			jStart = Math.max(1, i-w);
			jStop = Math.min(nc, i+w+1);
			
			for (j = jStart; j < jStop; j++) {
				cost = tools.distanceTo(Q[i], C[j]);
				temp = tools.Min3(D[i - 1][j - 1],D[i - 1][j],D[i][j - 1]);
				D[i][j] = cost + temp;
				distComputation++;
			}
			if (jStop < nc)
				D[i][jStop] = Double.POSITIVE_INFINITY;	
			if (i < nq-1)
				D[i+1][jStart] = Double.POSITIVE_INFINITY;
		}
		
		return D[nq-1][nc-1];
	}
	
	public final double compute(final double[] Q, final double[] C, final int W, final double[][] D) {
		/* Computes
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Compare time series
		 * 	W	: Warping window
		 * 	D	: Cost matrix
		 * Output:
		 * 	dtw	: DTW distance 
		 */
		final int nq = Q.length;					// length of query
		final int nc = C.length;					// length of compare
		final int w = Math.max(W, Math.abs(nq-nc));	// warping window
		final int mq = Math.min(nq, 1+w);			// boundary of the cost matrix
		final int mc = Math.min(nc, 1+w);				
		int i, j, jStart, jStop;					// counters
		double cost, temp;							// temporary variables
				
		D[0][0] = tools.distanceTo(Q[0], C[0]);
		distComputation = 1;
		
		for (i = 1; i < mq; i++) {
			D[i][0] = D[i-1][0] + tools.distanceTo(Q[i], C[0]);
		}
		if (i < nq)
			D[i][0] = Double.POSITIVE_INFINITY;
		
		for (j = 1; j < mc; j++) {
			D[0][j] = D[0][j-1] + tools.distanceTo(Q[0], C[j]);
		}
		if (j < nc)
			D[0][j] = Double.POSITIVE_INFINITY;
		
		distComputation += (mq + mc); 
		
		for (i = 1; i < nq; i++) {
			jStart = Math.max(1, i-w);
			jStop = Math.min(nc, i+w+1);
			
			for (j = jStart; j < jStop; j++) {
				cost = tools.distanceTo(Q[i], C[j]);
				temp = tools.Min3(D[i - 1][j - 1],D[i - 1][j],D[i][j - 1]);
				D[i][j] = cost + temp;
				distComputation++;
			}
			if (jStop < nc)
				D[i][jStop] = Double.POSITIVE_INFINITY;	
			if (i < nq-1)
				D[i+1][jStart] = Double.POSITIVE_INFINITY;
		}
		
		return D[nq-1][nc-1];
	}
	
	public final long getDistComputation() {
		return distComputation;
	}
}
