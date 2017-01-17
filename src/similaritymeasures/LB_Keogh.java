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
	 * 
	 * Last modified: 14/01/2017
	 */
	private final Tools tools = new Tools();	// tools for the functions
	private double[][] wedge;					// upper and lower envelope 
	private long distComputation = 0;			// distance computations
	
	public final double[][] envelope(final double[] query, final int w) {
		wedge = new double[2][query.length];		
		
		for (int i = 0; i < query.length; i++) {
			wedge[0][i] = Double.NEGATIVE_INFINITY;
			wedge[1][i] = Double.POSITIVE_INFINITY;
			final int jStart = (i - w < 0) ? 0 : i-w;
			final int jStop = (i + w + 1 > query.length) ? query.length : i+w+1;
			for (int j = jStart; j < jStop; j++) {
				wedge[0][i] = Math.max(wedge[0][i], query[j]);
				wedge[1][i] = Math.min(wedge[1][i], query[j]);
			}
			
		}
		distComputation = query.length;
		
		return wedge;
	}
	
	public final double[][] envelope(final ArrayList<double[]> candidates, final int w) {
		final int nc = candidates.get(0).length;	// length of query
		wedge = new double[2][nc];		
		
		for (int i = 0; i < candidates.size(); i++) {
			wedge[0][i] = Double.NEGATIVE_INFINITY;
			wedge[1][i] = Double.POSITIVE_INFINITY;
			final int jStart = (i - w < 0) ? 0 : i-w;
			final int jStop = (i + w + 1 > nc) ? nc : i+w+1;
			for (int j = jStart; j < jStop; j++) {
				for (int k = 0; k < candidates.size(); k++) {
					wedge[0][i] = Math.max(wedge[0][i], candidates.get(k)[j]);
					wedge[1][i] = Math.min(wedge[1][i], candidates.get(k)[j]);
				}
			}
			
		}
		distComputation = nc;
		
		return wedge;
	}
	
	public final double[][] envelope(final double[] query, final int w, final double[][] wedge) {		
		for (int i = 0; i < query.length; i++) {
			wedge[0][i] = Double.NEGATIVE_INFINITY;
			wedge[1][i] = Double.POSITIVE_INFINITY;
			final int jStart = (i - w < 0) ? 0 : i-w;
			final int jStop = (i + w + 1 > query.length) ? query.length : i+w+1;
			for (int j = jStart; j < jStop; j++) {
				wedge[0][i] = Math.max(wedge[0][i], query[j]);
				wedge[1][i] = Math.min(wedge[1][i], query[j]);
			}
			
		}
		distComputation = query.length;
		this.wedge = wedge;
		
		return wedge;
	}
	
	public final double[][] envelope(final ArrayList<double[]> candidate, final int w, final double[][] wedge) {
		final int nc = candidate.get(0).length;	// length of query
		
		for (int i = 0; i < nc; i++) {
			wedge[0][i] = Double.NEGATIVE_INFINITY;
			wedge[1][i] = Double.POSITIVE_INFINITY;
			final int jStart = (i - w < 0) ? 0 : i-w;
			final int jStop = (i + w + 1 > nc) ? nc : i+w+1;
			for (int j = jStart; j < jStop; j++) {
				for (int k = 0; k < candidate.size(); k++) {
					wedge[0][i] = Math.max(wedge[0][i], candidate.get(k)[j]);
					wedge[1][i] = Math.min(wedge[1][i], candidate.get(k)[j]);
				}
			}
			
		}
		distComputation = nc;
		this.wedge = wedge;
		
		return wedge;
	}
	
	public final double compute(final double[][] wedge, final double[] candidate) {
		final int minLength = Math.min(wedge[0].length, candidate.length);
		double res = 0;
		
		for (int i = 0; i < minLength; i++) {
			if (candidate[i] < wedge[1][i]) {
				res += tools.squaredEuclidean(wedge[1][i], candidate[i]);
			} else if (wedge[0][i] < candidate[i]) {
				res += tools.squaredEuclidean(wedge[0][i], candidate[i]);
			}
		}
		distComputation = minLength;
		
		return res;
	}
	
	public final double compute(final double[] candidate) {
		final int minLength = Math.min(wedge[1].length, candidate.length);
		double res = 0;
		
		for (int i = 0; i < minLength; i++) {
			if (candidate[i] < wedge[1][i]) {
				res += tools.squaredEuclidean(wedge[1][i], candidate[i]);
			} else if (wedge[0][i] < candidate[i]) {
				res += tools.squaredEuclidean(wedge[0][i], candidate[i]);
			}
		}
		distComputation = minLength;
		
		return res;
	}
	
	public final double computeEA(final double[][] wedge, final double[] query, final double r) {
		final int minLength = Math.min(wedge[1].length, query.length);
		double res = 0;
		
		for (int i = 0; i < minLength; i++) {
			if (query[i] < wedge[1][i]) {
				res += tools.squaredEuclidean(wedge[1][i], query[i]);
			} else if (wedge[0][i] < query[i]) {
				res += tools.squaredEuclidean(wedge[0][i], query[i]);
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
