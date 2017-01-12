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

public class ED {
	/* ED
	 * This is a class for Euclidean Distance
	 * 
	 * Last modified: 12/01/2017
	 */
	private static long distComputation = 0;
	
	public final double compute(final double[] Q, final double[] C) {
		/* Computes Euclidean distance
		 * Inputs: 
		 * 	Q 	: Query time series
		 * 	C	: Compare time series
		 * Output:
		 * 	d	: Euclidean distance 
		 */
		double distance;
		double temp = 0;
		int i;
		int sequenceLength = Q.length;
		
		distComputation = 0;
		
		for (i = 0; i < sequenceLength; i++) {
			distance = (Q[i]-C[i]);
			distance = distance * distance;
			temp += distance;
			distComputation++;
		}
		return Math.sqrt(temp);
	}
	
	public final long getDistComputation() {
		return distComputation;
	}

}
