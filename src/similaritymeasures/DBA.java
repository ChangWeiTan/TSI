/*******************************************************************************
 * Copyright (C) 2017 Francois Petitjean, Chang Wei Tan
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

public class DBA {
	/* DBA
	 * This is a class DTW Barycenter Averaging (DBA)
	 * Original Code from Francois Petitjean, https://github.com/fpetitjean/DBA
	 * 
	 * Last modified: 12/01/2017
	 */
	private final static Tools tools = new Tools();	// tools for the functions
	private final static int NIL = -1;				// constant for not available
	private final static int DIAGONAL = 0;			
	private final static int LEFT = 1;
	private final static int UP = 2;
	private final static int MAX_SEQ_LENGTH = 3000;
	private final static double[][] costMatrix = new double[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	private final static int[][] pathMatrix = new int[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	private final static int[][] optimalPathLength = new int[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	private static long distComputation = 0;
	
	public final static double[] update(double[] C, final ArrayList<double[]> sequences, final int w) {
		/* Updates the average sequence with a warping window
		 * Inputs:
		 * 	C			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	w			: warping window
		 * Outputs:
		 * 	C			: average sequence
		 */
		@SuppressWarnings("unchecked")
		final ArrayList<Double>[] tupleAssociation = new ArrayList[C.length];
		for (int i = 0; i < tupleAssociation.length; i++) {
			tupleAssociation[i] = new ArrayList<Double>(sequences.size());
		}
		double res = 0.0;
		int nbTuplesAverageSeq, i, j, indiceRes;
		int centerLength = C.length;
		int seqLength, jStart, jStop;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = tools.distanceTo(C[0], T[0]);
			pathMatrix[0][0] = NIL;
			optimalPathLength[0][0] = 0;
			distComputation = 1;
			
			for (i = 1; i < Math.min(centerLength, 1+w); i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + tools.distanceTo(C[i], T[0]);
				pathMatrix[i][0] = UP;
				optimalPathLength[i][0] = i;
			}
			if (i < centerLength)
				costMatrix[i][0] = Double.POSITIVE_INFINITY;
			
			for (j = 1; j < Math.min(seqLength, 1+w); j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + tools.distanceTo(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
				optimalPathLength[0][j] = j;
			}
			if (j < seqLength)
				costMatrix[0][j] = Double.POSITIVE_INFINITY;

			distComputation += (Math.min(centerLength, 1+w) + Math.min(seqLength, 1+w));
			
			for (i = 1; i < centerLength; i++) {
				jStart = Math.max(1, i-w);
				jStop = Math.min(seqLength, i+w+1);
				
				for (j = jStart; j < jStop; j++) {
					indiceRes = tools.ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + tools.distanceTo(C[i], T[j]);
					distComputation++;
				}
				if (jStop < seqLength)
					costMatrix[i][jStop] = Double.POSITIVE_INFINITY;
				if (i < centerLength - 1)
					costMatrix[i+1][jStart] = Double.POSITIVE_INFINITY;
			}

			nbTuplesAverageSeq = optimalPathLength[centerLength-1][seqLength-1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;
			
			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].add(T[j]);
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						break;
					case LEFT:
						j = j - 1;
						break;
					case UP:
						i = i - 1;
						break;
				}
			}
		}
		
		for (int t = 0; t < centerLength; t++) {
			C[t] = barycenter((tupleAssociation[t].toArray()));
		}
		return C;
	}
	
	public final static double[] update(double[] C, final double[][] sequences, final int w) {
		/* Updates the average sequence with a warping window
		 * Inputs:
		 * 	C			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	w			: warping window
		 * Outputs:
		 * 	C			: average sequence
		 */
		@SuppressWarnings("unchecked")
		final ArrayList<Double>[] tupleAssociation = new ArrayList[C.length];
		for (int i = 0; i < tupleAssociation.length; i++) {
			tupleAssociation[i] = new ArrayList<Double>(sequences.length);
		}
		double res = 0.0;
		int nbTuplesAverageSeq, i, j, indiceRes;
		int centerLength = C.length;
		int seqLength, jStart, jStop;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = tools.distanceTo(C[0], T[0]);
			pathMatrix[0][0] = NIL;
			optimalPathLength[0][0] = 0;
			distComputation = 1;
			
			for (i = 1; i < Math.min(centerLength, 1+w); i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + tools.distanceTo(C[i], T[0]);
				pathMatrix[i][0] = UP;
				optimalPathLength[i][0] = i;
			}
			if (i < centerLength)
				costMatrix[i][0] = Double.POSITIVE_INFINITY;
			
			for (j = 1; j < Math.min(seqLength, 1+w); j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + tools.distanceTo(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
				optimalPathLength[0][j] = j;
			}
			if (j < seqLength)
				costMatrix[0][j] = Double.POSITIVE_INFINITY;

			distComputation += (Math.min(centerLength, 1+w) + Math.min(seqLength, 1+w));
			
			for (i = 1; i < centerLength; i++) {
				jStart = Math.max(1, i-w);
				jStop = Math.min(seqLength, i+w+1);
				
				for (j = jStart; j < jStop; j++) {
					indiceRes = tools.ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + tools.distanceTo(C[i], T[j]);
					distComputation++;
				}
				if (jStop < seqLength)
					costMatrix[i][jStop] = Double.POSITIVE_INFINITY;
				if (i < centerLength - 1)
					costMatrix[i+1][jStart] = Double.POSITIVE_INFINITY;
			}

			nbTuplesAverageSeq = optimalPathLength[centerLength-1][seqLength-1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;
			
			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].add(T[j]);
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						break;
					case LEFT:
						j = j - 1;
						break;
					case UP:
						i = i - 1;
						break;
				}
			}
		}
		
		for (int t = 0; t < centerLength; t++) {
			C[t] = barycenter((tupleAssociation[t].toArray()));
		}
		return C;
	}
	
	public final static double[] update(double[] C, final ArrayList<double[]> sequences) {
		/* Updates the average sequence
		 * Inputs:
		 * 	C			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * Outputs:
		 * 	C			: average sequence
		 */
		@SuppressWarnings("unchecked")
		final ArrayList<Double>[] tupleAssociation = new ArrayList[C.length];
		for (int i = 0; i < tupleAssociation.length; i++) {
			tupleAssociation[i] = new ArrayList<Double>(sequences.size());
		}
		int nbTuplesAverageSeq, i, j, indiceRes;
		double res = 0.0;
		int centerLength = C.length;
		int seqLength;
		
		distComputation = 0;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = tools.distanceTo(C[0], T[0]);
			pathMatrix[0][0] = NIL;
			optimalPathLength[0][0] = 0;
			distComputation++;

			for (i = 1; i < centerLength; i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + tools.distanceTo(C[i], T[0]);
				pathMatrix[i][0] = UP;
				optimalPathLength[i][0] = i;
				distComputation++;
			}
			for (j = 1; j < seqLength; j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + tools.distanceTo(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
				optimalPathLength[0][j] = j;
				distComputation++;
			}

			for (i = 1; i < centerLength; i++) {
				for (j = 1; j < seqLength; j++) {
					indiceRes = tools.ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + tools.distanceTo(C[i], T[j]);
					distComputation++;
				}
			}

			nbTuplesAverageSeq = optimalPathLength[centerLength - 1][seqLength - 1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;

			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].add(T[j]);
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						break;
					case LEFT:
						j = j - 1;
						break;
					case UP:
						i = i - 1;
						break;
				}

			}

		}

		for (int t = 0; t < centerLength; t++) {
			C[t] = barycenter((tupleAssociation[t].toArray()));
		}
		return C;

	}
	
	public final static double[] update(double[] C, final double[][] sequences) {
		/* Updates the average sequence
		 * Inputs:
		 * 	C			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * Outputs:
		 * 	C			: average sequence
		 */
		@SuppressWarnings("unchecked")
		final ArrayList<Double>[] tupleAssociation = new ArrayList[C.length];
		for (int i = 0; i < tupleAssociation.length; i++) {
			tupleAssociation[i] = new ArrayList<Double>(sequences.length);
		}
		int nbTuplesAverageSeq, i, j, indiceRes;
		double res = 0.0;
		int centerLength = C.length;
		int seqLength;
		
		distComputation = 0;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = tools.distanceTo(C[0], T[0]);
			pathMatrix[0][0] = NIL;
			optimalPathLength[0][0] = 0;
			distComputation++;

			for (i = 1; i < centerLength; i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + tools.distanceTo(C[i], T[0]);
				pathMatrix[i][0] = UP;
				optimalPathLength[i][0] = i;
				distComputation++;
			}
			for (j = 1; j < seqLength; j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + tools.distanceTo(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
				optimalPathLength[0][j] = j;
				distComputation++;
			}

			for (i = 1; i < centerLength; i++) {
				for (j = 1; j < seqLength; j++) {
					indiceRes = tools.ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + tools.distanceTo(C[i], T[j]);
					distComputation++;
				}
			}

			nbTuplesAverageSeq = optimalPathLength[centerLength - 1][seqLength - 1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;

			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].add(T[j]);
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						break;
					case LEFT:
						j = j - 1;
						break;
					case UP:
						i = i - 1;
						break;
				}

			}

		}

		for (int t = 0; t < centerLength; t++) {
			C[t] = barycenter((tupleAssociation[t].toArray()));
		}
		return C;
	}
	
	public final double[] compute(double[] T, final double[][] sequences, final int Imax, final int w) {
		/* Compute
		 * Inputs:
		 * 	T			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * 	w			: warping window
		 * Outputs:
		 * 	T			: average sequence
		 */
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences, w);
		}
		
		return T;
	}
	
	public final double[] compute(double[] T, final ArrayList<double[]> sequences, final int Imax, final int w) {
		/* Compute
		 * Inputs:
		 * 	T			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * 	w			: warping window
		 * Outputs:
		 * 	T			: average sequence
		 */
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences, w);
		}
		
		return T;
	}
	
	public final double[] compute(double[] T, final double[][] sequences, final int Imax) {
		/* Compute
		 * Inputs:
		 * 	T			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * Outputs:
		 * 	T			: average sequence
		 */
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences);
		}
		
		return T;
	}
	
	public final double[] compute(double[] T, final ArrayList<double[]> sequences, final int Imax) {
		/* Compute
		 * Inputs:
		 * 	T			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * Outputs:
		 * 	T			: average sequence
		 */
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences);
		}
		
		return T;
	}
	
	public final double[] compute(final double[][] sequences, final int Imax, final int w) {
		/* Compute
		 * Inputs:
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * 	w			: warping window
		 * Outputs:
		 * 	T			: average sequence
		 */
		double randIndex = Math.random() * sequences.length;
		double[] T = sequences[(int) randIndex];
		
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences, w);
		}
		
		return T;
	}
	
	public final double[] compute(final ArrayList<double[]> sequences, final int Imax, final int w) {
		/* Compute
		 * Inputs:
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * 	w			: warping window
		 * Outputs:
		 * 	T			: average sequence
		 */
		double randIndex = Math.random() * sequences.size();
		double[] T = sequences.get((int) randIndex);
		
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences, w);
		}
		
		return T;
	}
	
	public final double[] compute(final double[][] sequences, final int Imax){
		/* Compute
		 * Inputs:
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * Outputs:
		 * 	T			: average sequence
		 */
		double randIndex = Math.random() * sequences.length;
		double[] T = sequences[(int) randIndex];
		
		for (int i=0; i < Imax; i++)
		{
			T = update(T, sequences);
		}
		
		return T;
	}
	
	public final double[] compute(final ArrayList<double[]> sequences, final int Imax) {
		/* Compute
		 * Inputs:
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * Outputs:
		 * 	T			: average sequence
		 */
		double randIndex = Math.random() * sequences.size();
		double[] T = sequences.get((int) randIndex);
		
		for (int i=0; i < Imax; i++) {
			T = update(T, sequences);
		}
		
		return T;
	}
	
	private final static double barycenter(final Object... tab) {
		/* Barycenter
		 * Inputs:
		 * 	T			: average sequence
		 * 	sequences	: set of time series to be averaged
		 * 	Imax		: maximum iteration
		 * Outputs:
		 * 	T			: average sequence
		 */
		if (tab.length < 1) {
			throw new RuntimeException("empty double tab");
		}
		double sum = 0.0;
		sum = 0.0;
		for (Object o : tab) {
			sum += ((Double) o);
		}
		return sum / tab.length;
	}
	
	/* Get number of distance computation */
	public final long getnbDistComputation() {
		return distComputation;
	}
}
