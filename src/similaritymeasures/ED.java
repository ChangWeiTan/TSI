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
