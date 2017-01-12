package classifiers;

import java.util.ArrayList;
import similaritymeasures.DistanceMeasures;

public class ClassifierNNDTW {
	/* LB_Keogh-NN-DTW Classifier
	 * This is a class for LB_Keogh-NN-DTW classifier
	 * 
	 * Last modified: 12/01/2017
	 */
	private final DistanceMeasures distComputer = new DistanceMeasures();	// computer for distance measures
	private final int NIL = -100;											// constant for not available								
	private int printInterval;												// interval for printing 
	private int steps;														// steps to record error (steps=1, means record error for each training example)
	private int nbSteps;													// number of steps
	private long startTime, stopTime, saveTime;								// start, stop, save time per query
	private double totalQueryTime, elapsedTime;								// total query time, elapsed time = time per query
	
	private ArrayList<double[]> TRAIN;		// training dataset
	private ArrayList<Integer> TrainClass;	// classes for the training dataset
	private ArrayList<Integer> TrainIndex;	// index for the training dataset
	private int TrainSize;					// size of training dataset
	private int ActualClass;				// actual class of the query
	private int predictClass;
	
	private double errorRate;				// error rate 
	private double distComputation = 0;		// distance computations
	private double averageQueryTime = 0;	// average query time
	private double[] distance1NN;			// distance to NN for each query
	private int[] index1NN;					// nearest neighbour for each query
	
	private int bsfIndex;					// best so far nearest neighbour index
	private double bsfDist;					// best so far distance
	
	private double[] averageErrorPerQuery;	// average error per query
	private double[] averageDistPerQuery;	// average distance computations per query
	private double[] averageTimePerQuery;	// average time per query
	
	private double[] errorPerQuery;			// error per query at different time intervals
	private double[] timePerQuery;			// time per query at different time intervals
	private double[] distPerQuery;			// distance computations per query at different time intervals
	private int[] seenSoFarPerQuery;		// number of time series seen so far
	
	public ClassifierNNDTW(){}
	
	public ClassifierNNDTW(final ArrayList<double[]> D, final ArrayList<Integer> dataClass, final ArrayList<Integer> index, final int step, final int stepSize) {
		TRAIN = D;
		TrainClass = dataClass;
		TrainSize = D.size();
		TrainIndex = index;
		steps = step;
		nbSteps = stepSize;
		printInterval = (int) 1;
		distComputation = 0;
	}
	
	/* Methods */	
	public final double performance(final ArrayList<double[]> testSet, final ArrayList<Integer> testClass, final int w) {
		/* Performance
		 * This is used to validate the performance of the classifier using a testing dataset.
		 * We record the error rate at different time interval
		 * Inputs: 
		 * 	testSet		: Testing dataset
		 * 	testClass	: Class for testing dataset
		 * 	w			: Warping window  
		 * Output: 
		 * 	errorRate 	: Error rate of the validation
		 */
		final int testSize = testClass.size();	// size of testing dataset
		int predictClass;						// predict class of the query
		int errorCount = 0;						// counter for errors
		int i;									// for loop counter
		double errorSoFar;						// error rate so far
		
		errorRate = NIL;
		totalQueryTime = 0;
		distComputation = 0;
		distance1NN = new double[testSize];
		index1NN = new int[testSize];
		
		// initialize the results for different time intervals
		seenSoFarPerQuery = new int[nbSteps];
		errorPerQuery = new double[nbSteps];
		timePerQuery = new double[nbSteps];
		distPerQuery = new double[nbSteps];
		averageErrorPerQuery = new double[nbSteps];
		averageTimePerQuery = new double[nbSteps];
		averageDistPerQuery = new double[nbSteps];
		
		for (i = 0; i < testSize; i++){
			// print status after certain intervals
			if (i % printInterval == 0) {
				System.out.print("Classifying Test " + (i+1) + ", ");
			}
			
			ActualClass = testClass.get(i);					// actual class of the query
			predictClass = classify(testSet.get(i), w, i); 	// predicted class of the query
			totalQueryTime += elapsedTime;					// update total query time
			
			// update error count
			if (predictClass != ActualClass) {
				errorCount++;
			}
			
			// update error rate so far & store the best so far distance & nearest neighbour index
			errorSoFar = (double) errorCount/(i+1);
			distance1NN[i] = bsfDist;
			index1NN[i] = bsfIndex;
			
			// print status after certain intervals
			if (i % printInterval == 0) {
				System.out.println("took " + elapsedTime + "s, NN Index: " + bsfIndex + ", Error so far: " + errorSoFar);
			}
		}
		
		// compute error rate, number of distance computations & average query time per query at different time intervals
		for (i = 0; i < nbSteps; i++) {
			averageErrorPerQuery[i] = (double) errorPerQuery[i]/testSize;
			averageTimePerQuery[i] = (double) timePerQuery[i]/testSize;
			averageDistPerQuery[i] = (double) distPerQuery[i]/testSize;
		}
		
		// compute error rate, number of distance computations & average query time
		errorRate = (double) errorCount/testSize;
		distComputation = (double) distComputation/testSize;
		averageQueryTime = (double) totalQueryTime/testSize;
		
		return errorRate;
	}
		
	private final int classify(final double[] Q, final int w, final int testNum) {
		/* Classify
		 * This is used to classify a query time series using LB_Keogh-NN-DTW 
		 * We record the error rate at different time interval
		 * Inputs: 
		 * 	Q				: Query time series
		 * 	w				: Warping window  
		 * 	testNum			: Query number, use for recording the error rates, only record when we are classifying for the first time
		 * Output: 
		 * 	predictClass 	: Predicted class for that query
		 */
		int i;				// for loop counter
		int count = 0; 		// counter for storing the results at different time intervals
		long n = 0;			// distance computations for the experiment
		double lbDistance, dtwDistance;	// dtw and lb distances
		double[] C;			// candidate time series
		double[][] W;		// upper and lower envelope sequences for LB_Keogh

		predictClass = NIL;
		bsfIndex = NIL;
		bsfDist = Double.POSITIVE_INFINITY;
		elapsedTime = 0;
		startTime = System.nanoTime();
		saveTime = 0;
		
		// compute envelope for query time series and use it to compare to the training dataset
		W = distComputer.envelope(Q, w);
		distComputation += distComputer.getLBDistComputation();
		n += distComputer.getLBDistComputation();
		
		// go through the whole training dataset
		for (i = 0; i < TrainSize; i++) {
			C = TRAIN.get(i);
			// compute lower-bound distance
			lbDistance = distComputer.lb_keogh(W, C);
			distComputation += distComputer.getLBDistComputation();
			n += distComputer.getLBDistComputation();
			
			if (lbDistance < bsfDist) {
				// if lower-bound is smaller than best so far, compute dtw distance
				dtwDistance = distComputer.dtw(Q, C, w);
				distComputation += distComputer.getDTWDistComputation();
				n += distComputer.getDTWDistComputation();
				
				if (dtwDistance < bsfDist) {
					// update nearest candidate if dtw is smaller than best so far
					bsfIndex = TrainIndex.get(i);;
					predictClass = TrainClass.get(i);
					bsfDist = dtwDistance;
				}
			}
			stopTime = System.nanoTime();
			
			if (i == 0 || i == (TrainSize-1) || Math.floorMod(i+1, steps) == 0) {
				// record error rate at fixed time intervals/number of time series seen
				elapsedTime = (double) (stopTime - startTime - saveTime)/1e9;
					
				if (predictClass != ActualClass) {
					errorPerQuery[count]++;
				}
					
				timePerQuery[count] += elapsedTime;
				distPerQuery[count] += n;
				if (testNum == 0) {
					seenSoFarPerQuery[count] = i+1;
				}
				count++;
			}
			saveTime = System.nanoTime() - stopTime;
		}
		
		// compute time per query
		stopTime = System.nanoTime();
		elapsedTime = (double) (stopTime - startTime - saveTime)/1e9;
		
		return predictClass;
	}
	
	public final double getErrorRate() {
		return errorRate;
	}
	
	public final double getDistComputation() {
		return distComputation;
	}
		
	public final double getAverageQueryTime() {
		return averageQueryTime;
	}
	
	public final int[] getIndex1NN() {
		return index1NN;
	}
	
	public final double[] getAverageErrorPerQuery(){
		return averageErrorPerQuery;
	}
	
	public final double[] getAverageTimePerQuery() {
		return averageTimePerQuery;
	}
	
	public final double[] getAverageDistPerQuery() {
		return averageDistPerQuery;
	}
	
	public final int[] getSeenSoFarPerQuery() {
		return seenSoFarPerQuery;
	}
}
