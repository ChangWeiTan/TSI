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
package classifiers;

import java.util.ArrayList;
import similaritymeasures.DistanceMeasures;

public class ClassifierNNED {
	/* NNED Classifier
	 * 
	 * Last modified: 01/10/2016
	 */
	private final DistanceMeasures distComputer = new DistanceMeasures();	// computer for distance measures
	private final int nil = -100;					// constant for not available								
	private int printInterval;						// interval for printing 
	private int steps;								// steps to record error (steps=1, means record error for each training example)
	private int intervals;							// time intervals
	private long startTime, stopTime, saveTime;		// start stop and save time per query
	private double totalQueryTime, elapsedTimePerQuery;		// total query time, elapsed time = time per query
	
	private ArrayList<double[]> trainingDataset;		// training dataset
	private ArrayList<Integer> trainingDatasetClass;	// classes for the training dataset
	private ArrayList<Integer> trainingDatasetIndex;	// index for the training dataset
	private int trainingDatasetSize;					// size of training dataset
	private int actualClass;							// actual class of the query
	
	private double errorRate;						// error rate 
	private double distComputation = 0;				// distance computations
	private double averageQueryTime = 0;				// average query time
	private double[] distance1NN;							// distance recorded for each query
	private int[] index1NN;									// nearest neighbour index for each query
	
	private int bsfIndex;									// best so far nearest neighbour index
	private double bsfDist;									// best so far distance
	
	private double[] averageErrorPerQuery;					// average error per query
	private double[] averageDistPerQuery;					// average distance computations per query
	private double[] averageTimePerQuery;					// average time per query
	
	private double[] errorPerQuery;							// error per query at different time intervals
	private double[] timePerQuery;							// time per query at different time intervals
	private double[] distPerQuery;							// distance computations per query at different time intervals
	private int[] seenSoFarPerQuery;						// number of time series seen so far
	
	public ClassifierNNED(){}
	
	public ClassifierNNED(final ArrayList<double[]> D, final ArrayList<Integer> dataClass, final ArrayList<Integer> index, final int step, final int stepSize) {
		trainingDataset = D;
		trainingDatasetClass = dataClass;
		trainingDatasetSize = D.size();
		trainingDatasetIndex = index;
		steps = step;
		intervals = stepSize;
		printInterval = (int) 1;
		distComputation = 0;
	}
		
	public final double performance(final ArrayList<double[]> testSet, final ArrayList<Integer> testClass) {
		int errorCount = 0;						
		
		// initialize
		totalQueryTime = 0;
		distComputation = 0;
		distance1NN = new double[testSet.size()];
		index1NN = new int[testSet.size()];
		
		// initialize the results for different time intervals
		seenSoFarPerQuery = new int[intervals];
		errorPerQuery = new double[intervals];
		timePerQuery = new double[intervals];
		distPerQuery = new double[intervals];
		averageErrorPerQuery = new double[intervals];
		averageTimePerQuery = new double[intervals];
		averageDistPerQuery = new double[intervals];
		
		for (int i = 0; i < testSet.size(); i++){
			// print status after certain intervals
			if (i % printInterval == 0) {
				System.out.print("Classifying Test " + (i+1) + ", ");
			}
			
			actualClass = testClass.get(i);							// actual class of the query
			final int predictClass = classify(testSet.get(i), i); 	// predicted class of the query
			totalQueryTime += elapsedTimePerQuery;					// update total query time
			
			// update error count
			if (predictClass != actualClass) {
				errorCount++;
			}
			
			// update error rate so far & store the best so far distance & nearest neighbour index
			final double errorSoFar = (double) errorCount/(i+1);
			distance1NN[i] = bsfDist;
			index1NN[i] = bsfIndex;
			
			// print status after certain intervals
			if (i % printInterval == 0) {
				System.out.println("took " + elapsedTimePerQuery + "s, NN Index: " + bsfIndex + ", Error so far: " + errorSoFar);
			}
		}
		
		// compute error rate, number of distance computations & average query time per query at different time intervals
		for (int i = 0; i < intervals; i++) {
			averageErrorPerQuery[i] = (double) errorPerQuery[i]/testSet.size();
			averageTimePerQuery[i] = (double) timePerQuery[i]/testSet.size();
			averageDistPerQuery[i] = (double) distPerQuery[i]/testSet.size();
		}
		
		// compute error rate, number of distance computations & average query time
		errorRate = (double) errorCount/testSet.size();
		distComputation = (double) distComputation/testSet.size();
		averageQueryTime = (double) totalQueryTime/testSet.size();
		
		return errorRate;
	}
	
	private final int classify(final double[] query, final int testNum) {
		int predictClass = nil;
		int count = 0; 					
		long distComputationPerQuery = 0;		

		// initialize
		bsfIndex = nil;
		bsfDist = Double.POSITIVE_INFINITY;
		elapsedTimePerQuery = 0;
		startTime = System.nanoTime();
		saveTime = 0;
		
		for (int i = 0; i < trainingDatasetSize; i++) {
			final double edDistance = distComputer.ed(query, trainingDataset.get(i));
			distComputationPerQuery += distComputer.getLBDistComputation();
			
			if (edDistance < bsfDist) {
				predictClass = trainingDatasetClass.get(i);
				bsfIndex = trainingDatasetIndex.get(i);
				bsfDist = edDistance;
			}
			stopTime = System.nanoTime();
			
			if (i == 0 || i == (trainingDatasetSize-1) || Math.floorMod(i+1, steps) == 0) {
				// record error rate at fixed time intervals/number of time series seen
				elapsedTimePerQuery = (double) (stopTime - startTime - saveTime)/1e9;
				
				if (predictClass != actualClass) {
					errorPerQuery[count]++;
				}
				
				timePerQuery[count] += elapsedTimePerQuery;
				distPerQuery[count] += distComputationPerQuery;
				if (testNum == 0) {
					seenSoFarPerQuery[count] = i+1;
				}
				count++;
			}
			saveTime = System.nanoTime() - stopTime;
		}
		
		// compute time per query
		stopTime = System.nanoTime();
		elapsedTimePerQuery = (double) (stopTime - startTime - saveTime)/1e9;
		distComputation += distComputationPerQuery;
		
		return predictClass;
	}
	
	public double getErrorRate() {
		return errorRate;
	}
	
	public double getDistComputation() {
		return distComputation;
	}
		
	public double getAverageQueryTime() {
		return averageQueryTime;
	}
	
	public int[] getIndex1NN() {
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
