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
import java.util.Arrays;
import pqueues.BranchPriorityQueue;
import pqueues.MaxPriorityQueue;
import pqueues.PriorityQueue;
import similaritymeasures.DistanceMeasures;
import tools.Tools;
import tree.TreeNode;

public class ClassifierTSI {
	/* Time Series Indexing (TSI) Classifier
	 * This is a class for Time Series Indexing (TSI) classifier
	 * 
	 * Last modified: 12/01/2017
	 */
	private final DistanceMeasures distComputer = new DistanceMeasures();	// computer for distance measures
	private final Tools tools = new Tools();								// tools for the functions
	private final int MaxPQSize = 10000;									// maximum size for priority queue
	private final int NIL = -100;											// constant for not available
	private int printInterval;												// interval for printing
	private int steps;														// steps to record error (steps=1, means record error for each training example)
	private int nbSteps;													// number of steps
	private int levels;														// tree level
	private long startTime, stopTime, saveTime;								// start and stop time per query
	private double totalQueryTime, elapsedTime;								// total query time, elapsed time = time per query
	
	private int branchFactor = 0;			// number of clusters
	private int timeseriesToExamine = 0;	// number of time series to examine
	private int kNearestNeighbour = 1;		// k nearest neighbour
	
	private ArrayList<double[]> TRAIN;		// training dataset
	private ArrayList<Integer> TrainClass;	// classes for the training dataset
	private int TrainSize;					// index for the training dataset
	private int ActualClass;				// size of training dataset
	private int ActualNNIndex;				// actual class of the query
	
	private double distComputation = 0;				
	private double distComputationPerQuery = 0;
	
	// search tree variables
	private BranchPriorityQueue dtwBranchPQ = null;
	private BranchPriorityQueue lbBranchPQ = null;
	
	private static int numTimeseriesSeen = 0;
	
	private int nbTimeseriesToActual = NIL;
	private double timeToActual = NIL;
	
	// outputs
	private MaxPriorityQueue kNNIndex = null;	// a maximum priority queue that store the results
	private int nnIndex = NIL;					// the nearest neighbor index
	private double averageQueryTime = NIL;		// average query time 
	private double precision = NIL;				// precision in searching for nearest neighbor
	private int[] index1NN;						// predicted nn index for each test time series
	private double precisionSoFar;				// precision so far for each test time series
	private double errorSoFar;					// error rate so far for each test time series
	
	// experiments parameters
	private double[] preError;
	private double[] preTime;
	private double[] preDist;
	private ArrayList<Double> preErrorAll = new ArrayList<Double>();
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
import pqueues.BranchPriorityQueue;
import pqueues.MaxPriorityQueue;
import pqueues.MinPriorityQueue;
import similaritymeasures.DistanceMeasures;
import tools.Tools;
import tree.TreeNode;

public class ClassifierTSI {
	/* Time Series Indexing (TSI) Classifier
	 * This is a class for Time Series Indexing (TSI) classifier
	 * 
	 * Last modified: 16/01/2017
	 */
	private final DistanceMeasures distComputer = new DistanceMeasures();	// computer for distance measures
	private final Tools tools = new Tools();		// tools for the functions
	private final int maxPQSize = 10000;			// maximum size for priority queue
	private final int nil = -100;					// constant for not available
	private int printInterval;						// interval for printing
	private int steps;								// steps to record error (steps=1, means record error for each training example)
	private int intervals;							// number of steps
	private int levels;								// tree level
	private long startTime, stopTime, saveTime;		// start and stop time per query
	private double totalQueryTime, elapsedTime;		// total query time, elapsed time = time per query
	
	private int branchFactor = 0;					// number of clusters
	private int nbTimeseriesToExamine = 0;			// number of time series to examine
	private int nearestNeighbour = 1;				// k nearest neighbour
	
	private ArrayList<double[]> trainingDataset;		// training dataset
	private ArrayList<Integer> trainingDatasetClass;	// classes for the training dataset
	private int trainingDatasetSize;					// index for the training dataset
	private int actualClass;							// size of training dataset
	private int actualNNIndex;							// actual class of the query
	
	private double distComputation = 0;				
	private double distComputationPerQuery = 0;
	
	// search tree variables
	private BranchPriorityQueue dtwBranchPQ = null;
	private BranchPriorityQueue lbBranchPQ = null;
	
	private int nbTimeseriesSeen = 0;
	
	private int nbTimeseriesToActualNN = nil;
	private double timeToActualNN = nil;
	
	// outputs
	private MaxPriorityQueue kNNIndex = null;		// a maximum priority queue that store the results
	private int nnIndex = nil;						// the nearest neighbor index
	private double averageQueryTime = nil;			// average query time 
	private double precision = nil;					// precision in searching for nearest neighbor
	private int[] index1NN;							// predicted nn index for each test time series
	private double precisionSoFar;					// precision so far for each test time series
	private double errorSoFar;						// error rate so far for each test time series
	
	// L experiments parameters
	private double[] preError;
	private double[] preTime;
	private double[] preDist;
	private ArrayList<Double> preErrorAll = new ArrayList<Double>();
	private ArrayList<Double> preTimeAll = new ArrayList<Double>();
	private ArrayList<Double> preDistAll = new ArrayList<Double>();
	private ArrayList<Integer> preLevels = new ArrayList<Integer>();
	private ArrayList<Integer> seenCounter = new ArrayList<Integer>();
	private double[] averagePrecisionPerQuery;		// average precision at different L for all the experiment runs
	private double[] averageErrorPerQuery;			// average error rate at different L
	private double[] averageDistPerQuery;			// average distance computations at different L for all the experiment runs
	private double[] averageTimePerQuery;			// average query time at different L for all the experiment runs
	private double[] precisionPerQuery;				// precision at different L for all the experiment runs
	private double[] errorPerQuery;					// error rate at different L for all the experiment runs
	private double[] timePerQuery;					// elapsed time at different L for all the experiment runs
	private double[] distPerQuery;					// distance computations at different L for all the experiment runs
	private int[] seenSoFarPerQuery;				// keep track of how many has seen so far for all the experiment runs
	private int arrayCount;							// counter for seen so far
	private int maxCount;
	
	public ClassifierTSI(final ArrayList<double[]> data, final ArrayList<Integer> dataClass, final int K, 
			final int k, final int L, final int step, final int stepSize, final int nbpre) {
		trainingDataset = data;
		trainingDatasetClass = dataClass;
		trainingDatasetSize = data.size();
		branchFactor = K;
		nbTimeseriesToExamine = L;
		nearestNeighbour = k;
		steps = step;
		intervals = stepSize;
		printInterval = (int) (0.1 * trainingDatasetSize);
		nbTimeseriesSeen = 0;
		distComputation = 0;
	}
	
	/*  Build Tree Algorithms */
	public final TreeNode buildTree(TreeNode parent, final ArrayList<double[]> data, final ArrayList<Integer> dataClass, final ArrayList<Integer> dataIndex, 
			final int Imax, final int w) {		
		if (data.size() <= branchFactor) {
			if (data.size() > 0) {
				parent.createLeaf(dataIndex, dataClass);
				nbTimeseriesSeen += data.size();
			}
		} else {
			// print status after certain intervals
			if (nbTimeseriesSeen%printInterval == 0 || nbTimeseriesSeen == trainingDatasetSize) {
				System.out.println("Completing: " + (double) nbTimeseriesSeen/ (double) trainingDatasetSize);
			}
			// initialize
			final int[] nbDataInClusters = new int[branchFactor]; 	// number of timeseries in each cluster
			final int[] clusterIndex = new int[data.size()];		// index of each cluster wrt dataset [1:N(1), N(1)+1:N(2), ... N(K-1):N(K)]
			int i, j;												// for loops counters
			
			final KMeansDTW kmeans = new KMeansDTW();
			final double[][] centroids = kmeans.compute(branchFactor, data, dataClass, dataIndex, w, Imax, nbDataInClusters, clusterIndex);
			
			int jStart = 0;				// counters to find the data in each cluster

			// go through all the clusters to build the tree
			for (i = 0; i < branchFactor; i++) {
				// initialize the loop
				final ArrayList<double[]> cluster = new ArrayList<double[]>(nbDataInClusters[i]);
				final ArrayList<Integer> clusterClass = new ArrayList<Integer>(nbDataInClusters[i]);
				final ArrayList<Integer> cIndex = new ArrayList<Integer>(nbDataInClusters[i]);
				
				// clustering around the centroids using the clusterIndex wrt to dataset found earlier
				// add them into arraylist
				int jEnd = jStart + nbDataInClusters[i]; 				
				for (j = jStart; j < jEnd; j++) {
					cluster.add(data.get(clusterIndex[j]));				// add all the timeseries in that cluster
					clusterClass.add(dataClass.get(clusterIndex[j]));	// class of each time series in that cluster
					cIndex.add(dataIndex.get(clusterIndex[j]));			// actual index of the time series wrt to original training dataset
				}
				jStart = jEnd;											// set the next start point
				
				// create a child, non leaf using the centroids 
				TreeNode childNode = new TreeNode();
				childNode.createNonLeaf(centroids[i], clusterClass, cIndex);
				
				// recursively build the tree
				childNode = buildTree(childNode, cluster, clusterClass, cIndex, Imax, w);
				
				// add child to the parent
				parent.addChild(childNode);
				
				// print status after certain intervals
				if (nbTimeseriesSeen%printInterval == 0 || nbTimeseriesSeen == trainingDatasetSize) {
					System.out.println("Completed: " + (double) nbTimeseriesSeen/ (double) trainingDatasetSize);
				}
			}
		}
		
		return parent;
	}
	
	/* Search Tree */	
	public final double performance(final TreeNode tree, final ArrayList<double[]> testingDataset, 
			final ArrayList<Integer> testingDatasetClass, final int w, final int[] nnDTWIndex) {
		int wrongCount = 0;						// counting the number of times predicted class is wrong
		int precisionCount = 0;					// counting the number of times predicted NN is correct
		double cumTimeToNN = 0, avgTimeToNN, avgNbTimeToNN;
		int cumNbToNN = 0;
		
		// initialize
		totalQueryTime = 0;		// total time to search
		index1NN = new int[testingDataset.size()];
		distComputation = 0;
		
		// initialize the results for different time intervals
		maxCount = 0;
		precisionPerQuery = new double[intervals];
		errorPerQuery = new double[intervals];
		timePerQuery = new double[intervals];
		distPerQuery = new double[intervals];
		seenSoFarPerQuery = new int[intervals];
		averagePrecisionPerQuery = new double[intervals];
		averageErrorPerQuery = new double[intervals];
		averageTimePerQuery = new double[intervals];
		averageDistPerQuery = new double[intervals];
		double[] tmpTime = new double[testingDataset.size()];
		double[] tmpError = new double[testingDataset.size()];
		double[] tmpDist = new double[testingDataset.size()];
		int[] counts = new int[testingDataset.size()];
		
		// go through the whole testing set
		for (int i = 0; i < testingDataset.size(); i++) {
			// results before seeing the first time series in the dataset
			preErrorAll = new ArrayList<Double>();
			preTimeAll = new ArrayList<Double>();
			preDistAll = new ArrayList<Double>();
			
			actualNNIndex = nnDTWIndex[i];		// correct nearest neighbour
			actualClass = testingDatasetClass.get(i);		// correct class
			nbTimeseriesToActualNN = nil;
			timeToActualNN = nil;
			final int predictClass = searchTree(tree, testingDataset.get(i), nbTimeseriesToExamine, w);	// predicted class of the query
			totalQueryTime += elapsedTime;		// update total query time
			final int predictNNIndex = getNN();	// predicted nearest neighbour index of the query
			index1NN[i] = predictNNIndex;		// store it away
			distComputation += distComputationPerQuery;
			
			// update error count
			if (predictClass != actualClass) {
				wrongCount++;
			}
			// update precision
			if (predictNNIndex == actualNNIndex) {
				precisionCount++;
			}
			
			// update so far
			precisionSoFar = (double) precisionCount/(i+1);
			errorSoFar = (double) wrongCount/(i+1);
			
			// get results before seeing the first time series in the dataset
			int count = 0;
			for (int j = 0; j < preTimeAll.size(); j++) {
				tmpError[count] += preErrorAll.get(j);
				tmpTime[count] += preTimeAll.get(j);
				tmpDist[count] += preDistAll.get(j);
				counts[count]++;
				count++;
			}
			maxCount = Math.max(maxCount, count);
			
			// print status after certain intervals
			if (i%printInterval == 0) {
				System.out.println(i + ", KNNIndex: " + predictNNIndex + "(" + actualNNIndex + "), Predicted: " +  
						predictClass + "(" + actualClass + "), Time: " + timeToActualNN + ", Nb: " + nbTimeseriesToActualNN + 
						", Error Rate: " + errorSoFar + ", Seen " + nbTimeseriesSeen + ", Precision: " + precisionSoFar);
			}
			cumTimeToNN += timeToActualNN;
			cumNbToNN += nbTimeseriesToActualNN;
		}
		
		// compute error rate, precision, number of distance computations & average query time per query at different time intervals
		for (int i = 0; i < intervals; i++) {
			averagePrecisionPerQuery[i] = (double) precisionPerQuery[i]/testingDataset.size();
			averageErrorPerQuery[i] = (double) errorPerQuery[i]/testingDataset.size();
			averageTimePerQuery[i] = (double) timePerQuery[i]/testingDataset.size();
			averageDistPerQuery[i] = (double) distPerQuery[i]/testingDataset.size();
		}
		
		// results before seeing the first time series in the dataset
		for (int i = 0; i < testingDataset.size(); i++) {
			if (counts[i] >= 0.5*testingDataset.size()) {
				tmpError[i] /= counts[i];
				tmpTime[i] /= counts[i];
				tmpDist[i] /= counts[i];
			} else {
				maxCount = i-1;
				break;
			}
		}
		preError = new double[maxCount];
		preTime = new double[maxCount];
		preDist = new double[maxCount];
		System.arraycopy(tmpError, 0, preError, 0, maxCount);
		System.arraycopy(tmpDist, 0, preDist, 0, maxCount);
		System.arraycopy(tmpTime, 0, preTime, 0, maxCount);
		
		// compute error rate, precision, number of distance computations & average query time
		final double errorRate = (double) wrongCount/testingDataset.size();
		averageQueryTime = (double) totalQueryTime/testingDataset.size();
		precision = (double) precisionCount/testingDataset.size();
		distComputation = (double) distComputation/testingDataset.size();
		
		avgTimeToNN = (double) cumTimeToNN / testingDataset.size();
		avgNbTimeToNN = (double) cumNbToNN / testingDataset.size();
		// print final results
		System.out.println("Error rate: " + errorRate + ", Avg Query Time: " + averageQueryTime + 
				", Avg Distance Computations: " + distComputation + ", Precision: " + precision);
		System.out.println("Average Time To Actual: " + avgTimeToNN + ", Average Nb To Actual: " + avgNbTimeToNN);
		return errorRate;
	}
	
	public final int searchTree(final TreeNode tree, final double[] query, final int L, final int w) {
		kNNIndex = new MaxPriorityQueue(nearestNeighbour);	// nn prirority queue (max PQ)
		dtwBranchPQ = new BranchPriorityQueue(maxPQSize);	// priority queue for DTW branches (min PQ)
		lbBranchPQ = new BranchPriorityQueue(maxPQSize);	// priority queue for LB branches (min PQ)
		distComputationPerQuery = 0;						// distance computation per query for experiment
		arrayCount = 0;										// counter for experimental results
		nbTimeseriesSeen = 0;								// number of time series seen
		
		startTime = System.nanoTime();						// classification start time 
		saveTime = 0;
		
		// compute envelope for query time series and use it to compare to the training dataset
		final double[][] wedge = distComputer.envelope(query, w);
		distComputationPerQuery += distComputer.getLBDistComputation();
		
		// traverse tree to leaf
		traverseTree(tree, query, wedge, w);
		while ((!dtwBranchPQ.isEmpty() || !lbBranchPQ.isEmpty()) && nbTimeseriesSeen < L) {
			// while priority queues are not empty and haven't seen L time series
			double minLB = Double.POSITIVE_INFINITY;	// minimum lower bound distance
			double minD = Double.POSITIVE_INFINITY;		// minimum dtw distance
			
			// update minimum lower-bound and dtw distance
			if (!lbBranchPQ.isEmpty()) {
				minLB = lbBranchPQ.firstDistance();
			}
			if (!dtwBranchPQ.isEmpty()) {
				minD = dtwBranchPQ.firstDistance();
			}
			
			while (minLB < minD) {
				// if minimum lower-bound distance is smaller than minimum dtw distance, compute dtw  
				final TreeNode topLBBranch = lbBranchPQ.pop();	// minimum branch from lower-bound PQ
				final double dtwDistance = distComputer.dtw(query, topLBBranch.getCentroid(), w);
				distComputationPerQuery += distComputer.getDTWDistComputation();
				
				dtwBranchPQ.insert(dtwDistance, topLBBranch);	// add that branch to dtw PQ 
				
				if (!lbBranchPQ.isEmpty()) {
					minLB = lbBranchPQ.firstDistance();
				} else {
					minLB = Double.POSITIVE_INFINITY;
				}
				if (!dtwBranchPQ.isEmpty()) {
					minD = dtwBranchPQ.firstDistance();
				}
			}
			
			if (!dtwBranchPQ.isEmpty()) {
				// if dtw queue is not empty, dequeue the first one as the most promising branch to traverse
				final TreeNode N = dtwBranchPQ.pop();
				traverseTree(N, query, wedge, w);
			}
		}
		
		final int predictClass = trainingDatasetClass.get(getNN());	// predicted class of the query
		
		// compute time per query
		stopTime = System.nanoTime();
		elapsedTime = (double) (stopTime - startTime)/1e9;
		
		return predictClass;
	}
		
	public final void traverseTree(final TreeNode tree, final double[] query, final double[][] wedge, final int w){	
		if (tree.isLeaf()){ 
			// if branch is leaf add to the nn priority queue
			final ArrayList<Integer> nodeIndex = tree.getDataIndex();	// data index of training dataset stored in the leaf 		
			
			// order the time series using a minimum priority queue based on lb distance
			final int[] index = new int[nodeIndex.size()];
			final double[] lbDistances = new double[nodeIndex.size()];	
			
			for (int j = 0; j < nodeIndex.size(); j++) {
				index[j] = nodeIndex.get(j);
				lbDistances[j] = distComputer.lb_keogh(wedge, trainingDataset.get(index[j])); 
				distComputationPerQuery += distComputer.getLBDistComputation();
			}
			tools.sort(lbDistances, index, 0, lbDistances.length-1);
			
			// apply LB_Keogh-NN-DTW with LB_Keogh to all the time series in this node
			for (int j = 0; j < nodeIndex.size(); j++) {
				// get the worst distance from result queue
				final double worstSoFarDist = kNNIndex.firstDistance();
								
				if (lbDistances[j] < worstSoFarDist) {
					final double dtwDistance = distComputer.dtw(query, trainingDataset.get(index[j]), w);
					distComputationPerQuery += distComputer.getDTWDistComputation();
					
					if (dtwDistance < worstSoFarDist) {
						// if dtw distance is better than the worst knn distance, add it to the results queue
						// if result queue is full, pop it and insert 
						if (kNNIndex.isFull())
							kNNIndex.pop();
						kNNIndex.insert(dtwDistance, index[j]);
					}
				}
				stopTime = System.nanoTime();
				nbTimeseriesSeen++;
				
				recordResults();
				saveTime = System.nanoTime() - stopTime;
				
				if (nbTimeseriesSeen == nbTimeseriesToExamine) {
					// stop searching if seen enough time series 
					break;
				}
			}
		} else {
			final ArrayList<TreeNode> children = tree.getChildren();	// list of childeren node
			final int numCluster = children.size();						// number of cluster 
			
			if (numCluster > 0) {
				// find the closest centroid and go there
				
				// order the centroids using a minimum priority queue based on lb distance
				final int[] index = new int[numCluster];
				final double[] lbDistances = new double[numCluster];
				for (int i = 0; i < numCluster; i++) {
					index[i] = i;
					lbDistances[i] = distComputer.lb_keogh(wedge, children.get(i).getCentroid()); 
					distComputation += distComputer.getLBDistComputation();
					distComputationPerQuery += distComputer.getLBDistComputation();
				}
				tools.sort(lbDistances, index, 0, lbDistances.length-1);
				
				double bestSoFar = Double.POSITIVE_INFINITY;		// best so far distance
				final double[] distances = new double[numCluster];	// vector for dtw distance
				final boolean[] dtwFlag = new boolean[numCluster];	// indicator for branches where dtw has been computed
				int nearestIndex = nil;								// nearest centroid index
				int bestClass = nil;
				
				// apply LB_Keogh-NN-DTW to all the centroids in this node
				for (int i = 0; i < numCluster; i++) {
					if (lbDistances[i] < bestSoFar) {
						// if lower-bound is smaller than best so far, compute dtw distance
						final double dtwDistance = distComputer.dtw(query, children.get(index[i]).getCentroid(), w);
						distComputationPerQuery += distComputer.getDTWDistComputation();
						distances[index[i]] = dtwDistance;
						dtwFlag[index[i]] = true;
						
						if (dtwDistance < bestSoFar) {
							// update nearest candidate if dtw is smaller than best so far
							bestSoFar = dtwDistance;
							nearestIndex = index[i];
							bestClass = children.get(index[i]).getCentroidLabel();
						}
						// record results before seeing the first time series
						recordPreResults(bestClass);
					} 
				}
				
				if (nbTimeseriesSeen < nbTimeseriesToExamine) {
					if (children.get(nearestIndex).leafSize() > 0) {
						// call traverse tree
						traverseTree(children.get(nearestIndex), query, wedge, w);
					} else {
						if (distances[nearestIndex] < kNNIndex.firstDistance()) {
							if (kNNIndex.isFull())
								kNNIndex.pop();
							kNNIndex.insert(distances[nearestIndex], children.get(nearestIndex).getDataIndex().get(0));
						}
						stopTime = System.nanoTime();
						nbTimeseriesSeen++;
						
						recordResults();
						saveTime = System.nanoTime() - stopTime;
					}
				}
				
				// get the nearest centroid and enqueue the rest to the priority queue
				for (int i = 0; i < numCluster; i++) {
					if (i != nearestIndex) {
						if (dtwFlag[i]) {
							if (children.get(i).leafSize() > 0) {
								dtwBranchPQ.insert(distances[i], children.get(i));
							} else {
								if (distances[i] < kNNIndex.firstDistance()) {
									if (kNNIndex.isFull())
										kNNIndex.pop();
									kNNIndex.insert(distances[i], children.get(i).getDataIndex().get(0));
								}
								stopTime = System.nanoTime();
								nbTimeseriesSeen++;
								recordResults();
								saveTime = System.nanoTime() - stopTime;
							}
						} else {
							if (children.get(i).leafSize() > 0) {
								lbBranchPQ.insert(distances[i], children.get(i));
							} else {
								if (distances[i] < kNNIndex.firstDistance()) {
									distances[i] = distComputer.dtw(query, children.get(i).getCentroid(), w);
									if (distances[i] < kNNIndex.firstDistance()) {
										if (kNNIndex.isFull())
											kNNIndex.pop();
										kNNIndex.insert(distances[i], children.get(i).getDataIndex().get(0));
									}
								}
								stopTime = System.nanoTime();
								nbTimeseriesSeen++;
								recordResults();
								saveTime = System.nanoTime() - stopTime;
							}
						}
						
						if (nbTimeseriesSeen == nbTimeseriesToExamine) {
							// stop searching if seen enough time series 
							break;
						}
					}
				}
			}
		}
	}
	
	private final void recordResults() {
		if (nbTimeseriesSeen == 1 || nbTimeseriesSeen == (trainingDatasetSize) || Math.floorMod(nbTimeseriesSeen, steps) == 0) {
			// record results at different time intervals
			elapsedTime = (double) (stopTime - startTime - saveTime)/1e9;
			if (trainingDatasetClass.get(kNNIndex.firstDataIndex()) != actualClass){
				errorPerQuery[arrayCount]++;
			}
			if (kNNIndex.firstDataIndex() == actualNNIndex) {
				precisionPerQuery[arrayCount]++;
				if (timeToActualNN == nil) {
					timeToActualNN = elapsedTime;
					nbTimeseriesToActualNN = nbTimeseriesSeen;
				}
			}
			
			timePerQuery[arrayCount] += elapsedTime;
			distPerQuery[arrayCount] += distComputationPerQuery;
			if (seenSoFarPerQuery[arrayCount] == 0) {
				seenSoFarPerQuery[arrayCount] = nbTimeseriesSeen;
			}
			arrayCount++;
		}
	}
	
	private final void recordPreResults(int bestClass) {
		if (nbTimeseriesSeen == 0) {
			stopTime = System.nanoTime();
			elapsedTime = (double) (stopTime - startTime)/1e9;
			
			seenCounter.add(1);
			preLevels.add(levels);
			if (bestClass != actualClass)
				preErrorAll.add(1.0);
			else
				preErrorAll.add(0.0);
			preTimeAll.add(elapsedTime);
			preDistAll.add(distComputationPerQuery);
		}
	}

	public final int getNN() {
		if (kNNIndex != null && !kNNIndex.isEmpty()) {
			int size = kNNIndex.sizeOf();
			MinPriorityQueue results = new MinPriorityQueue(size);
			for (int i = 0; i < size; i++) {
				kNNIndex.pop();
				results.insert(kNNIndex.popDistance(), kNNIndex.popDataIndex());
			}
			nnIndex = results.firstDataIndex();
		}
		return nnIndex;
	}
	
	public double getDistComputation() {
		return distComputation;
	}
		
	public double getAverageQueryTime() {
		return averageQueryTime;
	}
	
	public double getPrecision() {
		return precision;
	}
	
	public int[] getIndex1NN() {
		return index1NN;
	}
		
	public final double[] getAveragePrecisionPerQuery() {
		return averagePrecisionPerQuery;
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
	
	public final double[] getPreError() {
		return preError;
	}
	
	public final double[] getPreTime() {
		return preTime;
	}
	
	public final double[] getPreDist() {
		return preDist;
	}
}

			final int[] clusterIndex = new int[sizeDataset];		// index of each cluster wrt dataset [1:N(1), N(1)+1:N(2), ... N(K-1):N(K)]
			int i, j;												// for loops counters
			double[][] newCentroids;								// new centroids
			int iterations = 0;										// k-means iteration
			boolean converged = false;								// k-means convergence
			
			// initialize centroids using kmeans++
			double[][] centroids = kmeanspp(branchFactor, D, w, centroidIndex, clusterSeed); 	
			
			// do k-means clustering
			while (!converged && iterations < Imax) {
				// find the new centroids using the mean (DBA) of each cluster
				newCentroids = clusterMean(branchFactor, centroids, clusterSeed, D, w, nbDataInClusters, clusterIndex);
				
				// check if kmeans converges, when two centroids are equal
				converged = Arrays.deepEquals(centroids,newCentroids);
				
				// assign the new centroids to the old centroids
				if (!converged) {
					centroids = newCentroids;
				}
				
				// increment iteration 
				iterations++;
				
				// cluster around centroids
				if (!converged && iterations < Imax) {
					clusterAroundCentroid(D, centroids, w, clusterSeed);
				}
			}
			
			int jStart = 0, jEnd;				// counters to find the data in each cluster
			ArrayList<double[]> cluster;		// list of clusters
			ArrayList<Integer> clusterClass;	// class of data in cluster 
			ArrayList<Integer> cIndex;			// index data in cluster
			TreeNode childNode;					// new child
			
			// go through all the clusters to build the tree
			for (i = 0; i < branchFactor; i++) {
				// initialize the loop
				cluster = new ArrayList<double[]>(nbDataInClusters[i]);
				clusterClass = new ArrayList<Integer>(nbDataInClusters[i]);
				cIndex = new ArrayList<Integer>(nbDataInClusters[i]);
				
				// clustering around the centroids using the clusterIndex wrt to dataset found earlier
				// add them into arraylist
				jEnd = jStart + nbDataInClusters[i]; 					// set end point
				for (j = jStart; j < jEnd; j++) {
					cluster.add(D.get(clusterIndex[j]));				// add all the timeseries in that cluster
					clusterClass.add(dataClass.get(clusterIndex[j]));	// class of each time series in that cluster
					cIndex.add(dataIndex.get(clusterIndex[j]));			// actual index of the time series wrt to original training dataset
				}
				jStart = jEnd;											// set the next start point
				
				// create a child, non leaf using the centroids 
				childNode = new TreeNode();
				childNode.createNonLeaf(centroids[i], clusterClass, cIndex);
				
				// recursively build the tree
				childNode = buildTree(childNode, cluster, clusterClass, cIndex, Imax, w);
				
				// add child to the parent
				parent.addChild(childNode);
				
				// print status after certain intervals
				if (numTimeseriesSeen%printInterval == 0 || numTimeseriesSeen == TrainSize) {
					System.out.println("Completed: " + (double) numTimeseriesSeen/ (double) TrainSize);
				}
			}
		}
		
		return parent;
	}
	
	private final void clusterAroundCentroid(final ArrayList<double[]> D, final double[][] centroids, final int w, final int[] clusterSeed) {
		/* Cluster Around Centroid
		 * Clusters time series around the nearest cluster centroids
		 * It updates clusterSeed
		 * Inputs:
		 * 	D			: Training dataset
		 * 	centroids	: Cluster centroids
		 * 	w			: Warping window  
		 * 	clusterSeed	: Seed for data in each cluster 
		 */
		final int sizeDataset = D.size();			// size of the dataset
		int i, j;									// for loops counters
		double bestSoFar, lbDistance, dtwDistance;	// best so far, lower-bounds, dtw distances
		double[] C;									// candidate time series
		double[][] W;								// upper and lower envelope sequences for LB_Keogh
		
		// go through the whole training dataset
		for (i = 0; i < sizeDataset; i++) {
			C = D.get(i);
			
			// compute envelope for each time series in training set and use it to compare to the centroids
			W = distComputer.envelope(C, w);
			distComputation += distComputer.getLBDistComputation();
			
			// LB_Keogh-NN-DTW
			bestSoFar = Double.POSITIVE_INFINITY;
			for (j = 0; j < centroids.length; j++) {
				// compute lower-bound first
				lbDistance = distComputer.lb_keogh(W, centroids[j]);
				distComputation += distComputer.getLBDistComputation();
				if (lbDistance < bestSoFar) {
					// if lower-bound is smaller than best so far, compute dtw distance
					dtwDistance = distComputer.dtw(centroids[j], C, w);
					distComputation += distComputer.getDTWDistComputation();
					if (dtwDistance < bestSoFar) {
						// update nearest centroid if dtw is smaller than best so far
						bestSoFar = dtwDistance;
						clusterSeed[i] = j+1;
					}
				}
			}
		}
	}
		
	private final double[][] kmeanspp(final int K, final ArrayList<double[]> D, final int w, final int[] centroidIndex, final int[] clusterSeed) {
		/* K-Means++ Initialization
		 * Initialize k-means using k-means++
		 * Inputs:
		 * 	K 				: Branching factor, number of clusters
		 * 	D				: Training dataset
		 * 	w				: Warping window  
		 * 	centroidIndex	: Store the index of centroids in the dataset
		 * 	clusterSeed		: Seed for data in each cluster
		 * Outputs:
		 * 	centroid		: cluster centroids 
		 */
		final int sizeDataset = D.size();								// size of the dataset	
		final int lengthTimeseries = D.get(0).length;					// length of time series
		final double[][] centroids = new double[K][lengthTimeseries];	// centroids of each cluster
		final double[] dist2NearCentroid = new double[sizeDataset];		// distance to nearest centroids
		final boolean[] taken = new boolean[sizeDataset];				// indicator of the timeseries that is being assigned as centroid

		int cluster = 1;		// cluster label
		int centroidCount = 1;	// counting the number of centroids found so far
		int i, j;				// for loop counters
		int nextCentroidIndex;	// index of the next centroid
		double[][] W;			// upper and lower envelope sequences for LB_Keogh
		double[] distSum;		// cummulative distance, 
		double[] C; 			// compare time series
		double[] p; 			// temporary centroid
		double r, lbDistance, dtwDistance;	// r = random number, lower-bound and dtw distance

		// choose initial centroid at random
		final int firstCentroidIndex = tools.randInt(sizeDataset);			
		final double[] firstCentroid = D.get(firstCentroidIndex);
		centroids[centroidCount-1] = firstCentroid;
		centroidIndex[centroidCount-1] = firstCentroidIndex;
		clusterSeed[firstCentroidIndex] = cluster;
		taken[firstCentroidIndex] = true;

		// cluster the data around first centroid, store the distances
		for (i = 0; i < sizeDataset; i++) {
		    if (i != firstCentroidIndex) { 
			C = D.get(i);
			dist2NearCentroid[i] = distComputer.dtw(firstCentroid, C, w);
			distComputation += distComputer.getDTWDistComputation();
			clusterSeed[i] = cluster;	// assign that timeseries to the first cluster
		    }
		}

		// continue to find the remaining centroids
		while (centroidCount < K) {
			nextCentroidIndex = NIL; 			// index to the next centroid
			distSum = new double[sizeDataset];	// cumulative distance

			// calculating the cumulative sum of the distance to nearest centroid
			distSum[0] = dist2NearCentroid[0]; 		
		   	for (i = 1; i < sizeDataset; i++) {
				if (!taken[i]) {
					distSum[i] = distSum[i-1] + dist2NearCentroid[i];
				}
		   	}

	    	// pick a number in between 0 to cumulative distance randomly
	    	r = tools.randDouble() * distSum[sizeDataset-1];

	    	// find the next centroid index when sum >= r
	    	for (i = 0; i < sizeDataset; i++) {
				if (!taken[i] && distSum[i] >= r) {
					nextCentroidIndex = i;
					break;
				}
	    	}

		    // If it's not set to >= 0, the point wasn't found in the previous
		   	// for loop, probably because distances are extremely small. pick
		    // the last available point.
		    if (nextCentroidIndex == NIL) {
				for (i = sizeDataset - 1; i >= 0; i--) {
		    		if (!taken[i]) {
						nextCentroidIndex = i;
						break;
		    		}
				}
		    }

		    if (nextCentroidIndex != NIL) {
		    	// pick the next centroid from the dataset
		    	p = D.get(nextCentroidIndex);
				centroids[centroidCount] = p;
				centroidIndex[centroidCount] = nextCentroidIndex;
				clusterSeed[nextCentroidIndex] = ++cluster;
				taken[nextCentroidIndex] = true;

				if (centroidCount < K) {
				    // re-cluster the data around nearest centroids using NN-DTW with LB_Keogh
					// compute the envelop of the next centroid
					W = distComputer.envelope(p, w);
					distComputation += distComputer.getLBDistComputation();
					for (j = 0; j < sizeDataset; j++) {
						if (!taken[j]) {
							// compute lower bound distance to the candidates
							C = D.get(j);
							lbDistance = distComputer.lb_keogh(W, C);
							distComputation += distComputer.getLBDistComputation();
							if (lbDistance < dist2NearCentroid[j]){
								// if the lower bound distance is less than the current nearest distance, 
								// compute dtw distance
								dtwDistance = distComputer.dtw(p, C, w);
								distComputation += distComputer.getDTWDistComputation();
							    if (dtwDistance < dist2NearCentroid[j]) {
									// if it is the nearest, update the distance and cluster seed
									dist2NearCentroid[j] = dtwDistance;
									clusterSeed[j] = cluster;
							    } 
							}
						}
				    }
				}
				centroidCount++;
			} else {
				break;
		    }
		}

		return centroids;
	}
	
	private final double[][] clusterMean(final int K, final double[][] oldCentroids, final int[] clusterSeed, 
			final ArrayList<double[]> X, final int w, final int[] nbDataInClusters, final int[] clusterIndex) {
		/* Computes the mean of each cluster using DBA
		 * Inputs:
		 * 	oldCentroids		: Old centroids of the clusters
		 * 	clusterSeed			: Seed for data in each cluster
		 * 	X					: Training dataset
		 * 	w					: Warping window  
		 * 	nbDataInClusters	: Stores the number of data in each cluster
		 * 	clusterIndex		: Index of each data in each cluster
		 * Outputs:
		 * 	newCentroids		: New cluster centroids 
		 */
		final int lengthTimeseries = X.get(0).length;						// length of time series
		final double[][] newCentroids = new double[K][lengthTimeseries];	// the new centroids
		int dataCount = 0;													// iterates over the whole dataset
		int jStart = 0;														// determines the starting value when finding the clusters around centroids
		int i, j, clusterCount, clusterCount2;								// counters 
		double[][] cluster;													// list of clusters
		
		// go through all the centroids
		for (i = 0; i < K; i++) {
			newCentroids[i] = oldCentroids[i].clone(); 	// copy old to new
			clusterCount = 0;							// counting the number of data in the cluster
			
			// go through the dataset, count the number of data in the cluster and 
			// store the index of the data in the cluster
			for (j = 0; j < X.size(); j++) {
				if (clusterSeed[j] == (i+1)) {
					clusterIndex[dataCount++] = j;
					clusterCount++;
				}
			}
			
			// if there is data in that cluster
			if (clusterCount > 0) {
				// store the data count, for references later
				nbDataInClusters[i] = clusterCount;	
				
				clusterCount2 = 0;
				cluster = new double[clusterCount][lengthTimeseries];
				
				for (j = jStart; j < (jStart + clusterCount); j++) {
					cluster[clusterCount2++] = X.get(clusterIndex[j]);
				}
				
				jStart += clusterCount; // update starting index
				
				// compute mean using DBA
				newCentroids[i] = distComputer.dba(newCentroids[i], cluster, 10, w);
				distComputation += distComputer.getDBADistComputation();
			} else {
				nbDataInClusters[i] = 0; 
			}
		}
		
		return newCentroids;
	}
	
	/* Search Tree */	
	public final double performance(final TreeNode T, final ArrayList<double[]> TestSet, 
			final ArrayList<Integer> TestClass, final int w, final int[] nnDTWIndex) {
		/* Performance
		 * This is used to validate the performance of the classifier using a testing dataset.
		 * We record the error rate at different time interval
		 * Inputs: 
		 * 	T 			: Hierarchical K-Means tree
		 * 	TestSet		: Testing dataset
		 * 	TestClass	: Class for testing dataset
		 * 	w			: Warping window  
		 * 	nnDTWIndex	: Nearest neighbour index from LB_Keogh-NN-DTW 
		 * Output: 
		 * 	errorRate 	: Error rate of the validation
		 */
		final int TestSize = TestSet.size();	// size of testing dataset
		int wrongCount = 0;						// counting the number of times predicted class is wrong
		int precisionCount = 0;					// counting the number of times predicted NN is correct
		int PredictClass = NIL;					// predicted class for each query
		int PredictNNIndex = NIL;				// predicted nearest neighbour index
		int i;									// for loop counter
		double[] Q;								// query time series
		double cumTime = 0, avgTime, avgNb;
		int cumNb = 0;
		
		// initialize
		totalQueryTime = 0;		// total time to search
		index1NN = new int[TestSize];
		distComputation = 0;
		
		// initialize the results for different time intervals
		precisionPerQuery = new double[nbSteps];
		errorPerQuery = new double[nbSteps];
		timePerQuery = new double[nbSteps];
		distPerQuery = new double[nbSteps];
		seenSoFarPerQuery = new int[nbSteps];
		averagePrecisionPerQuery = new double[nbSteps];
		averageErrorPerQuery = new double[nbSteps];
		averageTimePerQuery = new double[nbSteps];
		averageDistPerQuery = new double[nbSteps];
		double[] tmpTime = new double[TestSize];
		double[] tmpError = new double[TestSize];
		double[] tmpDist = new double[TestSize];
		int[] counts = new int[TestSize];
		maxCount = 0;
		
		// go through the whole testing set
		for (i = 0; i < TestSize; i++) {
			// results before seeing the first time series in the dataset
			preErrorAll = new ArrayList<Double>();
			preTimeAll = new ArrayList<Double>();
			preDistAll = new ArrayList<Double>();
			
			ActualNNIndex = nnDTWIndex[i];		// correct nearest neighbour
			ActualClass = TestClass.get(i);		// correct class
			nbTimeseriesToActual = NIL;
			timeToActual = NIL;
			Q = TestSet.get(i);					// query time series
			PredictClass = searchTree(T, Q, timeseriesToExamine, w);	// predicted class of the query
			totalQueryTime += elapsedTime;		// update total query time
			PredictNNIndex = getNN();			// predicted nearest neighbour index of the query
			index1NN[i] = PredictNNIndex;		// store it away
			
			// update error count
			if (PredictClass != ActualClass) {
				wrongCount++;
			}
			// update precision
			if (PredictNNIndex == ActualNNIndex) {
				precisionCount++;
			}
			
			// update so far
			precisionSoFar = (double) precisionCount/(i+1);
			errorSoFar = (double) wrongCount/(i+1);
			
			// get results before seeing the first time series in the dataset
			int count = 0;
			for (int j = 0; j < preTimeAll.size(); j++) {
				tmpError[count] += preErrorAll.get(j);
				tmpTime[count] += preTimeAll.get(j);
				tmpDist[count] += preDistAll.get(j);
				counts[count]++;
				count++;
			}
			maxCount = Math.max(maxCount, count);
			
			// print status after certain intervals
			if (i%printInterval == 0) {
				System.out.println(i + ", KNNIndex: " + PredictNNIndex + "(" + ActualNNIndex + "), Predicted: " +  
						PredictClass + "(" + ActualClass + "), Time: " + timeToActual + ", Nb: " + nbTimeseriesToActual + 
						", Error Rate: " + errorSoFar + ", Seen " + numTimeseriesSeen + ", Precision: " + precisionSoFar);
			}
			cumTime += timeToActual;
			cumNb += nbTimeseriesToActual;
		}
		
		// compute error rate, precision, number of distance computations & average query time per query at different time intervals
		for (i = 0; i < nbSteps; i++) {
			averagePrecisionPerQuery[i] = (double) precisionPerQuery[i]/TestSize;
			averageErrorPerQuery[i] = (double) errorPerQuery[i]/TestSize;
			averageTimePerQuery[i] = (double) timePerQuery[i]/TestSize;
			averageDistPerQuery[i] = (double) distPerQuery[i]/TestSize;
		}
		
		// results before seeing the first time series in the dataset
		for (i = 0; i < TestSize; i++) {
			if (counts[i] >= 0.5*TestSize) {
				tmpError[i] /= counts[i];
				tmpTime[i] /= counts[i];
				tmpDist[i] /= counts[i];
			} else {
				maxCount = i-1;
				break;
			}
		}
		preError = new double[maxCount];
		preTime = new double[maxCount];
		preDist = new double[maxCount];
		System.arraycopy(tmpError, 0, preError, 0, maxCount);
		System.arraycopy(tmpDist, 0, preDist, 0, maxCount);
		System.arraycopy(tmpTime, 0, preTime, 0, maxCount);
		
		// compute error rate, precision, number of distance computations & average query time
		final double errorRate = (double) wrongCount/TestSize;
		averageQueryTime = (double) totalQueryTime/TestSize;
		precision = (double) precisionCount/TestSize;
		distComputation = (double) distComputation/TestSize;
		
		avgTime = (double) cumTime / TestSize;
		avgNb = (double) cumNb / TestSize;
		// print final results
		System.out.println("Error rate: " + errorRate + ", Avg Query Time: " + averageQueryTime + 
				", Avg Distance Computations: " + distComputation + ", Precision: " + precision);
		System.out.println("Average Time To Actual: " + avgTime + ", Average Nb To Actual: " + avgNb);
		return errorRate;
	}
	
	public final int searchTree(final TreeNode T, final double[] Q, final int L, final int w) {
		/* Search Tree 
		 * This searches the k-means tree and returns the nearest neighbour 
		 * We record the error rate at different time interval
		 * Inputs: 
		 * 	T 	: Hierarchical K-Means tree
		 * 	Q	: Query time series
		 * 	L	: Maximum number of time series to examine
		 * 	w	: Warping window  
		 * Output: 
		 * 	errorRate 	: Error rate of the validation
		 */
		kNNIndex = new MaxPriorityQueue(kNearestNeighbour);	// nn prirority queue (max PQ)
		dtwBranchPQ = new BranchPriorityQueue(MaxPQSize);	// priority queue for DTW branches (min PQ)
		lbBranchPQ = new BranchPriorityQueue(MaxPQSize);	// priority queue for LB branches (min PQ)
		distComputationPerQuery = 0;						// distance computation per query for experiment
		arrayCount = 0;										// counter for experimental results
		numTimeseriesSeen = 0;								// number of time series seen
		
		startTime = System.nanoTime();						// classification start time 
		saveTime = 0;
		
		// compute envelope for query time series and use it to compare to the training dataset
		final double[][] W = distComputer.envelope(Q, w);
		distComputation += distComputer.getLBDistComputation();
		distComputationPerQuery += distComputer.getLBDistComputation();
		
		// traverse tree to leaf
		traverseTree(T, Q, W, w);
		while ((!dtwBranchPQ.isEmpty() || !lbBranchPQ.isEmpty()) && numTimeseriesSeen < L) {
			// while priority queues are not empty and haven't seen L time series
			double minLB = Double.POSITIVE_INFINITY;	// minimum lower bound distance
			double minD = Double.POSITIVE_INFINITY;		// minimum dtw distance
			
			// update minimum lower-bound and dtw distance
			if (!lbBranchPQ.isEmpty()) {
				minLB = lbBranchPQ.firstDistance();
			}
			if (!dtwBranchPQ.isEmpty()) {
				minD = dtwBranchPQ.firstDistance();
			}
			
			while (minLB < minD) {
				// if minimum lower-bound distance is smaller than minimum dtw distance, compute dtw  
				lbBranchPQ.pop();
				final TreeNode topLBBranch = lbBranchPQ.popData();	// minimum branch from lower-bound PQ
				final double dtwDistance = distComputer.dtw(Q, topLBBranch.getCentroid(), w);
				distComputation += distComputer.getDTWDistComputation();
				distComputationPerQuery += distComputer.getDTWDistComputation();
				
				dtwBranchPQ.insert(dtwDistance, 0, topLBBranch);	// add that branch to dtw PQ 
				
				if (!lbBranchPQ.isEmpty()) {
					minLB = lbBranchPQ.firstDistance();
				} else {
					minLB = Double.POSITIVE_INFINITY;
				}
				if (!dtwBranchPQ.isEmpty()) {
					minD = dtwBranchPQ.firstDistance();
				}
			}
			
			if (!dtwBranchPQ.isEmpty()) {
				// if dtw queue is not empty, dequeue the first one as the most promising branch to traverse
				dtwBranchPQ.pop();
				final TreeNode N = dtwBranchPQ.popData();
				traverseTree(N, Q, W, w);
			}
		}
		
		final int predictClass = TrainClass.get(getNN());	// predicted class of the query
		
		// compute time per query
		stopTime = System.nanoTime();
		elapsedTime = (double) (stopTime - startTime)/1e9;
		
		return predictClass;
	}
		
	public final void traverseTree(final TreeNode N, final double[] Q, final double[][] W, final int w){
		/* Traverse tree 
		 * This searches the k-means tree and returns the nearest neighbour 
		 * We record the error rate at different time interval
		 * Inputs: 
		 * 	N 			: Hierarchical K-Means tree
		 * 	Q			: Query time series
		 * 	W			: Upper and lower envelope of the query
		 * 	w			: Warping window  
		 * Output: 
		 * 	errorRate 	: Error rate of the validation
		 */
		int bestClass = NIL;							// best so far class
		double bestSoFar = Double.POSITIVE_INFINITY;	// best so far distance 
		double dtwDistance, lbDistance;					// dtw and lower-bound distance
		double[] C;										// compare time series
		
		if (N.isLeaf()){ 
			// if branch is leaf add to the nn priority queue
			final ArrayList<Integer> nodeIndex = N.getDataIndex();	// data index of training dataset stored in the leaf 
			final int numTimeseries = nodeIndex.size();				// number of data in the leaf
			
			// order the time series using a minimum priority queue based on lb distance
			final int[] index = new int[numTimeseries];
			final double[] lbDistances = new double[numTimeseries];	
			
			for (int j = 0; j < numTimeseries; j++) {
				index[j] = nodeIndex.get(j);
				C = TRAIN.get(index[j]);
				lbDistances[j] = distComputer.lb_keogh(W, C); 
				distComputation += distComputer.getLBDistComputation();
				distComputationPerQuery += distComputer.getLBDistComputation();
			}
			tools.sort(lbDistances, index, 0, lbDistances.length-1);
			
			// apply LB_Keogh-NN-DTW with LB_Keogh to all the time series in this node
			for (int j = 0; j < numTimeseries; j++) {
				// get the worst distance from result queue
				final double worstSoFarDist = kNNIndex.firstDistance();
								
				if (lbDistances[j] < worstSoFarDist) {
					// if lb distance is better than the worst knn distance, compute dtw distance
					C = TRAIN.get(index[j]);
					dtwDistance = distComputer.dtw(Q, C, w);
					distComputation += distComputer.getDTWDistComputation();
					
					if (dtwDistance < worstSoFarDist) {
						// if dtw distance is better than the worst knn distance, add it to the results queue
						// if result queue is full, pop it and insert 
						if (kNNIndex.isFull())
							kNNIndex.pop();
						kNNIndex.insert(dtwDistance, index[j]);
					}
				}
				stopTime = System.nanoTime();
				numTimeseriesSeen++;
				
				recordResults();
				saveTime = System.nanoTime() - stopTime;
				
				if (numTimeseriesSeen == timeseriesToExamine) {
					// stop searching if seen enough time series 
					break;
				}
			}
		} else {
			// get all the children in this node
			final ArrayList<TreeNode> children = N.getChildren();	// list of childeren node
			final int numCluster = children.size();						// number of cluster 
			
			if (numCluster > 0) {
				// find the closest centroid and go there
				
				// order the centroids using a minimum priority queue based on lb distance
				final int[] index = new int[numCluster];
				final double[] lbDistances = new double[numCluster];
				for (int i = 0; i < numCluster; i++) {
					index[i] = i;
					C = children.get(i).getCentroid();
					lbDistances[i] = distComputer.lb_keogh(W, C); 
					distComputation += distComputer.getLBDistComputation();
					distComputationPerQuery += distComputer.getLBDistComputation();
				}
				tools.sort(lbDistances, index, 0, lbDistances.length-1);
				
				bestSoFar = Double.POSITIVE_INFINITY;				// best so far distance
				final double[] distances = new double[numCluster];	// vector for dtw distance
				final boolean[] dtwFlag = new boolean[numCluster];	// indicator for branches where dtw has been computed
				int nearestIndex = NIL;								// nearest centroid index
				int lbIndex;										// lower-bound index
				
				// apply LB_Keogh-NN-DTW to all the centroids in this node
				for (int i = 0; i < numCluster; i++) {
					lbDistance = lbDistances[i];
					lbIndex = index[i];
					distances[lbIndex] = lbDistance;
					
					if (lbDistance < bestSoFar) {
						// if lower-bound is smaller than best so far, compute dtw distance
						dtwDistance = distComputer.dtw(Q, children.get(lbIndex).getCentroid(), w);
						distComputation += distComputer.getDTWDistComputation();
						distComputationPerQuery += distComputer.getDTWDistComputation();
						distances[lbIndex] = dtwDistance;
						dtwFlag[lbIndex] = true;
						
						if (dtwDistance < bestSoFar) {
							// update nearest candidate if dtw is smaller than best so far
							bestSoFar = dtwDistance;
							nearestIndex = lbIndex;
							bestClass = children.get(lbIndex).getCentroidLabel();
						}
						// record results before seeing the first time series
						recordPre(bestClass);
					} 
				}
				
				if (numTimeseriesSeen < timeseriesToExamine) {
					if (children.get(nearestIndex).leafSize() > 0) {
						// call traverse tree
						traverseTree(children.get(nearestIndex), Q, W, w);
					} else {
						if (distances[nearestIndex] < kNNIndex.firstDistance()) {
							if (kNNIndex.isFull())
								kNNIndex.pop();
							kNNIndex.insert(distances[nearestIndex], children.get(nearestIndex).getDataIndex().get(0));
						}
						stopTime = System.nanoTime();
						numTimeseriesSeen++;
						
						recordResults();
						saveTime = System.nanoTime() - stopTime;
					}
				}
				
				// get the nearest centroid and enqueue the rest to the priority queue
				for (int i = 0; i < numCluster; i++) {
					if (i != nearestIndex) {
						if (dtwFlag[i]) {
							if (children.get(i).leafSize() > 0) {
								dtwBranchPQ.insert(distances[i], 0, children.get(i));
							} else {
								if (distances[i] < kNNIndex.firstDistance()) {
									if (kNNIndex.isFull())
										kNNIndex.pop();
									kNNIndex.insert(distances[i], children.get(i).getDataIndex().get(0));
								}
								stopTime = System.nanoTime();
								numTimeseriesSeen++;
								recordResults();
								saveTime = System.nanoTime() - stopTime;
							}
						} else {
							if (children.get(i).leafSize() > 0) {
								lbBranchPQ.insert(distances[i], 0, children.get(i));
							} else {
								if (distances[i] < kNNIndex.firstDistance()) {
									distances[i] = distComputer.dtw(Q, children.get(i).getCentroid(), w);
									if (distances[i] < kNNIndex.firstDistance()) {
										if (kNNIndex.isFull())
											kNNIndex.pop();
										kNNIndex.insert(distances[i], children.get(i).getDataIndex().get(0));
									}
								}
								stopTime = System.nanoTime();
								numTimeseriesSeen++;
								recordResults();
								saveTime = System.nanoTime() - stopTime;
							}
						}
						
						if (numTimeseriesSeen == timeseriesToExamine) {
							// stop searching if seen enough time series 
							break;
						}
					}
				}
			}
		}
	}
	
	private final void recordResults() {
		/* Record results
		 * This record results after seeing a candidate from the tree
		 * Input: 
		 * 	bestClass 	: Best so far class
		 */
		if (numTimeseriesSeen == 1 || numTimeseriesSeen == (TrainSize) || Math.floorMod(numTimeseriesSeen, steps) == 0) {
			// record results at different time intervals
			elapsedTime = (double) (stopTime - startTime - saveTime)/1e9;
			if (TrainClass.get(kNNIndex.firstDataIndex()) != ActualClass){
				errorPerQuery[arrayCount]++;
			}
			if (kNNIndex.firstDataIndex() == ActualNNIndex) {
				precisionPerQuery[arrayCount]++;
				if (timeToActual == NIL) {
					timeToActual = elapsedTime;
					nbTimeseriesToActual = numTimeseriesSeen;
				}
			}
			
			timePerQuery[arrayCount] += elapsedTime;
			distPerQuery[arrayCount] += distComputationPerQuery;
			if (seenSoFarPerQuery[arrayCount] == 0) {
				seenSoFarPerQuery[arrayCount] = numTimeseriesSeen;
			}
			arrayCount++;
		}
	}
	
	private final void recordPre(int bestClass) {
		/* Record pre results
		 * This record results before seeing first time series in the training dataset
		 * Input: 
		 * 	bestClass 	: Best so far class
		 */
		if (numTimeseriesSeen == 0) {
			stopTime = System.nanoTime();
			elapsedTime = (double) (stopTime - startTime)/1e9;
			
			seenCounter.add(1);
			preLevels.add(levels);
			if (bestClass != ActualClass)
				preErrorAll.add(1.0);
			else
				preErrorAll.add(0.0);
			preTimeAll.add(elapsedTime);
			preDistAll.add(distComputationPerQuery);
		}
	}

	public final int getNN() {
		/* Get Nearest Neighbor
		 * This returns the nearest neighbor index
		 * Output: 
		 * 	nnIndex	: nearest neighbour index
		 */
		if (kNNIndex != null && !kNNIndex.isEmpty()){
			int size = kNNIndex.sizeOf();
			PriorityQueue results = new PriorityQueue(size);
			for (int i = 0; i < size; i++) {
				kNNIndex.pop();
				results.insert(kNNIndex.popDistance(), kNNIndex.popDataIndex());
			}
			nnIndex = results.firstDataIndex();
		}
		return nnIndex;
	}
	
	public final double getDistComputation() {
		return distComputation;
	}
		
	public final double getAverageQueryTime() {
		return averageQueryTime;
	}
	
	public final double getPrecision() {
		return precision;
	}
	
	public final int[] getIndex1NN() {
		return index1NN;
	}
		
	public final double[] getAveragePrecisionPerQuery() {
		return averagePrecisionPerQuery;
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
	
	public final double[] getPreError() {
		return preError;
	}
	
	public final double[] getPreTime() {
		return preTime;
	}
	
	public final double[] getPreDist() {
		return preDist;
	}
}
