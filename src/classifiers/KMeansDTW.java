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

import similaritymeasures.DistanceMeasures;
import tools.Tools;

public class KMeansDTW {
	private final DistanceMeasures distComputer = new DistanceMeasures();
	private final Tools tools = new Tools();
	private final int nil = -100;
	private int[] clusterSeed;
	
	public final double[][] compute(final int K, final ArrayList<double[]> data, final ArrayList<Integer> dataClass, final ArrayList<Integer> dataIndex, 
			final int w, final int Imax, final int[] nbDataInClusters, final int[] clusterIndex) {
		final int[] centroidIndex = new int[K];
		final int[] clusterSeed = new int[data.size()];
		
		// initialize kmeans using kmeans++
		double[][] centroids = kmeanspp(K, data, w, centroidIndex, clusterSeed);
		
		for (int i = 0; i < Imax; i++) {
			System.out.println("Iterations " + (i+1));
			// recalculate new center
			double[][] newCentroids = clusterMean(K, centroids, clusterSeed, data, w, nbDataInClusters, clusterIndex);
			
			if (Arrays.deepEquals(centroids,newCentroids)) {
				break;
			} else {
				centroids = newCentroids;
				clusterAroundCentroid(data, centroids, w, clusterSeed);
			}
		}
		
		this.clusterSeed = clusterSeed;
		
		return centroids;
	}
	
	private final double[][] kmeanspp(final int K, final ArrayList<double[]> data, final int w, final int[] centroidIndex, final int[] clusterSeed) {
        final double[][] centroids = new double[K][data.get(0).length];	// centroids of each cluster
        final double[] dist2NearCentroid = new double[data.size()];		// distance to nearest centroids
        final boolean[] taken = new boolean[data.size()];				// indicator of the timeseries that is being assigned as centroid
        
        int cluster = 0;					// cluster label
        int centroidCount = 1;				// counting the number of centroids found so far
        int nextCentroidIndex;				// index of the next centroid
        double[] distSum;					// cummulative distance
        
        // choose initial centroid at random
        final int firstCentroidIndex = tools.randInt(data.size());			
        final double[] firstCentroid = data.get(firstCentroidIndex);
        centroids[centroidCount-1] = firstCentroid;
        centroidIndex[centroidCount-1] = firstCentroidIndex;
        clusterSeed[firstCentroidIndex] = cluster;
        taken[firstCentroidIndex] = true;
        
        // cluster the data around first centroid, store the distances
        for (int i = 0; i < data.size(); i++) {
            if (i != firstCentroidIndex) { 
            	dist2NearCentroid[i] = distComputer.dtw(firstCentroid, data.get(i), w);
                clusterSeed[i] = cluster;	// assign that timeseries to the first cluster
            }
        }
        
        // continue to find the remaining centroids
        while (centroidCount < K) {
        	nextCentroidIndex = nil; 				// index to the next centroid
        	distSum = new double[data.size()];		// cumulative distance
        	        	
        	// calculating the cumulative sum of the distance to nearest centroid
        	distSum[0] = dist2NearCentroid[0]; 		
            for (int i = 1; i < data.size(); i++) {
                if (!taken[i]) {
                	distSum[i] = distSum[i-1] + dist2NearCentroid[i];
                }
            }
            
            // pick a number in between 0 to cumulative distance randomly
            final double r = tools.randDouble() * distSum[data.size()-1];
            
            // find the next centroid index when sum >= r
            for (int i = 0; i < data.size(); i++) {
                if (!taken[i]) {
                    if (distSum[i] >= r) {
                    	nextCentroidIndex = i;
                        break;
                    }
                }
            }
            
            // If it's not set to >= 0, the point wasn't found in the previous
            // for loop, probably because distances are extremely small. pick
            // the last available point.
            if (nextCentroidIndex == nil) {
                for (int i = data.size() - 1; i >= 0; i--) {
                    if (!taken[i]) {
                    	nextCentroidIndex = i;
                        break;
                    }
                }
            }
            
            if (nextCentroidIndex != nil) {
            	// pick the next centroid from the dataset
                final double[] nextCentroid = data.get(nextCentroidIndex);
                centroids[centroidCount] = nextCentroid;
                centroidIndex[centroidCount] = nextCentroidIndex;
                clusterSeed[nextCentroidIndex] = ++cluster;
                taken[nextCentroidIndex] = true;
                                
                if (centroidCount < K) {
                	final double[][] wedge = distComputer.envelope(nextCentroid, w);
                	for (int j = 0; j < data.size(); j++) {
                        if (!taken[j]) {
                        	final double lbDistance = distComputer.lb_keogh(wedge, data.get(j));
                        	if (lbDistance < dist2NearCentroid[j]){
	                        	final double dtwDistance = distComputer.dtw(nextCentroid, data.get(j), w);
	                            if (dtwDistance < dist2NearCentroid[j]) {
	                            	dist2NearCentroid[j] = dtwDistance;
	                                clusterSeed[j] = cluster;
	                            } 
                        	}
                        }
                    }
                }
                centroidCount++;
            } 
            else {
                break;
            }
        }
        
        return centroids;
	}
	
	private final double[][] clusterMean(final int K, final double[][] oldCentroids, final int[] clusterSeed, 
			final ArrayList<double[]> data, final int w, final int[] nbDataInClusters, final int[] clusterIndex) {
		final double[][] newCentroids = new double[K][data.get(0).length];	
		int dataCount = 0;	// iterates over the whole dataset
		int jStart = 0;		// determines the starting value when finding the clusters around centroids
		
		// go through all the centroids
		for (int i = 0; i < K; i++) {
			int clusterSize = 0;						
			
			for (int j = 0; j < data.size(); j++) {
				if (clusterSeed[j] == i) {
					clusterIndex[dataCount++] = j;
					clusterSize++;
				}
			}
			
			// if there is data in that cluster
			if (clusterSize > 0) {
				// store the data count, for references later
				nbDataInClusters[i] = clusterSize;	
				
				int clusterCount2 = 0;
				final double[][] cluster = new double[clusterSize][data.get(0).length];
				
				for (int j = jStart; j < (jStart + clusterSize); j++) {
					cluster[clusterCount2++] = data.get(clusterIndex[j]);
				}
				
				jStart += clusterSize;
				
				// compute mean using DBA
				newCentroids[i] = distComputer.dba(oldCentroids[i].clone(), cluster, 10, w);
			} else {
				nbDataInClusters[i] = 0; 
				newCentroids[i] = null; 
			}
		}
		
		return newCentroids;
	}
	
	private final void clusterAroundCentroid(final ArrayList<double[]> data, final double[][] centroids, final int w, final int[] clusterSeed) {
		for (int i = 0; i < data.size(); i++) {
			final double[][] wedge = distComputer.envelope(data.get(i), w);
			
			// LB_Keogh-NN-DTW
			double bestSoFar = Double.POSITIVE_INFINITY;
			for (int j = 0; j < centroids.length; j++) {
				final double lbDistance = distComputer.lb_keogh(wedge, centroids[j]);
				if (lbDistance < bestSoFar) {
					final double dtwDistance = distComputer.dtw(centroids[j], data.get(i), w);
					if (dtwDistance < bestSoFar) {
						bestSoFar = dtwDistance;
						clusterSeed[i] = j;
					}
				}
			}
		}
	}
	
	public final int[] getClusterSeed() {
		return clusterSeed;
	}
}
