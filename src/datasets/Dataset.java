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
package datasets;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import tools.Tools;

public class Dataset {
	/* Dataset
	 * This is a class for time series dataset
	 * 
	 * Last modified: 12/01/2017
	 */
	private static final Tools tools = new Tools();	// tools for the functions
	private static String datasetName;				// dataset name
	private static int nbClass;						// number of class
	private static int lengthTimeseries;			// length of time series
	private static int sizeTrain;					// size of training dataset
	private static int sizeTest;					// size of testing dataset
	private static int warpingWindow;				// warping window size
										
	private static ArrayList<double[]> TRAIN;		// training dataset
	private static ArrayList<double[]> TEST;		// testing dataset
	private static ArrayList<Integer> TrainClass;	// class of training dataset
	private static ArrayList<Integer> TestClass;	// class of testing dataset
	private static ArrayList<Integer> TrainIndex;	// index of training dataset
	private static ArrayList<Integer> TestIndex;	// index of testing dataset
		
	private ArrayList<double[]> randDataset;		// random dataset
	private int[] randIndex;						// index of random dataset
	private ArrayList<Integer> randClass;			// class of random dataset
	
	private ArrayList<double[]> subDataset;			// sub dataset 
	private ArrayList<Integer> subClass;			// class of sub dataset
	
	/* Constructors */
	public Dataset(final String datasetname) {
		datasetName = datasetname;
	}
	
	/* Method */
	public final void loadTest(final String path) {
		/* Load testing dataset 
		 * Inputs: 
		 * 	path	: project path 
		 */
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// testing dataset file
		final String testFile = path + datasetName + "/" + datasetName + "_TEST";
		
		// read dataset
		tools.readDatasetCSV(testFile, sizeTest);
		TestIndex = tools.getDataIndexAftRead();
		TestClass = tools.getDataLabelAftRead();
		TEST = tools.getTimeseriesAftRead();
	}
	
	public final void loadTestcsv(final String path) {
		/* Load testing dataset in csv
		 * Inputs: 
		 * 	path	: project path 
		 */
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// testing dataset file
		final String testCSV = path + datasetName + "/" + datasetName + "_TEST.csv";
		
		// read dataset from csv
		tools.readDatasetCSV(testCSV, sizeTest);
		TestIndex = tools.getDataIndexAftRead();
		TestClass = tools.getDataLabelAftRead();
		TEST = tools.getTimeseriesAftRead();
	}
	
	public final void loadTrain(final String path) {
		/* Load training dataset 
		 * Inputs: 
		 * 	path	: project path 
		 */
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// training dataset file
		final String trainFile = path + datasetName + "/" + datasetName + "_TRAIN";
		
		// read dataset from csv
		tools.readDatasetCSV(trainFile, sizeTest);
		TrainIndex = tools.getDataIndexAftRead();
		TrainClass = tools.getDataLabelAftRead();
		TRAIN = tools.getTimeseriesAftRead();
	}
	
	public final void loadTraincsv(final String path) {
		/* Load training dataset in csv 
		 * Inputs: 
		 * 	path	: project path 
		 */
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// training dataset file
		final String trainCSV = path + datasetName + "/" + datasetName + "_TRAIN.csv";
		
		// read dataset from csv
		tools.readDatasetCSV(trainCSV, sizeTest);
		TrainIndex = tools.getDataIndexAftRead();
		TrainClass = tools.getDataLabelAftRead();
		TRAIN = tools.getTimeseriesAftRead();
	}
	
	public final void saveIndex1NN(final int[] index1NN, final String projectPath) {
		/* Save Index 1NN
		 * Inputs: 
		 * 	index1NN	: 1NN index
		 * 	projectPath	: project path 
		 */
		final String csvFile = projectPath + datasetName + "_1NN_LB_index1NN.csv";	// file name
		BufferedWriter br = null;
		
		try {
			br = new BufferedWriter(new FileWriter(csvFile));
			for (int i = 0; i < index1NN.length; i++) {
				br.append((char) index1NN[i]); 
				br.append('\n');
			}
		} catch (Exception e) {
		    e.printStackTrace();
		} 
	}
	
	public final int[] read1NNIndex(final String path) {
		/* Read the actual NN Index
		 * Inputs: 
		 * 	path		: project path 
		 * Output:
		 * 	index1NN	: 1NN index
		 */
		final int[] index1NN = new int[sizeTest];		// 1NN index
		final String csvSplitby = ",";					// delimiter
		final String csvFile = path + datasetName + "_1NN_LB_index1NN.csv";	// file name
		String line = "";
		BufferedReader br = null;
		
		System.out.println("Reading " + csvFile);
		try {
			int count = 0;
			br = new BufferedReader(new FileReader(csvFile));
			while ((line = br.readLine()) != null) 
			{
				final String[] data = line.split(csvSplitby);
				index1NN[count++] = Integer.parseInt(data[0]);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
				
		return index1NN;
	}
	
	private final void getDatasetProperties(final String path) {
		/* Get dataset properties
		 * Inputs: 
		 * 	path		: project path 
		 */
		final String filename = path + datasetName + "/" + datasetName + "_properties.csv";	// file name
		final String csvSplitby = ",";	// delimiter
		String line = "";	
		BufferedReader br = null;
		
		try {
			br = new BufferedReader(new FileReader(filename));
			while ((line = br.readLine()) != null) {
				String[] data = line.split(csvSplitby);
				nbClass = Integer.parseInt(data[0]);			// number of class
				sizeTrain = Integer.parseInt(data[1]);			// size of training dataset
				sizeTest = Integer.parseInt(data[2]);			// size of testing dataset
				lengthTimeseries = Integer.parseInt(data[3]);	// length of time series
				warpingWindow = Integer.parseInt(data[4]);		// warping window size
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}	
	
	public final int[] randomize(final ArrayList<double[]> dataset, final ArrayList<Integer> label, final ArrayList<Integer> index) {
		/* Randomize the given dataset 
		 * it returns an integer array
		 * Inputs: 
		 * 	dataset		: Dataset to be randomized
		 * 	label		: Class of the dataset
		 * 	index 		: Index of the dataset
		 * Output:
		 * 	randIndex	: Randomized index
		 */
		System.out.println("Randomizing dataset");
		
		final int NumTimeseries = dataset.size();	// size of dataset
		
		// temporary data index
		final ArrayList<Integer> list = new ArrayList<Integer>(NumTimeseries);
		list.addAll(index);
		
		// initialize
		randDataset = new ArrayList<double[]>(NumTimeseries);
		randClass = new ArrayList<Integer>(NumTimeseries);
		randIndex = new int[NumTimeseries]; 
		int count = 0;
		
		while (list.size() > 0) 
		{
			final Random rand = new Random();
			final int idx = rand.nextInt(list.size());
			randIndex[count] = list.get(idx);
			randDataset.add(dataset.get(randIndex[count]));
			randClass.add(label.get(randIndex[count]));
			count++;
			list.remove(idx);
		}
		
		return randIndex;
	}
	
	public final ArrayList<Integer> randomize2(final ArrayList<double[]> dataset, final ArrayList<Integer> label, final ArrayList<Integer> index) {
		/* Randomize the given dataset and returns arraylist
		 * it returns an arraylist
		 * Inputs: 
		 * 	dataset		: Dataset to be randomized
		 * 	label		: Class of the dataset
		 * 	index 		: Index of the dataset
		 * Output:
		 * 	randIndex	: Randomized index
		 */
		System.out.println("Randomizing dataset");
		
		final int NumTimeseries = index.size();	// size of dataset
		
		// temporary data index
		final ArrayList<Integer> list = new ArrayList<Integer>(NumTimeseries);
		list.addAll(index);
		
		// initialize
		randDataset = new ArrayList<double[]>(NumTimeseries);
		randClass = new ArrayList<Integer>(NumTimeseries);
		ArrayList<Integer> randIndex = new ArrayList<Integer>(NumTimeseries); 
		int count = 0;
		
		while (list.size() > 0) {
			final Random rand = new Random();
			final int idx = rand.nextInt(list.size());
			randIndex.add(list.get(idx));
			randDataset.add(dataset.get(randIndex.get(count)));
			randClass.add(label.get(randIndex.get(count)));
			count++;
			list.remove(idx);
		}
		
		return randIndex;
	}
	
	public final ArrayList<double[]> extractDataset(final ArrayList<double[]> dataset, final int[] index) {
		/* Extract dataset from the given index
		 * Inputs: 
		 * 	dataset		: Dataset to be randomized
		 * 	index 		: Index to be selected
		 * Output:
		 * 	newDataset	: Extracted dataset
		 */
		final int datasetSize = index.length;
		final ArrayList<double[]> newDataset = new ArrayList<double[]>(datasetSize);
		for (int i = 0; i < datasetSize; i++) {
			newDataset.add(dataset.get(index[i]));
		}
		
		return newDataset;
	}
	
	public final ArrayList<Integer> extractLabel(final ArrayList<Integer> dataset, final int[] index) {
		/* Extract dataset label from the given index
		 * Inputs: 
		 * 	dataset		: Dataset to be randomized
		 * 	index 		: Index to be selected
		 * Output:
		 * 	newLabel	: Extracted dataset label
		 */
		final int datasetSize = index.length;
		final ArrayList<Integer> newLabel = new ArrayList<Integer>(datasetSize);
		for (int i = 0; i < datasetSize; i++) {
			newLabel.add(dataset.get(index[i]));
		}
		
		return newLabel;
	}
	
	public final ArrayList<double[]> getSubDataset() {
		return subDataset;
	}
	
	public final ArrayList<Integer> getSubLabel() {
		return subClass;
	}
	
	public final ArrayList<double[]> getRandDataset() {
		return randDataset;
	}
	
	public final ArrayList<Integer> getRandLabel() {
		return randClass;
	}
	
	public final String datasetName(){
		return datasetName;
	}
	
	public final int classNb() {
		return nbClass;
	}
	
	public final int lengthTS() {
		return lengthTimeseries;
	}
	
	public final int trainSize() {
		return sizeTrain;
	}
	
	public final int testSize() {
		return sizeTest;
	}
	
	public final int window() {
		return warpingWindow;
	}
	
	public final ArrayList<double[]> TrainSet() {
		return TRAIN;
	}
	
	public final ArrayList<double[]> TestSet() {
		return TEST;
	}
	
	public final ArrayList<Integer> ClassTrain() {
		return TrainClass;
	}
	
	public final ArrayList<Integer> ClassTest() {
		return TestClass;
	}
	
	public final ArrayList<Integer> IndexTrain() {
		return TrainIndex;
	}
	
	public final ArrayList<Integer> IndexTest() {
		return TestIndex;
	}
}
