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
	 * 
	 * Last modified: 14/01/2017
	 */
	private static final Tools tools = new Tools();	// tools for the functions
	private static String datasetName;				// dataset name
	private static int nbClass;						// number of class
	private static int lengthTimeseries;			// length of time series
	private static int sizeTrain;					// size of training dataset
	private static int sizeTest;					// size of testing dataset
	private static int warpingWindow;				// warping window size
										
	private static ArrayList<double[]> trainingDataset;		// training dataset
	private static ArrayList<double[]> testingDataset;		// testing dataset
	private static ArrayList<Integer> trainingDatasetClass;	// class of training dataset
	private static ArrayList<Integer> testingDatasetClass;	// class of testing dataset
	private static ArrayList<Integer> trainingDatasetIndex;	// index of training dataset
	private static ArrayList<Integer> testingDatasetIndex;	// index of testing dataset
		
	private ArrayList<double[]> randDataset;		// random dataset
	private int[] randIndex;						// index of random dataset
	private ArrayList<Integer> randClass;			// class of random dataset
	
	/* Constructors */
	public Dataset(final String datasetname) {
		datasetName = datasetname;
	}
	
	public final void loadTestSet(final String path) {
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// testing dataset file
		final String testFile = path + datasetName + "/" + datasetName + "_TEST";
		
		// read dataset
		tools.readDatasetCSV(testFile, sizeTest);
		testingDatasetIndex = tools.getDataIndexAftRead();
		testingDatasetClass = tools.getDataLabelAftRead();
		testingDataset = tools.getTimeseriesAftRead();
	}
	
	public final void loadTestSetCSV(final String path) {
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// testing dataset file
		final String testCSV = path + datasetName + "/" + datasetName + "_TEST.csv";
		
		// read dataset from csv
		tools.readDatasetCSV(testCSV, sizeTest);
		testingDatasetIndex = tools.getDataIndexAftRead();
		testingDatasetClass = tools.getDataLabelAftRead();
		testingDataset = tools.getTimeseriesAftRead();
	}
	
	public final void loadTrainSet(final String path) {
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// training dataset file
		final String trainFile = path + datasetName + "/" + datasetName + "_TRAIN";
		
		// read dataset from csv
		tools.readDatasetCSV(trainFile, sizeTest);
		trainingDatasetIndex = tools.getDataIndexAftRead();
		trainingDatasetClass = tools.getDataLabelAftRead();
		trainingDataset = tools.getTimeseriesAftRead();
	}
	
	public final void loadTrainSetCSV(final String path) {
		getDatasetProperties(path);	// get dataset properties, number of class, warping window etc
		
		// training dataset file
		final String trainCSV = path + datasetName + "/" + datasetName + "_TRAIN.csv";
		
		// read dataset from csv
		tools.readDatasetCSV(trainCSV, sizeTest);
		trainingDatasetIndex = tools.getDataIndexAftRead();
		trainingDatasetClass = tools.getDataLabelAftRead();
		trainingDataset = tools.getTimeseriesAftRead();
	}
	
	public final void save1NNIndex(final int[] index1NN, final String projectPath) {
		final String csvFile = projectPath + datasetName + "_1NN_LB_index1NN.csv";	
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
		final int[] index1NN = new int[sizeTest];		
		final String delimiter = ",";					
		final String csvFile = path + datasetName + "_1NN_LB_index1NN.csv";	
		String line = "";
		BufferedReader br = null;
		
		System.out.println("Reading " + csvFile);
		try {
			int count = 0;
			br = new BufferedReader(new FileReader(csvFile));
			while ((line = br.readLine()) != null) 
			{
				final String[] data = line.split(delimiter);
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
		final String filename = path + datasetName + "/" + datasetName + "_properties.csv";	
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
	
	public final int[] randomize(final ArrayList<double[]> dataset, final ArrayList<Integer> dataClass, final ArrayList<Integer> dataIndex) {
		System.out.println("Randomizing dataset");
		
		// temporary data index
		final ArrayList<Integer> tempData = new ArrayList<Integer>(dataset.size());
		tempData.addAll(dataIndex);
		
		// initialize
		randDataset = new ArrayList<double[]>(dataset.size());
		randClass = new ArrayList<Integer>(dataset.size());
		randIndex = new int[dataset.size()]; 
		int dataCount = 0;
		
		while (tempData.size() > 0) {
			final Random rand = new Random();
			final int idx = rand.nextInt(tempData.size());
			randIndex[dataCount] = tempData.get(idx);
			randDataset.add(dataset.get(randIndex[dataCount]));
			randClass.add(dataClass.get(randIndex[dataCount]));
			dataCount++;
			tempData.remove(idx);
		}
		
		return randIndex;
	}
	
	public final ArrayList<Integer> randomize2(final ArrayList<double[]> dataset, final ArrayList<Integer> dataClass, final ArrayList<Integer> dataIndex) {
		System.out.println("Randomizing dataset");
		
		// temporary data index
		final ArrayList<Integer> tempData = new ArrayList<Integer>(dataset.size());
		tempData.addAll(dataIndex);
		
		// initialize
		randDataset = new ArrayList<double[]>(dataset.size());
		randClass = new ArrayList<Integer>(dataset.size());
		ArrayList<Integer> randIndex = new ArrayList<Integer>(dataset.size()); 
		int dataCount = 0;
		
		while (tempData.size() > 0) {
			final Random rand = new Random();
			final int idx = rand.nextInt(tempData.size());
			randIndex.add(tempData.get(idx));
			randDataset.add(dataset.get(randIndex.get(dataCount)));
			randClass.add(dataClass.get(randIndex.get(dataCount)));
			dataCount++;
			tempData.remove(idx);
		}
		
		return randIndex;
	}
	
	public final ArrayList<double[]> extractDataset(final ArrayList<double[]> dataset, final int[] dataIndex) {
		final ArrayList<double[]> newDataset = new ArrayList<double[]>(dataIndex.length);
		for (int i = 0; i < dataIndex.length; i++) {
			newDataset.add(dataset.get(dataIndex[i]));
		}
		
		return newDataset;
	}
	
	public final ArrayList<Integer> extractLabel(final ArrayList<Integer> dataset, final int[] dataIndex) {
		final ArrayList<Integer> newLabel = new ArrayList<Integer>(dataIndex.length);
		for (int i = 0; i < dataIndex.length; i++) {
			newLabel.add(dataset.get(dataIndex[i]));
		}
		
		return newLabel;
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
		return trainingDataset;
	}
	
	public final ArrayList<double[]> TestSet() {
		return testingDataset;
	}
	
	public final ArrayList<Integer> ClassTrain() {
		return trainingDatasetClass;
	}
	
	public final ArrayList<Integer> ClassTest() {
		return testingDatasetClass;
	}
	
	public final ArrayList<Integer> IndexTrain() {
		return trainingDatasetIndex;
	}
	
	public final ArrayList<Integer> IndexTest() {
		return testingDatasetIndex;
	}
}
