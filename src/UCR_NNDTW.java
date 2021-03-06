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
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import classifiers.ClassifierNNDTW;
import datasets.Dataset;

public class UCR_NNDTW {
	/* This is for LB_Keogh-NN-DTW Experiment on UCR datasets
	 * 
	 * Last modified: 14/01/2017
	 */
	public static void main (String[] args) {
		String projectPath = "";
		String datasetName = "ArrowHead";	
		int w = 0;	
		int testNum = 11;
		int step = 1;
		
		if (args.length > 0) {projectPath = args[0] + "/";}
		if (args.length > 1) {datasetName = args[1];}
		if (args.length > 2) {testNum = Integer.parseInt(args[2]);}
		if (args.length > 3) {step = Integer.parseInt(args[3]);}
		
		final Dataset UCR = new Dataset(datasetName);		// tools to read and load UCR dataset
		UCR.loadTestSet(projectPath + "dataset/UCR_Time_Series_Archive/");
		UCR.loadTrainSet(projectPath + "dataset/UCR_Time_Series_Archive/");
		w = UCR.window();
		
		if (args.length > 4) {w = Integer.parseInt(args[4]);}
		
		int stepSize = 0;
		if (step != 1) {stepSize = (int) Math.ceil((double) UCR.trainSize()/step) + 1;}
		else {stepSize = (int) Math.ceil((double) UCR.trainSize()/step);}
		
		final double[][] averageErrorPerQuery = new double[testNum][stepSize];
		final double[][] averageTimePerQuery = new double[testNum][stepSize];
		final double[][] averageDistPerQuery = new double[testNum][stepSize];
		int[] seenSoFar = new int[stepSize];
		
		for (int i = 0; i < testNum; i++) {
			ArrayList<double[]> trainingDataset = UCR.TrainSet();
			ArrayList<Integer> trainingDatasetClass = UCR.ClassTrain();
			ArrayList<Integer> trainingDatasetIndex = UCR.IndexTrain();
			trainingDatasetIndex = UCR.randomize2(trainingDataset, trainingDatasetClass, trainingDatasetIndex);
			trainingDataset = UCR.getRandDataset();
			trainingDatasetClass = UCR.getRandLabel();

			ArrayList<double[]> testingDataset = UCR.TestSet();
			ArrayList<Integer> testingDatasetClass = UCR.ClassTest();
			ArrayList<Integer> testingDatasetIndex = UCR.IndexTest();
			testingDatasetIndex = UCR.randomize2(testingDataset, testingDatasetClass, testingDatasetIndex);
			testingDataset = UCR.getRandDataset();
			testingDatasetClass = UCR.getRandLabel();
			
			final ClassifierNNDTW classifier = new ClassifierNNDTW(trainingDataset, trainingDatasetClass, trainingDatasetIndex, step, stepSize);
			final double errorRate = classifier.performance(testingDataset, testingDatasetClass, w);
			final double averageTime = classifier.getAverageQueryTime();
			final double distComputations = classifier.getDistComputation();
			
			averageErrorPerQuery[i] = classifier.getAverageErrorPerQuery();
			averageTimePerQuery[i] = classifier.getAverageTimePerQuery();
			averageDistPerQuery[i] = classifier.getAverageDistPerQuery();
			seenSoFar = classifier.getSeenSoFarPerQuery();
			
			System.out.println((i+1) + ", Error rate: " + errorRate + ", Average Query Time is: " + averageTime + ", Average Distance Computations: " + distComputations);
		}
		final String csvFile = projectPath + "outputs/experiments/" + datasetName + "_NNDTW.csv";
		writeToCsv(csvFile, datasetName, averageErrorPerQuery, averageTimePerQuery, averageDistPerQuery, seenSoFar);
	}
	
	private final static void writeToCsv(final String csvFile, final String datasetName, 
			final double[][] error, final double[][] time, 
			final double[][] dist, final int[] seen) {
		try {
			FileWriter writer = new FileWriter(csvFile);
			writer.append(datasetName + '\n');
			writer.append("Runs,Error rate" + '\n');
			writer.append("Seen,");
			writer.append(seen[0] + ",");
			for (int i = 1; i < seen.length; i++)
				writer.append(Integer.toString(seen[i]) + ",");
			writer.append('\n');
			for (int i = 1; i < error.length; i++) {
				writer.append((i-1) + "," + error[i][0]);
				for (int j = 1; j < error[0].length; j++)
					writer.append("," + Double.toString(error[i][j]));
				writer.append('\n');
			}
			writer.append('\n');
			
			writer.append("Runs, Time" + '\n');
			writer.append("Seen,");
			writer.append(seen[0] + ",");
			for (int i = 1; i < seen.length; i++)
				writer.append(Integer.toString(seen[i]) + ",");
			writer.append('\n');
			for (int i = 1; i < time.length; i++) {
				writer.append((i-1) + "," + time[i][0]);
				for (int j = 1; j < time[0].length; j++)
					writer.append("," + Double.toString(time[i][j]));
				writer.append('\n');
			}
			writer.append('\n');
			
			writer.append("Runs,Distance" + '\n');
			writer.append("Seen,");
			writer.append(seen[0] + ",");
			for (int i = 1; i < seen.length; i++)
				writer.append(Integer.toString(seen[i]) + ",");
			writer.append('\n');
			for (int i = 1; i < dist.length; i++) {
				writer.append((i-1) + "," + dist[i][0]);
				for (int j = 1; j < dist[0].length; j++)
					writer.append("," + Double.toString(dist[i][j]));
				writer.append('\n');
			}
			writer.append('\n');
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
