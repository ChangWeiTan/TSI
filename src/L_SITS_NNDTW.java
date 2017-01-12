import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import classifiers.ClassifierNNDTW;
import datasets.Dataset;

public class L_SITS_NNDTW {
	/* L_SITS_NNDTW
	 * This is a class for LB_Keogh-NN-DTW Experiment on SITS dataset
	 * 
	 * Last modified: 01/10/2016
	 */
	public static void main (String[] args) {
		String projectPath = "";				// project path
		String datasetName = "SITS1M_fold1";	// dataset name
		
		if (args.length > 0) {projectPath = args[0] + "/";}
		if (args.length > 1) {datasetName = args[1];}
		int w = 0;	
		int testNum = 1;
		int step = 10;
		if (args.length > 2)
		{
			testNum = Integer.parseInt(args[2]);
			step = Integer.parseInt(args[3]);
		}
		
		final Dataset SITS = new Dataset(datasetName);		// tools to read and load SITS dataset
		SITS.loadTestcsv(projectPath + "dataset/SITS_2006_NDVI_C/");
		SITS.loadTraincsv(projectPath + "dataset/SITS_2006_NDVI_C/");
		w = SITS.window();
		
		if (args.length > 4) {w = Integer.parseInt(args[4]);}
		
		int stepSize = 0;
		if (step != 1) {stepSize = (int) Math.ceil((double) SITS.trainSize()/step) + 1;}
		else {stepSize = (int) Math.ceil((double) SITS.trainSize()/step);}
		
		final double[][] averageErrorPerQuery = new double[testNum][stepSize];
		final double[][] averageTimePerQuery = new double[testNum][stepSize];
		final double[][] averageDistPerQuery = new double[testNum][stepSize];
		int[] seenSoFar = new int[stepSize];
		
		for (int i = 0; i < testNum; i++)
		{
			ArrayList<double[]> TRAIN = SITS.TrainSet();
			ArrayList<Integer> TrainClass = SITS.ClassTrain();
			ArrayList<Integer> TrainIndex = SITS.IndexTrain();
			TrainIndex = SITS.randomize2(TRAIN, TrainClass, TrainIndex);
			TRAIN = SITS.getRandDataset();
			TrainClass = SITS.getRandLabel();
			ArrayList<double[]> TEST = SITS.TestSet();
			ArrayList<Integer> TestClass = SITS.ClassTest();
			
			final ClassifierNNDTW classifier = new ClassifierNNDTW(TRAIN, TrainClass, TrainIndex, step, stepSize);
			final double errorRate = classifier.performance(TEST, TestClass, w);
			final double averageTime = classifier.getAverageQueryTime();
			final double distComputations = classifier.getDistComputation();
			
			averageErrorPerQuery[i] = classifier.getAverageErrorPerQuery();
			averageTimePerQuery[i] = classifier.getAverageTimePerQuery();
			averageDistPerQuery[i] = classifier.getAverageDistPerQuery();
			seenSoFar = classifier.getSeenSoFarPerQuery();
			
			System.out.println((i+1) + ", Error rate: " + errorRate + ", Average Query Time is: " + averageTime + ", Average Distance Computations: " + distComputations);
		}
		final String csvFile = projectPath + "outputs/L experiment/" + datasetName + "_NNDTW.csv";
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
			for (int i = 0; i < error.length; i++) {
				writer.append("," + error[i][0]);
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
			for (int i = 0; i < time.length; i++) {
				writer.append("," + time[i][0]);
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
			for (int i = 0; i < dist.length; i++) {
				writer.append("," + dist[i][0]);
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
