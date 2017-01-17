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
package tools;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class Tools {
	/* Tools
	 * This is a class for miscellaneous functions 
	 * 
	 * Last modified: 14/01/2017
	 */
	private ArrayList<double[]> timeseries;
	private ArrayList<Integer> dataLabel;
	private ArrayList<Integer> dataIndex;
	
	public final double squaredEuclidean(final double a, final double b) {
		return (a - b) * (a - b);
	}
	
	public final double min3(final double a, final double b, final double c) {
		if (a < b) {
			if (a < c) 
				return a;
			else 
				return c;
		} else {
			if (b < c) 
				return b;
			else 
				return c;
		}
	}
	
	public final int argMin3(final double a, final double b, final double c) {
		if (a < b) 	{
			if (a < c) 
				return 0;
			else 
				return 2;
		} else {
			if (b < c) 
				return 1;
			else 
				return 2;
		}
	}
	
	public final int mode(final ArrayList<Integer> n) {
    	final Object[] N = n.toArray();
		Arrays.sort(N);
		int count1 = 0;
		int count2 = 0;
		int popular1 = 0;
		int popular2 = 0;
		int i, j;
		
		for (i = 0; i < N.length; i++) {
			popular1 = (int) N[i];
			count1 = 1;
			for (j = i + 1; j < N.length; j++) {
				if (popular1 == (int) N[j])
					count1++;
			}
			if (count1 > count2) {
				popular2 = popular1;
				count2 = count1;
			} else if (count1 == count2) {
				popular2 = Math.min(popular1, popular2);
			}
		}
		
		return popular2;
	}
	
	public final int randInt(final int maxNum) {
		final Random randGen = new Random();
		return randGen.nextInt(maxNum);
	}
	
	public final double randDouble(){
		final Random randGen = new Random();
		return randGen.nextDouble();
	}
	
	public final void sort(double[] numbers, int[] index, final int low, final int high) {
		QuickSort.sort(numbers, index, low, high);
	}
	
	public final int[] mapRandArray(final int[] array, final ArrayList<Integer> index, final int[] mapIndex) {
		final int[] tempArray = new int[array.length];
		int i, j;
		
		// go through all the elements in the array
		for (i = 0; i < array.length; i++){
			// find the map index
			for (j = 0; j < mapIndex.length; j++) {
				if (array[index.get(i)] == mapIndex[j]) {
					tempArray[i] = j;
					break;
				}
			}
		}
		return tempArray;
	}
	
	public final int[] mapRandArray(final int[] array, final int[] index, final int[] mapIndex) {
		final int[] tempArray = new int[array.length];
		int i, j;
		
		// go through all the elements in the array
		for (i = 0; i < array.length; i++){
			// find the map index
			for (j = 0; j < mapIndex.length; j++) {
				if (array[index[i]] == mapIndex[j])	{
					tempArray[i] = j;
					break;
				}
			}
		}
		return tempArray;
	}
	
	public final void readDatasetCSV(final String csvFile, final int sizeDataset) {
		timeseries = new ArrayList<double[]>(); 
		dataLabel = new ArrayList<Integer>();
		dataIndex = new ArrayList<Integer>();
		
		BufferedReader br = null;
		String line = "";
		String delimiter = ",";
		
		int nbTimeseries = 0;
		double maxlen = Double.NEGATIVE_INFINITY;
		double minlen = Double.POSITIVE_INFINITY;
		
		System.out.println("Reading " + csvFile);
		try {
			br = new BufferedReader(new FileReader(csvFile),1024*1024*100);
			while ((line = br.readLine()) != null) {
				final String[] data = line.split(delimiter);
				final double[] ts = new double[data.length - 1];
				maxlen = Math.max(maxlen, data.length-1);
				minlen = Math.min(minlen, data.length-1);
				
				dataIndex.add(nbTimeseries);
				dataLabel.add(Integer.parseInt(data[0]));
		
				for (int i = 1; i < data.length; i++) 
					ts [i-1] = Double.parseDouble(data[i]);
				
				this.timeseries.add(ts);
				nbTimeseries++;
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
		
		timeseries.trimToSize();
		dataIndex.trimToSize();
		dataLabel.trimToSize();
	}
	
	public final ArrayList<double[]> getTimeseriesAftRead() {
		return timeseries;
	}
	
	public final ArrayList<Integer> getDataLabelAftRead() {
		return dataLabel;
	}
	
	public final ArrayList<Integer> getDataIndexAftRead() {
		return dataIndex;
	}
}
