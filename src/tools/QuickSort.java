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

public class QuickSort {
	/* Quick sort
	 * This is a class for Quick Sort
	 * 
	 * Last modified: 12/01/2017
	 */
	public final static void sort(double[] numbers, int[] index, final int low, final int high) {
		int i = low, j = high;
		
		final double pivot = numbers[(low + high)/2];
		
		while (i <= j) {
			while (numbers[i] < pivot)
				i++;
			while (numbers[j] > pivot)
				j--;
			if (i <= j) {
				swap(numbers, index, i, j);
				i++;
				j--;
			}
		}
		if (low < j)
			sort(numbers, index, low, j);
		if (i < high)
			sort(numbers, index, i, high);
	}
	
	private final static void swap(double[] numbers, int[] index, final int i, final int j)	{
		final double tempNum = numbers[i];
		final int tempIndex = index[i];
		numbers[i] = numbers[j];
		index[i] = index[j];
		numbers[j] = tempNum;
		index[j] = tempIndex;
	}
}
