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
