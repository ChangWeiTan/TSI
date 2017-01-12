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
package pqueues;

public class PriorityQueue {
	/* Priority Queue
	 * This is a min priority queue class for storing nn results
	 * 
	 * Last modified: 12/01/2017
	 */
	private final int NIL = -1;		// constant for not available
	private int lastElement;		// last element of the queue
	private int capacity;			// capacity of the queue
	private double[] distance;		// distance, used to sort the queue
	private int[] dataIndex;		// data index
	
	private double topDistance;		// head of the queue
	private int topDataIndex;		// first data index
	
	/* Constructor */
	public PriorityQueue(final int c) {
		capacity = c+1;						// capacity of the queue
		// initialize
		distance = new double[capacity];
		dataIndex = new int[capacity];
		lastElement = 0;
	}
	
	public final void insert(final double d, final int i) {
		/* Insert
		 * This is the enqueue function of the priority queue
		 * Inputs: 
		 * 	d	: distance to query
		 * 	i 	: index of the branch
		 */
		
		// increment counter
		if (lastElement < capacity-1) {
			lastElement++;
		}
		
		// insert to the last element
		distance[lastElement] = d;
		dataIndex[lastElement] = i;
		
		// go through the queue to make sure that the top is minimum
		int B = lastElement;
		int A = B/2;
		while (B > 1 && distance[A] > distance[B]) {
			this.swap(B, A);
			B = A;
			A = B/2;
		}
	}
	
	public final void pop() {	
		/* Pop
		 * This is the dequeue function of the priority queue
		 * The dequeued data is stored in global variable
		 */
		
		// get dequeued data
		topDistance = distance[1];
		topDataIndex = dataIndex[1];

		// swap last and first data
		swap(1, lastElement);
		
		// remove last data
		distance[lastElement] = -1;
		dataIndex[lastElement] = -1;
		
		lastElement--;
		
		// go through the queue to make sure that the top is minimum
		int A = 1;
		int B = 2*A;
		while(B <= lastElement) {
			if (B+1 <= lastElement && distance[B+1] < distance[B]) {
				B++;
			}
			
			if (distance[A] <= distance[B]) {
				break;
			}
			
			swap(A, B);
			A = B;
			B = 2*A;
		}

	}
	
	public final double firstDistance() {
		// get the distance of the first element without dequeueing
		if (isEmpty())
			return Double.POSITIVE_INFINITY;
		else
			return distance[1];
	} 
	
	public final double popDistance() {
		// get the distance after dequeueing
		return topDistance;
	}
	
	public final int firstDataIndex() {
		if (isEmpty())
			return NIL;
		else
			return dataIndex[1];
	}
	
	public final int popDataIndex() {
		return topDataIndex;
	}
	
	public final boolean isFull() {
		return lastElement == capacity-1;
	}
	
	public final boolean isEmpty() {
		return lastElement == 0;
	}
	
	public final int sizeOf() {
		return lastElement;
	}
	
	private final void swap(final int A, final int B) {
		final double tempDist = distance[A];
		final int tempIndex = dataIndex[A];
		
		distance[A] = distance[B];
		dataIndex[A] = dataIndex[B];
		
		distance[B] = tempDist;
		dataIndex[B] = tempIndex;
	}
}
