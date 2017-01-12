package pqueues;
import tree.TreeNode;

public class BranchPriorityQueue {
	/* Branch Priority Queue
	 * This is a min priority queue class for storing branches
	 * 
	 * Last modified: 12/01/2017
	 */
	// properties of the queue
	private final int NIL = -1;		// constant for not available
	private int lastElement;		// last element of the queue
	private int capacity;			// capacity of the queue
	private double[] distance;		// distance, used to sort the queue
	private int[] dataIndex;		// data index
	private TreeNode[] data;		// data to be stored
	
	private double topDistance;		// head of the queue
	private int topDataIndex;		// first data index
	private TreeNode topData;		// first data of the queue
	
	public BranchPriorityQueue(final int c) {
		capacity = c+1;
		distance = new double[capacity];	
		dataIndex = new int[capacity];
		data = new TreeNode[capacity];
		lastElement = 0;
	}
	
	public final void insert(final double d, final int i, final TreeNode N) {
		/* Insert
		 * This is the enqueue function of the priority queue
		 * Inputs: 
		 * 	d	: distance to query
		 * 	i 	: index of the branch
		 * 	N	: branch to be stored
		 */
		
		// increment counter
		if (lastElement < capacity-1) {
			lastElement++;
		}
		
		// insert to the last element
		distance[lastElement] = d;
		dataIndex[lastElement] = i;
		data[lastElement] = N;
		
		// go through the queue to make sure that the top is minimum
		int B = lastElement;
		int A = B/2;
		while (B > 1 && distance[A] > distance[B]) {
			swap(B, A);
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
		topData = data[1];
		
		// swap last and first data
		swap(1, lastElement);
		
		// remove last data
		distance[lastElement] = -1;
		dataIndex[lastElement] = -1;
		data[lastElement] = null;
		
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
	
	public final TreeNode popData() {
		return topData;
	}
	
	public final TreeNode firstData() {
		return data[1];
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
		final TreeNode tempData = data[A];
		
		distance[A] = distance[B];
		dataIndex[A] = dataIndex[B];
		data[A] = data[B];
		
		distance[B] = tempDist;
		dataIndex[B] = tempIndex;
		data[B] = tempData;
	}
}
