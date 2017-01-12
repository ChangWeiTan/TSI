package similaritymeasures;

import java.util.ArrayList;

public class DistanceMeasures {
	/* Distance Measures
	 * This is a class for all the distance measures used
	 * 
	 * Last modified: 12/01/2017
	 */
	private final static DTW dtwComputer = new DTW();
	private final static LB_Keogh lbKeoghComputer = new LB_Keogh();
	private final static DBA dbaComputer = new DBA();
	private final static ED edComputer = new ED();
	
	/**************** ED ********************/
	public final double ed(final double[] Q, final double[] C) {
		return edComputer.compute(Q, C);
	}
	
	public final long getEDDistComputation() {
		return edComputer.getDistComputation();
	}
	
	/**************** DTW *******************/
	public final double dtw(final double[] Q, final double[] C) {
		return dtwComputer.compute(Q, C);
	}
	
	public final double dtw(final double[] Q, final double[] C, final double[][] D) {
		return dtwComputer.compute(Q, C, D);
	}
	
	public final double dtw(final double[] Q, final double[] C, final int w) {
		return dtwComputer.compute(Q, C, w);
	}
	
	public final double dtw(final double[] Q, final double[] C, final int w, final double[][] D) {
		return dtwComputer.compute(Q, C, w, D);
	}
	
	public final long getDTWDistComputation() {
		return dtwComputer.getDistComputation();
	}
	
	/************** LB_Keogh ******************/
	public final double lb_keogh(final double[][] W, final double[] C) {
		return lbKeoghComputer.compute(W, C);
	}
	
	public final double lb_keogh_EA(final double[][] W, final double[] Q, final double r) {
		return lbKeoghComputer.computeEA(W, Q, r);
	}
	
	public final double lb_keogh(final double[] C) {
		return lbKeoghComputer.compute(C);
	}
	
	public final double[][] envelope(final double[] Q, final double[] C, final int w) {
		return lbKeoghComputer.envelope(Q, C, w);
	}
	
	public final double[][] envelope(final double[] Q, final ArrayList<double[]> C, final int w) {
		return lbKeoghComputer.envelope(Q, C, w);
	}
	
	public final double[][] envelope(final double[] Q, final double[] C, final int w, final double[][] W) {
		return lbKeoghComputer.envelope(Q, C, w, W);
	}
	
	public final double[][] envelope(final double[] Q, final  ArrayList<double[]> C, final int w, final double[][] W) {
		return lbKeoghComputer.envelope(Q, C, w, W);
	}
	
	public final double[][] envelope(final double[] Q, final int w) {
		return lbKeoghComputer.envelope(Q, w);
	}
	
	public final double[][] envelope(final ArrayList<double[]> C, final int w) {
		return lbKeoghComputer.envelope(C, w);
	}
	
	public final double[][] envelope(final double[] Q, final int w, final double[][] W) {
		return lbKeoghComputer.envelope(Q, w, W);
	}
	
	public final double[][] envelope(final ArrayList<double[]> C, final int w, final double[][] W) {
		return lbKeoghComputer.envelope(C, w, W);
	}
	
	public final long getLBDistComputation() {
		return lbKeoghComputer.getDistComputation();
	}
	
	/**************** DBA *******************/
	public final double[] dba(double[] T, final double[][] sequences,  final int Imax, final int w) {
		 return dbaComputer.compute(T, sequences, Imax, w);
	}
	
	public final double[] dba(double[] T, final ArrayList<double[]> sequences,  final int Imax, final int w) {
		 return dbaComputer.compute(T, sequences, Imax, w);
	}
	
	public final double[] dba(double[] T, final double[][] sequences,  final int Imax){
		return dbaComputer.compute(T, sequences, Imax);
	}
	
	public final double[] dba(double[] T, final ArrayList<double[]> sequences,  final int Imax){
		return dbaComputer.compute(T, sequences, Imax);
	}
	
	public final double[] dba(final double[][] sequences, final int Imax, final int w) {
		return dbaComputer.compute(sequences, Imax, w);
	}
	
	public final double[] dba(final ArrayList<double[]> sequences, final int Imax, final int w) {
		return dbaComputer.compute(sequences, Imax, w);
	}
	
	public final double[] dba(final double[][] sequences, final int Imax) {
		return dbaComputer.compute(sequences, Imax);
	}
	
	public final double[] dba(final ArrayList<double[]> sequences, final int Imax) {
		return dbaComputer.compute(sequences, Imax);
	}
	
	public final long getDBADistComputation() {
		return lbKeoghComputer.getDistComputation();
	}
}
