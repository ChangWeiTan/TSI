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
package tree;

import java.util.ArrayList;
import tools.Tools;

public class TreeNode implements java.io.Serializable {
	/* Tree Node
	 * This is a class for a node of the Hierarchical K-Means Tree
	 * It is a general tree structure
	 * 
	 * Last modified: 12/01/2017
	 */
	private static final Tools tools = new Tools();		// tools for the functions
	private static final long serialVersionUID = 1L;	
	private static final int NA = -1;					// constant for not available
	
	private boolean thisIsRoot = false; 				// indication if this node is root
	private boolean thisIsLeaf = false; 				// indication if this node is leaf
	private double[] centroid = null;					// stores centroid of that cluster
	private int centroidLabel = NA;						// stores class for each centroid
	private int centroidIndex = NA;
    private ArrayList<Integer> dataLabel = null;		// stores class for each data points
    private ArrayList<Integer> dataIndex = null;		// stores data index in training set
    private ArrayList<TreeNode> children = new ArrayList<TreeNode>();	// children nodes of that cluster
    private double[][] wedges = null;
    private int nbTimeseries = 0;
    
    /* Constructor */
    public TreeNode() {	
    	thisIsRoot = false;
    	thisIsLeaf = false;
    }
    
    /* Methods */
    public final void addChild(final TreeNode child) {
		children.add(child);
    }
    
    public final void addWedge(final double[][] wedge) {
    	wedges = wedge;
    }
    
    public final void createNonLeaf(final double[] Centroid, final ArrayList<Integer> Label) {
    	// initialize array then copy new centroid
    	centroid = new double[Centroid.length];
    	centroid = Centroid;
    	
    	centroidLabel = tools.mode(Label);
    }
    
    public final void createNonLeaf(final double[] Centroid, final ArrayList<Integer> Label, final ArrayList<Integer> dataIndex) {
    	// initialize array then copy new centroid
    	this.dataIndex = dataIndex;
    	centroid = new double[Centroid.length];
    	centroid = Centroid;
    	
    	centroidLabel = tools.mode(Label);
    }
    
    public final void createNonLeaf(final double[] Centroid, final int Label, final int dataIndex) {
    	// initialize array then copy new centroid
    	this.centroidIndex = dataIndex;
    	centroid = new double[Centroid.length];
    	centroid = Centroid;
    	
    	centroidLabel = Label;
    }
    
    public final void createNonLeaf(final ArrayList<Integer> Label, final ArrayList<Integer> dataIndex, final int num) {
    	this.dataIndex = dataIndex;
    	
    	centroidLabel = tools.mode(Label);
    	
    	nbTimeseries = num;
    }
    
    public final void createNonLeaf(final int centroidLabel, final int centroidIndex, final int num, final ArrayList<Integer> dataLabel, final ArrayList<Integer> dataIndex) {
    	this.centroidIndex = centroidIndex;
    	this.centroidLabel = centroidLabel;
    	this.dataIndex = dataIndex;
    	this.dataLabel = dataLabel;
    	
    	nbTimeseries = num;
    }
    
    public final void createNonLeaf(final double[] Centroid) {
    	// initialize array then copy new centroid
    	centroid = new double[Centroid.length];
    	centroid = Centroid;
    }
    
    public final void createLeaf(final ArrayList<Integer> dataIndex, final ArrayList<Integer> dataClass) {
    	this.thisIsLeaf = true;
    	this.dataIndex = new ArrayList<Integer>();
    	this.dataLabel = new ArrayList<Integer>();
    	
    	// initialize array then copy index
    	this.dataIndex.addAll(dataIndex);				
    			
    	// initialize array then copy data one by one
    	this.dataLabel.addAll(dataClass);
    }
    
    public final void createLeaf(final int dataIndex, final int dataClass) {
    	this.thisIsLeaf = true;
    	this.dataIndex = new ArrayList<Integer>();
    	this.dataLabel = new ArrayList<Integer>();
    	
    	// initialize array then copy index
    	this.dataIndex.add(dataIndex);				
    			
    	// initialize array then copy data one by one
    	this.dataLabel.add(dataClass);
    }
    
    public final void createLeaf(final int indexOfData, final int TrainClass, final int Label) {
    	this.thisIsLeaf = true;
    	
    	// initialize array then copy index
    	this.dataIndex = new ArrayList<Integer>();		
    	this.dataIndex.add(indexOfData);
    	
    	this.centroidLabel = Label;
    	
    	this.dataLabel = new ArrayList<Integer>();			
    	this.dataLabel.add(TrainClass);
    }
    
    public final double[] getCentroid() {
    	return this.centroid;
    }
    
    public final double[][] getWedges() {
    	return this.wedges;
    }

    public final int getCentroidLabel() {
    	return this.centroidLabel;
    }

    public final ArrayList<Integer> getDataLabel() {
    	return this.dataLabel;
    }
    
    public final ArrayList<Integer> getDataIndex() {
    	return dataIndex;
    }

    public final ArrayList<TreeNode> getChildren() {
    	return this.children;
    }
    
    public final TreeNode getChild(int i) {
    	if (i > 0 && i <= children.size())
    		return this.children.get(i);
    	else
    		return null;
    }
    
    public final int numTimeseries() {
    	return this.nbTimeseries;
    }
    
    public final int numChildren() {
    	return this.children.size();
    }
    
    public final int leafSize() {
    	return this.dataIndex.size();
    }
    
    public final int getCentroidIndex() {
    	return centroidIndex;
    }
    
    public final void setRoot(final int Label) {
    	this.thisIsRoot = true;
    	this.centroidLabel = Label;
    }
    
    public final void setRoot() {
    	this.thisIsRoot = true;
    }

    public final boolean isLeaf() {
    	return this.thisIsLeaf;
    }
    
    public final boolean isMediodLeaf() {
    	return this.dataIndex != null;
    }

    public boolean isRoot() {
    	return this.thisIsRoot;
    }
}
