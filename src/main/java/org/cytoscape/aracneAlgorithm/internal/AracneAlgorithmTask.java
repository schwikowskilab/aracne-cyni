/*
  File: EqualDiscretizationTask.java

  Copyright (c) 2006, 2010, The Cytoscape Consortium (www.cytoscape.org)

  This library is free software; you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published
  by the Free Software Foundation; either version 2.1 of the License, or
  any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
  documentation provided hereunder is on an "as is" basis, and the
  Institute for Systems Biology and the Whitehead Institute
  have no obligations to provide maintenance, support,
  updates, enhancements or modifications.  In no event shall the
  Institute for Systems Biology and the Whitehead Institute
  be liable to any party for direct, indirect, special,
  incidental or consequential damages, including lost profits, arising
  out of the use of this software and its documentation, even if the
  Institute for Systems Biology and the Whitehead Institute
  have been advised of the possibility of such damage.  See
  the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/
package org.cytoscape.aracneAlgorithm.internal;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.HashMap;
import java.util.Map;
import java.util.Collections;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JOptionPane;

import org.cytoscape.model.CyColumn;
import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNetworkTableManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric.Mutual_Info;
import org.cytoscape.cyni.*;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.vizmap.VisualMappingManager;
import java.util.Comparator;

import org.apache.commons.math3.distribution.*;





/**
 * The CyniSampleAlgorithmTask provides a simple example on how to create a cyni task
 */
public class AracneAlgorithmTask extends AbstractCyniTask {
	
	private final CyTable mytable;
	private CyCyniMetric selectedMetric;
	private final List<String> attributeArray;
	private CyLayoutAlgorithmManager layoutManager;
	private CyRootNetworkManager rootMgr;
	private CyniNetworkUtils netUtils;
	private boolean edgePresence[][];
	private double MIScore[][];
	private int miSteps;
	private double bandwith;
	private double dpiTol;
	private boolean manualKernel;
	private Mutual_Info.ALGORITHM algorithm;
	private boolean usePValue;
	private double pvalue;
	private double miTh;
	public static double KERNEL_ALPHA = 0.52477;
    public static double KERNEL_BETA = -0.24;
    public static double THRESHOLD_ALPHA = 1.062;
    public static double THRESHOLD_BETA = -48.7;
    public static double THRESHOLD_GAMMA = -0.634;
	

	/**
	 * Creates a new object.
	 */
	public AracneAlgorithmTask(final String name, final AracneAlgorithmContext context, CyNetworkFactory networkFactory, CyNetworkViewFactory networkViewFactory,
			CyNetworkManager networkManager,CyNetworkTableManager netTableMgr, CyRootNetworkManager rootNetMgr, VisualMappingManager vmMgr,
			CyNetworkViewManager networkViewManager,CyLayoutAlgorithmManager layoutManager, 
			CyCyniMetricsManager metricsManager, CyTable selectedTable)
	{
		super(name, context,networkFactory,networkViewFactory,networkManager, networkViewManager,netTableMgr,rootNetMgr, vmMgr);
		
		this.mytable = selectedTable;
		this.rootMgr = rootNetMgr;
		this.layoutManager = layoutManager;
		this.netUtils = new CyniNetworkUtils(networkViewFactory,networkManager,networkViewManager,netTableMgr,rootNetMgr,vmMgr);
		this.attributeArray = context.attributeList.getSelectedValues();
		this.bandwith = context.kernelWidth;
		this.dpiTol = context.dpiTol;
		this.manualKernel = context.manualKernel;
		this.miSteps = context.miSteps;
		selectedMetric = metricsManager.getCyniMetric("MIAracneMetric");
		this.miTh = context.miTh;
		this.pvalue = context.pvalue;
		if(context.algoChooser.getSelectedValue().matches("Naive Bayes"))
			this.algorithm = Mutual_Info.ALGORITHM.NAIVE_BAYES;
		if(context.algoChooser.getSelectedValue().matches("Adaptive Partitioning"))
			this.algorithm = Mutual_Info.ALGORITHM.ADAPTIVE_PARTITIONING;
		if(context.algoChooser.getSelectedValue().matches("Fixed Bandwith"))
			this.algorithm = Mutual_Info.ALGORITHM.FIXED_BANDWIDTH;
		if(context.algoChooser.getSelectedValue().matches("Variable Bandwith"))
			this.algorithm = Mutual_Info.ALGORITHM.VARIABLE_BANDWIDTH;
		if(context.thresholdChooser.getSelectedValue().matches("P-Value Threshold"))
			this.usePValue = true;
		else
			this.usePValue = false;
	}

	/**
	 *  Perform actualtask.
	 */
	@Override
	final protected void doCyniTask(final TaskMonitor taskMonitor) {
		
		Double progress = 0.0d;
		CyNetwork networkSelected = null;
		String networkName;
		CyNode node1,node2;
		Integer numNodes = 1;
		double chosenTh;
		CyLayoutAlgorithm layout;
		CyNetworkView newNetworkView ;
		CyTable nodeTable, edgeTable, netTable;
		Double step;
		CyEdge edge;
		ArrayList<Integer> index = new ArrayList<Integer>(1);
		int nRows,threadNumber;
		double threadResults[] = new double[nThreads];
		double result;
		CyNode mapRowNodes[];
		int threadIndex[] = new int[nThreads];
		threadNumber=0;
		List<CyEdge> edgeList;
		Arrays.fill(threadResults, 0.0);
		
		index.add(1);
		
        //step = 1.0 /  attributeArray.size();
        
        taskMonitor.setStatusMessage("ARACNE Algorithm running ...");
		taskMonitor.setProgress(progress);
		
		//Create new network
		CyNetwork newNetwork = netFactory.createNetwork();
		
		
		//Check if a network is associated to the selected table
		networkSelected = netUtils.getNetworkAssociatedToTable(mytable);
		
		// Create the CyniTable
		AracneCyniTable data = new AracneCyniTable(mytable,attributeArray.toArray(new String[0]), false, false, selectedOnly);
		
		data.computeMarkerVariance();
        data.computeBandwidth();
        data.computeMarkerRanks();
		
		//Set the name of the network, another name could be chosen
		networkName = "ARACNE Inference " + newNetwork.getSUID();
		if (newNetwork != null && networkName != null) {
			CyRow netRow = newNetwork.getRow(newNetwork);
			netRow.set(CyNetwork.NAME, networkName);
		}
		
		nodeTable = newNetwork.getDefaultNodeTable();
		edgeTable = newNetwork.getDefaultEdgeTable();
		netTable = newNetwork.getDefaultNetworkTable();
		netUtils.addColumns(networkSelected,newNetwork,mytable,CyNode.class, CyNetwork.LOCAL_ATTRS);
	
		edgeTable.createColumn("Metric", String.class, false);
		edgeTable.createColumn("Mutual Information", Double.class, false);	
		netTable.createColumn("Mutual Information Algorithm", String.class, false);
		netTable.getRow(newNetwork.getSUID()).set("Mutual Information Algorithm", algorithm.toString());
		
		nRows = data.nRows();
		step = 1.0 / nRows;
		
		mapRowNodes = new CyNode[nRows];
		Arrays.fill(mapRowNodes, null);
		threadResults = new double[nRows];
		threadIndex = new int[nRows];
		
		if(usePValue)
			chosenTh = findThreshold(data.nColumns(), pvalue);
		else
			chosenTh = miTh;
		
		if(!manualKernel)
			bandwith = KERNEL_ALPHA * Math.pow(data.nColumns(), KERNEL_BETA);
		
		Map<String,Object> params = new HashMap<String,Object>();
		params.put("KernelWidth", bandwith);
		params.put("MiSteps", miSteps);
		params.put("Size", data.nColumns());
		params.put("Type", algorithm);
		selectedMetric.setParameters(params);
		
		taskMonitor.setStatusMessage("ARACNE Algorithm running: Statistics calculated ...");
		
		edgePresence = new boolean[nRows][nRows];
		MIScore = new double [nRows][nRows];
		// Create the thread pools
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);

		for (int i = 0; i < nRows; i++) 
		{
					threadNumber = 0;
					if (cancelled)
						break;

					for (int j = i+1; j < nRows; j++) 
					{
						if (cancelled)
							break;
						
						index.set(0, j);
						executor.execute(new ThreadedGetMetric(data,i,index,threadNumber));
						threadIndex[threadNumber] = j;
						threadNumber++;
					}
					executor.shutdown();
					// Wait until all threads are finish
					try {
			         	executor.awaitTermination(7, TimeUnit.DAYS);
			        } catch (Exception e) {}
					
					for (int j = i+1; j < nRows; j++) 
					{
						if(MIScore[i][j]> chosenTh)
						{
							edgePresence[i][j]= true;
							edgePresence[j][i]= true;
							if(mapRowNodes[i] == null)
							{
								node1 = newNetwork.addNode();
								netUtils.cloneRow(newNetwork,CyNode.class,mytable.getRow(data.getRowLabel(i)), newNetwork.getRow(node1));
								if(newNetwork.getRow(node1).get(CyNetwork.NAME,String.class ) == null || newNetwork.getRow(node1).get(CyNetwork.NAME,String.class ).isEmpty() == true)
									newNetwork.getRow(node1).set(CyNetwork.NAME, "Node " + numNodes);
								if(newNetwork.getRow(node1).get(CyNetwork.SELECTED,Boolean.class ) == true)
									newNetwork.getRow(node1).set(CyNetwork.SELECTED, false);
								mapRowNodes[i] =node1;
								numNodes++;
							}
							if(mapRowNodes[j] == null)
							{
								node2 = newNetwork.addNode();
								netUtils.cloneRow(newNetwork,CyNode.class,mytable.getRow(data.getRowLabel(j)), newNetwork.getRow(node2));
								if(newNetwork.getRow(node2).get(CyNetwork.NAME,String.class ) == null || newNetwork.getRow(node2).get(CyNetwork.NAME,String.class ).isEmpty() == true)
									newNetwork.getRow(node2).set(CyNetwork.NAME, "Node " + numNodes);
								if(newNetwork.getRow(node2).get(CyNetwork.SELECTED,Boolean.class ) == true)
									newNetwork.getRow(node2).set(CyNetwork.SELECTED, false);
								mapRowNodes[j] = node2;
								numNodes++;
							}
									
							if(!newNetwork.containsEdge(mapRowNodes[i], mapRowNodes[j]))
							{
								edge = newNetwork.addEdge(mapRowNodes[i], mapRowNodes[j], false);
								newNetwork.getRow(edge).set("Mutual Information", MIScore[i][j]);
								newNetwork.getRow(edge).set("Metric",selectedMetric.toString());
								newNetwork.getRow(edge).set("name", newNetwork.getRow(mapRowNodes[i]).get("name", String.class)
										+ " (Aracne) " + newNetwork.getRow( mapRowNodes[j]).get("name", String.class));
							}
						}
						else
						{
							edgePresence[i][j]= false;
							edgePresence[j][i]= false;
							
						}
							
						
					}
					progress = progress + step;
					taskMonitor.setProgress(progress);
					
					threadNumber = 0;
					executor = Executors.newFixedThreadPool(nThreads);
		}
		
		if(dpiTol < 1)
		{
			taskMonitor.setStatusMessage("ARACNE Algorithm running: Applying DPI ...");
		
			Vector<ArrayValuePair> miVector = new Vector<ArrayValuePair>();
	        // put all the values of the i-th row in the matrix on a vector for sorting
	        // the ArrayValuePair structure is borrowed from the Microarray_Set class,
	        // it should actually be called EdgeMIPair here. They share the same
	        // structure consisting of two fields: a int and a double. The assoiciated
	        // sorting class will sort by the double fields.
	        for (int i = 0; i < nRows; i++) {
	        	miVector.clear();
	        	if (cancelled)
					break;
	        	 
	        	for (int j = 0; j < nRows; j++) 
				{
	        		if(edgePresence[i][j])
	        			miVector.add(new ArrayValuePair(j,MIScore[i][j]));
				}
	        	
	        	 // sort the vector
	            Collections.sort(miVector, new SortDecreasing_ArrayValuePair());
	
	            // for each value in the vector (these are genes directly connected to row_idx.
	            // i.e. geneId1 is connected to i)
	            for (int j = 0; j < miVector.size(); j++) {
	                // (geneId1, valueAB) are the (id, value) of the selected gene, starting with the largest mi
	                int geneId1 = miVector.get(j).getId();
	                double valueAB = miVector.get(j).getValue();
	
	                // Set the limits
	                double minMI = valueAB / (1.0 - dpiTol);
	                
	                // Loop on all the genes that have an equal or higher mutual information
	                // (these genes are also connected to i) bool connected = true;
	                for (int k = 0; k < j; k++) {
	                    int geneId2 = miVector.get(k).getId();
	
	                    double valueAC = miVector.get(k).getValue();
	                    // We know that valueAB is smaller than valueAC by default; then tolerance is always positive
	                    if (valueAC <= minMI) {
	                        // we can exit from this loop because all future values are smaller
	                        break;
	                    }
	                    // We know that valueAB is smaller than valueAC by default.
	                    // Therefore, if it is also smaller than valueBC, then
	                    // geneId1 is connected to gene i via geneId2
	                    double valueBC = MIScore[ geneId1][ geneId2];
	                    if (valueBC > minMI && edgePresence[geneId1][geneId2]) {
	                    	if(newNetwork.containsEdge(mapRowNodes[i], mapRowNodes[geneId1]))
	                    	{
	                    		edgeList = newNetwork.getConnectingEdgeList(mapRowNodes[i], mapRowNodes[geneId1], CyEdge.Type.UNDIRECTED);
	                    		rootMgr.getRootNetwork(newNetwork).removeEdges(edgeList);
	                    		//newNetwork.getDefaultEdgeTable().deleteRows(Collections.singletonList(edgeList.get(0).getSUID()));
	                    		//newNetwork.removeEdges(edgeList);
	                    		//System.out.println("removed " + i + " " + geneId1);
	                    	}
	                        break;
	                    }
	                }
	            }
	
	        }

		}
		/*for (int i = 0; i < nRows; i++) 
		{
			if (cancelled)
				break;
			
			for (int j = i+1; j < nRows; j++) 
			{
				if(!edgePresence[i][j])
					continue;
				if (cancelled)
					break;
				for (int k = j+1; k < nRows; k++) 
				{
					if(!edgePresence[j][k] || !edgePresence[i][k])
						continue;
					
					applyDPI(i,j,k);
					
				}
				
			}
		}*/
		
		
		//Display the new network
		if (!cancelled)
		{
			newNetworkView = netUtils.displayNewNetwork(newNetwork, networkSelected,false);
			taskMonitor.setProgress(1.0d);
			layout = layoutManager.getDefaultLayout();
			Object context = layout.getDefaultLayoutContext();
			insertTasksAfterCurrentTask(layout.createTaskIterator(newNetworkView, context, CyLayoutAlgorithm.ALL_NODE_VIEWS,""));
		}
		
		taskMonitor.setProgress(1.0d);
	}
	
	private double findThreshold(int n, double pvalue) {
       
        return (THRESHOLD_ALPHA - Math.log(pvalue)) / ((-THRESHOLD_BETA) + (-THRESHOLD_GAMMA) * n);
    }

	private void applyDPI(int index1, int index2, int index3)
	{
		
	}
	
	private class ThreadedGetMetric implements Runnable {
		private ArrayList<Integer> index2;
		private int index1;
		private AracneCyniTable tableData;
		
		ThreadedGetMetric(AracneCyniTable data,int index1, ArrayList<Integer> parentsToIndex,int pos)
		{
			this.index2 = new ArrayList<Integer>( parentsToIndex);
			this.index1 = index1;
			this.tableData = data;
			
		}
		
		public void run() {
			MIScore[index1][index2.get(0)] = selectedMetric.getMetric(tableData, tableData, index1, index2);

		}
		

	}
	
	static class ArrayValuePair {
        int arrayId;
        double value;

        ArrayValuePair(int Id, double v) {
            set(Id, v);
        }

        void set(int iId, double v) {
            arrayId = iId;
            value = v;
        }

        int getId() {
            return arrayId;
        }

        double getValue() {
            return value;
        }
        
        @Override
        public boolean equals(Object o) {
            ArrayValuePair avp = (ArrayValuePair)o;
            return (avp.arrayId == arrayId && avp.value == value);
        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 71 * hash + this.arrayId;
            hash = 71 * hash + (int) (Double.doubleToLongBits(this.value) ^ (Double.doubleToLongBits(this.value) >>> 32));
            return hash;
        }
    }

    static class SortIncreasing_ArrayValuePair implements Comparator<ArrayValuePair> {
        public int compare(ArrayValuePair o1, ArrayValuePair o2) {
            return (int) Math.signum(o1.getValue() - o2.getValue());
        }
    }

    static class SortDecreasing_ArrayValuePair implements Comparator<ArrayValuePair> {
        public int compare(ArrayValuePair o1, ArrayValuePair o2) {
            return (int) Math.signum(o2.getValue() - o1.getValue());
        }
    }

	
	
}
