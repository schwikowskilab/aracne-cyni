package fr.systemsbiology.aracneAlgorithm.internal;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.HashMap;
import java.util.Map;
import java.util.Collections;
import java.util.Vector;
import java.util.concurrent.*;
import java.io.*;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.cytoscape.model.CyColumn;
import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNetworkTableManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import fr.systemsbiology.aracneAlgorithm.internal.mutualInfoMetric.Mutual_Info;
import fr.systemsbiology.cyni.*;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.property.CyProperty;
import java.util.Comparator;



public class AracneAlgorithmTask extends AbstractCyniTask {
	
	private final CyTable mytable;
	private CyCyniMetric selectedMetric;
	private final List<String> attributeArray;
	private CyLayoutAlgorithmManager layoutManager;
	private CyRootNetworkManager rootMgr;
	private CyniNetworkUtils netUtils;
	private MatrixEdgePair matrix;
	private int miSteps;
	private double bandwith;
	private double dpiTol;
	private boolean manualKernel;
	private Mutual_Info.ALGORITHM algorithm;
	private boolean usePValue;
	private File fileHub;
	private File fileTFList;
	private String mapCol;
	private double pvalue;
	private double miTh;
	private String mode;
	private double chosenTh;
	private Thresholds thParams;
	double KERNEL_ALPHA = 0.52477;
    double KERNEL_BETA = -0.24;
    double THRESHOLD_ALPHA = 1.062;
    double THRESHOLD_BETA = -48.7;
    double THRESHOLD_GAMMA = -0.634;
    /** Default configuration directory used for all Cytoscape configuration files */
	public static final String DEFAULT_CONFIG_DIR = CyProperty.DEFAULT_PROPS_CONFIG_DIR ;
	
	private static final String DEF_USER_DIR = System.getProperty("user.home");
	
	private String kernel_file = join(File.separator, DEF_USER_DIR, DEFAULT_CONFIG_DIR, "3", "aracne", "config_kernel.txt");
	
	private String threshold_file = join(File.separator, DEF_USER_DIR, DEFAULT_CONFIG_DIR, "3", "aracne","config_threshold.txt");
	private String aracne_path = join(File.separator, DEF_USER_DIR, DEFAULT_CONFIG_DIR, "3", "aracne");
	private static int iteration = 0;

	

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
		iteration++;
		fileHub = null;
		fileTFList = null;
		thParams = null;
		mapCol = context.colMapping.getSelectedValue();
		
		mode = context.mode.getSelectedValue();
		
		if(context.hubFile != null && context.hubFile.exists())
			fileHub = context.hubFile;
		if(context.TFFile != null && context.TFFile.exists())
			fileTFList = context.TFFile;
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
		CyLayoutAlgorithm layout;
		CyNetworkView newNetworkView ;
		CyTable nodeTable, edgeTable, netTable;
		Double step;
		CyEdge edge;
		ArrayList<Integer> index = new ArrayList<Integer>(1);
		int nRows,iterSize,id;
		CyNode mapRowNodes[];
		Vector <Integer>ids;
		Map<Integer, Integer> transfac = new HashMap<Integer, Integer>();
		ArrayList<NodePair> removedNodes = new ArrayList<NodePair>();
		ArrayList<String> geneNames;
		
		index.add(1);
		
		
		geneNames = new ArrayList<String>();
		
        //step = 1.0 /  attributeArray.size();
		
		taskMonitor.setTitle("ARACNE Algorithm");
        
        taskMonitor.setStatusMessage("ARACNE Algorithm: Initializing Aracne data ...");
		taskMonitor.setProgress(progress);
		
		//Create new network
		CyNetwork newNetwork = netFactory.createNetwork();
		
		
		//Check if a network is associated to the selected table
		networkSelected = netUtils.getNetworkAssociatedToTable(mytable);
		
		Map<String,Object> params = new HashMap<String,Object>();
		params.put("Type", algorithm);
		selectedMetric.setParameters(params);
		
		// Create the CyniTable
		CyniTable data = selectedMetric.getCyniTable(mytable,attributeArray.toArray(new String[0]), false, false, selectedOnly);
		
		params.put("KernelWidth", bandwith);
		params.put("MiSteps", miSteps);
		params.put("Size", data.nColumns());
		selectedMetric.setParameters(params);
		
		ids = new Vector<Integer>();
		
		if(data.hasAnyMissingValue())
		{
			SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					JOptionPane.showMessageDialog(null, "The data selected contains missing values.\n " +
							"Therefore, this algorithm can not proceed with these conditions.\n" +
							"Please, use one of the imputation data algorithms to estimate the missing values.", "Warning", JOptionPane.WARNING_MESSAGE);
				}
			});
			newNetwork.dispose();
			return;
		}
		
		if(!mode.matches(AracneAlgorithmContext.MODE_PREPROCESSING))
		{
			if(fileHub != null)
			{
				taskMonitor.setStatusMessage("ARACNE Algorithm: Reading Hub List ...");	
				CyRow selectedRow = null;
				if(readProbeList(fileHub.toString(),geneNames) != -1)
				{
					for(int i=0;i<geneNames.size();i++)
					{
						for(CyRow row : mytable.getAllRows())
						{
							if(row.get(mapCol, String.class).matches(geneNames.get(i)))
							{
								selectedRow = row;
								break;
							}
						}
						if(selectedRow != null)
						{
							id = data.getRowIndex(selectedRow.getRaw(mytable.getPrimaryKey().getName()));
							selectedRow = null;
							if(id != -1)
								ids.add(id);
						}
					}
				}
			}
			
			if(fileTFList != null)
			{
				taskMonitor.setStatusMessage("ARACNE Algorithm: Reading Transcription Factor List ...");	
				CyRow selectedRow = null;
				if(readProbeList(fileTFList.toString(),geneNames) != -1)
				{
					for(int i=0;i<geneNames.size();i++)
					{
						for(CyRow row : mytable.getAllRows())
						{
							if(row.get(mapCol, String.class).matches(geneNames.get(i)))
							{
								selectedRow = row;
								break;
							}
						}
						if(selectedRow != null)
						{
							id = data.getRowIndex(selectedRow.getRaw(mytable.getPrimaryKey().getName()));
							selectedRow = null;
							if(id != -1)
								transfac.put(id, 1);
						}
					}
				}
			}
		}
		
		//Set the name of the network, another name could be chosen
		networkName = "ARACNE Inference " + iteration;
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
		
		
		mapRowNodes = new CyNode[nRows];
		Arrays.fill(mapRowNodes, null);
		
		taskMonitor.setStatusMessage("ARACNE Algorithm: Statistics calculated ...");	
		
		
		if(mode.matches(AracneAlgorithmContext.MODE_PREPROCESSING) || mode.matches(AracneAlgorithmContext.MODE_COMPLETE))
        {
			taskMonitor.setStatusMessage("ARACNE Algorithm: Calculating Thresholds ...");
        	AracneCyniTable data2 = new AracneCyniTable(mytable,attributeArray.toArray(new String[0]), false, false, selectedOnly,algorithm);
        	/*data2.computeMarkerVariance();
            data2.computeBandwidth();
            data2.computeMarkerRanks();*/
            
           /* if(algorithm.equals(Mutual_Info.ALGORITHM.FIXED_BANDWIDTH) || algorithm.equals(Mutual_Info.ALGORITHM.ADAPTIVE_PARTITIONING))
            	data2.addNoise();*/
            
            thParams = new Thresholds (data2,selectedMetric,threshold_file,kernel_file,aracne_path ,taskMonitor);
            progress = 0.1;
            taskMonitor.setProgress(0.1);
        }
	
		if (cancelled)
			return;
		if (mode.matches(AracneAlgorithmContext.MODE_PREPROCESSING)) {
			thParams.generateKernelWidthConfiguration(algorithm);
			if (cancelled)
				return;
            findKernelWidth(data.nColumns());
            params.put("KernelWidth", bandwith);
			selectedMetric.setParameters(params);
			selectedMetric.initMetric();
            thParams.generateMutualInformationThresholdConfiguration(algorithm);
            taskMonitor.setStatusMessage("ARACNE Algorithm: Saving Thresholds ...");
            newNetwork.dispose();
            return;
        } else if (mode.matches(AracneAlgorithmContext.MODE_COMPLETE)) {
        	selectedMetric.initMetric();
        	thParams.generateKernelWidthConfiguration(algorithm);
        	if (cancelled)
				return;
        	thParams.generateMutualInformationThresholdConfiguration(algorithm);
        	progress = 0.3;
        	taskMonitor.setProgress(progress);
        }
        else
        {
        	if(!manualKernel)
    			findKernelWidth(data.nColumns());
        }
		
		if(usePValue)
			chosenTh = findThreshold(data.nColumns());
		else
			chosenTh = miTh;
		
		params.clear();
		params.put("KernelWidth", bandwith);
		selectedMetric.setParameters(params);
		selectedMetric.initMetric();
		
		if(ids.isEmpty())
		{
			iterSize = nRows;
			step = 1.0 / nRows;
			matrix = new MatrixEdgePair(nRows,chosenTh);
		}
		else
		{
			iterSize = ids.size();
			step = 1.0 / iterSize;
			matrix = new MatrixEdgePair(chosenTh,ids);
		}
		
		// Create the thread pools
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		
		taskMonitor.setStatusMessage("ARACNE Algorithm: Calculating MI for all possible pairs ...");

		for (int i = 0; i < iterSize; i++) 
		{
			int j;
			
			if (cancelled)
				break;
					
			if(ids.isEmpty())
			{
				id = i;
				j = i+1;
			}
			else
			{
				id = ids.get(i);
				j = 0;
			}

			for (; j < nRows; j++) 
			{
						
				if (cancelled)
					break;
						
				index.set(0, j);
				executor.execute(new ThreadedGetMetric(data,id,index));
						
						
			}
			executor.shutdown();
			// Wait until all threads are finish
			try {
			   	executor.awaitTermination(7, TimeUnit.DAYS);
			} catch (Exception e) {}
					
					
			progress = progress + step;
			taskMonitor.setProgress(progress);
			executor = Executors.newFixedThreadPool(nThreads);
		}
		
		if(dpiTol < 1)
		{
			taskMonitor.setStatusMessage("ARACNE Algorithm: Applying DPI ...");
		
			Vector<ArrayValuePair> miVector = new Vector<ArrayValuePair>();
	        // put all the values of the i-th row in the matrix on a vector for sorting
	        // the ArrayValuePair structure is borrowed from the Microarray_Set class,
	        // it should actually be called EdgeMIPair here. They share the same
	        // structure consisting of two fields: a int and a double. The assoiciated
	        // sorting class will sort by the double fields.
	        for (int i = 0; i < iterSize; i++) {
	        	miVector.clear();
	        	if (cancelled)
					break;
	        	if(ids.isEmpty())
					id = i;
				else
					id = ids.get(i);
					
	        	 
	        	for (int j = 0; j < nRows; j++) 
				{
	        		//if(edgePresence[i][j])
	        		if(matrix.getPresence(id, j))
	        			miVector.add(new ArrayValuePair(j,matrix.getScore(id, j)));
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
	                    double valueBC = matrix.getScore(geneId1, geneId2);//MIScore[ geneId1][ geneId2];
	                    if (valueBC > minMI && matrix.getPresence(geneId1, geneId2)) {
	                    	// if TF annotation information is provided, the triangle will be
	                        // broken only if certein logic is met
	                        if (!(transfac.isEmpty())) {
	                            if (protectedByTFLogic(transfac, id, geneId1, geneId2)) {
	                                continue;
	                            }
	                        }

	                        removedNodes.add(new NodePair(id, geneId1));
	                        	
	                        break;
	                    }
	                }
	            }
	
	        }
	        
	        if(removedNodes.size() > 0)
	        {
	        	for(NodePair entry : removedNodes)
	        		matrix.removeScore(entry.node1,entry.node2);
	        	
	        	removedNodes.clear();
	        }

		}
		
		if (!cancelled)
		{
			taskMonitor.setStatusMessage("ARACNE Algorithm: Building network view ...");
			for (int i = 0; i < iterSize; i++) 
			{
				int j;
				
				if (cancelled)
					break;
						
				if(ids.isEmpty())
				{
					id = i;
					j = i+1;
				}
				else
				{
					id = ids.get(i);
					j = 0;
				}

				for( ; j < nRows; j++) 
				{
					
					if (cancelled)
						break;
							
					if(matrix.getPresence(id, j))
					{
						if(mapRowNodes[id] == null)
						{
							node1 = newNetwork.addNode();
							netUtils.cloneNodeRow(newNetwork,mytable.getRow(data.getRowLabel(id)),node1);
							if(newNetwork.getRow(node1).get(CyNetwork.NAME,String.class ) == null || newNetwork.getRow(node1).get(CyNetwork.NAME,String.class ).isEmpty() == true)
							{
								if(mytable.getPrimaryKey().getType().equals(String.class) && networkSelected == null)
									newNetwork.getRow(node1).set(CyNetwork.NAME,mytable.getRow(data.getRowLabel(id)).get(mytable.getPrimaryKey().getName(),String.class));
								else
									newNetwork.getRow(node1).set(CyNetwork.NAME, "Node " + numNodes);
							}
							if(newNetwork.getRow(node1).get(CyNetwork.SELECTED,Boolean.class ) == true)
								newNetwork.getRow(node1).set(CyNetwork.SELECTED, false);
							mapRowNodes[id] =node1;
							numNodes++;
						}
						if(mapRowNodes[j] == null)
						{
							node2 = newNetwork.addNode();
							netUtils.cloneNodeRow(newNetwork,mytable.getRow(data.getRowLabel(j)), node2);
							if(newNetwork.getRow(node2).get(CyNetwork.NAME,String.class ) == null || newNetwork.getRow(node2).get(CyNetwork.NAME,String.class ).isEmpty() == true)
							{
								if(mytable.getPrimaryKey().getType().equals(String.class) && networkSelected == null)
									newNetwork.getRow(node2).set(CyNetwork.NAME,mytable.getRow(data.getRowLabel(j)).get(mytable.getPrimaryKey().getName(),String.class));
								else
									newNetwork.getRow(node2).set(CyNetwork.NAME, "Node " + numNodes);
							}
							if(newNetwork.getRow(node2).get(CyNetwork.SELECTED,Boolean.class ) == true)
								newNetwork.getRow(node2).set(CyNetwork.SELECTED, false);
							mapRowNodes[j] = node2;
							numNodes++;
						}
										
						if(!newNetwork.containsEdge(mapRowNodes[id], mapRowNodes[j]))
						{
							edge = newNetwork.addEdge(mapRowNodes[id], mapRowNodes[j], false);
							newNetwork.getRow(edge).set("Mutual Information",matrix.getScore(id, j));
							newNetwork.getRow(edge).set(CyEdge.INTERACTION,((Double)matrix.getScore(id, j)).toString());
							newNetwork.getRow(edge).set("Metric",selectedMetric.toString());
							newNetwork.getRow(edge).set("name", newNetwork.getRow(mapRowNodes[id]).get("name", String.class)
									+ " (Aracne) " + newNetwork.getRow( mapRowNodes[j]).get("name", String.class));
						}
					}
				}
			}
			
		}
		
		
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
    void findKernelWidth(int n) {
        File f = new File(kernel_file);
        //f.deleteOnExit();
        if (f.exists()) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(f));
                String line = null;
                while ((line = br.readLine()) != null) {
                    if (!line.startsWith(">")) {
                        String[] tokens = line.split("\t");
                        KERNEL_ALPHA = Double.parseDouble(tokens[0].trim());
                        KERNEL_BETA = Double.parseDouble(tokens[1].trim());
                    }
                }
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }
        bandwith = KERNEL_ALPHA * Math.pow(n, KERNEL_BETA);
    }
    
    double findThreshold(int n) {
        File f = new File(threshold_file);
        //f.deleteOnExit();
        if (f.exists()) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(f));
                String line = null;
                while ((line = br.readLine()) != null) {
                    if (!line.startsWith(">")) {
                        String[] tokens = line.split("\t");
                        THRESHOLD_ALPHA = Double.parseDouble(tokens[0].trim());
                        THRESHOLD_BETA = Double.parseDouble(tokens[1].trim());
                        THRESHOLD_GAMMA = Double.parseDouble(tokens[2].trim());
                    }
                }
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }
       return ((THRESHOLD_ALPHA - Math.log(pvalue)) / ((-THRESHOLD_BETA) + (-THRESHOLD_GAMMA) * n));
    }
	
	
	static int readProbeList(String infile, List<String> probeList) {
        int count = 0;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(infile));

            String line = reader.readLine();
            while (line != null) {
                probeList.add(line);
                count++;
                line = reader.readLine();
            }
        } catch (IOException e) {
            System.err.println("Problem reading probe list file " + infile);
            e.printStackTrace();
            return -1;
        }
        return count;
    }
	
	boolean protectedByTFLogic(Map<Integer, Integer> transfac, int geneId1, int geneId2, int geneId3) {
        boolean isA = (transfac.get(geneId1) != null);
        boolean isB = (transfac.get(geneId2) != null);
        boolean isC = (transfac.get(geneId3) != null);
        
        if ( isC ) {
            return false;
        } else if (!(isA || isB)) {
            return false;
        }
        return true;
    }
	
	private static String join(String separator, String... parts) {
		StringBuilder builder = new StringBuilder();
		boolean isFirst = true;
		for (String part : parts) {
			if (!isFirst) {
				builder.append(separator);
			} else {
				isFirst = false;
			}
			builder.append(part);
		}
		return builder.toString();
	}
	
	private class ThreadedGetMetric implements Runnable {
		private ArrayList<Integer> index2;
		private int index1;
		private CyniTable tableData;
		
		ThreadedGetMetric(CyniTable data,int index1, ArrayList<Integer> parentsToIndex)
		{
			this.index2 = new ArrayList<Integer>( parentsToIndex);
			this.index1 = index1;
			this.tableData = data;
			
		}
		
		public void run() {
			matrix.setScore(index1,index2.get(0), selectedMetric.getMetric(tableData, tableData, index1, index2));

		}
		

	}
	
	class MatrixEdgePair{
		Map<Integer, Map<Integer,Double>> mi;
		double th;
		int size;
		boolean simetric;
		
		MatrixEdgePair(int size, double threshold)
		{
			mi = new ConcurrentHashMap<Integer,  Map<Integer,Double>>();
			th = threshold;
			simetric = true;
			this.size = size;
			for(int i=0; i<size;i++)
				mi.put(i, new ConcurrentHashMap<Integer,Double>());
		}
		MatrixEdgePair( double threshold, Vector<Integer> hubIds)
		{
			mi = new ConcurrentHashMap<Integer,  Map<Integer,Double>>();
			th = threshold;
			simetric = false;
			size = Collections.max(hubIds);
			for(int i=0; i<hubIds.size();i++)
				mi.put(hubIds.get(i), new ConcurrentHashMap<Integer,Double>());
		}
		
		boolean getPresence(int row, int col)
		{
			if(row == col)
				return false;
			
			if(getScore(row,col)> th)
				return true;
			else
				return false;
		}
		
		double getScore(int row, int col)
		{
			int i,j;
			double score = -1;
			Map<Integer,Double> temp;
			if(row == col)
				return -1.0;
			if(simetric)
			{
				i= Math.min(row, col);
				j = Math.max(row, col);
			}
			else
			{
				i = row;
				j = col;
			}
			temp = mi.get(i);
			
			if(temp != null)
			{
				if(temp.get(j) != null)
					score = temp.get(j);
			}
			
			if(!simetric && score == -1)
			{
				i = col;
				j = row;
				
				temp = mi.get(i);
				
				if(temp != null)
				{
					if(temp.get(j) != null)
						score = temp.get(j);
				}
				
			}
			
			return score;
		}
		
		void setScore(int row, int col, double score)
		{
			int i,j;
			Map<Integer,Double> temp;
			if(row == col)
				return;
			
			if(score < th)
				return;
			if(simetric)
			{
				i= Math.min(row, col);
				j = Math.max(row, col);
				if(i > size || j > size)
					return;
			}
			else
			{
				i = row;
				j = col;
				if(i > size )
					return;
			}
			
			temp = mi.get(i);
			temp.put(j,score);
			
			
		}
		void removeScore(int row, int col)
		{
			int i,j;
			Map<Integer,Double> temp;
			if(row == col)
				return;
			
			if(simetric)
			{
				i= Math.min(row, col);
				j = Math.max(row, col);
				if(i > size || j > size)
					return;
			}
			else
			{
				i = row;
				j = col;
				if(i > size )
					return;
			}
			temp = mi.get(i);
					
			if(temp.get(j) != null)
				temp.remove(j);

		}
	}
	
	static class NodePair {
	
		public int node1,node2;
		NodePair(int node1,int node2)
		{
			this.node1 = node1;
			this.node2 = node2;
			
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

    @Override
	public void cancel() {
		cancelled = true;
		if(thParams != null)
			thParams.setCancel();
		
	}
	
}
