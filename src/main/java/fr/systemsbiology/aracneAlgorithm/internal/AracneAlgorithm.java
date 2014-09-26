package fr.systemsbiology.aracneAlgorithm.internal;


import java.util.*;

import fr.systemsbiology.cyni.*;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkViewFactory;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TunableSetter;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNetworkTableManager;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.vizmap.VisualMappingManager;



public class AracneAlgorithm extends AbstractCyniAlgorithm {
	
	
	private CyTable selectedTable;
	/**
	 * Creates a new Inference Algorithm object.
	 */
	public AracneAlgorithm() {
		super("aracne","ARACNE Algorithm",true,CyniCategory.INDUCTION);
	
	}

	public TaskIterator createTaskIterator(CyniAlgorithmContext context, CyTable table,CyNetworkFactory networkFactory, CyNetworkViewFactory networkViewFactory,
			CyNetworkManager networkManager,CyNetworkTableManager netTableMgr, CyRootNetworkManager rootNetMgr,VisualMappingManager vmMgr,
			CyNetworkViewManager networkViewManager, CyLayoutAlgorithmManager layoutManager, CyCyniMetricsManager metricsManager) {
			selectedTable = table;
			return new TaskIterator(new AracneAlgorithmTask(getName(),(AracneAlgorithmContext) context,networkFactory,networkViewFactory,
					networkManager,netTableMgr,rootNetMgr,vmMgr,networkViewManager,layoutManager,metricsManager, selectedTable));
	}
	
	public CyniAlgorithmContext createCyniContext(CyTable table, CyCyniMetricsManager metricsManager, TunableSetter tunableSetter,Map<String, Object> mparams) {
		CyniAlgorithmContext context;
		selectedTable = table;
		context = new AracneAlgorithmContext(selectedTable);
		if(mparams != null && !mparams.isEmpty())
			tunableSetter.applyTunables(context, mparams);
		return context;
	}
	
}
