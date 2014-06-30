/*
  File: CyniSampleAlgorithm.java

  Copyright (c) 2010, The Cytoscape Consortium (www.cytoscape.org)

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



/**
 * The BasicInduction provides a very simple Induction, suitable as
 * the default Induction for Cytoscape data readers.
 */
public class AracneAlgorithm extends AbstractCyniAlgorithm {
	
	
	private CyTable selectedTable;
	/**
	 * Creates a new EqualDiscretization object.
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
