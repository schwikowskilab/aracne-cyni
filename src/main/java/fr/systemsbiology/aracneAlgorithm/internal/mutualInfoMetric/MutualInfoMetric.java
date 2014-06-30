/*
  File: CyniSampleMetric.java

  Copyright (c) 2010-2012, The Cytoscape Consortium (www.cytoscape.org)

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
package fr.systemsbiology.aracneAlgorithm.internal.mutualInfoMetric;

import fr.systemsbiology.aracneAlgorithm.internal.AracneCyniTable;
import fr.systemsbiology.cyni.*;
import org.cytoscape.model.CyTable;

import java.util.*;


/**
 * The BasicInduction provides a very simple Induction, suitable as
 * the default Induction for Cytoscape data readers.
 */
public class MutualInfoMetric extends AbstractCyniMetric {
	/**
	 * Creates a new  object.
	 */
	public MutualInfoMetric() {
		super("MIAracneMetric","Continous Mutual Information Metric");
		addTag(CyniMetricTags.INPUT_NUMBERS.toString());
		mi = null;
		type =  Mutual_Info.ALGORITHM.FIXED_BANDWIDTH;
		size = 1;
		kernelWidth = 0;
	}
	
    private  int miSteps = 6;
    private double kernelWidth;
    private int size;
    private Mutual_Info.ALGORITHM type;
    private Mutual_Info mi;    
	
	public Double getMetric(CyniTable table1, CyniTable table2, int indexBase, List<Integer> indexToCompare) { 
		double result = 0.0;
		AracneCyniTable aracneTable1 = null;
		AracneCyniTable aracneTable2 = null;
		
		int index2 = indexToCompare.get(0);
		
		
		size = Math.min(table1.nColumns(), table2.nColumns());
        
        
		if(table1 instanceof AracneCyniTable)
			aracneTable1 = (AracneCyniTable) table1;
		if(table2 instanceof AracneCyniTable)
			aracneTable2 = (AracneCyniTable) table2;
       
		if(aracneTable2 != null && aracneTable2 != null && mi != null)
			result = mi.Compute_Pairwise_MI(aracneTable1.getRankedValues(indexBase), aracneTable2.getRankedValues(index2), aracneTable1.getBandwith(indexBase), aracneTable2.getBandwith(index2));
		
		return  result;
	}

	public void setParameters(Map<String,Object> params){
		
		if(params.containsKey("KernelWidth"))
			kernelWidth = (Double) params.get("KernelWidth");
		
		if(params.containsKey("MiSteps"))
			miSteps = (Integer) params.get("MiSteps");
		
		if(params.containsKey("Size"))
			size = (Integer) params.get("Size");
		
		if(params.containsKey("Type"))
		{
			type =  Enum.valueOf(Mutual_Info.ALGORITHM.class,params.get("Type").toString());
			
		}
			//type = (Mutual_Info.ALGORITHM) params.get("Type");
		
		
	}
	
	@Override
	public  CyniTable getCyniTable( CyTable table, String[] attributes, boolean transpose, boolean ignoreMissing, boolean selectedOnly)
	{
		return  new AracneCyniTable(table,attributes, false, false, selectedOnly,type);
	}
	
	public void initMetric()
	{
		mi = new Mutual_Info(size,miSteps,kernelWidth,type);
	}
	

	
}
