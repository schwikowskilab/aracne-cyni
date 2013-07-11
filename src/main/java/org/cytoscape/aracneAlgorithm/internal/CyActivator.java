package org.cytoscape.aracneAlgorithm.internal;

import org.cytoscape.application.swing.CySwingApplication;

import org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric.*;

import org.cytoscape.application.swing.CyAction;
import org.cytoscape.cyni.*;

import org.osgi.framework.BundleContext;

import org.cytoscape.service.util.AbstractCyActivator;

import java.util.Properties;


public class CyActivator extends AbstractCyActivator {
	public CyActivator() {
		super();
	}


	public void start(BundleContext bc) {

		//Define new Cyni Algorithm
		AracneAlgorithm test = new AracneAlgorithm();
		//Register new Cyni Algorithm
		registerService(bc,test,CyCyniAlgorithm.class, new Properties());

		//Define new Cyni Metric
		MutualInfoMetric miMetric = new MutualInfoMetric();
		//Register new Cyni Algorithm
		registerService(bc,miMetric,CyCyniMetric.class, new Properties());

	}
}

