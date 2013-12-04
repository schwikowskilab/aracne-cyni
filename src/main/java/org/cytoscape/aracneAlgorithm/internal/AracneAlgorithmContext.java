package org.cytoscape.aracneAlgorithm.internal;

import java.io.*;
import java.util.List;

import org.cytoscape.cyni.CyniAlgorithmContext;
import org.cytoscape.model.CyTable;
import org.cytoscape.work.TunableValidator.ValidationState;
import org.cytoscape.work.util.*;

import org.cytoscape.work.Tunable;
import org.cytoscape.work.TunableValidator;

public class AracneAlgorithmContext extends CyniAlgorithmContext implements TunableValidator {
	
	
	@Tunable(description="Mutual Information Algorithm Type:",groups="Algorithm Definition")
	public ListSingleSelection<String> algoChooser = new ListSingleSelection<String>("Naive Bayes","Adaptive Partitioning","Fixed Bandwith","Variable Bandwith");
	
	
	@Tunable(description="Manual Kernel Width Definition",dependsOn="algoChooser=Fixed Bandwith",groups="Algorithm Definition")
	public boolean manualKernel = false;
	
	@Tunable(description="Kernel Width(0,1): ",dependsOn="manualKernel=true",groups="Algorithm Definition")
	public double kernelWidth = 0.0;
	
	@Tunable(description="DPI Tolerance[0,1]: ",groups="Algorithm Definition")
	public double dpiTol = 0.0;
	
	@Tunable(description="Mutual Information Steps: ",groups="Algorithm Definition")
	public int miSteps = 6;
	
	@Tunable(description="Hub Genes File ",groups="Hub/Transcription Factor Definition",params="input=true")
	public File hubFile ;
	
	@Tunable(description="Transcription Factor List: ",groups="Hub/Transcription Factor Definition",params="input=true")
	public File TFFile ;
	
	@Tunable(description="Gene/TF column name mapping:",groups="Hub/Transcription Factor Definition")
	public ListSingleSelection<String> colMapping ;
	
	@Tunable(description="Which threshold to use:",groups="Threshold Definition", xorChildren=true)
	public ListSingleSelection<String> thresholdChooser = new ListSingleSelection<String>("MI Threshold","P-Value Threshold");
	
	@Tunable(description="Mutual Information Threshold:",groups={"Threshold Definition","MI Threshold Definition"},xorKey="MI Threshold")
	public double miTh = 0.5;
	
	@Tunable(description="P-Value Threshold(0,1]:",groups={"Threshold Definition","P-Value Threshold Definition"},xorKey="P-Value Threshold")
	public double pvalue = 0.5;
	
	@Tunable(description="Data Attributes", groups="Sources for Network Inference")
	public ListMultipleSelection<String> attributeList;

	
	private List<String> attributes;
	

	public AracneAlgorithmContext(CyTable table ) {
		super(true);
		attributes = getAllAttributesNumbers(table);
		if(attributes.size() > 0)
		{
			attributeList = new  ListMultipleSelection<String>(attributes);
		}
		else
		{
			attributeList = new  ListMultipleSelection<String>("No sources available");
		}
		algoChooser.setSelectedValue("Naive Bayes");
		colMapping =  new  ListSingleSelection<String>(getAllAttributesStrings(table));
		if(table.getPrimaryKey().getType() == String.class)
			colMapping.setSelectedValue(table.getPrimaryKey().getName());
	}
	
	@Override
	public ValidationState getValidationState(final Appendable errMsg) {
		if(attributeList.getPossibleValues().get(0).matches("No sources available") || attributeList.getSelectedValues().size() == 0) {
			try {
				errMsg.append("No sources selected to apply the algorithm or there are no available. Please, select sources from the list if available.");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
			
		}
		
		if ((kernelWidth <= 0 || kernelWidth >= 1) && manualKernel)
		{
			try {
				errMsg.append("Kernel Width needs to be between 0 and 1 (0,1).");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
		}
		if (dpiTol < 0 ||dpiTol > 1) 
		{
			try {
				errMsg.append("DPI tolerance must be within [0,1]!");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
		}
		if (miSteps <= 0) 
		{
			try {
				errMsg.append("Number of mutual information steps must be postive!");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
		}
		if ((pvalue <= 0 || pvalue > 1)  && thresholdChooser.getSelectedValue().matches("P-Value Threshold"))
		{
			try {
				errMsg.append("P-value must be in the range (0,1]!");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
		}
		if (  miTh < 0 && thresholdChooser.getSelectedValue().matches("MI Threshold"))
		{
			try {
				errMsg.append("MI threshold must be non-negative!");
			} catch (IOException e) {
				e.printStackTrace();
				return ValidationState.INVALID;
			}
			return ValidationState.INVALID;
		}
		return ValidationState.OK;
	}
}
