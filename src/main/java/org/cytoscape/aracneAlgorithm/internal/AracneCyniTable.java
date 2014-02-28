package org.cytoscape.aracneAlgorithm.internal;

import org.cytoscape.cyni.*;
import org.cytoscape.model.CyTable;

import java.util.*;
import org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric.*;
import org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric.Mutual_Info.ALGORITHM;



public class AracneCyniTable extends CyniTable {
	
	private double bandwith[];
	private Vector<Vector<Gene>> rankedValues;
	private double variance[];


	public AracneCyniTable(  CyTable table, String[] attributes, boolean transpose, boolean ignoreMissing, boolean selectedOnly,ALGORITHM type) {
		super(table,attributes,transpose,ignoreMissing,  selectedOnly);
		rankedValues = new Vector<Vector<Gene>>(this.nRows());
		if(!hasAnyMissingValue())
		{
			computeMarkerVariance();
			computeBandwidth();
			computeMarkerRanks();
			if(type.equals(Mutual_Info.ALGORITHM.FIXED_BANDWIDTH) || type.equals(Mutual_Info.ALGORITHM.ADAPTIVE_PARTITIONING))
	        	addNoise();
		}
	}
	
	public AracneCyniTable(CyniTable table)
	{
		super(table);
		rankedValues = new Vector<Vector<Gene>>(this.nRows());
		computeMarkerVariance();
		computeBandwidth();
		computeMarkerRanks();
	}
	
	public void computeBandwidth( ) {
		
		int n = this.nColumns();
		bandwith = new double[this.nRows()];
        Vector<Double> data = new Vector<Double>(n);
        for (int i = 0; i < this.nRows(); i++) {

            double prop = 1.06; // Gaussian
            int dim = 1; // Dimension of data

            double std_dev = Math.sqrt(variance[i]);
            data.clear();
            for (int j = 0; j < n; j++) {
                data.add((double) doubleValue(i, j));
            }

            Collections.sort(data);

            double iqr = Util.interQuartileRange(data, n);
            double iqrSig = .7413 * iqr; // find interquartile range sigma est.

            if (iqrSig == 0) {
                iqrSig = std_dev;
            }
            double sig = Math.min(std_dev, iqrSig);

            // Computing bandwidth
            double h = prop * sig * Math.pow((double) n, ((double) -1) / ((double) (4 + dim)));
            bandwith[i] = h;
        }

	}
	
	public double getBandwith(int idx)
	{
		return bandwith[idx];
	}
	
	public double getVariance(int idx)
	{
		return variance[idx];
	}
	
	public void computeMarkerRanks() {
		int n = this.nColumns();
		
        for (int i = 0; i < this.nRows(); i++) {

            Vector<Gene> data = new Vector<Gene>(n);
            for (int j = 0; j < n; j++) {
                Gene g = new Gene(j,  doubleValue(i, j));
                data.add(g);
            }
            Collections.sort(data, new Sort_Gene());
            
            for (int j = 0; j < n; j++)
                data.get(j).xi = j;
            rankedValues.add(i, data) ;
        }
    }
	
	public Vector<Gene> getRankedValues(int idx)
	{
		return rankedValues.get(idx);
	}
    
    public void computeMarkerVariance() {
      
    	variance = new double[this.nRows()];
        for (int i = 0; i < this.nRows(); i++) {
             variance[i] = calVariance(i);
            
        }
    }
    
    public void setVariance(double variance,int idx) {
        this.variance[idx] = variance;
    }

    public void setBandwidth(double bandwidth,int idx) {
        this.bandwith[idx] = bandwidth;
    }
    
    public void setRankedValues(Vector<Gene> r,int idx) {
        rankedValues.set(idx, r);
    }
    
    public void addNoise()
    {
    	Random rng = new Random(System.currentTimeMillis());
      
        for (int id = 0; id < nRows(); id++) {
            
            for (int mid = 0; mid < nColumns(); mid++) {
            	if(hasValue(id,mid))
            	{
            		double noise = rng.nextDouble() * 1e-10;
            		setValue(id,mid,(Double)(doubleValue(id,mid)+noise));
            	}
            }
        }
    }
    
    public class Sort_Gene implements Comparator<Gene> {
        
        public int compare(Gene a, Gene b) {
            if (a.x != b.x) {
                if (a.x < b.x) {
                    return -1;
                } else {
                    return 1;
                }
            } else {
                return (int) Math.signum(a.maId - b.maId);
            }
        }
    }
    
    double calVariance(int m) {

        double var;
        int n = this.nColumns();
        double s = 0;
        double ss = 0;
        for (int i = 0; i < n; i++) {
            double v = doubleValue(m, i);
            s += v;
            ss += v * v;
        }
        var = (ss - s * s / n) / (n - 1);
        return var;
    }

    

}


