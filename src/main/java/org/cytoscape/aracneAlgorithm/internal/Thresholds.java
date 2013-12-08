package org.cytoscape.aracneAlgorithm.internal;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import java.io.*;
import java.util.*;

import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.MathException;
import org.apache.commons.math.optimization.ConvergenceChecker;
import org.apache.commons.math.optimization.CostException;
import org.apache.commons.math.optimization.CostFunction;
import org.apache.commons.math.optimization.NelderMead;
import org.apache.commons.math.optimization.PointCostPair;
import org.apache.commons.math.random.NotPositiveDefiniteMatrixException;
import org.apache.commons.math.special.Erf;
import org.apache.commons.math.stat.regression.SimpleRegression;
import org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric.*;
import org.cytoscape.cyni.CyCyniMetric;

/**
 * @author manjunath at c2b2 dot columbia dot edu
 */
public class Thresholds {
	
	public static double KERNEL_ALPHA = 0.52477;
    public static double KERNEL_BETA = -0.24;
    public static double THRESHOLD_ALPHA = 1.062;
    public static double THRESHOLD_BETA = -48.7;
    public static double THRESHOLD_GAMMA = -0.634;

    static double a =KERNEL_ALPHA;
    static double b = KERNEL_BETA;
    static cern.jet.math.Functions F = cern.jet.math.Functions.functions;
    
    private AracneCyniTable cyniTable;
    private CyCyniMetric selectedMetric;
    
    public Thresholds(AracneCyniTable table,CyCyniMetric selectedMetric)
    {
    	cyniTable = table;
    	this.selectedMetric = selectedMetric;
    }

    public void generateMutualInformationThresholdConfiguration(Mutual_Info.ALGORITHM algorithm) {

        int ng = cyniTable.nRows();
        int ma = cyniTable.nColumns();
        DoubleMatrix2D data = new DenseDoubleMatrix2D(ng, ma);
        System.out.println();
        System.out.println("Generating Mutual Information threshold configuration...");
        System.out.println("ma: " + ma);
        for (int i = 0; i < ma; i++) {
            for (int j = 0; j < ng; j++) {
                data.set(j, i, cyniTable.doubleValue( j,i));
            }
        }

        double[][] p = {{0.3}, {0.35}, {0.4}, {0.45}, {0.5}, {0.55}, {0.6}, {0.65}, {0.7}, {0.75}, {0.8}, {0.85}, {0.9}, {0.95}, {1.0}};
        DoubleMatrix1D n = new DenseDoubleMatrix1D(p.length);
        for (int i = 0; i < p.length; i++) {
            n.set(i, Math.floor(ma * p[i][0]));
        }

        int rep = 3;
        int N = 100000;
        DoubleMatrix2D alpha = new DenseDoubleMatrix2D(n.size(), rep);
        DoubleMatrix2D beta = new DenseDoubleMatrix2D(n.size(), rep);

        for (int i = 0; i < n.size(); i++) {
            System.out.println("Sample size = " + n.get(i));
            double h = 0;
            if (algorithm.equals( Mutual_Info.ALGORITHM.FIXED_BANDWIDTH)) {
                h = a * Math.pow(Math.round(n.get(i)), b);
            }
            for (int j = 0; j < rep; j++) {
                System.out.println("Repeat " + j);
                int[] idx = Shuffle.sampleWithoutReplacement(ma, (int) Math.round(n.get(i)));
                int[] rows = new int[ng];
                for (int k = 0; k < ng; k++) {
                    rows[k] = k;
                }
                double[] coef = extrapolateMIThreshold(data.viewSelection(rows, idx), h, N, algorithm);
                beta.set(i, j, coef[1]);
                alpha.set(i, j, coef[0]);
            }
            System.out.println("End Sample size = " + n.get(i));
        }
        DoubleMatrix1D meanBeta = new DenseDoubleMatrix1D(n.size());
        double meanAlpha = 0;
        for (int i = 0; i < n.size(); i++) {
            double mb = 0;
            for (int j = 0; j < rep; j++) {
                mb += beta.get(i, j);
                meanAlpha += alpha.get(i, j);
            }
            meanBeta.set(i, (mb / (double) rep));
        }

        meanAlpha /= (double) (n.size() * rep);

        SimpleRegression regression = new SimpleRegression();
        DoubleMatrix2D temp = new DenseDoubleMatrix2D(n.size(), 2);
        temp.viewColumn(1).assign(meanBeta);
        temp.viewColumn(0).assign(n);
        regression.addData(temp.toArray());
        double slope = regression.getSlope();
        double intercept = regression.getIntercept();
        //System.out.println("meanAlpha " + meanAlpha);
        //System.out.println("intercept " + intercept);
        //System.out.println("slope " + slope);

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(AracneAlgorithmTask.threshold_file));
            bw.write(">" + "ConfigData");
            bw.newLine();
            bw.write(meanAlpha + "\t" + intercept + "\t" + slope);
            bw.newLine();
            bw.flush();
            bw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    double[] extrapolateMIThreshold(DoubleMatrix2D data, double h, final int N, Mutual_Info.ALGORITHM algorithm) {
        int ng = data.rows();
        int ma = data.columns();
        ArrayList<Integer> index = new ArrayList<Integer>(1);
        DoubleMatrix2D rdata = new DenseDoubleMatrix2D(ng, ma);
        AracneCyniTable tempTable = cyniTable;
        index.add(1);
       /* MarkerStats stats = new MarkerStats();
        for (Marker marker : Aracne.stats.markerStats.keySet()) {
            MarkerStats.MarkerStat markerStat = new MarkerStats.MarkerStat();
            markerStat.setVariance(Aracne.stats.markerStats.get(marker).getVariance());
            stats.markerStats.put(marker, markerStat);
        }*/
        for (int i = 0; i < ng; i++) {
            int[] randperm = Shuffle.permutation(ma);
            rdata.viewRow(i).assign(data.viewRow(i).viewSelection(randperm));
        }

        computeMarkerBandwidth(rdata, tempTable);
        computeMarkerRanks(rdata, tempTable);

        int it = 0;
        final DoubleMatrix1D nullDistribution = new DenseDoubleMatrix1D(N);

        while (it < N) {
            int i = Math.round(Shuffle.rand.nextFloat() * (ng - 1));
            int j = Math.round(Shuffle.rand.nextFloat() * (ng - 1));
            if (i == j) {
                while (i == j) {
                    j = Math.round(Shuffle.rand.nextFloat() * (ng - 1));
                }
            }

            index.set(0, j);
           /* double sigmaX = stats.getMarkerStat(Aracne.data.getMarkers().getMarker(i)).getBandwidth();
            double sigmaY = stats.getMarkerStat(Aracne.data.getMarkers().getMarker(j)).getBandwidth();

            Vector<Gene> g1 = stats.getMarkerStat(Aracne.data.getMarkers().getMarker(i)).getRankedValues();
            Vector<Gene> g2 = stats.getMarkerStat(Aracne.data.getMarkers().getMarker(j)).getRankedValues();*/
            
            double mi = selectedMetric.getMetric(tempTable, tempTable, i, index);
            nullDistribution.set(it, mi);

           /* if (algorithm.equals(Mutual_Info.ALGORITHM.FIXED_BANDWIDTH)) {
            	double mi = selectedMetric.getMetric(cyniTable, cyniTable, i, index);
                //double mi = Aracne.mi.Get_Mutual_Info_XY_FIXED_BANDWIDTH(g1, g2);
                nullDistribution.set(it, mi);
            } else if (algorithm.equals(Mutual_Info.ALGORITHM.VARIABLE_BANDWIDTH)) {
                double mi = Aracne.mi.Get_Mutual_Info_XY_VARIABLE_BANDWIDTH(g1, g2, sigmaX, sigmaY);
                nullDistribution.set(it, mi);
            } else if (algorithm.equals(Mutual_Info.ALGORITHM.ADAPTIVE_PARTITIONING)) {
                double mi = Aracne.mi.Get_Mutual_Info_XY_ADAPTIVE_PARTITIONING(g1, g2);
                nullDistribution.set(it, mi);
            }*/
            it++;
        }

        nullDistribution.assign(nullDistribution.viewSorted());

        final DoubleMatrix1D idx = new DenseDoubleMatrix1D(N / 100);
        double[] indices = new double[N / 100];
        for (int k = 100; k < N; k += 100) {
            idx.set((k / 100) - 1, k);
            indices[(k / 100) - 1] = (k / 100) - 1;
        }

        DoubleMatrix2D cumf = new DenseDoubleMatrix2D(idx.size(), 2);
        cumf.viewColumn(0).assign(indices);
        cumf.viewColumn(0).assign(new DoubleFunction() {

            public double apply(double x) {
                return nullDistribution.get((int) idx.get((int) x));
            }
        });

        cumf.viewColumn(1).assign(indices);
        cumf.viewColumn(1).assign(new DoubleFunction() {

            public double apply(double x) {
                return Math.log(((double) N - idx.get((int) x)) / (double) N);
            }
        });

        DoubleMatrix2D temp = new DenseDoubleMatrix2D(idx.size(), 2);
        temp.viewColumn(0).assign(cumf.viewColumn(0));
        temp.viewColumn(1).assign(cumf.viewColumn(1));

        SimpleRegression regression = new SimpleRegression();
        regression.addData(temp.toArray());
        double[] coefs = new double[2];
        coefs[0] = regression.getIntercept();
        coefs[1] = regression.getSlope();
        return coefs;
    }

    void computeMarkerBandwidth(DoubleMatrix2D set, AracneCyniTable table) {
        int n = set.columns();
        for (int i = 0; i < set.rows(); i++) {
           
            double prop = 1.06; // Gaussian
            int dim = 1; // Dimension of data

            double std_dev = Math.sqrt(table.getVariance(i));
            Vector<Double> data = new Vector<Double>(n);
            for (int j = 0; j < n; j++) {
                data.add((double) set.get(i, j));
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
            table.setBandwidth(h,i);
        }
    }

    void computeMarkerRanks(DoubleMatrix2D set, AracneCyniTable table) {
        int n = set.columns();
        for (int i = 0; i < set.rows(); i++) {
         
            Vector<Gene> data = new Vector<Gene>(n);
            for (int j = 0; j < n; j++) {
                Gene g = new Gene(j, (double) set.get(i, j));
                data.add(g);
            }
            Collections.sort(data, new Sort_Gene());
            for (int j = 0; j < n; j++) {
                data.get(j).xi = j;
            }
            table.setRankedValues(data,i);
        }
    }

    static DoubleMatrix2D kmeans(DoubleMatrix2D data, int k) {
        int n = data.rows();
        int m = data.columns();

        int[] ri = Shuffle.sampleWithoutReplacement(n, k);
        DoubleMatrix2D c = new DenseDoubleMatrix2D(k, m);
        int[] cin = new int[m];
        for (int i = 0; i < m; i++) {
            cin[i] = i;
        }
        c.assign(data.viewSelection(ri, cin));

        int iter = 0;
        int max_iter = 50;
        int changes = 1;
        DoubleMatrix2D distances = new DenseDoubleMatrix2D(n, k);
        DoubleMatrix1D costFunction = new DenseDoubleMatrix1D(max_iter + 1);
        DoubleMatrix1D dataLabels = new DenseDoubleMatrix1D(n);

        while (iter < max_iter && (changes > 0)) {
            iter++;
            DoubleMatrix1D oldLabels = new DenseDoubleMatrix1D(n);
            oldLabels.assign(dataLabels);

            Algebra algebra = new Algebra();
            for (int i = 0; i < k; i++) {
                DoubleMatrix2D c_k = new DenseDoubleMatrix2D(n, m);
                for (int j = 0; j < n; j++) {
                    c_k.viewRow(j).assign(c.viewRow(i));
                }
                DoubleMatrix1D ones = new DenseDoubleMatrix1D(m);
                ones.assign(1d);
                distances.viewColumn(i).assign(algebra.mult(data.copy().assign(c_k, F.minus).assign(F.square), ones));
            }

            DoubleMatrix1D minValues = new DenseDoubleMatrix1D(n);
            for (int i = 0; i < n; i++) {
                final double mv = distances.viewRow(i).viewSorted().get(0);
                for (int j = 0; j < distances.viewRow(i).size(); j++) {
                    if (distances.viewRow(i).get(j) == mv) {
                        dataLabels.set(i, j);
                    }
                }
            }

            costFunction.set(iter - 1, minValues.zSum());

            for (int i = 0; i < k; i++) {
                int[] indices = new int[n];

                for (int j = 0; j < n; j++) {
                    if (dataLabels.get(j) == i) {
                        indices[j] = j;
                    }
                }

                int[] mIndices = new int[m];
                for (int j = 0; j < m; j++) {
                    mIndices[j] = j;
                }
                DoubleMatrix2D mat = data.viewSelection(indices, mIndices);
                for (int j = 0; j < m; j++) {
                    c.set(i, j, (double) mat.viewColumn(j).zSum() / (double) n);
                }
            }
            changes = 0;
            for (int i = 0; i < n; i++) {
                if (oldLabels.get(i) != dataLabels.get(i)) {
                    changes++;
                }
            }
        }
        return c;
    }

    public void generateKernelWidthConfiguration(Mutual_Info.ALGORITHM algorithm) {
        if (algorithm.equals(Mutual_Info.ALGORITHM.FIXED_BANDWIDTH)) {
            int narray = cyniTable.nColumns();
            int ngene = cyniTable.nRows();
            double[] p = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
            int[] n = new int[p.length];
            System.out.println();
            System.out.println("Generating Kernel Width configuration...");
            System.out.println("narray: " + narray);
            for (int i = 0; i < p.length; i++) {
                n[i] = (int) Math.floor(narray * p[i]);
            }

            DoubleMatrix2D results = new DenseDoubleMatrix2D(n.length, 3);
            for (int i = 0; i < n.length; i++) {
                System.out.println("Sample size = " + n[i]);
                for (int j = 0; j < 3; j++) {
                    System.out.println("Repeat: " + j);
                    int[] idx = Util.Shuffle.sampleWithoutReplacement(narray, n[i]);
                    DoubleMatrix2D subset = new DenseDoubleMatrix2D(ngene, idx.length);
                    for (int k = 0; k < idx.length; k++) {
                        for (int l = 0; l < ngene; l++) {
                            subset.set(l, k, cyniTable.doubleValue(l,idx[k]));
                        }
                    }
                    double kw = findKernelWidth(subset, 1000, 0.1d);
                    results.set(i, j, kw);
                }
            }

            DoubleMatrix1D h_bar = new DenseDoubleMatrix1D(n.length);
            for (int i = 0; i < n.length; i++) {
                double meanH = 0d;
                for (int j = 0; j < 3; j++) {
                    meanH += results.get(i, j);
                }
                h_bar.set(i, Math.log(meanH / 3));
            }
            DoubleMatrix1D nMatrix = new DenseDoubleMatrix1D(n.length);
            for (int i = 0; i < n.length; i++) {
                nMatrix.set(i, Math.log(n[i]));
            }

            SimpleRegression regression = new SimpleRegression();
            DoubleMatrix2D temp = new DenseDoubleMatrix2D(n.length, 2);
            temp.viewColumn(1).assign(h_bar);
            temp.viewColumn(0).assign(nMatrix);
            regression.addData(temp.toArray());

            a = Math.exp(regression.getIntercept());
            b = regression.getSlope();
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(AracneAlgorithmTask.kernel_file));
                bw.write(">" + "ConfigData");
                bw.newLine();
                bw.write(a + "\t" + b);
                bw.newLine();
                bw.flush();
                bw.close();
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }
    }

    static double findKernelWidth(DoubleMatrix2D data, int n, double h0) {
        int ngene = data.rows();
        int narray = data.columns();
        DoubleMatrix2D rankTransformedData = new DenseDoubleMatrix2D(ngene, narray);
        DoubleMatrix2D copulaTransformedData = new DenseDoubleMatrix2D(ngene, narray);
        for (int i = 0; i < ngene; i++) {
            Vector<Gene> idata = new Vector<Gene>(narray);
            for (int j = 0; j < narray; j++) {
                Gene g = new Gene(j, data.get(i, j));
                g.maId = j;
                g.xi = j;
                idata.add(g);
            }
            Collections.sort(idata, new Sort_Gene());
            for (int j = 0; j < narray; j++) {
                rankTransformedData.set(i, idata.get(j).xi, j);
                copulaTransformedData.set(i, idata.get(j).xi, (j + 0.5) / (double) narray);
            }
        }

        int it = 1;

        DoubleMatrix2D h = new DenseDoubleMatrix2D(n, 1);

        //Create instance of NelderMead
        NelderMead min = new NelderMead();

        KernelWidthLOO kwloo = new KernelWidthLOO();
        kwloo.setIdata(rankTransformedData);
        kwloo.setTdata(copulaTransformedData);
        double maxH = 0;
        double sumH = 0;
        while (it <= n) {
            if (((it * 100) / (double) n) % 10 == 0) {
                //System.out.println(((it * 100) / n) + "% completed");
            }

            int g1 = Math.round(Shuffle.rand.nextInt(ngene - 1));
            int g2 = Math.round(Shuffle.rand.nextInt(ngene - 1));

            if (g1 == g2) {
                while (g1 == g2) {
                    g2 = Math.round(Shuffle.rand.nextInt(ngene - 1));
                }
            }

            kwloo.setG1(g1);
            kwloo.setG2(g2);

            KernelWidthConvergenceChecker kwcc = new KernelWidthConvergenceChecker();

            // simplex vertices
            double[][] vertices = {{-0.5}, {0.1}};

            //maximum number of iterations
            int nmax = 10000;

            PointCostPair pcp = null;
            try {
                // Nelder and Mead minimization procedure
                pcp = min.minimize(kwloo, nmax, kwcc, vertices, 3, System.currentTimeMillis());
            } catch (CostException ce) {
                ce.printStackTrace();
            } catch (ConvergenceException ce) {
                ce.printStackTrace();
            } catch (NotPositiveDefiniteMatrixException npdme) {
                npdme.printStackTrace();
            }

            // get the minimum value
            double minimum = pcp.getPoint()[0];
            maxH = (Math.exp(minimum) > maxH) ? Math.exp(minimum) : maxH;
            sumH += Math.exp(minimum);
            h.set(it - 1, 0, Math.exp(minimum));
            it++;
        }

        double h_opt = 0;
        if (maxH < 0.5) {
            h_opt = sumH / h.rows();
        } else {
            DoubleMatrix1D centers = kmeans(h, 5).viewColumn(0);
            h_opt = Double.MAX_VALUE;
            for (double d : centers.toArray()) {
                h_opt = (d < h_opt) ? d : h_opt;
            }
        }
        return h_opt;
    }

    public static class KernelWidthLOO implements CostFunction {

        DoubleMatrix2D idata;
        DoubleMatrix2D tdata;
        int g1, g2;

        public double cost(double[] args) {
            double logl = 0d;
            try {
                double h = (Math.abs(args[0]) == 0d) ? 0.0001 : Math.abs(args[0]);
                int m = tdata.columns();
                int[] c = new int[m];
                for (int i = 0; i < m; i++) {
                    c[i] = i;
                }
                DoubleMatrix2D id = idata.viewSelection(new int[]{g1, g2}, c);
                DoubleMatrix2D td = tdata.viewSelection(new int[]{g1, g2}, c);
                int n = id.rows();

                DoubleMatrix1D P_TABLE = new DenseDoubleMatrix1D(m);
                P_TABLE.assign(-1d);
                DoubleMatrix1D N_TABLE = new DenseDoubleMatrix1D(m);
                for (int i = 0; i < m; i++) {
                    double k = (i + 0.5) / (double) m;
                    double n_value = 0d;
                    double a = (1 - k) / (h * Util.M_SQRT2);
                    double b = -k / (h * Util.M_SQRT2);
                    if (Math.abs(a) > 6d || Math.abs(b) > 6d) {
                        n_value = 1d * Math.signum(a - b);
                    } else {
                        n_value = 0.5 * (Erf.erf(a) - Erf.erf(b));
                    }
                    N_TABLE.set(i, n_value);
                }

                int p = 0;

                while (p + 1 < n) {
                    double ll = 0d;
                    for (int i = 0; i < m; i++) {
                        double l = Double.MIN_VALUE;
                        for (int j = 0; j < m; j++) {
                            if (i == j) {
                                continue;
                            }
                            int dx = (int) Math.abs(id.get(p, i) - id.get(p, j));
                            int dy = (int) Math.abs(id.get(p + 1, i) - id.get(p + 1, j));

                            if (P_TABLE.get(dx) == -1d) {
                                double v = Math.exp(-Math.pow((td.get(p, i) - td.get(p, j)), 2d) / (2 * h * h));
                                P_TABLE.set(dx, v);
                            }
                            if (P_TABLE.get(dy) == -1d) {
                                double v = Math.exp(-Math.pow((td.get(p + 1, i) - td.get(p + 1, j)), 2d) / (2 * h * h));
                                P_TABLE.set(dy, v);
                            }
                            l += P_TABLE.get(dx) * P_TABLE.get(dy) / (N_TABLE.get((int) id.get(p, j)) * N_TABLE.get((int) id.get(p + 1, j)));
                        }
                        ll += Math.log(l) - (2 * Math.log(h)) - Math.log(2 * Util.M_PI * (m - 1));
                    }
                    p += 2;
                    logl += ll;
                }
                logl = -logl - Math.log(2 / (Util.M_PI * (1 + h * h)));
            } catch (MathException me) {
                me.printStackTrace();
            }
            return logl;
        }

        public void setG1(int g1) {
            this.g1 = g1;
        }

        public void setG2(int g2) {
            this.g2 = g2;
        }

        public void setIdata(DoubleMatrix2D i) {
            idata = i;
        }

        public void setTdata(DoubleMatrix2D t) {
            tdata = t;
        }
    }

    static class KernelWidthConvergenceChecker implements ConvergenceChecker {

        double tolerance = 1e-15;
        double previousMinValue = Double.MAX_VALUE;

        public boolean converged(PointCostPair[] pcps) {
            double minValue = Double.MAX_VALUE;
            for (PointCostPair pcp : pcps) {
                minValue = (minValue < pcp.getCost()) ? minValue : pcp.getCost();
            }
            if ((previousMinValue - minValue) <= tolerance) {
                return true;
            } else {
                previousMinValue = minValue;
            }
            return false;
        }
    }
    
    public static class Shuffle {
        
        static java.util.Random rand = new java.util.Random(System.currentTimeMillis());
        // construct a random permutation of the ints 0 .. n - 1
        // represented as an int array of length n
        //
        public static int[] permutation(int n) {
            
            assert n > 0;
            //intitial element order is irrelevant so long as each int 1..n occurs exactly once
            //inorder initialization assures that is the case
            
            int[] sample = new int[n];
            for (int k = 0; k < sample.length; k++) {
                sample[k] = k;
            }
            //loop invariant: the tail of the sample array is randomized.
            //Intitally the tail is empty; at each step move a random
            //element from front of array into the tail, then decrement boundary of tail
            int last = sample.length - 1;   //last is maximal index of elements not in the tail
            
            while (last > 0) {
                // Select random index in range 0..last, and swap its contents with those at last
                // The null swap is allowed; it should be possible that sample[k] does not change
                swap(rand.nextInt(last + 1), last, sample);
                last -= 1;
            }
            return sample;
        }
        
        // swap the elements at indices j and k
        // j and k need not be distinct, allowing for null swaps
        //
        private static void swap(int j, int k, int[] array) {
            int temp = array[k];
            array[k] = array[j];
            array[j] = temp;
        }
        
        public static int[] sampleWithoutReplacement(int N, int M) {
            
            // create permutation 0, 1, ..., N-1
            int[] perm = new int[N];
            for (int i = 0; i < N; i++) {
                perm[i] = i;            // create random sample in perm[0], perm[1], ..., perm[M-1]
            }
            for (int i = 0; i < M; i++) {
                
                // random integer between i and N - 1
                int r = i + (int) (rand.nextInt(N - i));
                
                // swap elements at indices i and r
                int t = perm[r];
                perm[r] = perm[i];
                perm[i] = t;
            }
            
            int[] subset = new int[M];
            for (int i = 0; i < M; i++) {
                subset[i] = perm[i];
            }
            return subset;
        }

        public static int[] sampleWithReplacement(int N, int M) {
            
            // create permutation 0, 1, ..., N-1
            int[] perm = new int[N];
            for (int i = 0; i < N; i++) {
                perm[i] = i;            // create random sample in perm[0], perm[1], ..., perm[M-1]
            }
            int[] subset = new int[M];
            for (int i = 0; i < M; i++) {
                
                // random integer between 0 and N - 1
                int r = (int) (rand.nextInt(N - 1));
                
                //sample from N
                subset[i] = perm[r];
            }            
            return subset;
        }
    }
    
    public static class Sort_Gene implements Comparator<Gene> {
        
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


}
