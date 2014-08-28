package fr.systemsbiology.aracneAlgorithm.internal.mutualInfoMetric;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import org.apache.commons.math3.special.Erf;

import java.util.Vector;
import java.util.HashMap;

/**
 * @author manjunath at c2b2 dot columbia dot edu
 */

public class Mutual_Info {
	
	public static enum ALGORITHM {NAIVE_BAYES, FIXED_BANDWIDTH, VARIABLE_BANDWIDTH, ADAPTIVE_PARTITIONING};

    static final int MIBLOCKS = 2;
    static final int BINS = 1000;
    int Microarray_Num;
    int miSteps;
    int[][] MI_Space;
    double[][] MI_Prob;
    double MA_Per_MI_Step;
    double[] histogram;
    double[] background;
    int Max_Histo_Bin;
    int Max_Background_Bin;    
    //type of calculation
    ALGORITHM type;    

    double variance2;
    double[][] prob_table;
    double[] norm_1D_table;
    double[][] norm_2D_table;
    double[] kernel_bandwidth;
    boolean Copula_Transform;

    cern.jet.math.Functions F = cern.jet.math.Functions.functions;

    public Mutual_Info(int maNo, int _miSteps, double sigma, ALGORITHM type) {
        this.Microarray_Num = maNo;
        this.variance2 = 2 * sigma * sigma;
        this.type = type;
        switch (type) {
            case NAIVE_BAYES:
                // initialize the lookup table for computing Gaussian kernels
                initProbTable(maNo, 0);

                norm_1D_table = new double[maNo];
                norm_2D_table = new double[maNo][maNo];

                miSteps = _miSteps;
                MA_Per_MI_Step = (double) Microarray_Num / (double) miSteps;
                MI_Space = new int[miSteps][];
                MI_Prob = new double[miSteps][];
                for (int i = 0; i < miSteps; i++) {
                    MI_Space[i] = new int[miSteps];
                    MI_Prob[i] = new double[miSteps];
                    for (int j = 0; j < miSteps; j++) {
                        MI_Space[i][j] = 0;
                        MI_Prob[i][j] = 0.0;
                    }
                }
                break;

            case ADAPTIVE_PARTITIONING:

                Set_Copula_Transform(false);

                break;

            case FIXED_BANDWIDTH:

                // set coupla transform
                Set_Copula_Transform(true);

                // initialize the lookup table for computing Gaussian kernels
                initProbTable(maNo, -1);

                // initialize the lookup table for computing normalization factors
                norm_1D_table = new double[maNo];

                norm_2D_table = new double[maNo][];
                for (int p = 0; p < maNo; p++) {
                    norm_2D_table[p] = new double[maNo];
                }

                Initialize_Norm_Table(sigma);

                break;

            case VARIABLE_BANDWIDTH:
                // do not use copula transform
                Set_Copula_Transform(false);
                break;

            default:
                System.err.println("Constructor: MI computation not supported");
                break;
        }
    }

    private void initProbTable(int maNo, int value) {
        prob_table = new double[maNo][];
        for (int i = 0; i < maNo; i++) {
            prob_table[i] = new double[maNo];
            for (int j = 0; j < maNo; j++) {
                prob_table[i][j] = value;
            }
        }
    }

    ALGORITHM Get_Type() {
        return type;
    }

    void Set_Copula_Transform(boolean t) {
        Copula_Transform = t;
    }

    boolean Is_Copula_Transform() {
        return Copula_Transform;
    }

    int Get_Microarray_Num() {
        return Microarray_Num;
    }

    void Initialize_Norm_Table(double sigma) {
        // initialize the 1D table
        for (int i = 0; i < Microarray_Num; i++) {
            // array id is from 0 to M-1, after copula transform, it becomes from 0 to 1
            // with 1/(2M) at both sides
            double x = Copula_Transform ? (((double) (i + 1) - 0.5) / Microarray_Num) : (double) (i);
           // try {
                norm_1D_table[i] = 0.5 * (Erf.erf((1 - x) / (sigma * Util.M_SQRT2)) - Erf.erf((0 - x) / (sigma * Util.M_SQRT2)));
           // } catch (MathException e) {
           //     e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
           // }
        }

        // initialize the 2D table
        for (int j = 0; j < Microarray_Num; j++) {
            for (int k = j; k < Microarray_Num; k++) {
                norm_2D_table[j][k] = norm_1D_table[j] * norm_1D_table[k];
                norm_2D_table[k][j] = norm_2D_table[j][k];
            }
        }
    }

    double Compute_Pairwise_MI(Vector<Gene> g1, Vector<Gene> g2, double sigmaX, double sigmaY) {
        double mi = 0.0;
        switch (type) {
            case NAIVE_BAYES:
                mi = Get_Mutual_Info_XY_NAIVE_BAYES(g1, g2);
                break;
            case FIXED_BANDWIDTH:
                mi = Get_Mutual_Info_XY_FIXED_BANDWIDTH(g1, g2);
                break;

            case ADAPTIVE_PARTITIONING:
                mi = Get_Mutual_Info_XY_ADAPTIVE_PARTITIONING(g1, g2);
                break;

            case VARIABLE_BANDWIDTH:
                mi = Get_Mutual_Info_XY_VARIABLE_BANDWIDTH(g1, g2, sigmaX, sigmaY);
                break;

            default:
                System.err.println("Pairwise: MI computation not supported");
                break;
        }
        return mi;
    }

    double Get_Mutual_Info_XY_FIXED_BANDWIDTH(Vector<Gene> g1, Vector<Gene> g2) {
        final int size = g1.size();

        Gene[] gene1 = g1.toArray(new Gene[]{});
        Gene[] gene2 = g2.toArray(new Gene[]{});
        
        for (int i = 0; i < size; i++) {
            //gene1[i].xi = i;
            if (Copula_Transform) {
                gene1[i].x = ((double) (i + 1) - 0.5) / size;
            }
        }

        for (int i = 0; i < size; i++) {
            //gene2[i].xi = i;
            if (Copula_Transform) {
                gene2[i].x = ((double) (i + 1) - 0.5) / size;
            }
        }

        double sum = 0.0;

        for (int i = 0; i < size; i++) {
            double v = Get_Kernel_FIXED_BANDWIDTH(gene1, gene2, i);
            sum += Math.log(v);
        }
        double mi = sum / (double) size;
        return (Math.max(mi, 0.0));

    }

    double Get_Mutual_Info_XY_VARIABLE_BANDWIDTH(Vector<Gene> g1, Vector<Gene> g2, double sigmaX, double sigmaY) {
        int size = g1.size();
        double sum = 0.0;
        for (int i = 0; i < size; i++) {
            double v = Get_Kernel_VARIABLE_BANDWIDTH(g1, g2, i, sigmaX, sigmaY);
            sum += Math.log(v);
        }
        double mi = sum / (double) size;
        return (Math.max(mi, 0.0));
    }

    double Get_Kernel_VARIABLE_BANDWIDTH(Vector<Gene> g1, Vector<Gene> g2, int i, double sigmaX, double sigmaY) {
        double fxy = 0.0;
        double fx = 0.0;
        double fy = 0.0;

        int size = g1.size();

        HashMap<Integer, Gene> yranks = new HashMap<Integer, Gene>();

        for (int j = 0; j < g2.size(); j++) {
            yranks.put(g2.get(j).maId, g2.get(j));
        }

        for (int j = 0; j < size; j++) {
            double mu_x = g1.get(j).x;
            double mu_y = yranks.get(g1.get(j).maId).x;

            double dx = Math.abs(g1.get(i).x - mu_x);
            double dy = Math.abs(yranks.get(g1.get(i).maId).x - mu_y);

            fx += Util.normalPDF(dx, sigmaX);
            fy += Util.normalPDF(dy, sigmaY);
            fxy += Util.multinormalPDF(dx, dy, sigmaX, sigmaY);
        }

        return ((fxy * size) / (fx * fy));
    }

    double Get_Mutual_Info_XY_NAIVE_BAYES(Vector<Gene> g1, Vector<Gene> g2) {
        double H_xy = 0;
        double H_x = -Math.log(1 / (double) miSteps);
        double H_y = H_x;

        final int size = g1.size();

        HashMap<Integer, Gene> yranks = new HashMap<Integer, Gene>();

        for (int j = 0; j < g2.size(); j++) {
            yranks.put(g2.get(j).maId, g2.get(j));
        }
        
        for (int i = 0; i < miSteps; i++) {
            for (int j = 0; j < miSteps; j++) {
                MI_Space[i][j] = 0;
            }
        }
        for (int i = 0; i < g1.size(); i++) {
            int x = (int) ((double) g1.get(i).xi / MA_Per_MI_Step);
            Gene g2i = yranks.get(g1.get(i).maId);
            int y = (int) ((double) g2i.xi / MA_Per_MI_Step);
            MI_Space[x][y]++;
        }

        //
        // Compute the value for the initial block
        //
        int count = 0;
        for (int di = 0; di < MIBLOCKS; di++) {
            for (int dj = 0; dj < MIBLOCKS; dj++) {
                count += MI_Space[di][dj];
            }
        }
        double sumP = 0;
        for (int i = 0; i < miSteps; i++) {
            int count1 = count;
            for (int j = 0; j < miSteps; j++) {
                double p = (double) count1 / (double) (g1.size() * MIBLOCKS * MIBLOCKS);
                sumP += p;
                if (p > 0) {
                    H_xy -= p * Math.log(p);
                }
                MI_Prob[i][j] = (double) count1;
                for (int off = 0; off < MIBLOCKS; off++) {
                    //
                    // Wrap around the torus
                    //
                    int x = (i + off) % miSteps;
                    int y = (j + MIBLOCKS) % miSteps;
                    count1 -= MI_Space[x][j];
                    count1 += MI_Space[x][y];
                }
            }
            for (int off = 0; off < MIBLOCKS; off++) {
                int x = (i + MIBLOCKS) % miSteps;
                count -= MI_Space[i][off];
                count += MI_Space[x][off];
            }
        }
        double mi = (H_x + H_y - H_xy) / (H_x + H_y);
        return mi;
    }

    double Get_Kernel_FIXED_BANDWIDTH(Gene[] g1, Gene[] g2, int i) {
        // mu_x, mu_y - the center of each gaussian kernel
        // ix, iy     - distance to the gaussian center in unit of index.
        // dx, dy     - the actual distance used to compute the kernel denstiy.
        //              if compula transformed, they are equal to ix and iy rescaled
        //              between 0 and 1;

        double fxy = 0.0;
        double fx = 0.0;
        double fy = 0.0;

        HashMap<Integer, Gene> yranks = new HashMap<Integer, Gene>();

        for (int j = 0; j < g2.length; j++) {
            yranks.put(g2[j].maId, g2[j]);
        }
        
        final int size = g1.length;
        Gene g1i = g1[i];
        Gene g2i = yranks.get(g1i.maId);

        for (int j = 0; j < size; j++) {
            Gene g1j = g1[j];
            Gene g2j = yranks.get(g1j.maId);


            int ix = g1i.xi - g1j.xi;
            if (ix < 0) {
                ix = ix * -1;
            }
            int iy = g2i.xi - g2j.xi;
            if (iy < 0) {
                iy = iy * -1;
            }
            double dx = g1i.x - g1j.x;
            double dy = g2i.x - g2j.x;

            // compute the kernel density as necessary
            if (prob_table[ix][iy] == -1.0) {
                // if (ix, iy) is not computed, at least one of ix and iy is not computed
                // otherwise, all three should have been computed
                // when either ix=0 or iy=0, the trctdGssnBi.getProbability(dx, dy) is
                // equal to trctdGssnUni.getProbability(dx) or trctdGssnUni.getProbability(dy)
                if (prob_table[ix][0] == -1.0) {
                    prob_table[ix][0] = Math.exp(-(dx * dx) / variance2);
                    prob_table[0][ix] = prob_table[ix][0];
                }

                if (prob_table[0][iy] == -1.0) {
                    prob_table[0][iy] = Math.exp(-(dy * dy) / variance2);
                    prob_table[iy][0] = prob_table[0][iy];
                }

                prob_table[ix][iy] = prob_table[ix][0] * prob_table[0][iy];
                prob_table[iy][ix] = prob_table[ix][iy];
            }

            fx += prob_table[ix][0] / norm_1D_table[g1j.xi];
            fy += prob_table[0][iy] / norm_1D_table[g2j.xi];
            fxy += prob_table[ix][iy] / norm_2D_table[g1j.xi][g2j.xi];
        }

        return ((fxy * size) / (fx * fy));
    }

    double Get_Mutual_Info_XY_ADAPTIVE_PARTITIONING(Vector<Gene> g1, Vector<Gene> g2) {

        DoubleMatrix1D xranks = new DenseDoubleMatrix1D(g1.size());
        for (int i = 0; i < g1.size(); i++) {
            xranks.set(g1.get(i).maId, i + 1);
        }

        DoubleMatrix1D yranks = new DenseDoubleMatrix1D(g2.size());
        for (int i = 0; i < g2.size(); i++) {
            yranks.set(g2.get(i).maId, i + 1);
        }

        double xcor = 0;
        int npar = 1;
        int run = 0;
        int N = g1.size();
        DoubleMatrix1D poc = new DenseDoubleMatrix1D(20);
        poc.assign(1);
        DoubleMatrix1D kon = new DenseDoubleMatrix1D(20);
        kon.assign(N);
        DoubleMatrix1D poradi = new DenseDoubleMatrix1D(N);
        DoubleMatrix1D NN = new DenseDoubleMatrix1D(4);
        DoubleMatrix2D marg = new DenseDoubleMatrix2D(4, 20);

        for (int i = 1; i <= N; i++) {
            poradi.set(i - 1, i);
        }

        int t[] = {1, 1, N, N};
        marg.set(0, 0, t[0]);
        marg.set(1, 0, t[1]);
        marg.set(2, 0, t[2]);
        marg.set(3, 0, t[3]);
        
        while (npar > 0) {
            run++;
            int apoc = (int)poc.get(npar - 1);
            int akon = (int)kon.get(npar - 1);
            
            DoubleMatrix1D apor = new DenseDoubleMatrix1D(akon - apoc + 1);

            apor.viewPart(0, akon - apoc + 1).assign(poradi.viewPart(apoc - 1, akon - apoc + 1));

            int Nex = akon - apoc + 1;
            
            DoubleMatrix2D rownpar = marg.viewPart(0, npar - 1, 4, 1);

            final int ave1 = (int) Math.floor((rownpar.get(0, 0) + rownpar.get(2, 0)) / 2);
            final int ave2 = (int) Math.floor((rownpar.get(1, 0) + rownpar.get(3, 0)) / 2);

            DoubleMatrix1D J1 = new DenseDoubleMatrix1D(apor.size());
            int[] is = new int[apor.size()];
            for (int j = 0; j < apor.size(); j++) is[j] = (int) apor.get(j) - 1;

            J1.assign(xranks.viewSelection(is), new DoubleDoubleFunction() {
                public double apply(double x, double y) {
                    return (y <= ave1) ? 1d: 0d;
                }                
            });

            DoubleMatrix1D J2 = new DenseDoubleMatrix1D(apor.size());
            J2.assign(yranks.viewSelection(is), new DoubleDoubleFunction() {
                public double apply(double x, double y) {
                    return (y <= ave2) ? 1d: 0d;
                }                
            });
                        
            
            DoubleMatrix2D I = new DenseDoubleMatrix2D(apor.size(), 4);

            I.viewColumn(0).assign(J1.copy().assign(J2.copy(), new DoubleDoubleFunction(){
                public double apply(double x, double y) {
                    return (x == 1d) && (y == 1d) ? 1 : 0;
                }
            }));

            I.viewColumn(1).assign(J1.copy().assign(J2.copy(), new DoubleDoubleFunction(){
                public double apply(double x, double y) {
                    return (x == 1d) && !(y == 1d) ? 1 : 0;
                }
            }));

            I.viewColumn(2).assign(J1.copy().assign(J2, new DoubleDoubleFunction(){
                public double apply(double x, double y) {
                    return !(x == 1d) && (y == 1d) ? 1 : 0;
                }
            }));

            I.viewColumn(3).assign(J1.copy().assign(J2, new DoubleDoubleFunction(){
                public double apply(double x, double y) {
                    return !(x == 1d) && !(y == 1d) ? 1 : 0;
                }
            }));

            NN.set(0, I.viewColumn(0).zSum());
            NN.set(1, I.viewColumn(1).zSum());
            NN.set(2, I.viewColumn(2).zSum());
            NN.set(3, I.viewColumn(3).zSum());
                                    
            DoubleMatrix2D amarg = new DenseDoubleMatrix2D(4, 4);

            amarg.set(0, 0, marg.get(0, npar - 1));
            amarg.set(1, 0, marg.get(1, npar - 1));
            amarg.set(2, 0, ave1);
            amarg.set(3, 0, ave2);

            amarg.set(0, 1, marg.get(0, npar - 1));
            amarg.set(1, 1, ave2 + 1);
            amarg.set(2, 1, ave1);
            amarg.set(3, 1, marg.get(3, npar - 1));

            amarg.set(0, 2, ave1 + 1);
            amarg.set(1, 2, marg.get(1, npar - 1));
            amarg.set(2, 2, marg.get(2, npar - 1));
            amarg.set(3, 2, ave2);

            amarg.set(0, 3, ave1 + 1);
            amarg.set(1, 3, ave2 + 1);
            amarg.set(2, 3, marg.get(2, npar - 1));
            amarg.set(3, 3, marg.get(3, npar - 1));
                        
            DoubleMatrix1D constant = new DenseDoubleMatrix1D(4);
            constant.assign(((double)Nex) / 4);
            double tst = NN.copy().assign(constant, F.minus).assign(F.square).zSum() * 4 / (double)Nex;
            if (tst > 7.8 || run == 1) {
                --npar;
                for (int i = 0; i < 4; i++) {
                    if (NN.get(i) > 2) {
                        ++npar;
                        akon = apoc + (int)NN.get(i) - 1;
                        poc.set(npar - 1, apoc);
                        kon.set(npar - 1, akon);
                        marg.viewColumn(npar - 1).assign(amarg.viewColumn(i));
                        
                        IntArrayList indices = new IntArrayList(); 
                        DoubleArrayList values = new DoubleArrayList();
                        I.viewColumn(i).getNonZeros(indices, values);
                        is = new int[indices.size()];
                        for (int j = 0; j < indices.size(); j++) is[j] = indices.get(j);
                        poradi.viewPart(apoc - 1, akon - apoc + 1).assign(apor.viewSelection(is));                        
                        apoc = akon + 1;
                    } else if (NN.get(i) > 0) {
                        DoubleMatrix1D t3 = amarg.viewColumn(i);
                        double Nx = t3.get(2) - t3.get(0) + 1;
                        double Ny = t3.get(3) - t3.get(1) + 1;
                        xcor += NN.get(i) * Math.log(NN.get(i) / (Nx * Ny));
                    }
                }
            } else {
                DoubleMatrix1D t3 = marg.viewColumn((npar - 1));
                double Nx = t3.get(2) - t3.get(0) + 1;
                double Ny = t3.get(3) - t3.get(1) + 1;
                xcor += Nex * Math.log(Nex / (Nx * Ny));
                --npar;
            }
        }
        xcor = xcor / N + Math.log(N);
        return xcor;
    }   
}