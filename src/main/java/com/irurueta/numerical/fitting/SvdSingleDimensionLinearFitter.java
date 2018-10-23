/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.SingularValueDecomposer;
import com.irurueta.numerical.NotReadyException;

/**
 * Fits provided data (x,y) to a function made of a linear combination of
 * functions used as a basis (i.e. f(x) = a * f0(x) + b * f1(x) + ...).
 * Where f0, f1, ... is the function basis which ideally should be formed by
 * orthogonal function.
 * This is an improvement on SimpleLinearFitter for cases where basis functions
 * cannot be clearly distinguished from provided data, to avoid matrix 
 * singularities and obtain better results.
 * This class is based on the implementation available at Numerical Recipes
 * 3rd Ed, page 795.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class SvdSingleDimensionLinearFitter extends SingleDimensionLinearFitter {

    /**
     * Default tolerance.
     */
    public static final double TOL = 1e-12;
    
    /**
     * Tolerance to define convergence threshold for SVD.
     */
    private double tol;
    
    /**
     * Constructor.
     */
    public SvdSingleDimensionLinearFitter() {
        super();
        tol = TOL;
    }
    
    /**
     * Constructor.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     * length.
     */
    public SvdSingleDimensionLinearFitter(double[] x, double[] y, double[] sig)
            throws IllegalArgumentException {
        super(x, y, sig);
        tol = TOL;
    }
    
    /**
     * Constructor.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length.
     */    
    public SvdSingleDimensionLinearFitter(double[] x, double[] y, double sig)
            throws IllegalArgumentException {
        super(x, y, sig);
        tol = TOL;
    }
    
    /**
     * Constructor.
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point.
     * @throws FittingException if evaluation fails.
     */
    public SvdSingleDimensionLinearFitter(LinearFitterSingleDimensionFunctionEvaluator evaluator)
            throws FittingException {
        super(evaluator);
        tol = TOL;
    }
    
    /**
     * Constructor.
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant.
     * @throws FittingException if evaluation fails.
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length.
     */    
    public SvdSingleDimensionLinearFitter(LinearFitterSingleDimensionFunctionEvaluator evaluator,
            double[] x, double[] y, double[] sig)
            throws FittingException, IllegalArgumentException {
        super(evaluator, x, y, sig); 
        tol = TOL;
    }
    
    /**
     * Constructor.
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant.
     * @throws FittingException if evaluation fails.
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length.
     */    
    public SvdSingleDimensionLinearFitter(LinearFitterSingleDimensionFunctionEvaluator evaluator,
            double[] x, double[] y, double sig)
            throws FittingException, IllegalArgumentException {
        super(evaluator, x, y, sig);
        tol = TOL;
    }

    /**
     * Returns tolerance to define convergence threshold for SVD.
     * @return tolerance to define convergence threshold for SVD.
     */
    public double getTol() {
        return tol;
    }
    
    /**
     * Sets tolerance to define convergence threshold for SVD.
     * @param tol tolerance to define convergence threshold for SVD.
     */
    public void setTol(double tol) {
        this.tol = tol;
    }    
    
    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value.
     * @throws FittingException if fitting fails.
     * @throws NotReadyException if enough input data has not yet been provided.
     */    
    @Override
    public void fit() throws FittingException, NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        try {
            resultAvailable = false;
            
            int i,j,k;
            double tmp,thresh,sum;
            Matrix aa = new Matrix(ndat, ma);
            double[] b = new double[ndat];
            for (i = 0; i < ndat; i++) {
                evaluator.evaluate(x[i], afunc);
		        tmp = 1.0 / sig[i];
		        for (j = 0; j < ma; j++) {
                    aa.setElementAt(i, j, afunc[j]*tmp);
                }
		        b[i] = y[i] * tmp;
            }
            
            SingularValueDecomposer svd = 
                    new SingularValueDecomposer(aa);
            svd.decompose();
            thresh = (tol > 0.0 ? tol * svd.getSingularValues()[0] : -1.0);
            svd.solve(b, thresh, a);
            chisq = 0.0;
            for (i = 0; i < ndat; i++) {
		        sum = 0.0;
		        for (j = 0; j < ma; j++) {
                    sum += aa.getElementAt(i, j)*a[j];
                }
		        chisq += Math.pow(sum-b[i], 2.0);
            }
            for (i = 0; i < ma; i++) {
                for (j = 0; j < i + 1; j++) {
                    sum = 0.0;
                    double[] w = svd.getSingularValues();
                    double tsh = svd.getNegligibleSingularValueThreshold();
                    Matrix v = svd.getV();
                    for (k=0;k<ma;k++) {
                        if (w[k] > tsh) {
                            sum += v.getElementAt(i, k) * v.getElementAt(j, k) / 
                                    Math.pow(w[k], 2.0);
                        }
                    }
                    covar.setElementAt(j, i, sum);
                    covar.setElementAt(i, j, sum);
		        }
            }
            
            resultAvailable = true;
        } catch (FittingException e) {
            throw e;
        } catch (Throwable t) {
            throw new FittingException(t);
        }        
    }    
}
