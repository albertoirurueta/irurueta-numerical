/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.optimization;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.*;

import java.util.Arrays;

/**
 * This class searches for a multi dimension function local minimum.
 * The local minimum is searched by starting the algorithm at a start point
 * and a given set of deltas to obtain surrounding points to the start point
 * where the algorithm will be started.
 * The Simplex algorithm will increase or decrease such deltas and move the 
 * start point around until the minimum is found.
 * The implementation of this class is based on Numerical Recipes 3rd ed. 
 * Section 10.5 page 502.
 */
@SuppressWarnings("WeakerAccess")
public class SimplexMultiOptimizer extends MultiOptimizer {
    /**
     * Maximum number of iterations.
     */
    public static final int NMAX = 5000;
    
    /**
     * Small value considered to be machine precision.
     */
    public static final double TINY = 1e-10;
    
    /**
     * Constant defining default tolerance or accuracy to be achieved on the
     * minimum being estimated by this class.
     */    
    public static final double DEFAULT_TOLERANCE = 3e-8;
    
    /**
     * Minimum allowed tolerance value.
     */    
    public static final double MIN_TOLERANCE = 0.0;
    
    /**
     * The fractional tolerance in the function value such that failure to
     * decrease by more than this amount on one iteration signals doneness.
     */    
    private double ftol;
    
    /**
     * Number of function evaluations.
     */
    private int nfunc;
    
    /**
     * Number of points in the simplex.
     */
    private int mpts;
    
    /**
     * Number of dimensions of current function being optimized.
     */
    private int ndim;
    
    /**
     * Function values at the vertices of the simples.
     */
    private double[] y;
    
    /**
     * Current simplex.
     */
    private Matrix p;
        
    /**
     * Empty constructor.
     */
    public SimplexMultiOptimizer() {
        super();
        nfunc = mpts = ndim = 0;
        y = null;
        p = null;
        ftol = DEFAULT_TOLERANCE;
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param startPoint Point where the algorithm is started.
     * @param delta Delta to find other points close to the start point so that
     * an initial simplex is defined.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public SimplexMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener, 
            double[] startPoint, double delta, double tolerance) {
        super(listener);
        nfunc = mpts = ndim = 0;
        y = null;
        p = null;
        internalSetTolerance(tolerance);
        internalSetSimplex(startPoint, delta);
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param startPoint Point where the algorithm is started.
     * @param deltas Set of deltas to find other points close to the start point
     * so that an initial simplex is defined.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative or if
     * deltas don't have the same length as provided start point.
     */
    public SimplexMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener,
            double[] startPoint, double[] deltas, double tolerance) {
        super(listener);
        nfunc = mpts = ndim = 0;
        y = null;
        p = null;
        internalSetTolerance(tolerance);
        internalSetSimplex(startPoint, deltas);        
    }
        
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param simplex Initial simplex to start the algorithm. The simplex is
     * a set of points around the minimum to be found.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public SimplexMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener, Matrix simplex, 
            double tolerance) {
        super(listener);
        nfunc = mpts = ndim = 0;
        y = null;
        p = null;
        internalSetTolerance(tolerance);
        internalSetSimplex(simplex);
    }
        
    /**
     * Returns current simplex.
     * A simplex is a set of points delimiting an n-dimensional region.
     * The simplex can be seen as the n-dimensional version of a bracket of 
     * values.
     * This class must be started at an initial simplex where a minimum will be
     * searched.
     * As the algorithm iterates, the simplex will be moved around and its size
     * will be reduced until a minimum is found.
     * @return Current simplex.
     * @throws NotAvailableException Raised if not provided or computed.
     */
    public Matrix getSimplex() throws NotAvailableException {
        if (!isSimplexAvailable()) {
            throw new NotAvailableException();
        }
        return p;
    }
    
    /**
     * Sets a simplex defined as a central point and a set of surrounding points
     * at distance delta. The simplex will be made of the computed surrounding 
     * points.
     * @param startPoint Central point.
     * @param delta Distance of surrounding points.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setSimplex(double[] startPoint, double delta) 
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetSimplex(startPoint, delta);
    }

    /**
     * Internal method to set a simplex as a central point and a set of 
     * surrounding points at distance delta. The simplex will be made of the 
     * computed surrounding points.
     * This method does not check whether this instance is locked.
     * @param startPoint Central point.
     * @param delta Distance of surrounding points.
     */
    private void internalSetSimplex(double[] startPoint, double delta) {
        double[] deltas = new double[startPoint.length];
        Arrays.fill(deltas, delta);
        internalSetSimplex(startPoint, deltas);
    }
    
    /**
     * Sets a simplex defined as a central point and a set of surrounding points
     * at their corresponding distance deltas[i], where i corresponds to one
     * position of provided array of distances. The simplex will be made of the
     * computed surrounding points.
     * @param startPoint Central point.
     * @param deltas Distances of surrounding points. Each surrounding point can
     * have a different distance than the others. The number of provided 
     * distances must be equal to the dimension or length of the start point 
     * array.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if startPoint and deltas don't
     * have the same length.
     */
    public void setSimplex(double[] startPoint, double[] deltas)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetSimplex(startPoint, deltas);
    }
    
    /**
     * Internal method to set a simplex defined as a central point and a set of
     * surrounding points at their corresponding distance deltas[i], where i
     * corresponds to one position of provided array of distances. The simplex
     * will be made of the computed surrounding points.
     * @param startPoint Central point.
     * @param deltas Distances of surrounding points. Each surrounding point can
     * have a different distance than the others. The number of provided
     * distances must be equal to the dimension or length of the start point
     * array.
     * @throws IllegalArgumentException Raised if startPoint and deltas don't
     * have the same length.
     */
    private void internalSetSimplex(double[] startPoint, double[] deltas) {
        int localndim = startPoint.length;
        
        if(deltas.length != localndim) {
            throw new IllegalArgumentException();
        }
        
        Matrix pp;
        try {
            pp = new Matrix(localndim + 1, localndim);
        } catch (WrongSizeException e) {
            throw new IllegalArgumentException(e);
        }
        for (int i = 0; i < localndim + 1; i++) {
            for (int j = 0; j < localndim; j++) {
                pp.setElementAt(i, j, startPoint[j]);
            }
            if (i != 0) {
                pp.setElementAt(i, i - 1, 
                    pp.getElementAt(i, i - 1) + deltas[i - 1]);
            }
        }
                
        internalSetSimplex(pp);
    }
    
    /**
     * Sets simplex as a matrix containing on each row a point of the simplex.
     * The number of columns defines the dimension of the points in the 
     * function, if not properly set, then minimize() method will fail when 
     * being called.
     * @param simplex Simplex of points.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setSimplex(Matrix simplex) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetSimplex(simplex);
    }

    /**
     * Internal method to Set simplex as a matrix containing on each row a point 
     * of the simplex.
     * The number of columns defines the dimension of the points in the 
     * function, if not properly set, then minimize() method will fail when 
     * being called.
     * This method does not check whether this instance is locked.
     * @param simplex Simplex of points.
     */    
    private void internalSetSimplex(Matrix simplex) {
        p = simplex;
    }
    
    /**
     * Returns boolean indicating whether a simplex has been provided and is
     * available for retrieval.
     * @return True if available, false otherwise.
     */
    public boolean isSimplexAvailable() {        
        return p != null;
    }
    
    /**
     * Returns function evaluations at simplex points or vertices.
     * @return Function evaluations at simplex vertices.
     * @throws NotAvailableException Raised if not available for retrieval.
     * This parameter will be available one minimization has been computed.
     */
    public double[] getEvaluationsAtSimplex() throws NotAvailableException {
        if (!areFunctionEvaluationsAvailable()) {
            throw new NotAvailableException();
        }
        if (!areFunctionEvaluationsAvailable()) {
            throw new NotAvailableException();
        }
        return y;
    }
    
    /**
     * Returns boolean indicating whether function evaluations at simplex 
     * vertices are available for retrieval.
     * Function evaluations at simplex vertices will be available one 
     * minimization has been computed.
     * @return True if available, false otherwise.
     */
    public boolean areFunctionEvaluationsAvailable() {
        return y != null;
    }
    
    /**
     * Returns tolerance or accuracy to be expected on estimated local minimum.
     * @return Tolerance or accuracy to be expected on estimated local minimum.
     */      
    public double getTolerance() {
        return ftol;
    }
    
    /**
     * Sets tolerance or accuracy to be expected on estimated local minimum.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided tolerance is 
     * negative.
     */       
    public void setTolerance(double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * Internal method to set tolerance or accuracy to be expected on estimated
     * local minimum.
     * This method does not check whether this instance is locked.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if provided tolerance is 
     * negative.
     */        
    private void internalSetTolerance(double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        ftol = tolerance;
    }
    
    /**
     * This function estimates a function minimum.
     * Implementations of this class will usually search a local minimum inside
     * a given simplex of points. The algorithm will iterate moving the simplex
     * around and reducing its size until a minimum is found.
     * @throws LockedException Raised if this instance is locked, because
     * estimation is being computed.
     * @throws NotReadyException Raised if this instance is not ready, because
     * a listener, a gradient listener and a start point haven't been provided.
     * @throws OptimizationException Raised if the algorithm failed because of
     * lack of convergence or because function couldn't be evaluated.
     */     
    @Override
    public void minimize() throws LockedException, NotReadyException,
            OptimizationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        
        //clean previous result and previous function evaluations
        try {
            int ihi;
            int ilo;
            int inhi;
            mpts = p.getRows();
            ndim = p.getColumns();
            
            double[] psum = new double[ndim];
            double[] pmin = new double[ndim];
            double[] x = new double[ndim];
            double[] v1 = new double[1];
            double[] v2 = new double[1];
            
            //instantiate new function evaluations
            y = new double[mpts];
            
            for (int i = 0; i < mpts; i++) {
                for (int j = 0; j < ndim; j++) {
                    x[j] = p.getElementAt(i, j);
                }
                y[i] = listener.evaluate(x);
            }
            
            nfunc = 0;
            getPsum(p, psum);
            
            for (;;) {
                ilo = 0;
                if (y[0] > y[1]) {
                    ihi = 0;
                    inhi = 1;
                } else {
                    ihi = 1;
                    inhi = 0;
                }
                for (int i = 0; i < mpts; i++) {
                    if (y[i] <= y[ilo]) {
                        ilo = i;
                    }
                    if (y[i] > y[ihi]) {
                        inhi = ihi;
                        ihi = i;
                    } else if (y[i] > y[inhi] && i != ihi) {
                        inhi = i;
                    }
                }
                double rtol = 2.0 * Math.abs(y[ihi] - y[ilo]) /
                        (Math.abs(y[ihi]) + Math.abs(y[ilo]) + TINY);
                
                if (rtol < ftol) {
                    v1[0] = y[0];
                    v2[0] = y[ilo];
                    swap(v1, v2);
                    y[0] = v1[0];
                    y[ilo] = v2[0];
                    for (int i = 0; i < ndim; i++) {
                        v1[0] = p.getElementAt(0, i);
                        v2[0] = p.getElementAt(ilo, i);
                        swap(v1, v2);
                        p.setElementAt(0, i, v1[0]);
                        p.setElementAt(ilo, i, v2[0]);
                        pmin[i] = p.getElementAt(0, i);
                    }
                    fmin = y[0];
                    
                    //save result
                    xmin = pmin;
                    resultAvailable = true;

                    return;
                }
                if (nfunc >= NMAX) {
                    throw new OptimizationException();
                }
                
                nfunc += 2;
                double ytry = amotry(p, y, psum, ihi, -1.0, listener);
                
                if(ytry <= y[ilo]) {
                    amotry(p, y, psum, ihi, 2.0, listener);
                } else if (ytry >= y[inhi]) {
                    double ysave = y[ihi];
                    ytry = amotry(p, y, psum, ihi, 0.5, listener);
                    if (ytry >= ysave) {
                        for (int i = 0; i < mpts; i++) {
                            if (i != ilo) {
                                for (int j = 0; j < ndim; j++) {
                                    psum[j] = 0.5 * (p.getElementAt(i, j) + 
                                            p.getElementAt(ilo, j));
                                    p.setElementAt(i, j, psum[j]);
                                }
                                y[i] = listener.evaluate(psum);
                            }
                        }
                        nfunc += ndim;
                        getPsum(p, psum);
                    }
                } else {
                    --nfunc;
                }
            }
        } catch (EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a local minimum.
     * An instance is ready once a listener and a simplex are provided.
     * @return True if this instance is ready, false otherwise.
     */    
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isSimplexAvailable();
    }
    
    /**
     * Computes the sum of the elements of each matrix column.
     * @param p Matrix where columns will be summed.
     * @param psum Array to contain computed sums. Provided array must have a
     * length equal to the number of matrix columns. This is an output 
     * parameter.
     */
    private void getPsum(Matrix p, double[] psum) {
        for (int j = 0; j < ndim; j++) {
            double sum = 0.0;
            for (int i = 0; i < mpts; i++) {
                sum += p.getElementAt(i, j);
            }
            psum[j] = sum;
        }
    }
    
    /**
     * Internal method to move simplex around.
     * @param p a matrix.
     * @param y an array.
     * @param psum another array.
     * @param ihi a value.
     * @param fac another value.
     * @param listener a listener to evaluate function to optimize.
     * @return a value returned by the evaluated function.
     * @throws EvaluationException Raised if function evaluation fails.
     */
    private double amotry(Matrix p, double[] y, double[] psum, int ihi, 
            double fac, MultiDimensionFunctionEvaluatorListener listener)
            throws EvaluationException {
        double[] ptry = new double[ndim];
        double fac1 = (1.0 - fac) / ndim;
        double fac2 = fac1 - fac;
        for (int j = 0; j < ndim; j++) {
            ptry[j] = psum[j] * fac1 - p.getElementAt(ihi, j) * fac2;
        }
        double ytry = listener.evaluate(ptry);
        if (ytry < y[ihi]) {
            y[ihi] = ytry;
            for (int j = 0; j < ndim; j++) {
                psum[j] += ptry[j] - p.getElementAt(ihi, j);
                p.setElementAt(ihi, j, ptry[j]);
            }
        }
        
        return ytry;
    }
    
    /**
     * Swaps a and b values.
     * @param a value.
     * @param b value.
     */
    private void swap(double[] a, double[] b) {
        double tmp = a[0];
        a[0] = b[0];
        b[0] = tmp;
    }
}
