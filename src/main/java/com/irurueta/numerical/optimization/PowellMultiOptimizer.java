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

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.*;

/**
 * This class searches for a multi dimension function local minimum.
 * The local minimum is searched by starting the algorithm at a start point
 * and a given direction which should be close and point to the local minimum to 
 * be found to achieve the best accuracy with the lowest number of iterations.
 * The implementation of this class is based on Numerical Recipes 3rd ed. 
 * Section 10.7 page 509.
 */
@SuppressWarnings("WeakerAccess")
public class PowellMultiOptimizer extends LineMultiOptimizer {
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
     * Maximum allowed iterations.
     */
    public static final double ITMAX = 200;
    
    /**
     * A small number.
     */
    public static final double TINY = 1e-25;
    
    /**
     * The fractional tolerance in the function value such that failure to
     * decrease by more than this amount on one iteration signals doneness.
     */
    private double tolerance;

    /**
     * Member contains number of iterations that were needed to estimate a 
     * minimum.
     */    
    private int iter;
    
    /**
     * Value of the function at the minimum.
     */
    private double fret;
    
    /*
     * set of directions.
     */
    private Matrix ximat;
    
    /**
     * Empty constructor.
     */
    public PowellMultiOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
        iter = 0;
        fret = 0.0;
        ximat = null;
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public PowellMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener, double tolerance) {
        super(listener);
        internalSetTolerance(tolerance);
        iter = 0;
        fret = 0.0;
        ximat = null;
    }
     
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param point Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param directions Set of directions to start looking for a minimum. 
     * Provided matrix must have the same rows as the number of dimensions of 
     * the function being evaluated. Each column will be a direction to search
     * for a minimum.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length or if provided tolerance is negative.
     */
    public PowellMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener, double[] point, 
            Matrix directions, double tolerance) {
        super(listener);
        internalSetPointAndDirections(point, directions);
        internalSetTolerance(tolerance);
        iter = 0;
        fret = 0.0;
    }
        
    /**
     * Returns set of directions to start looking for a minimum.
     * Returned matrix will have the same number of rows as the number of
     * dimensions of the function (or the start point).
     * Each column of the matrix will represent a vector containing a direction
     * to search for a minimum.
     * @return Set of directions.
     * @throws NotAvailableException Raised if this has not yet been provided or 
     * computed.
     */
    public Matrix getDirections() throws NotAvailableException {
        if (!areDirectionsAvailable()) {
            throw new NotAvailableException();
        }
        return ximat;
    }
    
    /**
     * Returns boolean indicating whether set of directions is available for
     * retrieval.
     * @return True if available, false otherwise.
     */
    public boolean areDirectionsAvailable() {
        return ximat != null;
    }

    /**
     * Sets start point and set of directions to start looking for minimum.
     * @param point Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param directions Set of directions to start looking for a minimum. 
     * Provided matrix must have the same rows as the number of dimensions of 
     * the function being evaluated. Each column will be a direction to search
     * for a minimum.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length.
     */    
    public void setPointAndDirections(double[] point, Matrix directions)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetPointAndDirections(point, directions);
    }
    
    /**
     * Sets start point where local minimum is searched nearby.
     * @param startPoint Start point to search for a local minimum
     * @throws LockedException Raised if this instance is locked.
     */    
    public void setStartPoint(double[] startPoint) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        p = startPoint;
    }
    
    /**
     * Returns tolerance or accuracy to be expected on estimated local minimum.
     * @return Tolerance or accuracy to be expected on estimated local minimum.
     */    
    public double getTolerance() {
        return tolerance;
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
        this.tolerance = tolerance;
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
     * This function estimates a function minimum.
     * Implementations of this class will usually search a local minimum close
     * to a start point and will start looking into provided start direction.
     * Minimization of a function f. Input consists of an initial starting
     * point p. The initial matrix ximat, whose columns contain the initial
     * set of directions, is set to the identity. Returned is the best point 
     * found, at which point fret is the minimum function value and iter is 
     * the number of iterations taken.
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
        
        buildDirections();
        
        //in case that set of directions have been directly provided, we check
        //for correctness
        int n = p.length;
        if (ximat.getRows() != ximat.getColumns() || ximat.getRows() != n) {
            locked = false;
            throw new OptimizationException();
        }
        
        double fptt;
        double[] pt = new double[n];
        double[] ptt = new double[n];
        
        //set vector of directions
        if (!isDirectionAvailable() || xi.length != n) {
            xi = new double[n];
        }
        
        try {
            fret = listener.evaluate(p);
            
            //save the initial point
            System.arraycopy(p, 0, pt, 0, n);
            
            for (iter = 0; ; ++iter) {
                double fp = fret;
                int ibig = 0;
                //Will be the biffest function decrease
                double del = 0.0;
                //In each iteration, loop over all directions in the set.
                for (int i = 0; i < n; i++) {
                    //Copy the direction contained in i-th column of ximat
                    ximat.getSubmatrixAsArray(0, i, n-1, i, xi);
                    fptt = fret;
                    
                    //minimize along it, and record it if it is the largest
                    //decrease so far.
                    fret = linmin();
                    if (fptt - fret > del) {
                        del = fptt - fret;
                        ibig = i + 1;
                    }
                }
                
                //Here comes the termination criterion.
                if (2.0 * (fp - fret) <= tolerance * (Math.abs(fp) + 
                        Math.abs(fret)) + TINY) {
                    break; //minimum has been found
                }
                if (iter == ITMAX) {
                    //to many iterations
                    locked = false;
                    throw new OptimizationException();
                }
                
                //Construct the extrapolated point and the average direction
                //moved. Save the old starting point
                for (int j = 0; j < n; j++) {
                    ptt[j] = 2.0 * p[j] - pt[j];
                    xi[j] = p[j] - pt[j];
                    pt[j] = p[j];
                }
                
                //Function value at extrapolated point
                fptt = listener.evaluate(ptt);
                
                if (fptt < fp) {
                    double t = 2.0 * (fp - 2.0 * fret + fptt) *
                            sqr(fp - fret - del) - del * sqr(fp - fptt);
                    
                    if (t < 0.0) {
                        //Move to the minimum of the new direction, and save the
                        //new direction
                        fret = linmin();
                        
                        ximat.setSubmatrix(0, ibig - 1, n - 1, ibig - 1, ximat, 
                                0, n - 1, n - 1, n - 1);
                        ximat.setSubmatrix(0, n - 1, n - 1, n - 1, xi);
                    }
                }
            }
        } catch (AlgebraException | EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }
        
        //set result
        xmin = p;
        resultAvailable = true;
        fmin = fret;
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a local minimum.
     * An instance is ready once a listener and a start point are provided.
     * @return True if this instance is ready, false otherwise.
     */    
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isStartPointAvailable();
    }

    /**
     * Internal method to set start point and set of directions to start looking
     * for minimum.
     * This method does not check whether this instance is locked.
     * @param point Start point where algorithm will be started. Start point
     * should be close to the local minimum to be found. Provided array must
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param directions Set of directions to start looking for a minimum.
     * Provided matrix must have the same rows as the number of dimensions of
     * the function being evaluated. Each column will be a direction to search
     * for a minimum.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length.
     */
    private void internalSetPointAndDirections(double[] point,
                                               Matrix directions) {
        if ((point.length != directions.getRows()) ||
                (point.length != directions.getColumns())) {
            throw new IllegalArgumentException();
        }
        p = point;
        ximat = directions;
        xi = null;
    }

    /**
     * Internal method to build or rebuild the set of directions if needed.
     * @throws NotReadyException Raised if no start point has yet been provided.
     */
    private void buildDirections() throws NotReadyException {
        if (!isStartPointAvailable()) {
            throw new NotReadyException();
        }
        
        int n = p.length;
        
        if (areDirectionsAvailable() && (ximat.getRows() == n) &&
                (ximat.getColumns() == n)) {
            return;
        }
        
        try {
            ximat = Matrix.identity(n, n);
        } catch (WrongSizeException ignore) {
            //never happens
        }
    }
    
    /**
     * Computes the squared value of provided double.
     * @param x Value to be squared.
     * @return Squared value.
     */
    private double sqr(double x) {
        return x * x;
    }
}
