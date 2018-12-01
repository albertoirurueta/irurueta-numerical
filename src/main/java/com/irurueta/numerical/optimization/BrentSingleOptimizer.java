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

import com.irurueta.numerical.*;

/**
 * This class uses Brent algorithm to determine a local function minimum for
 * single dimension functions.
 * Brent's algorithm will search for a local minimum inside the provided or
 * computed bracket of values.
 * Brent's algorithm is not the fastest among all optimization algorithms, but 
 * it is usually one that provides good convergence for most continuous 
 * functions.
 * It 's recommended to always set or compute a bracket of values, as the search
 * range is reduced and results usually become more accurate.
 * The implementation of this class is based on Numerical Recipes 3rd ed. 
 * Section 10.3. Page 496.
 */
@SuppressWarnings("WeakerAccess")
public class BrentSingleOptimizer extends BracketedSingleOptimizer {
    
    /**
     * Is the maximum allowed number of iterations.
     */
    public static final int ITMAX = 100;
    
    /**
     * Is the golden ratio.
     */
    public static final double CGOLD = 0.3819660;
    
    /**
     * Small number that protects against trying to achieve fractional accuracy
     * for a minimum that happens to be exactly zero.
     */
    public static final double ZEPS = 1e-10;
    
    /**
     * Constant defining the default accuracy of the estimated minimum.
     */
    public static final double DEFAULT_TOLERANCE = 3e-8;
    
    /**
     * Minimum allowed tolerance.
     */
    public static final double MIN_TOLERANCE = 0.0;
    
    /**
     * Tolerance value. The algorithm will iterate until the result converges
     * below this value of accuracy or until the maximum number of iterations is
     * achieved (and in such case, convergence will be assumed to have failed).
     */
    private double tolerance;

    /**
     * Empty constructor.
     */
    public BrentSingleOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
    }
        
    /**
     * Constructor. Creates an instance with provided bracket of values.
     * @param minEvalPoint Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint Maximum bracket evaluation point.
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws InvalidBracketRangeException Raised if the following condition is
     * not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public BrentSingleOptimizer(double minEvalPoint, double middleEvalPoint,
            double maxEvalPoint, double tolerance) 
            throws InvalidBracketRangeException {
        super(minEvalPoint, middleEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }
        
    /**
     * Constructor. Creates an instance with provided bracket of values and a
     * listener to get single dimension function evaluations.
     * @param listener Listener to evaluate a function.
     * @param minEvalPoint Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint Maximum bracket evaluation point.
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws InvalidBracketRangeException Raised if the following condition is
     * not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public BrentSingleOptimizer(
            SingleDimensionFunctionEvaluatorListener listener,
            double minEvalPoint, double middleEvalPoint, double maxEvalPoint,
            double tolerance) throws InvalidBracketRangeException {
        super(listener, minEvalPoint, middleEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }

    /**
     * Returns tolerance value, which is the accuracy to be obtained when a
     * minimum is estimated.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (and in 
     * such case, convergence will be assumed to have failed).
     * @return Tolerance value.
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * Sets algorithm's tolerance.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (an in 
     * such case, convergence will be assumed to have failed).
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws LockedException Raised if this instance is locked. This instance
     * will be locked while doing some operations. Attempting to change any 
     * parameter while being locked will raise this exception.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public void setTolerance(double tolerance)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * This function estimates a function minimum within provided or computed
     * bracket of values.
     * Given a function f, and given a bracketing triplet of abscissas ax, bx,
     * cx (such that bx is between ax and cx, and f(bx) is less than both f(ax)
     * and f(cx), this routine isolates the minimum to a fractional prevision of
     * about tolerance using Brent's method. The abscissa of the minimum is 
     * returned as xmin, and the function value of the minimum is returned as
     * fmin, the returned function value.
     * @throws LockedException Raised if this instance is locked, because
     * estimation is being computed.
     * @throws NotReadyException Raised if this instance is not ready because
     * either a listener or a bracket has not yet been provided or computed.
     * @throws OptimizationException Raised if the algorithm failed because of
     * lack of convergence or because function couldn't be evaluated.
     */
    @Override
    @SuppressWarnings("Duplicates")
    public void minimize() throws LockedException, NotReadyException,
            OptimizationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        double[] v1 = new double[1];
        double[] v2 = new double[2];
        double[] v3 = new double[3];        
        
        try {
            double a;
            double b;
            double d = 0.0;
            double etemp;
            double fu;
            double fv;
            double fw;
            double fx;
            double p;
            double q;
            double r;
            double tol1;
            double tol2;
            double u;
            double v;
            double w;
            double x;
            double xm;
            //This will be the distance moved on the step before last.
            double e = 0.0;
            
            //a and b must be in ascending order, but input abscissas need not
            //be.
            a = (ax < cx ? ax : cx);
            b = (bx > cx ? ax : cx);
            //Initializations...
            x = w = v = bx;
            fw = fv = fx = listener.evaluate(x);
            
            for (int iter = 0; iter < ITMAX; iter++) {
                //Main program loop
                xm = 0.5 * (a + b);
                tol1 = tolerance * Math.abs(x) + ZEPS;
                tol2 = 2.0 * tol1;
                
                if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                    //Test for done here.
                    fmin = fx;
                    xmin = x;
                    
                    resultAvailable = true;
                    locked = false;
                    return;
                }
                if (Math.abs(e) > tol1) {
                    //Construct a trial parabolic fit.
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    
                    if (q > 0.0) {
                        p = -p;
                    }
                    q = Math.abs(q);
                    etemp = e;
                    e = d;
                    
                    if (Math.abs(p) >= Math.abs(0.5 * q * etemp) || 
                            p <= q * (a - x) || p >= q * (b - x)) {
                        //noinspection all
                        e = x >= xm ? a - x : b - x;
                        d = CGOLD * (e);
                        //The above conditions determine the acceptability of
                        //the parabolic fit. Here we take the golden section 
                        //step into the larger of the two segments.
                    } else {
                        //Take the parabolic step
                        d = p / q;
                        u = x + d;
                        
                        if (u - a < tol2 || b - u < tol2) {
                            d = sign(tol1, xm - x);
                        }
                    }
                } else {
                    //noinspection all
                    e = x >= xm ? a - x : b - x;
                    d = CGOLD * (e);
                }
                
                u = Math.abs(d) >= tol1 ? x + d : x + sign(tol1, d);
                fu = listener.evaluate(u);
                
                //This is the one function evaluation per iteration
                if (fu <= fx) {
                    //No decide what to do with function evaluation
                    if (u >= x) {
                        a = x;
                    } else {
                        b = x;
                    }
                    //Housekeeping follows
                    v1[0] = v;
                    v2[0] = w;
                    v3[0] = x;
                    shft3(v1, v2, v3, u);
                    v = v1[0];
                    w = v2[0];
                    x = v3[0];
                    
                    v1[0] = fv;
                    v2[0] = fw;
                    v3[0] = fx;
                    shft3(v1, v2, v3, fu);
                    fv = v1[0];
                    fw = v2[0];
                    fx = v3[0];
                } else {
                    if (u < x) {
                        a = u;
                    } else {
                        b = u;
                    }
                    if (fu <= fw || w == x) {
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    } else if (fu <= fv || v == x || v == w) {
                        v = u;
                        fv = fu;
                    }
                }
                
                //Done with housekeeping. Back for another iteration.
            }
        } catch (EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }
        
        //Too many iterations in Brent!
        throw new OptimizationException();
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a minimum or not.
     * The instance is ready when both the listener and the bracket are 
     * available.
     * @return True if this instance is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable();
    }

    /**
     * Internal method to set algorithm tolerance. This method does not check
     * whether this instance is locked or not.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (and in
     * such case, convergence will be assumed to have failed).
     * @param tolerance Tolerance or accuracy to be obtained in estimated
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    private void internalSetTolerance(double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }
}
