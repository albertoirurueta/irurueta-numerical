/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.optimization.QuasiNewtonMultiOptimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 7, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.algebra.Matrix;
import com.irurueta.numerical.GradientFunctionEvaluatorListener;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;

/**
 * This class searches for a multi dimension function local minimum.
 * The local minimum is searched by starting the algorithm at a start point
 * The implementation of this class is based on Numerical Recipes 3rd ed. 
 * Section 10.9 page 521.
 */
public class QuasiNewtonMultiOptimizer extends MultiOptimizer {
    /**
     * Matimum number of iterations.
     */
    public static final int ITMAX = 200;
    
    /**
     * Machine precision.
     */
    public static final double EPS = 1e-12;
    
    /**
     * Convergence criterion on x values.
     */
    public static final double TOLX = 4.0 * EPS;
    
    /**
     * Scaled maximum step length allowed in line searches.
     */
    public static final double STPMX = 100.0;
    public static final double ALF = 1e-4;
    public static final double TOLX2 = 1e-12;
    
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
     * n-dimensional start point.
     */
    private double[] p;
    
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
     * Listener to obtain gradient values for the multi dimension function being
     * evaluated.
     * If the gradient is unknown (e.g. doesn't have a closed expression), the
     * provided listener could use a GradientEstimator to obtain one.
     */    
    GradientFunctionEvaluatorListener gradientListener;
    
    /**
     * value of the function at the minimum.
     */
    double fret;
    
    /**
     * Empty constructor.
     */
    public QuasiNewtonMultiOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
        p = null;
        iter = 0;
        gradientListener = null;
    }
     
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param gradientListener Listener to obtain gradient value for the 
     * multidimensional function being evaluated.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public QuasiNewtonMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener,
            GradientFunctionEvaluatorListener gradientListener, 
            double tolerance) throws IllegalArgumentException {
        super(listener);
        internalSetTolerance(tolerance);
        p = null;
        iter = 0;
        this.gradientListener = gradientListener;
    }
        
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param gradientListener Listener to obtain gradient value for the 
     * multidimensional function being evaluated.
     * @param startPoint Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public QuasiNewtonMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener,
            GradientFunctionEvaluatorListener gradientListener,
            double[] startPoint, double tolerance)
            throws IllegalArgumentException {
        super(listener);
        internalSetTolerance(tolerance);
        p = startPoint;
        iter = 0;
        this.gradientListener = gradientListener;
    }

    /**
     * This function estimates a function minimum.
     * Implementations of this class will usually search a local minimum close
     * to a start point.
     * Given a starting point p[0..n-1], the Broyden-Fletcher-Goldfarb-Sharno
     * variant of Davidon Fletcher-Powell minimization is performed on a 
     * function whose value and gradient are provided by respective listeners.
     * The convergence requirement on zeroing the gradient is provided by 
     * tolerance. This method estimates the location of the minimum, sets the
     * number of iterations that were required and the minimum value of the 
     * function at the minimum.
     * The internal method lnsrch is called within this method to perform
     * approximate line minimizations.
     * @throws LockedException Raised if this instance is locked, because
     * estimation is being computed.
     * @throws NotReadyException Raised if this instance is not ready, because
     * a listener, a gradient listener and a start point haven't been provided
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
        
        boolean validResult = false;
        int n = p.length;
        
        try {
            double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, 
                    temp, test;
            double[] dg = new double[n];
            double[] g = new double[n];
            double[] hdg = new double[n];
            double[] pnew = new double[n];
            double[] xi = new double[n];
            double[] fretArray = new double[1];
            boolean[] check = new boolean[1];
            
            Matrix hessin = new Matrix(n, n);
            fp = listener.evaluate(p);
            gradientListener.evaluateGradient(p, g);
            
            for (int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) hessin.setElementAt(i, j, 0.0);
                hessin.setElementAt(i, i, 1.0);
                xi[i] = -g[i];
                sum += p[i] * p[i];
            }
            stpmax = STPMX * Math.max(Math.sqrt(sum), (double)n);
            for (int its = 0; its < ITMAX; its++) {
                iter = its;
                lnsrch(p, fp, g, xi, pnew, fretArray, stpmax, check);
                fret = fretArray[0];
                fp = fret;
                for (int i = 0; i < n; i++) {
                    xi[i] = pnew[i] - p[i];
                    p[i] = pnew[i];
                }
                test = 0.0;
                for (int i = 0; i < n; i++) {
                    temp = Math.abs(xi[i]) / Math.max(Math.abs(p[i]), 1.0);
                    if(temp > test) test = temp;
                }
                if (test < TOLX) {
                    //minimum found
                    validResult = true;
                    break;
                }
                
                for (int i = 0; i < n; i++) {
                    dg[i] = g[i];
                }
                
                gradientListener.evaluateGradient(p, g);
                
                test = 0.0;
                den = Math.max(Math.abs(fret), 1.0);
                
                for (int i = 0; i < n; i++) {
                    temp = Math.abs(g[i]) * Math.max(Math.abs(p[i]), 1.0) / den;
                    if (temp > test) {
                        test = temp;
                    }
                }
                if (test < tolerance) {
                    //minimum found
                    validResult = true;
                    break;
                }
                
                for (int i = 0; i < n; i++) {
                    dg[i] = g[i] - dg[i];
                }
                
                for (int i = 0; i < n; i++) {
                    hdg[i] = 0.0;
                    for (int j = 0; j < n; j++) {
                        hdg[i] += hessin.getElementAt(i, j) * dg[j];
                    }
                }
                fac = fae = sumdg = sumxi = 0.0;
                for (int i = 0; i < n; i++) {
                    fac += dg[i] * xi[i];
                    fae += dg[i] * hdg[i];
                    sumdg += sqr(dg[i]);
                    sumxi += sqr(xi[i]);
                }
                
                if (fac > Math.sqrt(EPS * sumdg * sumxi)) {
                    fac = 1.0 / fac;
                    fad = 1.0 / fae;
                    for (int i = 0; i < n; i++) {
                        dg[i] = fac * xi[i] - fad * hdg[i];
                    }
                    
                    for (int i = 0; i < n; i++) {
                        for (int j = i; j < n; j++) {
                            hessin.setElementAt(i, j, 
                                    hessin.getElementAt(i, j) + 
                                    fac * xi[i] * xi[j] - 
                                    fad * hdg[i] * hdg[j] + 
                                    fae * dg[i] * dg[j]);
                            hessin.setElementAt(j, i, 
                                    hessin.getElementAt(i, j));
                        }
                    }
                }
                for (int i = 0; i < n; i++) {
                    xi[i] = 0.0;
                    for (int j = 0; j < n; j++) {
                        xi[i] -= hessin.getElementAt(i, j) * g[j];
                    }
                }
            }
            
            if (!validResult) {
                //too many iterations
                locked = false;
                throw new OptimizationException();
            }
        } catch (Throwable t) {
            locked = false;
            throw new OptimizationException(t);
        }
        
        //set result
        xmin = p;
        resultAvailable = true;
        fmin = fret;
        
        locked = false;
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a local minimum.
     * An instance is ready once a listener, a gradient listener and a start 
     * point are provided.
     * @return True if this instance is ready, false otherwise.
     */    
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isGradientListenerAvailable() &&
                isStartPointAvailable();
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
    private void internalSetTolerance(double tolerance) 
            throws IllegalArgumentException {
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
    public void setTolerance(double tolerance) throws LockedException,
            IllegalArgumentException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * Returns boolean indicating whether start point has already been provided
     * and is ready for retrieval.
     * @return True if available, false otherwise.
     */    
    public boolean isStartPointAvailable() {
        return p != null;
    }
    
    /**
     * Sets start point where algorithm will be started. Start point should be
     * close to the local minimum to be found.
     * @param startPoint Start point where algorithm will be started.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setStartPoint(double[] startPoint) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        p = startPoint;
    }
    
    /**
     * Returns start point where algorithm will be started. Start point should 
     * be close to the local minimum to be found.
     * @return Start point where algorithm will be started.
     * @throws NotAvailableException Raised if start point has not yet been
     * provided and is not available.
     */    
    public double[] getStartPoint() throws NotAvailableException {
        if (!isStartPointAvailable()) {
            throw new NotAvailableException();
        }
        return p;
    }
    
    /**
     * Returns gradient listener in charge of obtaining gradient values for the
     * function to be evaluated.
     * @return Gradient listener.
     * @throws NotAvailableException Raised if gradient listener has not yet 
     * been provided.
     */    
    public GradientFunctionEvaluatorListener getGradientListener()
            throws NotAvailableException {
        if (!isGradientListenerAvailable()) {
            throw new NotAvailableException();
        }
        return gradientListener;
    }
    
    /**
     * Sets gradient listener in charge of obtaining gradient values for the
     * function to be evaluated.
     * @param gradientListener Gradient listener.
     * @throws LockedException Raised if this instance is locked.
     */    
    public void setGradientListener(
            GradientFunctionEvaluatorListener gradientListener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.gradientListener = gradientListener;
    }
    
    /**
     * Returns boolean indicating whether a gradient listener has already been
     * provided and is available for retrieval.
     * @return True if available, false otherwise.
     */    
    public boolean isGradientListenerAvailable() {
        return gradientListener != null;
    }
    
    /**
     * Computes the squared value of provided double.
     * @param x Value to be squared.
     * @return Squared value.
     */    
    private double sqr(double x) {
        return x * x;
    }
    
    /**
     * Internal method to search for a minimum along a line.
     * NOTE: comments of params might be incorrect.
     * @param xold previous point.
     * @param fold previous function evaluation.
     * @param g gradient.
     * @param p point.
     * @param x point.
     * @param f function.
     * @param stpmax stores point maximum.
     * @param check checks line search.
     * @throws OptimizationException Raised if convergence was not achieved.
     */
    private void lnsrch(double[] xold, double fold, double[] g, double[] p, 
            double[] x, double[] f, double stpmax, boolean[] check) 
            throws OptimizationException {
        double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
        double rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
        int i, n = xold.length;
        check[0] = false;
        
        for (i = 0; i < n; i++) {
            sum += p[i] * p[i];
        }
        sum = Math.sqrt(sum);
        
        if (sum > stpmax) {
            for(i = 0; i < n; i++) {
                p[i] *= stpmax / sum;
            }
        }
        for(i = 0; i < n; i++) {
            slope += g[i] * p[i];
        }
        if (slope >= 0.0) {
            throw new OptimizationException();
        }
        test = 0.0;
        for (i = 0; i < n; i++) {
            temp = Math.abs(p[i]) / Math.max(Math.abs(xold[i]), 1.0);
            if(temp > test) test = temp;
        }
        alamin = TOLX2 / test;
        alam = 1.0;
        for (;;) {
            for (i = 0; i < n; i++) {
                x[i] = xold[i] + alam * p[i];
            }
            try {
                f[0] = listener.evaluate(x);
            } catch (Throwable t) {
                throw new OptimizationException(t);
            }
            
            if (alam < alamin) {
                for(i = 0; i < n; i++) x[i] = xold[i];
                check[0] = true;
                return;
            } else if (f[0] <= fold + ALF * alam * slope) {
                return;
            } else {
                if (alam == 1.0) {
                    tmplam = -slope / (2.0 * (f[0] - fold - slope));
                } else {
                    rhs1 = f[0] - fold - alam * slope;
                    rhs2 = f2 - fold - alam2 * slope;
                    a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) /
                            (alam - alam2);
                    b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 /
                            (alam2 * alam2)) / (alam - alam2);
                    if (a == 0.0) {
                        tmplam = -slope / (2.0 * b);
                    } else {
                        disc = b * b - 3.0 * a * slope;
                        if (disc < 0.0) {
                            tmplam = 0.5 * alam;
                        } else if (b <= 0.0) {
                            tmplam = (-b + Math.sqrt(disc)) /
                                (3.0 * a);
                        } else {
                            tmplam = -slope / (b + Math.sqrt(disc));
                        }
                    }
                    if (tmplam > 0.5 * alam) {
                        tmplam = 0.5 * alam;
                    }
                }
            }
            alam2 = alam;
            f2 = f[0];
            alam = Math.max(tmplam, 0.1 * alam);
        }
    }
}
