/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.roots.LaguerrePolynomialRootsEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 13, 2012
 */
package com.irurueta.numerical.roots;

import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import java.util.Arrays;

/**
 * This class estimates the roots of a polynomial of degree n.
 * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
 * then the array of parameters is [an, a(n-1), ... a1, a0]
 * This class supports polynomials having either real or complex parameters
 */
public class LaguerrePolynomialRootsEstimator extends PolynomialRootsEstimator{

    //In this implementation we have increased MR and MT to increase accuracy
    //by iterating a larger but finite number of times
    
    /**
     * Constant that affects the number of iterations
     */
    public static final int MR = 80;
    
    /**
     * Constant that affects the number of iterations
     */
    public static final int MT = 100;
    
    /**
     * Maximum number of iterations
     */
    public static final int MAXIT = MT * MR;
    
    /**
     * Constant considered as machine precision for Laguerre method
     */
    public static final double LAGUER_EPS = 1e-10;
    
    /**
     * Constant considered as machine precision
     */
    public static final double EPS = 1e-14;
    
    /**
     * Constant indicating whether roots will be refined
     */
    public static final boolean DEFAULT_POLISH_ROOTS = true;
    
    /**
     * Minimum allowed length in polynomial parameters
     */
    public static final int MIN_VALID_POLY_PARAMS_LENGTH = 2;
    
    /**
     * Array containing values for Laguerre method
     */
    private static final double[] frac = 
        {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
    
    
    /**
     * Indicates if roots should be refined
     */
    private boolean polishRoots;
    
    /**
     * Constructor
     * @param polishRoots Boolean to determine whether roots should be refined
     */
    public LaguerrePolynomialRootsEstimator(boolean polishRoots){
        super();
        this.polishRoots = polishRoots;
    }
    
    /**
     * Empty constructor
     */
    public LaguerrePolynomialRootsEstimator(){
        super();
        this.polishRoots = DEFAULT_POLISH_ROOTS;
    }
    
    /**
     * Constructor
     * @param polyParams Array containing polynomial parameters
     * @param polishRoots Boolean indicating whether roots will be refined
     * @throws IllegalArgumentException Raised if length of provided parameters
     * is not valid. It has to be greater or equal than 2
     */
    public LaguerrePolynomialRootsEstimator(Complex[] polyParams, 
            boolean polishRoots) throws IllegalArgumentException{
        super();
        this.polishRoots = polishRoots;
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Constructor
     * @param polyParams Array containing polynomial parameters
     * @throws IllegalArgumentException Raised if length of provided parameters
     * is not valid. It has to be greater or equal than 2
     */
    public LaguerrePolynomialRootsEstimator(Complex[] polyParams)
            throws IllegalArgumentException{
        super();
        this.polishRoots = DEFAULT_POLISH_ROOTS;
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Internal method to set parameters of a polynomial, taking into account 
     * that a polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * then the array of parameters is [an, a(n - 1), ... a1, a0]
     * Polynomial parameters can be either real or complex values
     * This method does not check if this class is locked.
     * @param polyParams Polynomial parameters
     * @throws IllegalArgumentException Raised if the length of the array is not
     * valid
     */    
    @Override
    protected final void internalSetPolynomialParameters(
            Complex[] polyParams) throws IllegalArgumentException {
        if(polyParams.length < MIN_VALID_POLY_PARAMS_LENGTH)
            throw new IllegalArgumentException();
        this.polyParams = polyParams;
    }
    
    /**
     * Estimates the roots of provided polynomial.
     * @throws LockedException Raised if this instance is locked estimating a
     * root.
     * @throws NotReadyException Raised if this instance is not ready because
     * polynomial parameters have not been provided
     * @throws RootEstimationException Raised if roots cannot be estimated for
     * some reason (lack of convergence, etc)
     */    
    @Override
    public void estimate() throws LockedException, NotReadyException,
            RootEstimationException{
        
        if(isLocked()) throw new LockedException();
        if(!isReady()) throw new NotReadyException();
        
        //polynomial must be at least degree 1
        if(polyParams.length < MIN_VALID_POLY_PARAMS_LENGTH)
            throw new RootEstimationException();
        
        locked = true;
        
        Complex[] a = polyParams;
        roots = new Complex[a.length - 1];
        
        int i;
        int[] its = new int[1];
        Complex x = new Complex();
        Complex b, c;
        int m = a.length - 1;
        
        Complex[] ad = Arrays.copyOf(a, a.length);
        for(int j = m - 1; j >= 0; j--){
            x.setRealAndImaginary(0.0, 0.0);
            Complex[] ad_v = Arrays.copyOf(ad, j + 2);
            internalLaguer(ad_v, x, its);
            if(Math.abs(x.getImaginary()) <= 2.0 * EPS * Math.abs(x.getReal()))
                x.clone().setRealAndImaginary(x.getReal(), 0.0);
            roots[j] = x.clone();
            b = ad[j + 1].clone();
            for(int jj = j; jj >= 0; jj--){
                c = ad[jj].clone();
                ad[jj] = b.clone();
                b.multiply(x);
                b.add(c);
            }
        }
        if(polishRoots)
            for(int j = 0; j < m; j++)
                internalLaguer(a, roots[j], its);
        for(int j = 1; j < m; j++){
            x = roots[j].clone();
            for(i = j - 1; i >= 0; i--){
                if(roots[i].getReal() <= x.getReal()) break;
                roots[i + 1] = roots[i].clone();
            }
            roots[i + 1] = x.clone();
        }
        
        locked = false;
    }
    
    /**
     * Returns boolean indicating whether roots are refined after an initial
     * estimation
     * @return True if roots are refined, false otherwise
     */
    public boolean areRootsPolished(){
        return polishRoots;
    }
    
    /**
     * Sets boolean indicating whether roots will be refined after an initial
     * estimation
     * @param enable True if roots will be refined, false otherwise
     * @throws LockedException Raised if this instance is locked
     */
    public void setPolishRootsEnabled(boolean enable) throws LockedException{
        if(isLocked()) throw new LockedException();
        polishRoots = enable;
    }
    
    /**
     * Internal method to compute a root after decomposing and decreasing the
     * degree of the polynomial.
     * @param a Remaining polynomial parameters (on 1st iteration, the whole
     * polynomial is provided, on subsequent iterations, the polynomial is
     * deflated and the degree is reduced)
     * @param x Estimated root
     * @param its number of iterations needed to achieve the estimation
     * @throws RootEstimationException Raised if root couldn't be estimated 
     * because of lack of convergence
     */
    private void internalLaguer(Complex[] a, Complex x, int[] its) 
            throws RootEstimationException{
        
        Complex x1, b, g, g2;
        Complex dx = new Complex();
        Complex d = new Complex();
        Complex f = new Complex();
        Complex h = new Complex();
        Complex sq = new Complex();
        Complex gp = new Complex();
        Complex gm = new Complex();
        int m = a.length - 1;
        for(int iter = 1; iter <= MAXIT; iter++){
            its[0] = iter;
            b = a[m].clone();
            double err = b.getModulus();
            d.setRealAndImaginary(0.0, 0.0);
            f.setRealAndImaginary(0.0, 0.0);
            double abx = x.getModulus();
            for(int j = m - 1; j >= 0; j--){
                f.multiply(x);
                f.add(d);
                
                d.multiply(x);
                d.add(b);
                
                b.multiply(x);
                b.add(a[j]);
                
                err = b.getModulus() + abx * err;
            }
            err *= LAGUER_EPS;
            if(b.getModulus() <= err) return;
            g = d.divideAndReturnNew(b);
            g2 = g.powAndReturnNew(2.0);
            f.divide(b, h);
            h.multiplyByScalar(-2.0);
            h.add(g2);
            
            h.multiplyByScalar(m, sq);
            sq.subtract(g2);
            sq.multiplyByScalar(m - 1);
            sq.sqrt();
            
            g.add(sq, gp);
            g.subtract(sq, gm);
            
            double abp = gp.getModulus();
            double abm = gm.getModulus();
            if(abp < abm) gp = gm;
            if(Math.max(abp, abm) > 0.0){
                dx.setRealAndImaginary(m, 0.0);
                dx.divide(gp);
            }else{
                dx.setModulusAndPhase(1.0 + abx, iter);
            }
            x1 = x.subtractAndReturnNew(dx);
            if(x.equals(x1)) return;
            if(iter % MT != 0) x.copyFrom(x1);
            else{
                int pos = Math.min(iter/MT, frac.length - 1);
                x.subtract(dx.multiplyByScalarAndReturnNew(frac[pos]));
            }
        }
        //too many iterations in Laguerre
        locked = false;
        throw new RootEstimationException();
    }    
}
