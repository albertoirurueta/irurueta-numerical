/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.roots.ThirdDegreePolynomialRootsEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 13, 2012
 */
package com.irurueta.numerical.roots;

import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;

/**
 * Class to estimate the roots of a third degree polynomial along with other
 * polynomial properties.
 * A second degree polynomial is defined by its parameters as p(x) = a * x^3 + 
 * b * x^2 + c * x + d, hence the polynomial can be simply be defined by an 
 * array of length 4 [d, c, b, a]
 * This class is based on:
 * http://en.wikipedia.org/wiki/Cubic_function
 */
public class ThirdDegreePolynomialRootsEstimator 
    extends PolynomialRootsEstimator{

    /**
     * Constant defining machine precision
     */        
    public static final double EPS = 1e-10;
    
    /**
     * Constant defining one third
     */
    public static final double THIRD = 1.0 / 3.0;
    
    /**
     * Constant defining the squared root of three
     */
    public static double ROOTTHREE = Math.sqrt(3.0);
    
    /**
     * Number of parameters valid for a third degree polynomial.
     * 
     */        
    public static final int VALID_POLY_PARAMS_LENGTH = 4;
    
    /**
     * Array containing parameters of a second degree polynomial
     */        
    private double[] realPolyParams;
    
    /**
     * Empty constructor
     */        
    public ThirdDegreePolynomialRootsEstimator(){
        super();
        realPolyParams = null;
    }
    
    /**
     * Constructor
     * @param polyParams Array containing polynomial parameters
     * @throws IllegalArgumentException  Raised if the length of the provided
     * array is not valid.
     */        
    public ThirdDegreePolynomialRootsEstimator(double[] polyParams) 
            throws IllegalArgumentException{
        super();      
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Set array of third degree polynomial parameters.
     * A third degree polynomial is defined by p(x) = a * x^3 + b * x^2 + c * x 
     * + d, and the array must be provided as [d, c, b, a].
     * Note: This class only supports real polynomial parameters
     * @param polyParams Array containing polynomial parameters
     * @throws LockedException Raised if this instance is locked
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */    
    public void setPolynomialParameters(double[] polyParams)
            throws LockedException, IllegalArgumentException{
        if(isLocked()) throw new LockedException();
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Internal method to set array of third degree polynomial parameters.
     * A third degree polynomial is defined by p(x) = a * x^3 + b * x^2 + c * d 
     * + d, and the array must be provided as [d, c, b, a].
     * Note: This class only supports real polynomial parameters
     * This method does not check if this instance is locked.
     * @param polyParams Array containing polynomial parameters
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */        
    private void internalSetPolynomialParameters(double[] polyParams) 
            throws IllegalArgumentException{
        if(polyParams.length < VALID_POLY_PARAMS_LENGTH)
            throw new IllegalArgumentException();
        if(!isThirdDegree(polyParams)) throw new IllegalArgumentException();
        
        this.realPolyParams = polyParams;
    }
    
    /**
     * Returns array of third degree polynomial parameters.
     * A third degree polynomial is defined by p(x) = a * x^3 + b * x^2 + c * x 
     * + d, and the array is returned as [d, c, b, a].
     * Note: This class only supports real polynomial parameters
     * @return Array of first degree polynomial parameters
     * @throws NotAvailableException Raised if polynomial parameter have not yet
     * been provided
     */        
    public double[] getRealPolynomialParameters() throws NotAvailableException{
        if(!arePolynomialParametersAvailable()) 
            throw new NotAvailableException();
        return realPolyParams;
    }
    
    /**
     * Returns boolean indicating whether REAL polynomial parameters have been 
     * provided and is available for retrieval.
     * Note: This class only supports real polynomial parameters
     * @return True if available, false otherwise
     */        
    @Override
    public boolean arePolynomialParametersAvailable(){
        return realPolyParams != null;
    }

    /**
     * This method will always raise a NotAvailableException because this class
     * only supports REAL polynomial parameters
     * @return this method always throws an exception
     * @throws com.irurueta.numerical.NotAvailableException always thrown
     * @deprecated
     */        
    @Override
    public Complex[] getPolynomialParameters() throws NotAvailableException{
        throw new NotAvailableException();
    }
    

    /**
     * This method will always raise an IllegalArgumentException because this
     * class only supports REAL polynomial parameters
     * @throws IllegalArgumentException always thrown
     * @deprecated 
     */        
    @Override
    protected void internalSetPolynomialParameters(Complex[] polyParams) 
            throws IllegalArgumentException {
        //complex values are not supported
        throw new IllegalArgumentException();
    }
    
    /**
     * Estimates the roots of provided polynomial.
     * @throws LockedException Raised if this instance is locked estimating 
     * roots.
     * @throws NotReadyException Raised if this instance is not ready because
     * polynomial parameters have not been provided
     * @throws RootEstimationException Raised if roots cannot be estimated for
     * some reason
     */        
    @Override
    public void estimate() throws LockedException, NotReadyException,
            RootEstimationException{
        if(isLocked()) throw new LockedException();
        if(!isReady()) throw new NotReadyException();
        
        locked = true;
        
        roots = new Complex[VALID_POLY_PARAMS_LENGTH - 1];
        
        double d = realPolyParams[0];
        double c = realPolyParams[1];
        double b = realPolyParams[2];
        double a = realPolyParams[3];
        
        Complex x1 = new Complex();
        Complex x2 = new Complex();
        Complex x3 = new Complex();
        solveCubic(a, b, c, d, x1, x2, x3);
        
        if(Double.isNaN(x1.getReal()) || Double.isNaN(x1.getImaginary()) ||
           Double.isNaN(x2.getReal()) || Double.isNaN(x2.getImaginary()) ||
           Double.isNaN(x3.getReal()) || Double.isNaN(x3.getImaginary())){

            locked = false;
            throw new RootEstimationException();
        }

        if(x1.getReal() < x2.getReal() && x1.getReal() < x3.getReal()){
            //x1 goes first
            roots[0] = x1;
            if(x2.getReal() < x3.getReal()){
                //x2 goes second and x3 goes third
                roots[1] = x2;
                roots[2] = x3;
            }else{
                //x3 goes second and x2 goes third
                roots[1] = x3;
                roots[2] = x2;
            }
        }else if(x2.getReal() < x1.getReal() && x2.getReal() < x3.getReal()){
            //x2 goes first
            roots[0] = x2;
            if(x1.getReal() < x3.getReal()){
                //x1 goes second and x3 goes third
                roots[1] = x1;
                roots[2] = x3;
            }else{
                //x3 goes second and x1 goes third
                roots[1] = x3;
                roots[2] = x1;
            }
        }else{
            //x3 goes first
            roots[0] = x3;
            if(x1.getReal() < x2.getReal()){
                //x1 goes second and x2 goes third
                roots[1] = x1;
                roots[2] = x2;
            }else{
                //x2 goes second and x1 goes second
                roots[1] = x2;
                roots[2] = x1;
            }
        }

        locked = false;
    }
    
    /**
     * Returns boolean indicating whether provided array of polynomial 
     * parameters correspond to a valid third degree polynomial.
     * A third degree polynomial is defined by p(x) = a * x^3 + b * x^2 + c *x + 
     * d, and the array is returned as [d, c, b, a].
     * Note: This class only supports real polynomial parameters
     * @param polyParams Array containing polynomial parameters
     * @return True if is a third degree polynomial, false otherwise
     */        
    public static boolean isThirdDegree(double[] polyParams){
        int length = polyParams.length;
        if(length >= VALID_POLY_PARAMS_LENGTH){
            if(Math.abs(polyParams[VALID_POLY_PARAMS_LENGTH - 1]) > EPS){
                for(int i = VALID_POLY_PARAMS_LENGTH; i < length; i++){
                    if(Math.abs(polyParams[i]) > EPS) return false;
                }
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether polynomial parameters provided to this
     * instance correspond to a valid third degree polynomial.
     * A third degree polynomial is defined by p(x) = a * x^3 + b * x^2 + c * x 
     * + d, and the array is returned as [d, c, b, a].
     * Note: This class only supports real polynomial parameters
     * @return True if is a second degree polynomial, false otherwise
     * @throws NotReadyException Raised if this instance is not ready because
     * an array of polynomial parameters has not yet been provided.
     */        
    public boolean isThirdDegree() throws NotReadyException{
        if(!isReady()) throw new NotReadyException();
        return isThirdDegree(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are three
     * distinct and real roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for third degree polynomials that have three distinct roots,
     * they must be either real or one real and 2 complex conjugate.
     * @param polyParams Array containing polynomial parameters
     * @return True if roots are distinct and real, false otherwise
     */    
    public static boolean hasThreeDistinctRealRoots(double[] polyParams){
        if(polyParams.length >= VALID_POLY_PARAMS_LENGTH)
            return getDiscriminant(polyParams) > EPS;
        return false;
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are three
     * distinct and real roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for third degree polynomials that have three distinct roots,
     * they must be either real or one real and 2 complex conjugate.
     * @return True if roots are distinct and real, false otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet 
     * been provided
     */    
    public boolean hasThreeDistinctRealRoots() throws NotReadyException{
        if(!isReady()) throw new NotReadyException();
        return hasThreeDistinctRealRoots(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether the polynomial has two real and equal 
     * roots and a third different one (multiplicity 2), or all three roots are 
     * real and equal (multiplicity 3).
     * @param polyParams Array containing polynomial parameters
     * @return True if there are roots with multiplicity greater than one, false 
     * otherwise
     */
    public static boolean hasMultipleRealRoot(double[] polyParams){
        if(polyParams.length >= VALID_POLY_PARAMS_LENGTH)
            return Math.abs(getDiscriminant(polyParams)) <= EPS;
        return false;
    }
    
    /**
     * Returns boolean indicating whether the polynomial has two real and equal
     * roots and a third different one (multiplicity 2), or all three roots are 
     * real and equal (multiplicity 3).
     * @return True if there are roots with multiplicity greater than one, false
     * otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet
     * been provided
     */
    public boolean hasMultipleRealRoot() throws NotReadyException{
        if(!isReady()) throw new NotReadyException();
        return hasMultipleRealRoot(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether the polynomial has one real root and
     * two complex conjugate roots.
     * @param polyParams Array containing polynomial parameters
     * @return True if polynomial has 1 real root and 2 complex conjugate roots,
     * false otherwise
     */
    public static boolean hasOneRealRootAndTwoComplexConjugateRoots(
            double[] polyParams){
        if(polyParams.length >= VALID_POLY_PARAMS_LENGTH)
            return getDiscriminant(polyParams) < -EPS;
        return false;
    }
    
    /**
     * Returns boolean indicating whether the polynomial has one real root and
     * two complex conjugate roots.
     * @return True if polynomial has 1 real root and 2 complex conjugate roots,
     * false otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet
     * been provided
     */
    public boolean hasOneRealRootAndTwoComplexConjugateRoots()
            throws NotReadyException{
        if(!isReady()) throw new NotReadyException();        
        return hasOneRealRootAndTwoComplexConjugateRoots(realPolyParams);
    }
    
    /**
     * Internal method to compute the discriminant of a 3rd degree polynomial.
     * Discriminants are helpful to determine properties of a 3rd degree 
     * polynomial
     * @param polyParams Array containing polynomial parameters
     * @return Value of discriminant
     */    
    private static double getDiscriminant(double[] polyParams){
        double d = polyParams[0];
        double c = polyParams[1];
        double b = polyParams[2];
        double a = polyParams[3];

        return 18.0 * a * b * c * d
            - 4.0 * b * b * b * d
            + b * b * c * c
            - 4.0 * a * c * c * c
            - 27.0 * a * a * d * d;
    }
    
    /**
     * Computes the cube root or x^(1/3) of provided value x
     * @param x Provided value
     * @return Cube root
     */
    private double cubeRoot(double x){
        if(x < 0.0) return -Math.pow(-x, THIRD);
        else return Math.pow(x, THIRD);
    }
    
    /**
     * Finds 3rd degree polynomial roots
     * @param a 1st parameter
     * @param b 2nd parameter
     * @param c 3rd parameter
     * @param d 4th parameter
     * @param x1 1st root (output parameter)
     * @param x2 2nd root (output parameter)
     * @param x3 3rd root (output parameter)
     */
    private void solveCubic(double a, double b, double c, double d,
            Complex x1, Complex x2, Complex x3){
        
        //find the discriminant
        double f, g, h;
        f = (3.0 * c / a - Math.pow(b, 2.0) / Math.pow(a, 2.0)) / 3.0;
        g = (2.0 * Math.pow(b, 3.0) / Math.pow(a, 3.0) - 9.0 * b * c / 
                Math.pow(a, 2.0) + 27.0 * d / a) / 27.0;
        h = Math.pow(g, 2.0) / 4.0 + Math.pow(f, 3.0) / 27.0;
        double absF = Math.abs(f);
        double absG = Math.abs(g);
        double absH = Math.abs(h);
        // evaluate discriminant
        if (absF <= EPS && absG <= EPS && absH <= EPS)
        {
            // 3 equal roots
            double x;
            // when f, g, and h all equal 0 the roots can be found by the following line
            x = -cubeRoot(d / a);
            x1.setRealAndImaginary(x, 0.0);
            x2.setRealAndImaginary(x, 0.0);
            x3.setRealAndImaginary(x, 0.0);
        }
        else if (h <= 0.0)
        {
            // 3 real roots
            double i, j, k, m, n, p;
            // complicated maths making use of the method
            i = Math.pow(Math.pow(g, 2.0) / 4 - h, 0.5);
            j = cubeRoot(i);
            k = Math.acos(-(g / (2.0 * i)));
            m = Math.cos(k / 3.0);
            n = ROOTTHREE * Math.sin(k / 3.0);
            p = -(b / (3.0 * a));
            // print solutions
            x1.setRealAndImaginary(2.0 * j * m + p, 0.0);
            x2.setRealAndImaginary(-j * (m + n) + p, 0.0);
            x3.setRealAndImaginary(-j * (m - n) + p, 0.0);
        }
        else if (h > 0)
        {
            // 1 real root and 2 complex roots
            double r, s, t, u, p;
            // complicated maths making use of the method
            r = -(g / 2) + Math.pow(h, 0.5);
            s = cubeRoot(r);
            t = -(g / 2) - Math.pow(h, 0.5);
            u = cubeRoot(t);
            p = -(b / (3 * a));
            // print solutions
            x1.setRealAndImaginary((s + u) + p, 0.0);
            double real = -(s + u) / 2 + p;
            x2.setReal(real);
            x3.setReal(real);
            double imag = (s - u) * ROOTTHREE / 2;
            x2.setImaginary(imag);
            x3.setImaginary(imag);
        }
    }
}
