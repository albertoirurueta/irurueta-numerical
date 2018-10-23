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
package com.irurueta.numerical.roots;

import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;

/**
 * Class to estimate the roots of a second degree polynomial along with other
 * polynomial properties.
 * A second degree polynomial is defined by its parameters as p(x) = a * x^2 + 
 * b * x + c, hence the polynomial can be simply be defined by an array of 
 * length 3 [c, b, a]
 * This class is based on:
 * http://en.wikipedia.org/wiki/Quadratic_function
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class SecondDegreePolynomialRootsEstimator 
    extends PolynomialRootsEstimator {
    
    /**
     * Constant defining machine precision.
     */    
    public static final double EPS = 1e-10;
    
    /**
     * Number of parameters valid for a second degree polynomial.
     * 
     */    
    public static final int VALID_POLY_PARAMS_LENGTH = 3;
    
    /**
     * Array containing parameters of a second degree polynomial.
     */    
    private double[] realPolyParams;    
    
    /**
     * Empty constructor.
     */    
    public SecondDegreePolynomialRootsEstimator() {
        super();
        realPolyParams = null;                 
    }
    
    /**
     * Constructor.
     * @param polyParams Array containing polynomial parameters.
     * @throws IllegalArgumentException  Raised if the length of the provided
     * array is not valid.
     */    
    public SecondDegreePolynomialRootsEstimator(double[] polyParams)
            throws IllegalArgumentException {
        super();
        internalSetPolynomialParameters(polyParams);
    }

    /**
     * Set array of second degree polynomial parameters.
     * A second degree polynomial is defined by p(x) = a * x^2 + b * x + c, and 
     * the array must be provided as [c, b, a].
     * Note: This class only supports real polynomial parameters
     * @param polyParams Array containing polynomial parameters.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */
    public void setPolynomialParameters(double[] polyParams)
            throws LockedException, IllegalArgumentException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Internal method to set array of second degree polynomial parameters.
     * A second degree polynomial is defined by p(x) = a * x^2 + b * x + c, and 
     * the array must be provided as [c, b, a].
     * Note: This class only supports real polynomial parameters
     * This method does not check if this instance is locked.
     * @param polyParams Array containing polynomial parameters.
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */    
    private void internalSetPolynomialParameters(double[] polyParams) 
            throws IllegalArgumentException {
        if (polyParams.length < VALID_POLY_PARAMS_LENGTH) {
            throw new IllegalArgumentException();
        }
        if (!isSecondDegree(polyParams)) {
            throw new IllegalArgumentException();
        }
        
        this.realPolyParams = polyParams;
    }
    
    /**
     * Returns array of second degree polynomial parameters.
     * A second degree polynomial is defined by p(x) = a * x^2 + b * x + c, and 
     * the array is returned as [c, b, a].
     * Note: This class only supports real polynomial parameters
     * @return Array of first degree polynomial parameters
     * @throws NotAvailableException Raised if polynomial parameter have not yet
     * been provided
     */    
    public double[] getRealPolynomialParameters() throws NotAvailableException {
        if (!arePolynomialParametersAvailable()) {
            throw new NotAvailableException();
        }
        return realPolyParams;
    }
    
    /**
     * Returns boolean indicating whether REAL polynomial parameters have been 
     * provided and is available for retrieval.
     * Note: This class only supports real polynomial parameters
     * @return True if available, false otherwise
     */    
    @Override
    public boolean arePolynomialParametersAvailable() {
        return realPolyParams != null;
    }

    /**
     * This method will always raise a NotAvailableException because this class
     * only supports REAL polynomial parameters
     * @return always throws NotAvailableException
     * @throws NotAvailableException always thrown
     * @deprecated
     */    
    @Override
    public Complex[] getPolynomialParameters() throws NotAvailableException {
        throw new NotAvailableException();
    }
    

    /**
     * This method will always raise an IllegalArgumentException because this
     * class only supports REAL polynomial parameters
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
            RootEstimationException {
        
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        
        roots = new Complex[VALID_POLY_PARAMS_LENGTH - 1];
        
        double c = realPolyParams[0];
        double b = realPolyParams[1];
        double a = realPolyParams[2];
        
        Complex x1 = new Complex();
        Complex x2 = new Complex();
        solveQuadratic(a, b, c, x1, x2);
        
        if (Double.isNaN(x1.getReal()) || Double.isNaN(x1.getImaginary()) ||
                Double.isNaN(x2.getReal()) || Double.isNaN(x2.getImaginary())) {
            locked = false;
            throw new RootEstimationException();
        }
        if (x1.getReal() < x2.getReal()) {
            //x1 goes first
            roots[0] = x1;
            roots[1] = x2;
        } else {
            //x2 goes first
            roots[0] = x2;
            roots[1] = x1;
        }
        
        locked = false;
    }    
    
    /**
     * Returns boolean indicating whether provided array of polynomial 
     * parameters correspond to a valid second degree polynomial.
     * A second degree polynomial is defined by p(x) = a * x^2 + b * x + c, and 
     * the array is returned as [a, b, a].
     * Note: This class only supports real polynomial parameters
     * @param polyParams Array containing polynomial parameters
     * @return True if is a second degree polynomial, false otherwise
     */    
    public static boolean isSecondDegree(double[] polyParams) {
        int length = polyParams.length;
        if (length >= VALID_POLY_PARAMS_LENGTH) {
            if (Math.abs(polyParams[VALID_POLY_PARAMS_LENGTH - 1]) > EPS) {
                for (int i = VALID_POLY_PARAMS_LENGTH; i < length; i++) {
                    if (Math.abs(polyParams[i]) > EPS) {
                        return false;
                    }
                }
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether polynomial parameters provided to this
     * instance correspond to a valid second degree polynomial.
     * A second degree polynomial is defined by p(x) = a * x^2 + b * x + c, and 
     * the array is returned as [c, b, a].
     * Note: This class only supports real polynomial parameters
     * @return True if is a second degree polynomial, false otherwise
     * @throws NotReadyException Raised if this instance is not ready because
     * an array of polynomial parameters has not yet been provided.
     */    
    public boolean isSecondDegree() throws NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }
        return isSecondDegree(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are two
     * distinct and real roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for second degree polynomials that have two distinct roots,
     * its roots must be either real or complex conjugate.
     * @param polyParams Array containing polynomial parameters
     * @return True if roots are distinct and real, false otherwise
     */
    public boolean hasTwoDistinctRealRoots(double[] polyParams) {
        if (polyParams.length >= VALID_POLY_PARAMS_LENGTH) {
            return getDiscriminant(polyParams) > EPS;
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are two
     * distinct and real roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for second degree polynomials that have two distinct roots,
     * its roots must be either real or complex conjugate.
     * @return True if roots are distinct and real, false otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet 
     * been provided
     */
    public boolean hasTwoDistinctRealRoots() throws NotReadyException{
        if (!isReady()) {
            throw new NotReadyException();
        }
        return hasTwoDistinctRealRoots(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether a second degree polynomial has
     * multiple roots (for the 2nd degree case this means 2 equal roots)
     * This is true for polynomials of the form (x - r)^2 = x^2 - 2 * r * x + 
     * r^2, where r is the double root
     * @param polyParams Array containing polynomial parameters
     * @return True if it has double root, false otherwise
     */
    public static boolean hasDoubleRoot(double[] polyParams){
        if (polyParams.length >= VALID_POLY_PARAMS_LENGTH) {
            return Math.abs(getDiscriminant(polyParams)) <= EPS;
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether this second degree polynomial has
     * multiple roots (for the 2nd degree case this means 2 equal roots)
     * This is true for polynomials of the form (x - r)^2 = x^2 - 2 * r * x + 
     * r^2, where r is the double root
     * @return True if it has double root, false otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet
     * been provided.
     */
    public boolean hasDoubleRoot() throws NotReadyException{
        if (!isReady()) {
            throw new NotReadyException();
        }
        return hasDoubleRoot(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are two
     * complex conjugate roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for second degree polynomials that have two distinct roots,
     * its roots must be either real or complex conjugate.
     * @param polyParams Array containing polynomial parameters
     * @return True if roots are complex conjugate, false otherwise
     */
    public static boolean hasTwoComplexConjugateRoots(double[] polyParams){
        if (polyParams.length >= VALID_POLY_PARAMS_LENGTH) {
            return getDiscriminant(polyParams) < -EPS;
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether the roots of the polynomial are two
     * complex conjugate roots or not.
     * Because this class only supports polynomials with real parameters, we
     * know that for second degree polynomials that have two distinct roots,
     * its roots must be either real or complex conjugate.
     * @return True if roots are complex conjugate, false otherwise
     * @throws NotReadyException Raised if polynomial parameters haven't yet 
     * been provided
     */    
    public boolean hasTwoComplexConjugateRoots() throws NotReadyException{
        if (!isReady()) {
            throw new NotReadyException();
        }
        return hasTwoComplexConjugateRoots(realPolyParams);
    }
    
    /**
     * Internal method to compute the discriminant of a 2nd degree polynomial.
     * Discriminants are helpful to determine properties of a 2nd degree 
     * polynomial
     * @param polyParams Array containing polynomial parameters
     * @return Value of discriminant
     */
    private static double getDiscriminant(double[] polyParams) {

        double c = polyParams[0];
        double b = polyParams[1];
        double a = polyParams[2];

        return b * b - 4.0 * a * c;
    }
    
    /**
     * Finds 2nd degree polynomial roots
     * @param a 1st parameter
     * @param b 2nd parameter
     * @param c 3rd parameter
     * @param x1 1st root (output parameter)
     * @param x2 2nd root (output parameter)
     */
    private void solveQuadratic(double a, double b, double c, Complex x1,
            Complex x2) {
        
        double discriminant = b * b - 4.0 * a * c;
        
        if (discriminant >= 0.0) {
            //real solutions (double or distinct)
            x1.setRealAndImaginary((-b + Math.sqrt(discriminant)) / (2.0 * a), 
                    0.0);
            x2.setRealAndImaginary((-b - Math.sqrt(discriminant)) / (2.0 * a), 
                    0.0);            
        } else {
            //complex conjugate solutions
            double real = -b / (2.0 * a);
            double imag = Math.sqrt(Math.abs(discriminant)) / (2.0 * a);
            x1.setRealAndImaginary(real, imag);
            x2.setRealAndImaginary(real, -imag);
        }
    }
}
