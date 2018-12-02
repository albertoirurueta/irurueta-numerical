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
package com.irurueta.numerical.roots;

import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;

/**
 * Class to estimate the root of a first degree polynomial along with other
 * polynomial properties.
 * A first degree polynomial is defined by its parameters as p(x) = a * x + b,
 * hence the polynomial can be simply be defined by an array of length 2 [b, a].
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class FirstDegreePolynomialRootsEstimator 
    extends PolynomialRootsEstimator {
    
    /**
     * Constant defining machine precision.
     */
    public static final double EPS = 1e-10;
    
    /**
     * Number of parameters valid for a first degree polynomial.
     * 
     */
    public static final int VALID_POLY_PARAMS_LENGTH = 2;
    
    /**
     * Array containing parameters of a first degree polynomial.
     */
    private double[] realPolyParams;    
    
    /**
     * Empty constructor.
     */
    public FirstDegreePolynomialRootsEstimator() {
        super();
        realPolyParams = null;                 
    }
    
    /**
     * Constructor.
     * @param polyParams Array containing polynomial parameters.
     * @throws IllegalArgumentException  Raised if the length of the provided
     * array is not valid.
     */
    public FirstDegreePolynomialRootsEstimator(double[] polyParams) {
        super();
        internalSetPolynomialParameters(polyParams);
    }

    /**
     * Set array of first degree polynomial parameters.
     * A first degree polynomial is defined by p(x) = a * x + b, and the array
     * must be provided as [b, a].
     * Note: This class only supports real polynomial parameters.
     * @param polyParams Array containing polynomial parameters.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */
    public void setPolynomialParameters(double[] polyParams)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Internal method to set array of first degree polynomial parameters.
     * A first degree polynomial is defined by p(x) = a * x + b, and the array
     * must be provided as [b, a].
     * Note: This class only supports real polynomial parameters.
     * This method does not check if this instance is locked.
     * @param polyParams Array containing polynomial parameters.
     * @throws IllegalArgumentException Raised if the length of the provided
     * array is not valid.
     */
    private void internalSetPolynomialParameters(double[] polyParams) {
        if (polyParams.length < VALID_POLY_PARAMS_LENGTH) {
            throw new IllegalArgumentException();
        }
        if (!isFirstDegree(polyParams)) {
            throw new IllegalArgumentException();
        }
        
        this.realPolyParams = polyParams;
    }
    
    /**
     * Returns array of first degree polynomial parameters.
     * A first degree polynomial is defined by p(x) = a * x + b, and the array
     * is returned as [b, a].
     * Note: This class only supports real polynomial parameters.
     * @return Array of first degree polynomial parameters.
     * @throws NotAvailableException if parameters are not available for retrieval.
     */
    public double[] getRealPolynomialParameters() throws NotAvailableException {
        if(!arePolynomialParametersAvailable()) {
            throw new NotAvailableException();
        }
        return realPolyParams;
    }
    
    /**
     * Returns boolean indicating whether REAL polynomial parameters have been 
     * provided and is available for retrieval.
     * Note: This class only supports real polynomial parameters.
     * @return True if available, false otherwise.
     */
    @Override
    public boolean arePolynomialParametersAvailable() {
        return realPolyParams != null;
    }

    /**
     * This method will always raise a NotAvailableException because this class
     * only supports REAL polynomial parameters.
     * @return throws NotAvailableException.
     * @throws NotAvailableException alwasy throws this exception.
     * @deprecated
     */
    @Override
    @Deprecated
    public Complex[] getPolynomialParameters() throws NotAvailableException {
        throw new NotAvailableException();
    }
    

    /**
     * This method will always raise an IllegalArgumentException because this
     * class only supports REAL polynomial parameters.
     * @deprecated 
     */
    @Override
    @Deprecated
    protected void internalSetPolynomialParameters(Complex[] polyParams) {
        //complex values are not supported
        throw new IllegalArgumentException();
    }
    
    /**
     * Estimates the root of provided polynomial.
     * @throws LockedException Raised if this instance is locked estimating a
     * root.
     * @throws NotReadyException Raised if this instance is not ready because
     * polynomial parameters have not been provided.
     */
    @Override
    public void estimate() throws LockedException, NotReadyException {
        
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        
        roots = new Complex[VALID_POLY_PARAMS_LENGTH - 1];
        
        double b = realPolyParams[0];
        double a = realPolyParams[1];
        
        double x = solveLinear(a, b);
        
        roots[0] = new Complex(x, 0.0);
        
        locked = false;
    }    
    
    /**
     * Returns boolean indicating whether provided array of polynomial 
     * parameters correspond to a valid first degree polynomial.
     * A first degree polynomial is defined by p(x) = a * x + b, and the array
     * is returned as [b, a].
     * Note: This class only supports real polynomial parameters
     * @param polyParams Array containing polynomial parameters
     * @return True if is a first degree polynomial, false otherwise
     */
    public static boolean isFirstDegree(double[] polyParams) {
        int length = polyParams.length;
        if (length >= VALID_POLY_PARAMS_LENGTH &&
                Math.abs(polyParams[VALID_POLY_PARAMS_LENGTH - 1]) > EPS) {
            for (int i = VALID_POLY_PARAMS_LENGTH; i < length; i++) {
                if (Math.abs(polyParams[i]) > EPS) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
    
    /**
     * Returns boolean indicating whether polynomial parameters provided to this
     * instance correspond to a valid first degree polynomial.
     * A first degree polynomial is defined by p(x) = a * x + b, and the array
     * is returned as [b, a].
     * Note: This class only supports real polynomial parameters.
     * @return True if is a first degree polynomial, false otherwise.
     * @throws NotReadyException Raised if this instance is not ready because
     * an array of polynomial parameters has not yet been provided.
     */
    public boolean isFirstDegree() throws NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }
        return isFirstDegree(realPolyParams);
    }
    
    /**
     * Returns boolean indicating whether estimated root is real.
     * Because this class only accepts real polynomial parameters, then the
     * estimated root will always be real, and consequently this method always
     * returns true.
     * @return True if estimated root is real, false otherwise.
     */
    public boolean isRealSolution() {
        return true;
    }
    
    /**
     * Internal method to estimate a root on a first degree polynomial.
     * @param a A parameter.
     * @param b B parameter.
     * @return Root.
     */
    private double solveLinear(double a, double b) {
        //a * x + b = 0
        return -b / a;
    }
}
