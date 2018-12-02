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

/**
 * Abstract class to estimate the roots of a polynomial.
 */
@SuppressWarnings("WeakerAccess")
public abstract class PolynomialRootsEstimator extends RootEstimator {
    
    /**
     * Array containing parameters of a polynomial, taking into account that
     * a polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * then the array of parameters is [a0, a1, ... a(n - 1), an]
     */
    protected Complex[] polyParams;
    
    /**
     * Array containing estimated roots.
     */
    protected Complex[] roots;
    
    /**
     * Empty constructor.
     */
    public PolynomialRootsEstimator() {
        super();
        polyParams = roots = null;
    }
    
    /**
     * Returns array containing polynomial parameters.
     * @return Array of polynomial parameters.
     * @throws NotAvailableException Raised if polynomial parameters have not
     * been provided and are not available for retrieval.
     */
    public Complex[] getPolynomialParameters() throws NotAvailableException {
        if (!arePolynomialParametersAvailable()) {
            throw new NotAvailableException();
        }
        return polyParams;
    }

    /**
     * Sets parameters of a polynomial, taking into account 
     * that a polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * then the array of parameters is [a0, a1, ... a(n - 1), an]
     * @param polyParams Polynomial parameters.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if the length of the array is not
     * valid depending on the subclass implementation.
     */
    public void setPolynomialParameters(Complex[] polyParams) 
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetPolynomialParameters(polyParams);
    }
    
    /**
     * Returns boolean indicating whether polynomial parameters have been
     * provided and are available for retrieval.
     * @return True if available, false otherwise.
     */
    public boolean arePolynomialParametersAvailable() {
        return polyParams != null;
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of the polynomial roots.
     * This instance is considered to be ready once polynomial parameters are
     * provided.
     * @return True if ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return arePolynomialParametersAvailable();
    }
    
    /**
     * Returns array of estimated polynomial roots.
     * The array will have a length equal to the polynomial degree.
     * If a polynomial has multiple roots, then such roots will be repeated.
     * @return Array of estimated polynomial roots
     * @throws NotAvailableException Raised if roots have not yet been estimated
     * and are not available for retrieval
     */
    public Complex[] getRoots() throws NotAvailableException {
        if (!areRootsAvailable()) {
            throw new NotAvailableException();
        }
        return roots;
    }
    
    /**
     * Returns boolean indicating whether roots have been estimated and are
     * available for retrieval.
     * @return True if available, false otherwise.
     */
    public boolean areRootsAvailable() {
        return roots != null;
    }

    /**
     * Internal method to set parameters of a polynomial, taking into account
     * that a polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * then the array of parameters is [a0, a1, ... a(n - 1), an]
     * This method does not check if this class is locked.
     * @param polyParams Polynomial parameters.
     * @throws IllegalArgumentException Raised if the length of the array is not
     * valid depending on the subclass implementation.
     */
    protected abstract void internalSetPolynomialParameters(
            Complex[] polyParams);
}
