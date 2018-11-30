/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

import com.irurueta.algebra.Complex;

/**
 * Utility class to evaluate complex polynomials.
 * This class is useful when the same complex polynomial needs to be evaluated
 * multiple times.
 */
@SuppressWarnings("WeakerAccess")
public class ComplexPolynomialEvaluator extends PolynomialEvaluator {
    
    /**
     * Polynomial coefficients.
     * A polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * Hence, the array of polynomial coefficients is [a0, a1, ... a(n-1), an]
     */
    private Complex[] mPolyParams;
    
    /**
     * Constructor.
     * @param polyParams polynomial coefficients
     * @throws IllegalArgumentException if provided array is null or has length 
     * 0.
     */
    public ComplexPolynomialEvaluator(Complex[] polyParams) {
        if (polyParams == null || polyParams.length == 0) {
            throw new IllegalArgumentException();
        }
        
        mPolyParams = polyParams;
    }
    
    /**
     * Gets polynomial parameters.
     * @return polynomial parameters.
     */
    public Complex[] getPolyParams() {
        return mPolyParams;
    }
    
    /**
     * Evaluates polynomial at provided point x.
     * @param x point where polynomial is evaluated.
     * @return result of evaluation.
     */
    public Complex evaluate(Complex x) {
        return evaluate(mPolyParams, x);
    }    
}
