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
package com.irurueta.numerical;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is more accurate although slower than 
 * DerivativeEstimator.
 */
@SuppressWarnings("WeakerAccess")
public class SymmetricDerivativeEstimator extends DerivativeEstimator {
    
    /**
     * Constructor.
     * @param listener listener to evaluate a single dimension function.
     */    
    public SymmetricDerivativeEstimator(
            SingleDimensionFunctionEvaluatorListener listener) {
        super(listener);
    }
    
    /**
     * Computes the function derivative at provided point x.
     * @param x Point where derivative is estimated.
     * @return Derivative of function at provided point.
     * @throws EvaluationException Raised if function cannot be properly 
     * evaluated.
     */    
    @Override
    public double derivative(double x) throws EvaluationException {
        try {
            double h = EPS * Math.abs(x);
            if (h == 0.0) {
                h = EPS; //Trick to reduce finite-precision error
            }
            
            double xh1 = x + h;
            double xh2 = x - h;
            //because of machine precision h could be different in both cases
            double h1 = xh1 - x;
            double h2 = x - xh2;
            
            double hh = h1 + h2; //this is more or less equal to 2.0 * h
            
            double fh1 = listener.evaluate(xh1);
            double fh2 = listener.evaluate(xh2);
            
            return (fh1 - fh2) / hh;
        } catch (Throwable t) {
            throw new EvaluationException(t);
        }
    }
}
