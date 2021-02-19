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

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Abstract class to find function roots of a single dimension function using
 * also its derivative information.
 * This class is meant to be extended by final implementations.
 */
@SuppressWarnings("WeakerAccess")
public abstract class DerivativeSingleRootEstimator 
    extends BracketedSingleRootEstimator {
    
    /**
     * Listener to evaluate a function's derivative. If the function's 
     * derivative is not known (e.g. a closed expression is not available), then
     * a DerivativeEstimator can be used inside the derivative listener 
     * implementation.
     */
    protected SingleDimensionFunctionEvaluatorListener derivativeListener;

    /**
     * Empty constructor.
     */    
    public DerivativeSingleRootEstimator() {
        super();
        derivativeListener = null;
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */    
    public DerivativeSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener,
            final double minEvalPoint, final double maxEvalPoint)
            throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        derivativeListener = null;
    }

    /**
     * Constructor
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param derivativeListener Listener to evaluate the function's derivative
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */        
    public DerivativeSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener,
            final SingleDimensionFunctionEvaluatorListener derivativeListener,
            final double minEvalPoint, final double maxEvalPoint)
            throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        this.derivativeListener = derivativeListener;
    }

    /**
     * Returns derivative listener to evaluate a function's derivative. 
     * If the function's derivative is not known (e.g. a closed expression is 
     * not available), then a DerivativeEstimator can be used inside the 
     * derivative listener implementation.
     * @return Derivative listener.
     * @throws NotAvailableException if listener is not available for retrieval.
     */
    public SingleDimensionFunctionEvaluatorListener getDerivativeListener()
            throws NotAvailableException {
        if (!isDerivativeListenerAvailable()) {
            throw new NotAvailableException();
        }
        return derivativeListener;
    }
    
    /**
     * Sets derivative listener to evaluate a function's derivative.
     * If the function's derivative is not known (e.g. a closed expression is 
     * not available), then a DerivativeEstimator can be used inside the 
     * derivative listener implementation.
     * @param derivativeListener Derivative listener to be set.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setDerivativeListener(
            final SingleDimensionFunctionEvaluatorListener derivativeListener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.derivativeListener = derivativeListener;
    }
    
    /**
     * Returns boolean indicating whether the derivative listener has been 
     * provided and is available for retrieval.
     * @return true if derivative listener is available, false otherwise
     */
    public boolean isDerivativeListenerAvailable() {
        return derivativeListener != null;
    }
    
    /**
     * Returns boolean indicating whether enough parameters have been provided
     * in order to start the estimation of the roots of a function.
     * An instance of this class is assumed to be ready when a listener, a
     * derivative listener and a bracket have been provided or computed.
     * @return True if this instance is ready to start the root estimation,
     * false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable() &&
                isDerivativeListenerAvailable();
    }
}
