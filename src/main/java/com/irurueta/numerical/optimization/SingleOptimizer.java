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

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Abstract class to find minima on single dimension functions.
 * Single dimension functions are functions having a single parameter and 
 * returning a single scalar value, such as f(x).
 * 
 * Subclasses of this class will implement specific methods to find function
 * minima.
 */
@SuppressWarnings("WeakerAccess")
public abstract class SingleOptimizer extends Optimizer {
    
    /**
     * Listener to evaluate single dimension functions.
     */
    protected SingleDimensionFunctionEvaluatorListener listener;
    
    /**
     * Value where minimum has been found.
     */
    protected double xmin;
    
    /**
     * Function evaluation at minimum that has been found.
     */
    protected double fmin;
    
    /**
     * Boolean indicating whether a minimum has been found and is available for
     * retrieval.
     */
    protected boolean resultAvailable;
    
    /**
     * Empty constructor.
     */
    public SingleOptimizer() {
        super();
        xmin = fmin = 0.0;
        resultAvailable = false;
    }
    
    /**
     * Constructor with listener.
     * @param listener Listener to evaluate a single dimension function where
     * minima is meant to be found.
     */
    public SingleOptimizer(SingleDimensionFunctionEvaluatorListener listener) {
        super();
        this.listener = listener;
        xmin = fmin = 0.0;
        resultAvailable = false;
    }
    
    /**
     * Returns listener to evaluate a single dimension function.
     * @return Listener to evaluate a single dimension function.
     * @throws NotAvailableException Raised if a listener has not yet been
     * provided.
     */
    public SingleDimensionFunctionEvaluatorListener getListener() 
            throws NotAvailableException {
        if (!isListenerAvailable()) {
            throw new NotAvailableException();
        }
        
        return listener;
    }
    
    /**
     * Sets listener.
     * @param listener Listener to evaluate a single dimension function.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setListener(SingleDimensionFunctionEvaluatorListener listener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        
        this.listener = listener;
    }
    
    /**
     * Returns boolean indicating whether a listener has been provided and is
     * available for retrieval.
     * @return True if listener is available, false otherwise.
     */
    public boolean isListenerAvailable() {
        return listener != null;
    }
    
    /**
     * Returns true if this instance is ready to start the minimum estimation,
     * false otherwise.
     * @return Boolean indicating whether this instance is ready to start the
     * minimum estimation.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable();
    }
    
    /**
     * Returns boolean indicating whether the estimated minimum is available for
     * retrieval.
     * @return True if result is available, false otherwise.
     */
    public boolean isResultAvailable() {
        return resultAvailable;
    }
    
    /**
     * Returns value of the minimum that has been found.
     * @return Value of the minimum that has been found.
     * @throws NotAvailableException Raised if minimum is not yet available for
     * retrieval.
     */
    public double getResult() throws NotAvailableException {
        if (!isResultAvailable()) {
            throw new NotAvailableException();
        }
        
        return xmin;
    }
    
    /**
     * Returns function evaluation at minimum that has been found.
     * @return Function evaluation at minimum that has been found.
     * @throws NotAvailableException raised if result is not yet available for retrieval.
     */
    public double getEvaluationAtResult() throws NotAvailableException {
        if (!isResultAvailable()) {
            throw new NotAvailableException();
        }
        
        return fmin;
    }
}
