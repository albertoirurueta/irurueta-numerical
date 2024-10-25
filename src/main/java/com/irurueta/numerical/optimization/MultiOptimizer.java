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
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;

/**
 * Abstract class to search for minima on multidimensional classes.
 * A multidimensional class is one having several input parameters (usually
 * provided as an array of values), and returning a single scalar value.
 */
public abstract class MultiOptimizer extends Optimizer {
    /**
     * Listener to evaluate a multidimensional function.
     */
    protected MultiDimensionFunctionEvaluatorListener listener;

    /**
     * Minimum that was estimated.
     */
    protected double[] xmin;

    /**
     * Function value at estimated minimum.
     */
    protected double fmin;

    /**
     * Boolean indicating whether a minimum has already been found or not.
     */
    protected boolean resultAvailable;

    /**
     * Empty constructor.
     */
    protected MultiOptimizer() {
        super();
        listener = null;
        xmin = null;
        fmin = 0.0;
    }

    /**
     * Constructor.
     *
     * @param listener Listener to evaluate a multidimensional function.
     */
    protected MultiOptimizer(final MultiDimensionFunctionEvaluatorListener listener) {
        super();
        this.listener = listener;
        xmin = null;
        fmin = 0.0;
    }

    /**
     * Returns listener to evaluate a multidimensional function
     *
     * @return Listener to evaluate a multidimensional function.
     * @throws NotAvailableException Raised if listener has not yet been
     *                               provided and is not available for retrieval.
     */
    public MultiDimensionFunctionEvaluatorListener getListener() throws NotAvailableException {
        if (!isListenerAvailable()) {
            throw new NotAvailableException();
        }
        return listener;
    }

    /**
     * Sets listener to evaluate a multidimensional function.
     *
     * @param listener Listener to evaluate a multidimensional function.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setListener(final MultiDimensionFunctionEvaluatorListener listener) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.listener = listener;
    }

    /**
     * Returns boolean indicating whether listener has been provided and is
     * available for retrieval.
     *
     * @return True if available, false otherwise.
     */
    public boolean isListenerAvailable() {
        return listener != null;
    }

    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a minimum.
     * Because this class is abstract, this method is meant to be overridden,
     * otherwise false will always be returned.
     *
     * @return Boolean indicating whether this instance is ready.
     */
    @Override
    public boolean isReady() {
        return false;
    }

    /**
     * Returns boolean indicating whether a minimum has been estimated and is
     * available for retrieval.
     *
     * @return True if result is available, false otherwise.
     */
    public boolean isResultAvailable() {
        return resultAvailable;
    }

    /**
     * Returns minimum point that was found.
     *
     * @return Minimum point
     * @throws NotAvailableException Raised if a minimum is not yet available
     *                               for retrieval.
     */
    public double[] getResult() throws NotAvailableException {
        if (!isResultAvailable()) {
            throw new NotAvailableException();
        }
        return xmin;
    }

    /**
     * Returns function evaluation at estimated minimum point.
     *
     * @return Function evaluation at minimum.
     * @throws NotAvailableException Raised if a minimum is not yet available
     *                               for retrieval.
     */
    public double getEvaluationAtResult() throws NotAvailableException {
        if (!isResultAvailable()) {
            throw new NotAvailableException();
        }
        return fmin;
    }
}
