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

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Abstract class to find roots of single dimension functions.
 * A root is the locus of points (set of points) where the value of a given
 * function equals to zero.
 * A single dimension function is one containing a single parameter and
 * returning a single value (i.e. f(x)).
 * Usually root estimators will only find a single root around an initial
 * coarsely estimated solution.
 */
public abstract class SingleRootEstimator extends RootEstimator {

    /**
     * Listener that evaluates a single dimension function in order to find its
     * root.
     */
    protected SingleDimensionFunctionEvaluatorListener listener;

    /**
     * Boolean indicating that a root has been computed and is available to be
     * retrieved.
     */
    protected boolean rootAvailable;

    /**
     * Root that has been found.
     */
    protected double root;

    /**
     * Empty constructor.
     */
    protected SingleRootEstimator() {
        super();
        listener = null;
        rootAvailable = false;
        root = 0.0;
    }

    /**
     * Constructor.
     *
     * @param listener Listener that evaluates a single dimension function in
     *                 order to find its root.
     */
    protected SingleRootEstimator(final SingleDimensionFunctionEvaluatorListener listener) {
        super();
        this.listener = listener;
        rootAvailable = false;
        root = 0.0;
    }

    /**
     * Returns listener that evaluates a single dimension function in order to
     * find its root.
     *
     * @return Listener that evaluates a single dimension function.
     * @throws NotAvailableException Raised if listener has not yet been
     *                               provided.
     */
    public SingleDimensionFunctionEvaluatorListener getListener() throws NotAvailableException {
        if (!isListenerAvailable()) {
            throw new NotAvailableException();
        }
        return listener;
    }

    /**
     * Sets listener that evaluates a single dimension function in order to find
     * its root.
     *
     * @param listener Listener that evaluates a single dimension function.
     * @throws LockedException Raised if this instance is already locked.
     */
    public void setListener(final SingleDimensionFunctionEvaluatorListener listener) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.listener = listener;
    }

    /**
     * Returns boolean indicating whether a listener has been provided.
     *
     * @return True if listener is available, false otherwise.
     */
    public boolean isListenerAvailable() {
        return listener != null;
    }

    /**
     * Returns boolean indicating whether enough parameters have been provided
     * in order to start the estimation of the roots of a function.
     *
     * @return True if this instance is ready to start the root estimation,
     * false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable();
    }

    /**
     * Returns boolean indicating whether a root has been estimated and is
     * available for retrieval.
     *
     * @return True if root is available, false otherwise.
     */
    public boolean isRootAvailable() {
        return rootAvailable;
    }

    /**
     * Returns estimated root for a single dimension function inside a given
     * bracket of values.
     *
     * @return Estimated root.
     * @throws NotAvailableException Raised if root has not yet been estimated.
     */
    public double getRoot() throws NotAvailableException {
        if (!isRootAvailable()) {
            throw new NotAvailableException();
        }
        return root;
    }
}
