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
import com.irurueta.numerical.NotReadyException;

/**
 * Abstract class to find function minima. Implementations will take into account whether the function is single or
 * multi-dimension, and will use different algorithms to find minima.
 */
public abstract class Optimizer {

    /**
     * Boolean indicating whether this instance is locked because computations are being done.
     */
    protected boolean locked;

    /**
     * Listener to handle minimization events.
     */
    protected OnIterationCompletedListener iterationCompletedListener;

    /**
     * Empty constructor.
     */
    protected Optimizer() {
        locked = false;
    }

    /**
     * Gets listener to handle minimization events.
     *
     * @return listener to handle minimization events.
     */
    public OnIterationCompletedListener getOnIterationCompletedListener() {
        return iterationCompletedListener;
    }

    /**
     * Sets listener to handle minimization events.
     *
     * @param iterationCompletedListener listener to handle minimization events.
     * @throws LockedException Raised if this instance is locked, because estimation is being computed.
     */
    public void setOnIterationCompletedListener(
            final OnIterationCompletedListener iterationCompletedListener)
            throws LockedException {
        if (locked) {
            throw new LockedException();
        }
        this.iterationCompletedListener = iterationCompletedListener;
    }

    /**
     * Returns boolean indicating whether this instance is locked. This instance will be locked while computations are being done. Attempting
     * to change any parameter while this instance is locked will raise a LockedException.
     *
     * @return True if this instance is locked, false otherwise.
     */
    public boolean isLocked() {
        return locked;
    }

    /**
     * This function estimates a function minimum. Implementations of this class will usually search a local minimum within a bracket of input
     * values. Because this is an abstract class, this method is meant to be overridden, otherwise a NotReadyException will always be thrown.
     *
     * @throws LockedException       Raised if this instance is locked, because estimation is being computed.
     * @throws NotReadyException     Raised if this instance is not ready, usually because listener has not yet been provided.
     * @throws OptimizationException Raised if the algorithm failed because of lack of convergence or because function couldn't be
     *                               evaluated.
     */
    public void minimize() throws LockedException, NotReadyException,
            OptimizationException {
        throw new NotReadyException();
    }

    /**
     * Returns boolean indicating whether this instance is ready. Usually an instance will be ready once its listener has been provided.
     * Because this is an abstract class, it will always return false;
     *
     * @return True if this instance is ready, false otherwise.
     */
    public boolean isReady() {
        return false;
    }
}
