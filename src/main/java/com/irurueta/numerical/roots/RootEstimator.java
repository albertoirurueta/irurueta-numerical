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
import com.irurueta.numerical.NotReadyException;

/**
 * Abstract class to find roots of functions.
 * A root is the locus of points (set of points) where the value of a given
 * function equals to zero.
 * Usually root estimators will only find a single root around an initial
 * coarsely estimated solution.
 */
@SuppressWarnings("WeakerAccess")
public abstract class RootEstimator {

    /**
     * Boolean indicating that this instance is locked because it is doing
     * computations.
     * Attempting to change any parameters while being locked will raise a
     * LockedException.
     */
    protected boolean locked;

    /**
     * Constructor.
     */
    public RootEstimator() {
        locked = false;
    }

    /**
     * Returns boolean indicating whether this instance is locked.
     * An instance is locked while it is doing computations. Attempting to
     * change any parameters while being locked will raise a LockedException.
     *
     * @return Boolean indicating whether this instance is locked.
     */
    public boolean isLocked() {
        return locked;
    }

    /**
     * Estimates the root or roots for a given function.
     *
     * @throws LockedException         Exception raised if this instance is already
     *                                 locked.
     * @throws NotReadyException       Exception raised if not enough parameters have
     *                                 been provided in order to start the estimation.
     * @throws RootEstimationException Raised if the root estimation failed for
     *                                 some other reason (usually inability to evaluate the function,
     *                                 numerical instability or convergence problems, or no roots are found).
     */
    public void estimate() throws LockedException, NotReadyException,
            RootEstimationException {
    }

    /**
     * Returns boolean indicating whether enough parameters have been provided
     * in order to start the estimation of the roots of a function.
     *
     * @return True if this instance is ready to start the root estimation,
     * false otherwise.
     */
    public boolean isReady() {
        return false;
    }
}
