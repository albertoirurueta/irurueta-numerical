/*
 * Copyright (C) 2013 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;

/**
 * Robust estimator to estimate some object in a robust manner
 *
 * @param <T> Object to be estimated (i.e. lines, cameras, etc)
 */
@SuppressWarnings("WeakerAccess")
public abstract class RobustEstimator<T> {

    /**
     * Default amount of progress variation before notifying a change in
     * estimation progress. By default this is set to 5%
     */
    public static final float DEFAULT_PROGRESS_DELTA = 0.05f;

    /**
     * Minimum allowed value for progress delta
     */
    public static final float MIN_PROGRESS_DELTA = 0.0f;

    /**
     * Maximum allowed value for progress delta.
     */
    public static final float MAX_PROGRESS_DELTA = 1.0f;

    /**
     * Listener to be notified of events such as when estimation starts, ends
     * or its progress significantly changes.
     */
    protected RobustEstimatorListener<T> mListener;

    /**
     * Indicates if this estimator is locked because an estimation is being
     * computed.
     */
    protected volatile boolean mLocked;

    /**
     * Amount of progress variation before notifying a progress change during
     * estimation.
     */
    protected float mProgressDelta;

    /**
     * Constructor.
     */
    public RobustEstimator() {
        mListener = null;
        mLocked = false;
        mProgressDelta = DEFAULT_PROGRESS_DELTA;
    }

    /**
     * Constructor.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     */
    public RobustEstimator(final RobustEstimatorListener<T> listener) {
        mListener = listener;
        mLocked = false;
        mProgressDelta = DEFAULT_PROGRESS_DELTA;
    }

    /**
     * Returns reference to listener to be notified of events such as when
     * estimation starts, ends or its progress significantly changes.
     *
     * @return listener to be notified of events.
     */
    public RobustEstimatorListener<T> getListener() {
        return mListener;
    }

    /**
     * Sets listener to be notified of events such as when estimation starts,
     * ends or its progress significantly changes.
     *
     * @param listener listener to be notified of events.
     * @throws LockedException if robust estimator is locked.
     */
    public void setListener(final RobustEstimatorListener<T> listener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        mListener = listener;
    }

    /**
     * Indicates whether listener has been provided and is available for
     * retrieval.
     *
     * @return true if available, false otherwise.
     */
    public boolean isListenerAvailable() {
        return mListener != null;
    }

    /**
     * Indicates if this instance is locked because estimation is being computed.
     *
     * @return true if locked, false otherwise.
     */
    public boolean isLocked() {
        return mLocked;
    }

    /**
     * Returns amount of progress variation before notifying a progress change
     * during estimation.
     *
     * @return amount of progress variation before notifying a progress change
     * during estimation.
     */
    public float getProgressDelta() {
        return mProgressDelta;
    }

    /**
     * Sets amount of progress variation before notifying a progress change
     * during estimation.
     *
     * @param progressDelta amount of progress variation before notifying a
     *                      progress change during estimation.
     * @throws IllegalArgumentException if progress delta is less than zero or
     *                                  greater than 1.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    @SuppressWarnings("Duplicates")
    public void setProgressDelta(final float progressDelta) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (progressDelta < MIN_PROGRESS_DELTA ||
                progressDelta > MAX_PROGRESS_DELTA) {
            throw new IllegalArgumentException();
        }
        mProgressDelta = progressDelta;
    }

    /**
     * Robustly estimates an instance of T.
     *
     * @return estimated object.
     * @throws LockedException          if robust estimator is locked.
     * @throws NotReadyException        if provided input data is not enough to start
     *                                  the estimation.
     * @throws RobustEstimatorException if estimation fails for any reason
     *                                  (i.e. numerical instability, no solution available, etc).
     */
    public abstract T estimate() throws LockedException, NotReadyException,
            RobustEstimatorException;

    /**
     * Returns data about inliers once estimation has been done.
     *
     * @return data about inliers or null if estimation has not been done.
     */
    public abstract InliersData getInliersData();

    /**
     * Returns method being used for robust estimation.
     *
     * @return method being used for robust estimation.
     */
    public abstract RobustEstimatorMethod getMethod();

    /**
     * Indicates if estimator is ready to start the estimation process.
     *
     * @return true if ready, false otherwise.
     */
    public boolean isReady() {
        if (mListener != null) {
            return mListener.isReady();
        } else {
            return false;
        }
    }
}
