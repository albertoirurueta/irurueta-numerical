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
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.numerical.robust.LMedSRobustEstimator;
import com.irurueta.numerical.robust.LMedSRobustEstimatorListener;
import com.irurueta.numerical.robust.RobustEstimator;
import com.irurueta.numerical.robust.RobustEstimatorException;
import com.irurueta.numerical.robust.RobustEstimatorMethod;

import java.util.ArrayList;
import java.util.List;

/**
 * Finds the best polynomial using LMedS algorithm.
 */
public class LMedSPolynomialRobustEstimator extends PolynomialRobustEstimator {

    /**
     * Default value to be used for stop threshold. Stop threshold can be used
     * to keep the algorithm iterating in case that best estimated threshold
     * using median of residuals is not small enough. Once a solutions is found
     * that generates a threshold below this value, the algorithm will stop.
     * Threshold will be used to compare either algebraic or geometric distance
     * of estimated polynomial respect each provided evaluation.
     */
    public static final double DEFAULT_STOP_THRESHOLD = 1e-6;

    /**
     * Minimum value that can be set as stop threshold.
     * Threshold must be strictly greater than 0.0.
     */
    public static final double MIN_STOP_THRESHOLD = 0.0;

    /**
     * Threshold to be used to keep the algorithm iterating in case that best
     * estimated threshold using median of residuals is not small enough. Once
     * a solution is found that generates a threshold below this value, the
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in case where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would
     * iterate for a long time trying to find the best solution when indeed
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the sopt threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * strill produce even smaller thresholds in estimated results.
     */
    private double mStopThreshold;

    /**
     * Constructor.
     */
    public LMedSPolynomialRobustEstimator() {
        super();
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMedSPolynomialRobustEstimator(final int degree) {
        super(degree);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param evaluations collections of polynomial evaluations.
     * @throws IllegalArgumentException if provided number of evaluations is
     *                                  less than the required minimum.
     */
    public LMedSPolynomialRobustEstimator(
            final List<PolynomialEvaluation> evaluations) {
        super(evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     */
    public LMedSPolynomialRobustEstimator(
            final PolynomialRobustEstimatorListener listener) {
        super(listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     *                                  provided number of evaluations is less than the required minimum for
     *                                  provided degree.
     */
    public LMedSPolynomialRobustEstimator(
            final int degree,
            final List<PolynomialEvaluation> evaluations) {
        super(degree, evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMedSPolynomialRobustEstimator(
            final int degree,
            final PolynomialRobustEstimatorListener listener) {
        super(degree, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided number of evaluations is
     *                                  less than the required minimum.
     */
    public LMedSPolynomialRobustEstimator(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener) {
        super(evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     *                                  provided number of evaluations is less than the required minimum for
     *                                  provided degree.
     */
    public LMedSPolynomialRobustEstimator(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener) {
        super(degree, evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Returns threshold to be used to keep the algorithm iterating in case that
     * best estimated threshold using median of residuals is not small enough.
     * Once a solution is found that generates a threshold below this value, the
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would
     * iterate for a long time trying to find the best solution when indeed
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results.
     *
     * @return stop threshold to stop the algorithm prematurely when a certain
     * accuracy has been reached.
     */
    public double getStopThreshold() {
        return mStopThreshold;
    }

    /**
     * Sets threshold to be used to keep the algorithm iterating in case that
     * best estimated threshold using median of residuals is not small enough.
     * Once a solution is found that generates a threshold below this value, the
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would
     * iterate for a long time trying to find the best solution when indeed
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results
     *
     * @param stopThreshold stop threshold to stop the algorithm prematurely
     *                      when a certain accuracy has been reached
     * @throws IllegalArgumentException if provided value is zero or negative
     * @throws LockedException          if robust estimator is locked because an
     *                                  estimation is already in progress
     */
    public void setStopThreshold(final double stopThreshold)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (stopThreshold <= MIN_STOP_THRESHOLD) {
            throw new IllegalArgumentException();
        }
        mStopThreshold = stopThreshold;
    }

    /**
     * Estimates polynomial.
     *
     * @return estimated polynomial.
     * @throws LockedException          if robust estimator is locked because an
     *                                  estimation is already in progress.
     * @throws NotReadyException        if provided input data is not enough to start
     *                                  the estimation.
     * @throws RobustEstimatorException if estimation fails for any other reason
     *                                  (i.e. numerical instability, no solution available, etc).
     */
    @Override
    public Polynomial estimate() throws LockedException, NotReadyException,
            RobustEstimatorException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        final LMedSRobustEstimator<Polynomial> innerEstimator =
                new LMedSRobustEstimator<>(
                        new LMedSRobustEstimatorListener<Polynomial>() {

                            // subset of evaluations picked on each iteration
                            private final List<PolynomialEvaluation> mSubsetEvaluations =
                                    new ArrayList<>();

                            @Override
                            public int getTotalSamples() {
                                return mEvaluations.size();
                            }

                            @Override
                            public int getSubsetSize() {
                                return mPolynomialEstimator.getMinNumberOfEvaluations();
                            }

                            @Override
                            public void estimatePreliminarSolutions(
                                    final int[] samplesIndices,
                                    final List<Polynomial> solutions) {
                                mSubsetEvaluations.clear();
                                for (final int samplesIndex : samplesIndices) {
                                    mSubsetEvaluations.add(mEvaluations.get(samplesIndex));
                                }

                                try {
                                    mPolynomialEstimator.setLMSESolutionAllowed(false);
                                    mPolynomialEstimator.setEvaluations(mSubsetEvaluations);

                                    final Polynomial polynomial = mPolynomialEstimator.estimate();
                                    solutions.add(polynomial);
                                } catch (final Exception e) {
                                    // if anything fails, no solution is added
                                }
                            }

                            @Override
                            public double computeResidual(
                                    final Polynomial currentEstimation, final int i) {
                                final PolynomialEvaluation eval = mEvaluations.get(i);
                                return getDistance(eval, currentEstimation);
                            }

                            @Override
                            public boolean isReady() {
                                return LMedSPolynomialRobustEstimator.this.isReady();
                            }

                            @Override
                            public void onEstimateStart(
                                    final RobustEstimator<Polynomial> estimator) {
                                if (mListener != null) {
                                    mListener.onEstimateStart(
                                            LMedSPolynomialRobustEstimator.this);
                                }
                            }

                            @Override
                            public void onEstimateEnd(
                                    final RobustEstimator<Polynomial> estimator) {
                                if (mListener != null) {
                                    mListener.onEstimateEnd(
                                            LMedSPolynomialRobustEstimator.this);
                                }
                            }

                            @Override
                            public void onEstimateNextIteration(
                                    final RobustEstimator<Polynomial> estimator,
                                    final int iteration) {
                                if (mListener != null) {
                                    mListener.onEstimateNextIteration(
                                            LMedSPolynomialRobustEstimator.this, iteration);
                                }
                            }

                            @Override
                            public void onEstimateProgressChange(
                                    final RobustEstimator<Polynomial> estimator,
                                    final float progress) {
                                if (mListener != null) {
                                    mListener.onEstimateProgressChange(
                                            LMedSPolynomialRobustEstimator.this, progress);
                                }
                            }
                        });

        try {
            mLocked = true;
            innerEstimator.setConfidence(mConfidence);
            innerEstimator.setMaxIterations(mMaxIterations);
            innerEstimator.setProgressDelta(mProgressDelta);
            innerEstimator.setStopThreshold(mStopThreshold);
            return innerEstimator.estimate();
        } finally {
            mLocked = false;
        }
    }

    /**
     * Returns method being used for robust estimation.
     *
     * @return method being used for robust estimation.
     */
    @Override
    public RobustEstimatorMethod getMethod() {
        return RobustEstimatorMethod.LMedS;
    }
}
