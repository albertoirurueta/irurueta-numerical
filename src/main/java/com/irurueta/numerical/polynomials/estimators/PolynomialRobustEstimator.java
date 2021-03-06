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
import com.irurueta.numerical.robust.RobustEstimatorException;
import com.irurueta.numerical.robust.RobustEstimatorMethod;

import java.util.List;

/**
 * This is an abstract class for algorithms to robustly find the best
 * Polynomial for provided collection of evaluations.
 * Implementations of this class should be able to detect and discard outliers
 * in order to find the best solution.
 */
@SuppressWarnings("Duplicates")
public abstract class PolynomialRobustEstimator {

    /**
     * Default robust estimator method when none is provided.
     * In general for Polynomial estimation is best to use PROSAC or RANSAC
     * than any other method, as it provides more robust methods.
     */
    public static final RobustEstimatorMethod DEFAULT_ROBUST_METHOD =
            RobustEstimatorMethod.PROSAC;

    /**
     * Default amount of progress variation before notifying a change in
     * estimation progress. By default this is set to 5%.
     */
    public static final float DEFAULT_PROGRESS_DELTA = 0.05f;

    /**
     * Minimum allowed value for progress delta.
     */
    public static final float MIN_PROGRESS_DELTA = 0.0f;

    /**
     * Maximum allowed value for progress delta.
     */
    public static final float MAX_PROGRESS_DELTA = 1.0f;

    /**
     * Constant defining default confidence of the estimated result, which is
     * 99%. This means that with a probability of 99% estimation will be
     * accurate because chosen sub-samples will be inliers.
     */
    public static final double DEFAULT_CONFIDENCE = 0.99;

    /**
     * Default maximum allowed number of iterations.
     */
    public static final int DEFAULT_MAX_ITERATIONS = 5000;

    /**
     * Minimum allowed confidence value.
     */
    public static final double MIN_CONFIDENCE = 0.0;

    /**
     * Maximum allowed confidence value.
     */
    public static final double MAX_CONFIDENCE = 1.0;

    /**
     * Minimum allowed number of iterations.
     */
    public static final int MIN_ITERATIONS = 1;

    /**
     * Flag indicating whether geometric distance is used by default or not
     * to find outliers.
     */
    public static final boolean DEFAULT_USE_GEOMETRIC_DISTANCE = false;

    /**
     * Collection of polynomial evaluations and their corresponding point of
     * evaluation used to determine a polynomial of required degree.
     */
    protected List<PolynomialEvaluation> mEvaluations;

    /**
     * Internal non robust estimator of polynomial estimator.
     */
    protected final LMSEPolynomialEstimator mPolynomialEstimator;

    /**
     * Listener to be notified of events such as when estimation starts, ends or
     * its progress significantly changes.
     */
    protected PolynomialRobustEstimatorListener mListener;

    /**
     * Indicates if this estimator is locked because an estimation is being
     * computed.
     */
    protected boolean mLocked;

    /**
     * Amount of progress variation before notifying a progress change during
     * estimation.
     */
    protected float mProgressDelta;

    /**
     * Amount of confidence expressed as a value between 0.0 and 1.0 (which is
     * equivalent to 100%9. The amount of confidence indicates the probability
     * that the estimated result is correct. Usually this value will be close
     * to 1.0, but not exactly 1.0.
     */
    protected double mConfidence;

    /**
     * Maximum allowed number of iterations. When the maximum number of
     * iterations is exceeded, result will not be available, however an
     * approximate result will be available for retrieval.
     */
    protected int mMaxIterations;

    /**
     * Indicates whether geometric distance will be used to find outliers or
     * algebraic distance will be used instead.
     */
    protected boolean mUseGeometricDistance;

    /**
     * Constructor.
     */
    protected PolynomialRobustEstimator() {
        mProgressDelta = DEFAULT_PROGRESS_DELTA;
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        mUseGeometricDistance = DEFAULT_USE_GEOMETRIC_DISTANCE;
        mPolynomialEstimator = new LMSEPolynomialEstimator();
    }

    /**
     * Constructor.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialRobustEstimator(final int degree) {
        mProgressDelta = DEFAULT_PROGRESS_DELTA;
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        mUseGeometricDistance = DEFAULT_USE_GEOMETRIC_DISTANCE;
        mPolynomialEstimator = new LMSEPolynomialEstimator(degree);
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided number of evaluations is
     *                                  less than the required minimum.
     */
    protected PolynomialRobustEstimator(final List<PolynomialEvaluation> evaluations) {
        this();
        internalSetEvaluations(evaluations);
    }

    /**
     * Constructor.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     */
    protected PolynomialRobustEstimator(
            final PolynomialRobustEstimatorListener listener) {
        this();
        mListener = listener;
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
    protected PolynomialRobustEstimator(final int degree,
                                        final List<PolynomialEvaluation> evaluations) {
        this(degree);
        internalSetEvaluations(evaluations);
    }

    /**
     * Constructor.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialRobustEstimator(final int degree,
                                        final PolynomialRobustEstimatorListener listener) {
        this(degree);
        mListener = listener;
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
    protected PolynomialRobustEstimator(final List<PolynomialEvaluation> evaluations,
                                        final PolynomialRobustEstimatorListener listener) {
        this(evaluations);
        mListener = listener;
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     *                                  provided number of evaluations is less than the required minimum for
     *                                  provided degree.
     */
    protected PolynomialRobustEstimator(final int degree,
                                        final List<PolynomialEvaluation> evaluations,
                                        final PolynomialRobustEstimatorListener listener) {
        this(degree, evaluations);
        mListener = listener;
    }

    /**
     * Gets collection of polynomial evaluations and their corresponding point
     * of evaluation used to determine a polynomial of required degree.
     *
     * @return collection of polynomial evaluations.
     */
    public List<PolynomialEvaluation> getEvaluations() {
        return mEvaluations;
    }

    /**
     * Sets collection of polynomial evaluations and their corresponding point
     * of evaluation used to determine a polynomial of required degree.
     *
     * @param evaluations collection of polynomial evaluations.
     * @throws LockedException          if estimator is locked.
     * @throws IllegalArgumentException if provided list of evaluations does
     *                                  not containing enough evaluations to estimate the polynomial using current
     *                                  settings.
     */
    public void setEvaluations(final List<PolynomialEvaluation> evaluations)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetEvaluations(evaluations);
    }

    /**
     * Gets minimum number of evaluations required to estimate a polynomial of
     * the specified degree.
     *
     * @return number of required evaluations.
     */
    public int getMinNumberOfEvaluations() {
        return mPolynomialEstimator.getMinNumberOfEvaluations();
    }

    /**
     * Gets listener to be notified of events such as when estimation starts,
     * ends or its progress significantly changes.
     *
     * @return listener to be notified of events.
     */
    public PolynomialRobustEstimatorListener getListener() {
        return mListener;
    }

    /**
     * Sets listener to be notified of events such as when estimation starts,
     * ends or its progress significantly changes.
     *
     * @param listener listener to be notified of events.
     */
    public void setListener(final PolynomialRobustEstimatorListener listener) {
        mListener = listener;
    }

    /**
     * Indicates if this estimator is locked because an estimation is being
     * computed.
     *
     * @return true if this estimator is locked, false otherwise.
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
    public void setProgressDelta(final float progressDelta)
            throws LockedException {
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
     * Returns amount of confidence expressed as a value between 0.0 and 1.0
     * (which is equivalent to 100%). The amount of confidence indicates the
     * probability that the estimated result is correct. Usually this value will
     * be close to 1.0, but not exactly 1.0.
     *
     * @return amount of confidence as a value between 0.0 and 1.0.
     */
    public double getConfidence() {
        return mConfidence;
    }

    /**
     * Sets amount of confidence expressed as a value between 0.0 and 1.0 (which
     * is equivalent to 100%). The amount of confidence indicates the
     * probability that the estimated result is correct. Usually this value will
     * be close to 1.0 but not exactly 1.0.
     *
     * @param confidence confidence to be set as a value between 0.0 and 1.0.
     * @throws IllegalArgumentException if provided value is not between 0.0 and
     *                                  1.0.
     * @throws LockedException          if this estimator is locked because an estimator
     *                                  is being computed.
     */
    public void setConfidence(final double confidence)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (confidence < MIN_CONFIDENCE || confidence > MAX_CONFIDENCE) {
            throw new IllegalArgumentException();
        }
        mConfidence = confidence;
    }

    /**
     * Returns maximum allowed number of iterations. If maximum allowed number
     * of iterations is achieved without converging to a result when calling
     * estimate(), a RobustEstimatorException will be raised.
     *
     * @return maximum allowed number of iterations.
     */
    public int getMaxIterations() {
        return mMaxIterations;
    }

    /**
     * Sets maximum allowed number of iterations. When the maximum number of
     * iterations is exceeded, result will not be available, however an
     * approximate result will be available for retrieval.
     *
     * @param maxIterations maximum allowed number of iterations to be set.
     * @throws IllegalArgumentException if provided value is less than 1.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setMaxIterations(final int maxIterations)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxIterations < MIN_ITERATIONS) {
            throw new IllegalArgumentException();
        }
        mMaxIterations = maxIterations;
    }

    /**
     * Indicates whether geometric distance will be used to find outliers or
     * algebraic distance will be used instead.
     *
     * @return true if geometric distance is used, false otherwise.
     */
    public boolean isGeometricDistanceUsed() {
        return mUseGeometricDistance;
    }

    /**
     * Specifies whether geometric distance will be used to find outliers or
     * algebraic distance will be used instead.
     *
     * @param geometricDistanceUsed true if geometric distance is used, false
     *                              otherwise.
     * @throws LockedException if this estimator is locked.
     */
    public void setGeometricDistanceUsed(final boolean geometricDistanceUsed)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        mUseGeometricDistance = geometricDistanceUsed;
    }

    /**
     * Gets degree of polynomial to be estimated.
     *
     * @return degree of polynomial to be estimated.
     */
    public int getDegree() {
        return mPolynomialEstimator.getDegree();
    }

    /**
     * Sets degree of polynomial to be estimated.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     * @throws LockedException          if this instance is locked.
     */
    public void setDegree(final int degree) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        mPolynomialEstimator.setDegree(degree);
    }

    /**
     * Determines whether estimation is ready to start with the given data and
     * required degree of polynomial to be estimated.
     *
     * @return true if estimator is ready, false otherwise.
     */
    public boolean isReady() {
        final int nParams = mPolynomialEstimator.getDegree() + 1;
        if (mEvaluations == null || mEvaluations.size() < nParams) {
            return false;
        }

        // also ensure that at least a direct or integral evaluation exists
        int count = 0;
        for (final PolynomialEvaluation eval : mEvaluations) {
            if (eval.getType() == PolynomialEvaluationType.DIRECT_EVALUATION ||
                    eval.getType() == PolynomialEvaluationType.INTEGRAL_EVALUATION ||
                    eval.getType() == PolynomialEvaluationType.INTEGRAL_INTERVAL) {
                count++;
            }
        }

        return count >= 1 && mEvaluations.size() >= nParams;
    }

    /**
     * Returns quality scores corresponding to each polynomial evaluation.
     * The larger the score value the better the quality of the evaluation.
     * This implementation always returns null.
     * Subclasses using quality scores must implement proper behaviour.
     *
     * @return quality scores corresponding to each evaluation.
     */
    public double[] getQualityScores() {
        // quality scores ignored
        return null;
    }

    /**
     * Sets quality scores corresponding to each polynomial evaluation.
     * The larger the score value the better the quality of the evaluation.
     * This implementation makes no action.
     * Subclasses using quality scores must implement proper behaviour.
     *
     * @param qualityScores quality scores corresponding to each evaluation.
     * @throws LockedException          if robust estimator is locked because an
     *                                  estimation is already in progress.
     * @throws IllegalArgumentException if provided quality scores length is
     *                                  smaller than minimum required number of evaluations.
     */
    public void setQualityScores(final double[] qualityScores) throws LockedException {
        // quality scores ignored
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
    public abstract Polynomial estimate() throws LockedException,
            NotReadyException, RobustEstimatorException;

    /**
     * Returns method being used for robust estimation.
     *
     * @return method being used for robust estimation.
     */
    public abstract RobustEstimatorMethod getMethod();

    /**
     * Creates a robust polynomial estimator using default method.
     *
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create() {
        return create(DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided degree and default
     * method.
     *
     * @param degree degree of polynomial to be estimated.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(final int degree) {
        return create(degree, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided evaluations and
     * default method.
     *
     * @param evaluations collection of polynomial evaluations.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final List<PolynomialEvaluation> evaluations) {
        return create(evaluations, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided listener and default
     * method.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final PolynomialRobustEstimatorListener listener) {
        return create(listener, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided degree, evaluations
     * and default method.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations) {
        return create(degree, evaluations, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided degree, listener
     * and default method.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final PolynomialRobustEstimatorListener listener) {
        return create(degree, listener, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided evaluations,
     * listener and default method.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener) {
        return create(evaluations, listener, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided degree, evaluations,
     * listener and default method.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener) {
        return create(degree, evaluations, listener, DEFAULT_ROBUST_METHOD);
    }

    /**
     * Creates a robust polynomial estimator using provided method.
     *
     * @param method method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator();
            case LMedS:
                return new LMedSPolynomialRobustEstimator();
            case MSAC:
                return new MSACPolynomialRobustEstimator();
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator();
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator();
        }
    }

    /**
     * Creates a robust polynomial estimator using provided degree and method.
     *
     * @param degree degree of polynomial to be estimated.
     * @param method method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(degree);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(degree);
            case MSAC:
                return new MSACPolynomialRobustEstimator(degree);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(degree);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(degree);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided evaluations and
     * method.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param method      method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(evaluations);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(evaluations);
            case MSAC:
                return new MSACPolynomialRobustEstimator(evaluations);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(evaluations);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(evaluations);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided listener and method.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @param method   method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final PolynomialRobustEstimatorListener listener,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(listener);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(listener);
            case MSAC:
                return new MSACPolynomialRobustEstimator(listener);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(listener);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(listener);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided degree, evaluations
     * and method.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param method      method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(degree, evaluations);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(degree, evaluations);
            case MSAC:
                return new MSACPolynomialRobustEstimator(degree, evaluations);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(degree,
                        evaluations);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(degree, evaluations);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided degree, listener
     * and method.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes.
     * @param method   method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final PolynomialRobustEstimatorListener listener,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(degree, listener);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(degree, listener);
            case MSAC:
                return new MSACPolynomialRobustEstimator(degree, listener);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(degree, listener);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(degree, listener);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided evaluations,
     * listener and method.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @param method      method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     */
    public static PolynomialRobustEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(evaluations,
                        listener);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(evaluations,
                        listener);
            case MSAC:
                return new MSACPolynomialRobustEstimator(evaluations, listener);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(evaluations,
                        listener);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(evaluations,
                        listener);
        }
    }

    /**
     * Creates a robust polynomial estimator using provided degree, evaluations,
     * listener and method.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events such as when estimation
     *                    starts, ends or its progress significantly changes.
     * @param method      method of a robust polynomial estimator.
     * @return an instance of a robust polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialRobustEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final PolynomialRobustEstimatorListener listener,
            final RobustEstimatorMethod method) {
        switch (method) {
            case RANSAC:
                return new RANSACPolynomialRobustEstimator(degree, evaluations,
                        listener);
            case LMedS:
                return new LMedSPolynomialRobustEstimator(degree, evaluations,
                        listener);
            case MSAC:
                return new MSACPolynomialRobustEstimator(degree, evaluations,
                        listener);
            case PROMedS:
                return new PROMedSPolynomialRobustEstimator(degree, evaluations,
                        listener);
            case PROSAC:
            default:
                return new PROSACPolynomialRobustEstimator(degree, evaluations,
                        listener);
        }
    }

    /**
     * Computes geometric or algebraic distance between provided polynomial
     * and evaluation.
     *
     * @param eval       polynomial evaluation.
     * @param polynomial polynomial.
     * @return distance.
     */
    protected double getDistance(final PolynomialEvaluation eval,
                                 final Polynomial polynomial) {
        if (mUseGeometricDistance) {
            return getGeometricOrAlgebraicDistance(eval, polynomial);
        } else {
            return getAlgebraicDistance(eval, polynomial);
        }
    }

    /**
     * Computes algebraic distance between provided polynomial and evaluation.
     *
     * @param eval       polynomial evaluation.
     * @param polynomial polynomial.
     * @return algebraic distance.
     */
    protected double getAlgebraicDistance(final PolynomialEvaluation eval,
                                          final Polynomial polynomial) {
        switch (eval.getType()) {
            case DIRECT_EVALUATION:
                return getAlgebraicDistance((DirectPolynomialEvaluation) eval,
                        polynomial);
            case DERIVATIVE_EVALUATION:
                return getAlgebraicDistance(
                        (DerivativePolynomialEvaluation) eval, polynomial);
            case INTEGRAL_EVALUATION:
                return getAlgebraicDistance((IntegralPolynomialEvaluation) eval,
                        polynomial);
            case INTEGRAL_INTERVAL:
                return getAlgebraicDistance(
                        (IntegralIntervalPolynomialEvaluation) eval, polynomial);
            default:
                return Double.MAX_VALUE;
        }
    }

    /**
     * Computes algebraic distance of between provided polynomial and direct
     * evaluation.
     *
     * @param eval       direct polynomial evaluation.
     * @param polynomial polynomial.
     * @return algebraic distance.
     */
    protected double getAlgebraicDistance(final DirectPolynomialEvaluation eval,
                                          final Polynomial polynomial) {
        final double x = eval.getX();
        final double y1 = eval.getEvaluation();
        final double y2 = polynomial.evaluate(x);
        return Math.abs(y2 - y1);
    }

    /**
     * Computes algebraic distance of a derivative between provided polynomial
     * and evaluation.
     *
     * @param eval       derivative polynomial evaluation.
     * @param polynomial polynomial.
     * @return algebraic distance.
     */
    protected double getAlgebraicDistance(final DerivativePolynomialEvaluation eval,
                                          final Polynomial polynomial) {
        final double x = eval.getX();
        final int order = eval.getDerivativeOrder();
        final double d1 = eval.getEvaluation();
        final double d2 = polynomial.evaluateNthDerivative(x, order);
        return Math.abs(d2 - d1);
    }

    /**
     * Computes algebraic distance of an integral between provided polynomial
     * and evaluation.
     *
     * @param eval       integration polynomial evaluation.
     * @param polynomial polynomial.
     * @return algebraic distance.
     */
    protected double getAlgebraicDistance(final IntegralPolynomialEvaluation eval,
                                          final Polynomial polynomial) {
        final double x = eval.getX();
        final int order = eval.getIntegralOrder();
        final double[] constants = eval.getConstants();
        final double i1 = eval.getEvaluation();
        final double i2 = polynomial.nthIntegrationAndReturnNew(order, constants).
                evaluate(x);
        return Math.abs(i2 - i1);
    }

    /**
     * Computes algebraic distance of an integration interval between provided
     * polynomial and evaluation.
     *
     * @param eval       integration interval polynomial evaluation.
     * @param polynomial polynomial.
     * @return algebraic distance.
     */
    protected double getAlgebraicDistance(
            final IntegralIntervalPolynomialEvaluation eval, final Polynomial polynomial) {
        final double startX = eval.getStartX();
        final double endX = eval.getEndX();
        final int order = eval.getIntegralOrder();
        final double[] constants = eval.getConstants();
        final double i1 = eval.getEvaluation();
        final double i2 = polynomial.nthOrderIntegrateInterval(startX, endX, order,
                constants);
        return Math.abs(i2 - i1);
    }

    /**
     * Commutes distance of evaluation respect to provided polynomial in
     * a geometric sense if evaluation is direct, otherwise returns algebraic
     * distance.
     *
     * @param eval       polynomial evaluation.
     * @param polynomial polynomial.
     * @return geometric distance for direct evaluation or algebraic distance
     * otherwise.
     */
    protected double getGeometricOrAlgebraicDistance(
            final PolynomialEvaluation eval, final Polynomial polynomial) {
        if (eval.getType() == PolynomialEvaluationType.DIRECT_EVALUATION) {
            return getGeometricDistance((DirectPolynomialEvaluation) eval,
                    polynomial);
        } else {
            return getAlgebraicDistance(eval, polynomial);
        }
    }

    /**
     * Computes distance of evaluation respect to provided polynomial in
     * a geometric sense by computing a tangent line to polynomial at point x
     * and comparing the distance of such line to provided evaluation point.
     *
     * @param eval       polynomial evaluation.
     * @param polynomial polynomial.
     * @return geometric distance.
     */
    protected double getGeometricDistance(final DirectPolynomialEvaluation eval,
                                          final Polynomial polynomial) {
        final double x = eval.getX();
        final double y1 = eval.getEvaluation();
        final double y2 = polynomial.evaluate(x);

        final double slope = polynomial.evaluateDerivative(x);
        final double a;
        final double b;
        final double c;
        if (Math.abs(slope) > 1.0) {
            a = 1.0;
            b = -1.0 / slope;
            c = -x + y2 / slope;
        } else {
            a = -slope;
            b = 1.0;
            c = slope * x - y2;
        }

        final double num = x * a + y1 * b + c;
        final double den = Math.sqrt(a * a + b * b);

        return Math.abs(num / den);
    }

    /**
     * Sets list of polynomial evaluations.
     * This method does not check whether estimator is locked.
     *
     * @param evaluations list of polynomial evaluations to estimate polynomial.
     * @throws IllegalArgumentException if provided list of polynomials is null
     *                                  or too small.
     */
    private void internalSetEvaluations(final List<PolynomialEvaluation> evaluations) {
        if (evaluations == null ||
                evaluations.size() < getMinNumberOfEvaluations()) {
            throw new IllegalArgumentException();
        }
        mEvaluations = evaluations;
    }
}
