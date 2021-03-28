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

import com.irurueta.algebra.Matrix;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;

import java.util.Arrays;
import java.util.List;

/**
 * This class defines the interface for an estimator of a polynomial of a given
 * degree using points where polynomials are evaluated.
 */
@SuppressWarnings("Duplicates")
public abstract class PolynomialEstimator {

    /**
     * Minimum allowed degree to be estimated
     */
    public static final int MIN_DEGREE = 1;

    /**
     * Default estimator type.
     */
    public static final PolynomialEstimatorType DEFAULT_ESTIMATOR_TYPE =
            PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR;

    /**
     * Degree of polynomial to be estimted.
     */
    protected int mDegree;

    /**
     * Collection of polynomial evaluations and their corresponding point of
     * evaluation used to determine a polynomial of required degree.
     */
    protected List<PolynomialEvaluation> mEvaluations;

    /**
     * True when estimator is estimating radial distorion.
     */
    protected boolean mLocked;

    /**
     * Listener to be notified of events such as when estimation starts, ends or
     * estimation progress changes.
     */
    protected PolynomialEstimatorListener mListener;

    /**
     * Constructor.
     */
    protected PolynomialEstimator() {
        mDegree = MIN_DEGREE;
    }

    /**
     * Constructor.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialEstimator(final int degree) {
        internalSetDegree(degree);
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     */
    protected PolynomialEstimator(final List<PolynomialEvaluation> evaluations) {
        this();
        mEvaluations = evaluations;
    }

    /**
     * Constructor.
     *
     * @param listener listener to be notified of events.
     */
    protected PolynomialEstimator(final PolynomialEstimatorListener listener) {
        this();
        mListener = listener;
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialEstimator(final int degree,
                                  final List<PolynomialEvaluation> evaluations) {
        this(degree);
        mEvaluations = evaluations;
    }

    /**
     * Constructor.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialEstimator(final int degree,
                                  final PolynomialEstimatorListener listener) {
        this(degree);
        mListener = listener;
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     */
    protected PolynomialEstimator(final List<PolynomialEvaluation> evaluations,
                                  final PolynomialEstimatorListener listener) {
        this(evaluations);
        mListener = listener;
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimted.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    protected PolynomialEstimator(final int degree,
                                  final List<PolynomialEvaluation> evaluations,
                                  final PolynomialEstimatorListener listener) {
        this(degree, evaluations);
        mListener = listener;
    }

    /**
     * Gets degree of polynomial to be estimated.
     *
     * @return degree of polynomial to be estimated.
     */
    public int getDegree() {
        return mDegree;
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

        internalSetDegree(degree);
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
     * @throws LockedException if this instance is locked.
     */
    public void setEvaluations(final List<PolynomialEvaluation> evaluations)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }

        mEvaluations = evaluations;
    }

    /**
     * Sets degree of polynomial to be estimated and collection of polynomial
     * evaluations and their corresponding point of evaluation used to determine
     * a polynomial of specified degree.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1.
     * @throws LockedException          if this instance is locked.
     */
    public void setDegreeAndEvaluations(
            final int degree, final List<PolynomialEvaluation> evaluations)
            throws LockedException {
        setDegree(degree);
        setEvaluations(evaluations);
    }

    /**
     * Determines whether estimation is ready to start with the given data
     * and required degree of polynomial to be estimated
     *
     * @return true if estimator is ready, false otherwise.
     */
    public boolean isReady() {

        final int nParams = mDegree + 1;
        if (mEvaluations == null || mEvaluations.size() < nParams) {
            return false;
        }

        //also ensure that at least a direct or integral evaluation exists
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
     * Gets minimum number of evaluations required to estimate a polynomial of
     * the specified degree.
     *
     * @param degree degree of polynomial to be estimated.
     * @return number of required evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static int getMinNumberOfEvaluations(final int degree) {
        if (degree < MIN_DEGREE) {
            throw new IllegalArgumentException();
        }

        return degree + 1;
    }

    /**
     * Gets minimum number of evaluations required to estimate a polynomial of
     * the specified degree.
     *
     * @return number of required evaluations.
     */
    public int getMinNumberOfEvaluations() {
        return getMinNumberOfEvaluations(mDegree);
    }

    /**
     * Indicates whether this instance is locked.
     *
     * @return true if this estimator is busy estimating a polynomial, false
     * otherwise.
     */
    public boolean isLocked() {
        return mLocked;
    }

    /**
     * Gets listener to be notified of events such as when estimation starts,
     * ends or estimation progress changes.
     *
     * @return listener to be notified of events.
     */
    public PolynomialEstimatorListener getListener() {
        return mListener;
    }

    /**
     * Sets listener to be notified of events such as when estimation starts,
     * ends or estimation progress changes.
     *
     * @param listener listener to be notified of events.
     * @throws LockedException if estimator is locked.
     */
    public void setListener(final PolynomialEstimatorListener listener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }

        mListener = listener;
    }

    /**
     * Estimates a polynomial based on provided evaluations.
     *
     * @return estimated polynomial.
     * @throws LockedException               if estimator is locked.
     * @throws NotReadyException             if estimator is not ready.
     * @throws PolynomialEstimationException if polynomial estimation fails.
     */
    public abstract Polynomial estimate() throws LockedException,
            NotReadyException, PolynomialEstimationException;

    /**
     * Returns type of polynomial estimator.
     *
     * @return type of polynomial estimator.
     */
    public abstract PolynomialEstimatorType getType();

    /**
     * Creates an instance of a polynomial estimator using default
     * type.
     *
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create() {
        return create(DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree and
     * default type.
     *
     * @param degree degree of polynomial to be estimated.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(final int degree) {
        return create(degree, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided evaluations,
     * and default type and degree.
     *
     * @param evaluations collection of polynomial evaluations.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final List<PolynomialEvaluation> evaluations) {
        return create(evaluations, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided listener
     * and default type and degree.
     *
     * @param listener listener to be notified of events.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final PolynomialEstimatorListener listener) {
        return create(listener, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * polynomial evaluations and default type.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations) {
        return create(degree, evaluations, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * listener and default type.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final PolynomialEstimatorListener listener) {
        return create(degree, listener, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided evaluations,
     * listener and default type.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorListener listener) {
        return create(evaluations, listener, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * evaluations, listener and default type.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorListener listener) {
        return create(degree, evaluations, listener, DEFAULT_ESTIMATOR_TYPE);
    }

    /**
     * Creates an instance of a polynomial estimator using provided type and
     * default degree.
     *
     * @param type type of polynomial estimator
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                return new WeightedPolynomialEstimator();
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator();
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree and
     * type.
     *
     * @param degree degree of polynomial to be estimated.
     * @param type   type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                return new WeightedPolynomialEstimator(degree);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(degree);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided evaluations
     * and type.
     *
     * @param evaluations collection of polynomial esvaluations.
     * @param type        type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                final double[] weights = new double[evaluations.size()];
                Arrays.fill(weights, 1.0);
                return new WeightedPolynomialEstimator(evaluations, weights);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(evaluations);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided listener
     * and type.
     *
     * @param listener listener to be notified of events.
     * @param type     type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final PolynomialEstimatorListener listener,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                return new WeightedPolynomialEstimator(listener);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(listener);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * evaluations and type.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param type        type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                final double[] weights = new double[evaluations.size()];
                Arrays.fill(weights, 1.0);
                return new WeightedPolynomialEstimator(degree, evaluations,
                        weights);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(degree, evaluations);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * listener and type.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @param type     type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree, final PolynomialEstimatorListener listener,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                return new WeightedPolynomialEstimator(degree, listener);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(degree, listener);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided evaluations,
     * listener and type.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @param type        type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     */
    public static PolynomialEstimator create(
            final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorListener listener,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                final double[] weights = new double[evaluations.size()];
                Arrays.fill(weights, 1.0);
                return new WeightedPolynomialEstimator(evaluations, weights,
                        listener);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(evaluations, listener);
        }
    }

    /**
     * Creates an instance of a polynomial estimator using provided degree,
     * evaluations, listener and type.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener    listener to be notified of events.
     * @param type        type of polynomial estimator.
     * @return an instance of a polynomial estimator.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public static PolynomialEstimator create(
            final int degree,
            final List<PolynomialEvaluation> evaluations,
            final PolynomialEstimatorListener listener,
            final PolynomialEstimatorType type) {
        switch (type) {
            case WEIGHTED_POLYNOMIAL_ESTIMATOR:
                final double[] weights = new double[evaluations.size()];
                Arrays.fill(weights, 1.0);
                return new WeightedPolynomialEstimator(degree, evaluations,
                        weights, listener);
            case LMSE_POLYNOMIAL_ESTIMATOR:
            default:
                return new LMSEPolynomialEstimator(degree, evaluations,
                        listener);
        }
    }

    /**
     * Fils row of system of equations for a direct polynomial evaluation.
     *
     * @param evaluation a direct polynomial evaluation.
     * @param a          system matrix.
     * @param b          values matrix.
     * @param row        row to be filled.
     */
    protected void fillDirectEvaluation(
            final DirectPolynomialEvaluation evaluation,
            final Matrix a, final Matrix b, final int row) {

        double powX = 1.0;
        final double x = evaluation.getX();
        for (int i = 0; i < a.getColumns(); i++) {
            a.setElementAt(row, i, powX);
            powX *= x;
        }

        b.setElementAtIndex(row, evaluation.getEvaluation());
    }

    /**
     * Fills row of system of equations for a derivative polynomial evaluation.
     *
     * @param evaluation a derivative polynomial evaluation.
     * @param a          system matrix.
     * @param b          values matrix.
     * @param row        row to be filled.
     */
    protected void fillDerivativeEvaluation(
            final DerivativePolynomialEvaluation evaluation, final Matrix a, final Matrix b,
            final int row) {

        final int order = evaluation.getDerivativeOrder();

        for (int i = 0; i < order; i++) {
            a.setElementAt(row, i, 0.0);
        }

        double powX = 1.0;
        final double x = evaluation.getX();
        for (int i = order; i < a.getColumns(); i++) {
            int param = i;
            for (int j = 1; j < order; j++) {
                param *= i - j;
            }
            a.setElementAt(row, i, param * powX);
            powX *= x;
        }

        b.setElementAtIndex(row, evaluation.getEvaluation());
    }

    /**
     * Fills row of system of equations for an integral polynomial evaluation.
     *
     * @param evaluation an integral polynomial evaluation.
     * @param a          system matrix.
     * @param b          values matrix.
     * @param row        row to be filled.
     * @throws PolynomialEstimationException if constant terms does not have
     *                                       proper size (it must be null or have order length).
     */
    protected void fillIntegralEvaluation(
            final IntegralPolynomialEvaluation evaluation,
            final Matrix a, final Matrix b, final int row)
            throws PolynomialEstimationException {

        final int order = evaluation.getIntegralOrder();
        final double[] constants = evaluation.getConstants();
        if (constants != null && constants.length != order) {
            throw new PolynomialEstimationException();
        }

        double accum = 0.0;
        double powX = 1.0;
        final double x = evaluation.getX();
        for (int i = 0; i < order; i++) {
            if (constants != null) {
                int param = 1;
                for (int k = 1; k <= i; k++) {
                    param *= k;
                }
                accum += constants[i] / (double) param * powX;
            }
            powX *= x;
        }

        for (int i = 0, j = order; i < a.getColumns(); i++, j++) {
            int param = j;
            for (int k = 1; k < order; k++) {
                param *= j - k;
            }
            a.setElementAt(row, i, powX / (double) param);
            powX *= x;
        }

        b.setElementAtIndex(row, evaluation.getEvaluation() - accum);
    }

    /**
     * Fills row of system of equations for a polynomial evaluation of an
     * interval integration.
     *
     * @param evaluation an interval integration of a polynomial evaluation.
     * @param a          system matrix.
     * @param b          values matrix.
     * @param row        row to be filled.
     * @throws PolynomialEstimationException if constant terms does not have
     *                                       proper size (it must be null or have order length).
     */
    protected void fillIntegralIntervalEvaluation(
            final IntegralIntervalPolynomialEvaluation evaluation, final Matrix a, final Matrix b,
            final int row) throws PolynomialEstimationException {

        final int order = evaluation.getIntegralOrder();
        final double[] constants = evaluation.getConstants();
        if (constants != null && constants.length != order) {
            throw new PolynomialEstimationException();
        }

        double accum = 0.0;
        double powStartX = 1.0;
        double powEndX = 1.0;
        final double startX = evaluation.getStartX();
        final double endX = evaluation.getEndX();
        for (int i = 0; i < order; i++) {
            if (constants != null) {
                int param = 1;
                for (int k = 1; k <= i; k++) {
                    param *= k;
                }
                accum += constants[i] / (double) param * (powEndX - powStartX);
            }
            powStartX *= startX;
            powEndX *= endX;
        }

        for (int i = 0, j = order; i < a.getColumns(); i++, j++) {
            int param = j;
            for (int k = 1; k < order; k++) {
                param *= j - k;
            }
            a.setElementAt(row, i, (powEndX - powStartX) / (double) param);
            powStartX *= startX;
            powEndX *= endX;
        }

        b.setElementAtIndex(row, evaluation.getEvaluation() - accum);
    }

    /**
     * Normalizes rows of system matrix and values matrix to increase accuracy
     * of linear system of equations to be solved.
     *
     * @param a   system matrix.
     * @param b   values matrix.
     * @param row row to normalize.
     */
    protected void normalize(final Matrix a, final Matrix b, final int row) {
        double sqrNorm = 0.0;
        for (int i = 0; i < a.getColumns(); i++) {
            sqrNorm += Math.pow(a.getElementAt(row, i), 2.0);
        }
        sqrNorm += Math.pow(b.getElementAtIndex(row), 2.0);

        final double norm = Math.sqrt(sqrNorm);

        for (int i = 0; i < a.getColumns(); i++) {
            a.setElementAt(row, i, a.getElementAt(row, i) / norm);
        }
        b.setElementAtIndex(row, b.getElementAtIndex(row) / norm);
    }

    /**
     * Internal method to set degree of polynomial to be estimated.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    private void internalSetDegree(final int degree) {
        if (degree < MIN_DEGREE) {
            throw new IllegalArgumentException("degree must be at least 1");
        }
        mDegree = degree;
    }
}
