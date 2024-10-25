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

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.Utils;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.numerical.robust.WeightSelection;
import com.irurueta.sorting.SortingException;

import java.util.List;

/**
 * This class implements a polynomial estimator using weighted evaluations.
 * Weights can be used so that evaluations assumed to have a better quality
 * (i.e. more precisely estimated) are considered to be more relevant.
 * It is discouraged to use a large number of evaluations, even if they are
 * correctly weighted, since as the number of evaluations increase so do the
 * rounding errors.
 */
@SuppressWarnings("DuplicatedCode")
public class WeightedPolynomialEstimator extends PolynomialEstimator {

    /**
     * Default number of evaluations to be weighted and taken into account.
     */
    public static final int DEFAULT_MAX_EVALUATIONS = 50;

    /**
     * Indicates if weights are sorted by default so that largest weighted
     * evaluations are used first.
     */
    public static final boolean DEFAULT_SORT_WEIGHTS = true;

    /**
     * Maximum number of evaluations to be weighted and taken into account.
     */
    private int maxEvaluations = DEFAULT_MAX_EVALUATIONS;

    /**
     * Indicates if weights are sorted by default so that largest weighted
     * evaluations are used first.
     */
    private boolean sortWeights = DEFAULT_SORT_WEIGHTS;

    /**
     * Array containing weights for all evaluations.
     */
    private double[] weights;

    /**
     * Constructor.
     */
    public WeightedPolynomialEstimator() {
        super();
    }

    /**
     * Constructor.
     *
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(final int degree) {
        super(degree);
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each evaluation. The
     *                    larger the value of a weight, the most significant the correspondence
     *                    will be.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size.
     */
    public WeightedPolynomialEstimator(
            final List<PolynomialEvaluation> evaluations, final double[] weights) {
        super();
        internalSetEvaluationsAndWeights(evaluations, weights);
    }

    /**
     * Constructor.
     *
     * @param listener listener to be notified of events.
     */
    public WeightedPolynomialEstimator(final PolynomialEstimatorListener listener) {
        super(listener);
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each evaluation. The
     *                    larger the value of a weight, the most significant the correspondence
     *                    will be.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size, or if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(
            final int degree, final List<PolynomialEvaluation> evaluations, final double[] weights) {
        super(degree);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }

    /**
     * Constructor.
     *
     * @param degree   degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(final int degree, final PolynomialEstimatorListener listener) {
        super(degree, listener);
    }

    /**
     * Constructor.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each evaluation. The
     *                    larger the value of a weight, the most significant the correspondence
     *                    will be.
     * @param listener    listener to be notified of events.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size.
     */
    public WeightedPolynomialEstimator(
            final List<PolynomialEvaluation> evaluations, final double[] weights,
            final PolynomialEstimatorListener listener) {
        super(listener);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }

    /**
     * Constructor.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each evaluation. The
     *                    larger the value of a weight, the most significant the correspondence
     *                    will be.
     * @param listener    listener to be notified of events.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size, or if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(
            final int degree, final List<PolynomialEvaluation> evaluations, final double[] weights,
            final PolynomialEstimatorListener listener) {
        super(degree, listener);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }

    /**
     * Sets collection of polynomial evaluations and their corresponding point
     * of evaluation used to determine a polynomial of required degree.
     * This method override always throws an IllegalArgumentException because it
     * is expected to provide both evaluations and their weights.
     *
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException always thrown.
     */
    @Override
    public void setEvaluations(final List<PolynomialEvaluation> evaluations) {
        throw new IllegalArgumentException("evaluations and weights must be provided at once");
    }

    /**
     * Sets collection of polynomial evaluations along with their corresponding
     * weights.
     *
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each polynomial
     *                    evaluation. The larger the value of a weight, the most significant the
     *                    evaluation will be.
     * @throws LockedException          if estimator is locked.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size.
     */
    public void setEvaluationsAndWeights(
            final List<PolynomialEvaluation> evaluations, final double[] weights) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetEvaluationsAndWeights(evaluations, weights);
    }

    /**
     * Sets degree of polynomial to be estimated and collection of polynomial
     * evaluations and their corresponding point of evaluation used to determine
     * a polynomial of specified degree.
     *
     * @param degree      degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights     array containing a weight amount for each polynomial
     *                    evaluation. The larger the value of a weight, the most significant the
     *                    evaluation will be.
     * @throws IllegalArgumentException if provided degree is less than 1 or
     *                                  if evaluations or weights are null or don't have the same size.
     * @throws LockedException          if this instance is locked.
     */
    public void setDegreeEvaluationsAndWeights(
            final int degree, final List<PolynomialEvaluation> evaluations, final double[] weights)
            throws LockedException {
        setDegree(degree);
        setEvaluationsAndWeights(evaluations, weights);
    }


    /**
     * Returns array containing a weight amount for each polynomial evaluation.
     * The larger the value of a weight, the most significant the correspondence
     * will be.
     *
     * @return array containing weights for each correspondence.
     */
    public double[] getWeights() {
        return weights;
    }

    /**
     * Returns boolean indicating whether weights have been provided and are
     * available for retrieval.
     *
     * @return true if weights are available, false otherwise.
     */
    public boolean areWeightsAvailable() {
        return weights != null;
    }

    /**
     * Returns maximum number of evaluations to be weighted and taken into
     * account.
     *
     * @return maximum number of evaluations to be weighted.
     */
    public int getMaxEvaluations() {
        return maxEvaluations;
    }

    /**
     * Sets maximum number of evaluations to be weighted and taken into account.
     * This method must be called after setting degree, because the minimum
     * number of required evaluations will be checked based on degree of
     * polynomial to be estimated.
     *
     * @param maxEvaluations maximum number of evaluations to be weighted.
     * @throws IllegalArgumentException if provided value is less than the
     *                                  minimum number of required evaluations.
     * @throws LockedException          if this instance is locked.
     */
    public void setMaxEvaluations(final int maxEvaluations) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxEvaluations < getMinNumberOfEvaluations()) {
            throw new IllegalArgumentException();
        }
        this.maxEvaluations = maxEvaluations;
    }

    /**
     * Indicates if weights are sorted by so that largest weighted evaluations
     * are used first.
     *
     * @return true if weights are sorted, false otherwise.
     */
    public boolean isSortWeightsEnabled() {
        return sortWeights;
    }

    /**
     * Specifies whether weights are sorted by so that largest weighted
     * evaluations are used first.
     *
     * @param sortWeights true if weights are sorted, false otherwise.
     * @throws LockedException if this instance is locked.
     */
    public void setSortWeightsEnabled(final boolean sortWeights) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }

        this.sortWeights = sortWeights;
    }

    /**
     * Indicates if this estimator is ready to start the estimation.
     * Estimator will be ready once enough evaluations and weights are provided.
     *
     * @return true if estimator is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return super.isReady() && areWeightsAvailable() && evaluations.size() == weights.length;
    }

    /**
     * Estimates a polynomial based on provided evaluations.
     *
     * @return estimated polynomial.
     * @throws LockedException               if estimator is locked.
     * @throws NotReadyException             if estimator is not ready.
     * @throws PolynomialEstimationException if polynomial estimation fails.
     */
    @Override
    public Polynomial estimate() throws LockedException, NotReadyException, PolynomialEstimationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        try {
            locked = true;
            if (listener != null) {
                listener.onEstimateStart(this);
            }

            final var selection = WeightSelection.selectWeights(weights, sortWeights, maxEvaluations);
            final var selected = selection.getSelected();
            final var nEvaluations = selection.getNumSelected();


            final var a = new Matrix(nEvaluations, degree + 1);
            final var b = new Matrix(nEvaluations, 1);

            var index = 0;
            var counter = 0;
            double weight;
            for (final var evaluation : evaluations) {
                if (selected[index]) {
                    weight = weights[index];

                    switch (evaluation.getType()) {
                        case DIRECT_EVALUATION:
                            fillDirectEvaluation((DirectPolynomialEvaluation) evaluation, a, b, counter);
                            break;
                        case DERIVATIVE_EVALUATION:
                            fillDerivativeEvaluation((DerivativePolynomialEvaluation) evaluation, a, b, counter);
                            break;
                        case INTEGRAL_EVALUATION:
                            fillIntegralEvaluation((IntegralPolynomialEvaluation) evaluation, a, b, counter);
                            break;
                        case INTEGRAL_INTERVAL:
                            fillIntegralIntervalEvaluation((IntegralIntervalPolynomialEvaluation) evaluation, a, b,
                                    counter);
                            break;
                        default:
                            continue;
                    }

                    normalize(a, b, counter, weight);
                    counter++;
                }

                index++;
            }

            final var params = Utils.solve(a, b);

            final var result = new Polynomial(params.toArray());

            if (listener != null) {
                listener.onEstimateEnd(this);
            }

            return result;
        } catch (final AlgebraException | SortingException e) {
            throw new PolynomialEstimationException(e);
        } finally {
            locked = false;
        }
    }

    /**
     * Returns type of polynomial estimator.
     *
     * @return type of polynomial estimator.
     */
    @Override
    public PolynomialEstimatorType getType() {
        return PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR;
    }

    /**
     * Normalizes rows of system matrix and values matrix to increase accuracy
     * of linear system of equations to be solved.
     *
     * @param a      system matrix.
     * @param b      values matrix.
     * @param row    row to normalize.
     * @param weight weight.
     */
    private static void normalize(final Matrix a, final Matrix b, final int row, final double weight) {
        var sqrNorm = 0.0;
        for (var i = 0; i < a.getColumns(); i++) {
            sqrNorm += Math.pow(a.getElementAt(row, i), 2.0);
        }
        sqrNorm += Math.pow(b.getElementAtIndex(row), 2.0);

        final var norm = Math.sqrt(sqrNorm);
        final var factor = weight / norm;

        for (var i = 0; i < a.getColumns(); i++) {
            a.setElementAt(row, i, a.getElementAt(row, i) * factor);
        }
        b.setElementAtIndex(row, b.getElementAtIndex(row) * factor);
    }

    /**
     * Internal method to set evaluations and weights.
     *
     * @param evaluations evaluations.
     * @param weights     weights.
     * @throws IllegalArgumentException if evaluations or weights are null or
     *                                  don't have the same size.
     */
    private void internalSetEvaluationsAndWeights(
            final List<PolynomialEvaluation> evaluations, final double[] weights) {
        if (weights == null || evaluations == null || weights.length != evaluations.size()) {
            throw new IllegalArgumentException();
        }
        this.evaluations = evaluations;
        this.weights = weights;
    }
}
