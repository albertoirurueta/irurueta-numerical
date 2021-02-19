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
package com.irurueta.numerical;

import com.irurueta.numerical.optimization.BrentSingleOptimizer;

/**
 * Class to estimate the most likely value from a series of samples assumed to
 * be normally distributed.
 * This implementation will first use an internal HistogramMaximumLikelihood
 * to estimate the most likely value using a histogram, and then it will
 * refine the solution by using a BrentSingleOptimizer in order to find the
 * maximum of the probability distribution function assuming that such function
 * is computed by aggregating small Gaussians (of size gaussianSigma) centered
 * at the location of each sample.
 */
@SuppressWarnings("WeakerAccess")
public class AccurateMaximumLikelihoodEstimator
        extends MaximumLikelihoodEstimator {

    /**
     * Boolean indicating if an initial solution should be obtained first by
     * using the Histogram method. It is suggested to always enable this option.
     */
    public static final boolean DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION = true;

    /**
     * Value to be considered as the machine precision.
     */
    public static final double EPS = 1e-9;

    /**
     * Boolean that indicates that an initial coarse solution will be computed
     * first by using an internal HistogramMaximumLikelihoodEstimator in order
     * to initialize the internal BrentSingleOptimizer to obtain a more accurate
     * solution.
     */
    private boolean useHistogramInitialSolution;

    /**
     * Internal maximum likelihood estimator based on the Histogram method.
     */
    private HistogramMaximumLikelihoodEstimator internalEstimator;

    /**
     * Internal optimizer to find the true maximum of the probability
     * distribution function. Because a BrentSingleOptimizer is only guaranteed
     * to obtain local minima/maxima, it is preferred to start the optimizer
     * near the true solution to be found, for that reason it is suggested to
     * always use the Histogram initial solution as coarse approximation to
     * start the optimizer and get a more accurate solution.
     */
    private BrentSingleOptimizer optimizer;

    /**
     * Constructor.
     *
     * @param gaussianSigma               Gaussian sigma to be used on each sample.
     * @param useHistogramInitialSolution Boolean indicating whether an internal
     *                                    HistogramMaximumLikelihoodEstimator will be used to obtain a coarse
     *                                    initial solution to initialize the BrentSingleOptimizer. It is suggested
     *                                    to set this value always to true.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero.
     */
    public AccurateMaximumLikelihoodEstimator(final double gaussianSigma,
                                              final boolean useHistogramInitialSolution) {
        super(gaussianSigma);
        this.useHistogramInitialSolution = useHistogramInitialSolution;
        internalEstimator = null;
        optimizer = null;
    }

    /**
     * Empty constructor.
     */
    public AccurateMaximumLikelihoodEstimator() {
        super();
        this.useHistogramInitialSolution =
                DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION;
        internalEstimator = null;
        optimizer = null;
    }

    /**
     * Constructor
     *
     * @param inputData                   Array containing input data where most likely value must
     *                                    be estimated from.
     * @param gaussianSigma               Gaussian sigma to be used on each sample.
     * @param useHistogramInitialSolution Boolean indicating whether an internal
     *                                    HistogramMaximumLikelihoodEstimator will be used to obtain a coarse
     *                                    initial solution to initialize the BrentSingleOptimizer. It is suggested
     *                                    to set this value always to true.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero.
     */
    public AccurateMaximumLikelihoodEstimator(
            final double[] inputData, final double gaussianSigma,
            final boolean useHistogramInitialSolution) {
        super(inputData, gaussianSigma);
        this.useHistogramInitialSolution = useHistogramInitialSolution;
        internalEstimator = null;
        optimizer = null;
    }

    /**
     * Constructor.
     *
     * @param minValue                    Minimum value assumed to be contained within input data
     *                                    array.
     * @param maxValue                    Maximum value assumed to be contained within input data
     *                                    array.
     * @param inputData                   Array containing input data where most likely value must
     *                                    be estimated from.
     * @param gaussianSigma               Gaussian sigma to be used on each sample.
     * @param useHistogramInitialSolution Boolean indicating whether an internal
     *                                    HistogramMaximumLikelihoodEstimator will be used to obtain a coarse
     *                                    initial solution to initialize the BrentSingleOptimizer. It is suggested
     *                                    to set this value always to true.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero, or if minValue &lt; maxValue.
     */
    public AccurateMaximumLikelihoodEstimator(
            final double minValue, final double maxValue,
            final double[] inputData, final double gaussianSigma,
            final boolean useHistogramInitialSolution) {
        super(minValue, maxValue, inputData, gaussianSigma);
        this.useHistogramInitialSolution = useHistogramInitialSolution;
        internalEstimator = null;
        optimizer = null;
    }

    /**
     * Returns method to be used for maximum likelihood estimation, which for
     * this class is MaximumLikelihoodEstimatorMethod.
     * ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR.
     *
     * @return Method for maximum likelihood estimation.
     */
    @Override
    public MaximumLikelihoodEstimatorMethod getMethod() {
        return MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR;
    }

    /**
     * Returns boolean that indicates that an initial coarse solution will be
     * computed first by using an internal HistogramMaximumLikelihoodEstimator
     * in order to initialize the internal BrentSingleOptimizer to obtain a more
     * accurate solution.
     *
     * @return True if an initial coarse solution if found by the Histogram
     * method, false otherwise.
     */
    public boolean isHistogramInitialSolutionUsed() {
        return useHistogramInitialSolution;
    }

    /**
     * Sets boolean that indicates that an initial coarse solution will be
     * computed first by using an internal HistogramMaximumLikelihoodEstimator
     * in order to initialize the internal BrentSingleOptimizer to obtain a more
     * accurate solution.
     *
     * @param used True if an initial coarse solution will be found by the
     *             Histogram method, false otherwise.
     * @throws LockedException Exception raised if this instance is locked.
     *                         This method can only be executed when computations finish and this
     *                         instance becomes unlocked.
     */
    public void setHistogramInitialSolutionUsed(final boolean used)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        useHistogramInitialSolution = used;
    }

    /**
     * Starts the estimation of the most likely value contained within provided
     * input data array.
     *
     * @return The most likely value.
     * @throws LockedException   Exception raised if this instance is locked.
     *                           This method can only be executed when computations finish and this
     *                           instance becomes unlocked.
     * @throws NotReadyException Exception raised if this instance is not yet
     *                           ready.
     * @see #isReady()
     */
    @Override
    public double estimate() throws LockedException, NotReadyException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        locked = true;

        final double minEvalPoint;
        final double middleEvalPoint;
        final double maxEvalPoint;

        if (useHistogramInitialSolution) {
            if (internalEstimator == null) {
                internalEstimator = new HistogramMaximumLikelihoodEstimator();
            }
            internalEstimator.setInputData(inputData);
            internalEstimator.setGaussianSigma(gaussianSigma);
            internalEstimator.computeMinMaxValues();

            middleEvalPoint = internalEstimator.estimate();

            double localMinValue = 0.0;
            double localMaxValue = 0.0;

            try {
                localMinValue = internalEstimator.getMinValue();
                localMaxValue = internalEstimator.getMaxValue();
            } catch (final NotAvailableException ignore) {
                // never happens
            }

            final int numberOfBins = internalEstimator.getNumberOfBins();

            final double delta = (localMaxValue - localMinValue) /
                    ((double) (numberOfBins - 1));

            // pick two values around initial coarse solution
            minEvalPoint = middleEvalPoint - delta;
            maxEvalPoint = middleEvalPoint + delta;

            if (!areMinMaxAvailable) {
                this.minValue = localMinValue;
                this.maxValue = localMaxValue;
                areMinMaxAvailable = true;
            }
        } else {
            if (!areMinMaxAvailable) {
                computeMinMaxValues();
            }

            // use min/max values as a bracket to obtain optimal solution
            minEvalPoint = minValue;
            maxEvalPoint = maxValue;
            middleEvalPoint = (minValue + maxValue) * 0.5;
        }

        if ((maxValue - minValue) < EPS) {
            // min-max limits are almost equal, so we return it as the solution
            locked = false;
            return middleEvalPoint;
        }

        double solution;
        try {
            // Use an optimizer to find maximum value on histogram (PDF)
            if (optimizer == null) {
                optimizer = new BrentSingleOptimizer(new EvaluatorListener(),
                        BrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_TOLERANCE);
            }

            optimizer.setBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);
            optimizer.minimize();
            solution = optimizer.getResult();
        } catch (final Exception ignore) {
            // if optimization fails, pick coarse solution if available
            if (useHistogramInitialSolution) {
                solution = middleEvalPoint;
            } else {
                // if coarse solution is not available, then compute it
                internalEstimator.setInputData(inputData);
                internalEstimator.setGaussianSigma(gaussianSigma);
                internalEstimator.computeMinMaxValues();
                solution = internalEstimator.estimate();
            }
        }

        locked = false;
        return solution;
    }

    /**
     * Internal class used by the BrentSingleOptimizer in order to evaluate
     * the aggregation of Gaussians for all the samples in input data array with
     * a high degree of precision.
     */
    private class EvaluatorListener
            implements SingleDimensionFunctionEvaluatorListener {


        /**
         * Evaluates the aggregation of Gaussians for all the samples in input
         * data array, by assuming that each sample has an associated small
         * Gaussian centered at the sample value and with a small sigma value.
         * The aggregation of Gaussians will generate the averaged PDF function
         * of all the input values.
         *
         * @param point Point where the aggregation of Gaussians will be
         *              evaluated
         * @return The value of the aggregation of samples at provided point.
         * @throws EvaluationException Raised if anything failed during the execution.
         *                             This exception should never be raised.
         */
        @Override
        public double evaluate(final double point) throws EvaluationException {
            double out = 0.0;
            double x;
            for (final double data : inputData) {
                x = point - data;
                out += Math.exp(-x * x / (2.0 * gaussianSigma * gaussianSigma)) /
                        (Math.sqrt(2.0 * Math.PI) * gaussianSigma);
            }

            // negate value because optimizers always attempt to
            // minimize function
            return -out;
        }

    }
}
