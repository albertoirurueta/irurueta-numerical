/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
import com.irurueta.sorting.Sorter;

import java.util.ArrayList;
import java.util.BitSet;

/**
 * This class implements MSAC (Median SAmple Consensus) algorithm to robustly
 * estimate a data model.
 * MSAC is a mixture between LMedS and RANSAC, where a fixed threshold is
 * used such as in RANSAC to determine the number of remaining iterations, and
 * the least median of residuals is used to pick the best solution, rather than
 * the one producing a higher number of inliers based on the fixed threshold,
 * such as in RANSAC.
 * This algorithm requires a threshold known beforehand such as RANSAC, but
 * might get better accuracy if the inlier samples are very accurate, since
 * the solution with the smallest median of residuals will be picked. In typical
 * situations however, this algorithm will produce similar results to RANSAC in
 * both terms of accuracy and computational cost, since typically inlier samples
 * tend to have certain error.
 *
 * @param <T> type of object to be estimated.
 */
public class MSACRobustEstimator<T> extends RobustEstimator<T> {

    /**
     * Constant defining default confidence of estimated result, which is 99%.
     * This means that with a probability of 99% estimation will be accurate
     * because chosen sub-samples will be inliers.
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
     * Minimum allowed threshold to determine inliers.
     */
    public static final double MIN_THRESHOLD = 0.0;

    /**
     * Amount of confidence expressed as a value between 0 and 1.0 (which is
     * equivalent to 100%). The amount of confidence indicates the probability
     * that the estimated result is correct. Usually this value will be close
     * to 1.0, but not exactly 1.0.
     */
    private double confidence;

    /**
     * Maximum allowed number of iterations. When the maximum number of
     * iterations is exceeded, result will not be available, however an
     * approximate result will be available for retrieval.
     */
    private int maxIterations;

    /**
     * Instance in charge of picking random subsets of samples.
     */
    private SubsetSelector subsetSelector;

    /**
     * Number of iterations to be done to obtain required confidence.
     */
    private int iters;

    /**
     * Best solution that has been found so far during an estimation.
     */
    private T bestResult;

    /**
     * Data related to inliers found for best result.
     */
    private MSACInliersData bestResultInliersData;

    /**
     * Data related to solution producing the largest number of inliers.
     */
    private MSACInliersData bestNumberInliersData;

    /**
     * Constructor.
     */
    public MSACRobustEstimator() {
        super();
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        iters = maxIterations;
        bestResult = null;
        bestResultInliersData = bestNumberInliersData = null;
    }

    /**
     * Constructor.
     *
     * @param listener listener to handle events raised by this estimator.
     */
    public MSACRobustEstimator(final MSACRobustEstimatorListener<T> listener) {
        super(listener);
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        iters = maxIterations;
        bestResult = null;
        bestResultInliersData = bestNumberInliersData = null;
    }

    /**
     * Returns amount of confidence expressed as a value between 0 and 1.0
     * (which is equivalent to 100%). The amount of confidence indicates the
     * probability that the estimated result is correct. Usually this value will
     * be close to 1.0, but not exactly 1.0.
     *
     * @return amount of confidence as a value between 0.0 and 1.0.
     */
    public double getConfidence() {
        return confidence;
    }

    /**
     * Sets amount of confidence expressed as a value between 0 and 1.0 (which
     * is equivalent to 100%). The amount of confidence indicates the
     * probability that the estimated result is correct. Usually this value will
     * be close to 1.0, but not exactly 1.0.
     *
     * @param confidence confidence to be set as a value between 0.0 and 1.0.
     * @throws IllegalArgumentException if provided value is not between 0.0 and
     *                                  1.0.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setConfidence(final double confidence) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (confidence < MIN_CONFIDENCE || confidence > MAX_CONFIDENCE) {
            throw new IllegalArgumentException();
        }
        this.confidence = confidence;
    }

    /**
     * Maximum allowed number of iterations. When the maximum number of
     * iterations is exceeded, result will not be available, however an
     * approximate result will be available for retrieval.
     *
     * @return maximum allowed number of iterations.
     */
    public int getMaxIterations() {
        return maxIterations;
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
    public void setMaxIterations(final int maxIterations) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxIterations < MIN_ITERATIONS) {
            throw new IllegalArgumentException();
        }
        this.maxIterations = maxIterations;
    }

    /**
     * Returns number of iterations to be done to obtain required confidence.
     *
     * @return number of iterations to be done to obtain required confidence.
     */
    public int getNIters() {
        return iters;
    }

    /**
     * Returns best solution that has been found so far during an estimation.
     *
     * @return best solution that has been found so far during an estimation.
     */
    public T getBestResult() {
        return bestResult;
    }

    /**
     * Returns data related to the best inliers found for best result.
     *
     * @return data related to inliers found for best result.
     */
    public MSACInliersData getBestResultInliersData() {
        return bestResultInliersData;
    }

    /**
     * Returns data related to solution producing the largest number of inliers.
     *
     * @return data related to solution producing the largest number of inliers.
     */
    public MSACInliersData getBestNumberInliersData() {
        return bestNumberInliersData;
    }

    /**
     * Indicates if estimator is ready to start the estimation process.
     *
     * @return true if ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        if (!super.isReady()) {
            return false;
        }
        return (listener instanceof MSACRobustEstimatorListener);
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
    @Override
    @SuppressWarnings("Duplicates")
    public T estimate() throws LockedException, NotReadyException, RobustEstimatorException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        try {
            final var listener = (MSACRobustEstimatorListener<T>) this.listener;

            locked = true;

            listener.onEstimateStart(this);

            final var totalSamples = listener.getTotalSamples();
            final var subsetSize = listener.getSubsetSize();
            final var threshold = listener.getThreshold();
            // only positive thresholds are allowed
            if (threshold < MIN_THRESHOLD) {
                throw new RobustEstimatorException();
            }

            var bestNumInliers = 0;
            var bestMedianResidual = Double.MAX_VALUE;
            iters = Integer.MAX_VALUE;
            int newNIters;
            var currentIter = 0;
            // reusable list that will contain preliminary solutions on each
            // iteration
            final var iterResults = new ArrayList<T>();
            bestResult = null; // best result found so far
            int currentInliers;
            // progress and previous progress to determine when progress
            // notification must occur
            var previousProgress = 0.0f;
            float progress;
            // indices of subset picked in one iteration
            final var subsetIndices = new int[subsetSize];
            final var residualsTemp = new double[totalSamples];

            if (subsetSelector == null) {
                // create new subset selector
                subsetSelector = SubsetSelector.create(totalSamples);
            } else {
                // set number of samples to current subset selector
                subsetSelector.setNumSamples(totalSamples);
            }

            // data related to inliers
            var inliersData = new MSACInliersData(totalSamples);
            // sorter to compute medians
            final var sorter = Sorter.<Double>create();

            while ((iters > currentIter) && (currentIter < maxIterations)) {
                // generate a random subset of samples
                subsetSelector.computeRandomSubsets(subsetSize, subsetIndices);

                // clear list of preliminary solutions before calling listener
                iterResults.clear();
                // compute solution for current iteration
                listener.estimatePreliminarSolutions(subsetIndices, iterResults);

                // iterate over all solutions that have been found
                for (final var iterResult : iterResults) {
                    // compute inliers
                    computeInliers(iterResult, threshold, residualsTemp, listener, sorter, inliersData);

                    // save solution that  minimizes the median residual
                    if (inliersData.isMedianResidualImproved()) {
                        // keep current solution
                        bestResult = iterResult;

                        // keep the best inliers data corresponding to best solution
                        // in case it can be useful along with the result
                        bestResultInliersData = inliersData;
                        bestMedianResidual = inliersData.getBestMedianResidual();
                    }

                    // if number of inliers have improved, update number of
                    // remaining iterations
                    currentInliers = inliersData.getNumInliers();
                    if (currentInliers > bestNumInliers) {
                        // update best number of inliers
                        bestNumInliers = currentInliers;

                        // keep inliers data corresponding to best number of
                        // inliers
                        bestNumberInliersData = inliersData;

                        // recompute number of times the algorithm needs to be
                        // executed depending on current number of inliers to
                        // achieve with probability mConfidence that we have
                        // inliers and probability 1 - mConfidence that we have
                        // outliers
                        final var probSubsetAllInliers = Math.pow((double) bestNumInliers / (double) totalSamples,
                                subsetSize);

                        if (Math.abs(probSubsetAllInliers) < Double.MIN_VALUE || Double.isNaN(probSubsetAllInliers)) {
                            newNIters = Integer.MAX_VALUE;
                        } else {
                            final var logProbSomeOutliers = Math.log(1.0 - probSubsetAllInliers);
                            if (Math.abs(logProbSomeOutliers) < Double.MIN_VALUE || Double.isNaN(logProbSomeOutliers)) {
                                newNIters = Integer.MAX_VALUE;
                            } else {
                                newNIters = (int) Math.ceil(Math.abs(Math.log(1.0 - confidence) / logProbSomeOutliers));
                            }
                        }
                        if (newNIters < iters) {
                            iters = newNIters;
                        }
                    }

                    // reset inliers data if either residual or number of inliers
                    // improved
                    if (inliersData.isMedianResidualImproved()) {
                        // create new inliers data instance until a new best solution
                        // is found
                        inliersData = new MSACInliersData(totalSamples);
                        // update the best median residual on new instance so that
                        // only better solutions that are found later can update
                        // inliers data
                        inliersData.update(bestMedianResidual, inliersData.getInliers(), inliersData.getResiduals(),
                                inliersData.getNumInliers(), false);
                    }
                }

                if (iters > 0) {
                    progress = Math.min((float) currentIter / (float) iters, 1.0f);
                } else {
                    progress = 1.0f;
                }
                if (progress - previousProgress > progressDelta) {
                    previousProgress = progress;
                    listener.onEstimateProgressChange(this, progress);
                }
                currentIter++;

                listener.onEstimateNextIteration(this, currentIter);
            }

            // no solution could be found after completing all iterations
            if (bestResult == null) {
                throw new RobustEstimatorException();
            }

            listener.onEstimateEnd(this);

            return bestResult;
        } catch (final SubsetSelectorException e) {
            throw new RobustEstimatorException(e);
        } finally {
            locked = false;
        }
    }

    /**
     * Returns data about inliers once estimation has been done.
     *
     * @return data about inliers or null if estimation has not been done.
     */
    @Override
    public InliersData getInliersData() {
        return getBestNumberInliersData();
    }

    /**
     * Returns method being used for robust estimation.
     *
     * @return method being used for robust estimation.
     */
    @Override
    public RobustEstimatorMethod getMethod() {
        return RobustEstimatorMethod.MSAC;
    }

    /**
     * Computes inliers data for current iteration.
     *
     * @param <T>           type of result to be estimated.
     * @param iterResult    result to be tested on current iteration.
     * @param threshold     threshold to determine whether samples are inliers or
     *                      not.
     * @param residualsTemp temporal array to store residuals, since median
     *                      computation requires modifying the original array.
     * @param listener      listener to obtain residuals for samples.
     * @param sorter        sorter instance to compute median of residuals.
     * @param inliersData   inliers data to be reused on each iteration
     */
    private static <T> void computeInliers(
            final T iterResult, final double threshold, final double[] residualsTemp,
            final LMedSRobustEstimatorListener<T> listener, final Sorter<Double> sorter,
            final MSACInliersData inliersData) {

        final var residuals = inliersData.getResiduals();
        final var inliers = inliersData.getInliers();
        var bestMedianResidual = inliersData.getBestMedianResidual();
        var medianResidualImproved = false;

        final var totalSamples = residuals.length;
        double residual;
        var numInliers = 0;
        // find residuals and inliers
        for (var i = 0; i < totalSamples; i++) {
            residual = Math.abs(listener.computeResidual(iterResult, i));
            if (residual < threshold) {
                residuals[i] = residual;
                numInliers++;
                inliers.set(i);
            } else {
                residuals[i] = threshold;
                inliers.clear(i);
            }
        }

        // compute median of residuals
        System.arraycopy(residuals, 0, residualsTemp, 0, residuals.length);
        final var medianResidual = sorter.median(residualsTemp);
        if (medianResidual < bestMedianResidual) {
            bestMedianResidual = medianResidual;
            medianResidualImproved = true;
        }


        // store values in inliers data, only if residuals improve
        if (medianResidualImproved) {
            inliersData.update(bestMedianResidual, inliers, residuals, numInliers, true);
        }
    }

    /**
     * Contains data related to inliers estimated in one iteration.
     */
    public static class MSACInliersData extends InliersData {
        /**
         * Best median of error found so far taking into account all provided
         * samples.
         */
        private double bestMedianResidual;

        /**
         * Efficiently stores which samples are considered inliers and which
         * ones aren't.
         */
        private BitSet inliers;

        /**
         * Indicates whether median residual computed in current iteration has
         * improved respect to previous iterations.
         */
        private boolean medianResidualImproved;

        /**
         * Constructor.
         *
         * @param totalSamples total number of samples.
         */
        protected MSACInliersData(final int totalSamples) {
            bestMedianResidual = Double.MAX_VALUE;
            inliers = new BitSet(totalSamples);
            residuals = new double[totalSamples];
            numInliers = 0;
            medianResidualImproved = false;
        }

        /**
         * Returns best median of error found so far taking into account all
         * provided samples.
         *
         * @return best median of error found so far taking into account all
         * provided samples.
         */
        public double getBestMedianResidual() {
            return bestMedianResidual;
        }

        /**
         * Returns efficient array indicating which samples are considered
         * inliers and which ones aren't.
         *
         * @return array indicating which samples are considered inliers and
         * which ones aren't.
         */
        @Override
        public BitSet getInliers() {
            return inliers;
        }

        /**
         * Returns boolean indicating whether median residual computed in
         * current iteration has improved respect to previous iterations.
         *
         * @return true if median residual improved, false otherwise.
         */
        public boolean isMedianResidualImproved() {
            return medianResidualImproved;
        }

        /**
         * Updates data contained in this instance.
         *
         * @param bestMedianResidual     best median of error found so far taking
         *                               into account all provided samples.
         * @param inliers                efficiently stores which samples are considered
         *                               inliers and which ones aren't.
         * @param residuals              residuals obtained for each sample of data.
         * @param numInliers             number of inliers found on current iteration.
         * @param medianResidualImproved indicates whether median residual
         *                               computed in current iteration has improved respect to previous
         *                               iteration.
         */
        protected void update(final double bestMedianResidual, final BitSet inliers,
                              final double[] residuals, final int numInliers, final boolean medianResidualImproved) {
            this.bestMedianResidual = bestMedianResidual;
            this.inliers = inliers;
            this.residuals = residuals;
            this.numInliers = numInliers;
            this.medianResidualImproved = medianResidualImproved;
        }
    }
}
