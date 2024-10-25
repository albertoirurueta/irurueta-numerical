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

import java.util.ArrayList;
import java.util.BitSet;

/**
 * This class implements RANSAC (RANdom SAmple Consensus) algorithm to robustly
 * estimate a data model.
 * RANSAC is based on the idea that given a proportion of outliers (that is
 * estimated while the algorithm is run), it is needed a certain number of
 * random sub-samples to obtain required data with a certain level of confidence.
 * To determine whether a sample is an outlier or not, provided constant
 * threshold is used. If a threshold cannot be determined beforehand, then it
 * is better to use LMedS algorithm, at the expense of slightly less accurate
 * results.
 *
 * @param <T> type of object to be estimated.
 */
@SuppressWarnings("DuplicatedCode")
public class RANSACRobustEstimator<T> extends RobustEstimator<T> {

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
     * Minimum allowed threshold to determine inliers.
     */
    public static final double MIN_THRESHOLD = 0.0;

    /**
     * Indicates that by default inliers will only be computed but not kept.
     */
    public static final boolean DEFAULT_COMPUTE_AND_KEEP_INLIERS = false;

    /**
     * Indicates that by default residuals will only be computed but not kept.
     */
    public static final boolean DEFAULT_COMPUTE_AND_KEEP_RESIDUALS = false;

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
    private int nIters;

    /**
     * Best solution that has been found so far during an estimation.
     */
    private T bestResult;

    /**
     * Data related to inliers found for best result.
     */
    private RANSACInliersData bestInliersData;

    /**
     * Indicates whether inliers must be computed and kept.
     */
    private boolean computeAndKeepInliers;

    /**
     * Indicates whether residuals must be computed and kept.
     */
    private boolean computeAndKeepResiduals;

    /**
     * Constructor.
     */
    public RANSACRobustEstimator() {
        super();
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        nIters = maxIterations;
        bestResult = null;
        bestInliersData = null;
        computeAndKeepInliers = DEFAULT_COMPUTE_AND_KEEP_INLIERS;
        computeAndKeepResiduals = DEFAULT_COMPUTE_AND_KEEP_RESIDUALS;
    }

    /**
     * Constructor with listener.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes, as well as in charge
     *                 of picking samples and doing per-iteration estimations.
     */
    public RANSACRobustEstimator(final RANSACRobustEstimatorListener<T> listener) {
        super(listener);
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        nIters = maxIterations;
        bestResult = null;
        bestInliersData = null;
        computeAndKeepInliers = DEFAULT_COMPUTE_AND_KEEP_INLIERS;
        computeAndKeepResiduals = DEFAULT_COMPUTE_AND_KEEP_RESIDUALS;
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
        return nIters;
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
     * Gets data related to inliers found for best result.
     *
     * @return data related to inliers found for best result.
     */
    public RANSACInliersData getBestInliersData() {
        return bestInliersData;
    }

    /**
     * Indicates whether inliers must be computed and kept.
     *
     * @return true if inliers must be computed and kept, false if inliers
     * only need to be computed but not kept.
     */
    public boolean isComputeAndKeepInliersEnabled() {
        return computeAndKeepInliers;
    }

    /**
     * Specifies whether inliers must be computed and kept.
     *
     * @param computeAndKeepInliers true if inliers must be computed and kept,
     *                              false if inliers only need to be computed but not kept.
     * @throws LockedException if estimator is locked.
     */
    public void setComputeAndKeepInliersEnabled(final boolean computeAndKeepInliers) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.computeAndKeepInliers = computeAndKeepInliers;
    }

    /**
     * Indicates whether residuals must be computed and kept.
     *
     * @return true if residuals must be computed and kept, false if residuals
     * only need to be computed but not kept.
     */
    public boolean isComputeAndKeepResidualsEnabled() {
        return computeAndKeepResiduals;
    }

    /**
     * Specifies whether residuals must be computed and kept.
     *
     * @param computeAndKeepResiduals true if residuals must be computed and
     *                                kept, false if residuals only need to be computed but not kept.
     * @throws LockedException if estimator is locked.
     */
    public void setComputeAndKeepResidualsEnabled(final boolean computeAndKeepResiduals) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.computeAndKeepResiduals = computeAndKeepResiduals;
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
        return (listener instanceof RANSACRobustEstimatorListener);
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
    public T estimate() throws LockedException, NotReadyException, RobustEstimatorException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        try {
            final var listener = (RANSACRobustEstimatorListener<T>) this.listener;

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
            nIters = Integer.MAX_VALUE;
            int newNIters;
            var currentIter = 0;
            // reusable list that will contain preliminary solutions on each
            // iteration
            final var iterResults = new ArrayList<T>();
            bestResult = null;
            int currentInliers;
            var previousProgress = 0.0f;
            float progress;
            final var subsetIndices = new int[subsetSize];

            BitSet inliers = null;
            double[] residuals = null;
            if (computeAndKeepInliers || computeAndKeepResiduals) {
                bestInliersData = new RANSACInliersData(totalSamples, computeAndKeepInliers, computeAndKeepResiduals);

                if (computeAndKeepInliers) {
                    inliers = new BitSet(totalSamples);
                }
                if (computeAndKeepResiduals) {
                    residuals = new double[totalSamples];
                }
            }

            if (subsetSelector == null) {
                // create new subset selector
                subsetSelector = SubsetSelector.create(totalSamples);
            } else {
                // set number of samples to current subset selector
                subsetSelector.setNumSamples(totalSamples);
            }

            while ((nIters > currentIter) && (currentIter < maxIterations)) {
                // generate a random subset of samples
                subsetSelector.computeRandomSubsets(subsetSize, subsetIndices);

                // clear list of preliminary solutions before calling listener
                iterResults.clear();
                // compute solution for current iteration
                listener.estimatePreliminarSolutions(subsetIndices, iterResults);

                for (final var iterResult : iterResults) {
                    // compute number of inliers
                    currentInliers = 0;
                    for (var i = 0; i < totalSamples; i++) {
                        final var error = listener.computeResidual(iterResult, i);
                        if (error <= threshold) {
                            currentInliers++;
                            // keep inlier data if needed
                            if (inliers != null) {
                                inliers.set(i);
                            }
                        } else {
                            // outlier
                            if (inliers != null) {
                                inliers.clear(i);
                            }
                        }

                        if (residuals != null) {
                            residuals[i] = error;
                        }
                    }

                    // save result that produces the largest number of inliers
                    // and update number of iterations
                    if (currentInliers > bestNumInliers) {
                        // update best number of inliers
                        bestNumInliers = currentInliers;
                        // keep current result
                        bestResult = iterResult;

                        // update best inlier data
                        if (bestInliersData != null) {
                            bestInliersData.update(inliers, residuals, bestNumInliers);
                        }

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
                        if (newNIters < nIters) {
                            nIters = newNIters;
                        }
                    }
                }

                if (nIters > 0) {
                    progress = Math.min((float) currentIter / (float) nIters, 1.0f);
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
        return getBestInliersData();
    }


    /**
     * Returns method being used for robust estimation.
     *
     * @return method being used for robust estimation.
     */
    @Override
    public RobustEstimatorMethod getMethod() {
        return RobustEstimatorMethod.RANSAC;
    }

    /**
     * Contains data related to estimated inliers.
     */
    public static class RANSACInliersData extends InliersData {

        /**
         * Efficiently stores which samples are considered inliers and which
         * ones aren't.
         */
        private BitSet inliers;

        /**
         * Constructor.
         *
         * @param totalSamples  total number of samples.
         * @param keepInliers   true to keep inliers, false otherwise.
         * @param keepResiduals true to keep residuals, false otherwise.
         */
        protected RANSACInliersData(final int totalSamples, final boolean keepInliers, final boolean keepResiduals) {
            if (keepInliers) {
                inliers = new BitSet(totalSamples);
            }
            if (keepResiduals) {
                residuals = new double[totalSamples];
            }
            numInliers = 0;
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
         * Updates data contained in this instance.
         *
         * @param inliers    efficiently stores which samples are considered
         *                   inliers and which ones aren't.
         * @param residuals  residuals obtained for each sample of data.
         * @param numInliers number of inliers found on current iteration.
         */
        protected void update(final BitSet inliers, final double[] residuals, final int numInliers) {
            var totalSamples = 0;
            if (inliers != null) {
                totalSamples = inliers.length();
            }
            if (residuals != null) {
                totalSamples = residuals.length;
            }

            if (this.inliers != null && inliers != null && this.residuals != null && residuals != null) {
                // update inliers and residuals
                for (var i = 0; i < totalSamples; i++) {
                    this.inliers.set(i, inliers.get(i));
                    this.residuals[i] = residuals[i];
                }
            } else if (this.inliers != null && inliers != null) {
                // update inliers but not the residuals
                for (var i = 0; i < totalSamples; i++) {
                    this.inliers.set(i, inliers.get(i));
                }
            } else if (this.residuals != null && residuals != null) {
                // update residuals but not inliers
                System.arraycopy(residuals, 0, this.residuals, 0, totalSamples);
            }
            this.numInliers = numInliers;
        }
    }
}
