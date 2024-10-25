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
import com.irurueta.sorting.SortingException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

/**
 * This class implements PROMedS (PROgressive least Median Sample) algorithm
 * to robustly estimate a data model.
 * This algorithm is a mixture between LMedS and PROSAC, taking the best of
 * both.
 * Firstly, it has the advantage that no threshold is required to be set
 * beforehand, the same as LMedS. Threshold to determine inliers is computed
 * dynamically, which helps for an easier setup that is problem independent and
 * depending on the accuracy of the inliers, results will be more accurate than
 * RANSAC or PROSAC, just the same as LMedS.
 * On the other hand, if certain information about the quality of the samples
 * is available, as in PROSAC, the algorithm takes advantage of this additional
 * information to prioritize the samples with higher quality in order to find
 * a solution much faster than RANSAC or LMedS.
 * Finally, if by any chance a threshold to determine inliers is also used, the
 * algorithm will try to get the solution that better fits in a pure median of
 * residuals model or in a threshold based one to determine inliers.
 * Hence, PROMedS can be as fast as PROSAC (which is typically about 100x faster
 * than RANSAC or LMedS), can obtain the same accuracy as LMedS (which can be
 * much better than RANSAC or PROSAC in certain scenarios), and has an easier
 * setup, which is problem independent because no threshold is required to be
 * known beforehand although one can be provided as well.
 *
 * @param <T> type of object to be estimated.
 */
@SuppressWarnings("Duplicates")
public class PROMedSRobustEstimator<T> extends RobustEstimator<T> {
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
     * Indicates whether the algorithm must stop prematurely when dynamically
     * computed threshold using median of residuals has a value lower than
     * provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     * By default, stop threshold is enabled, so that computational cost is
     * similar to RANSAC and only accuracy gets better if inliers are more
     * accurate.
     */
    public static final boolean DEFAULT_STOP_THRESHOLD_ENABLED = true;

    /**
     * Default maximum allowed outliers proportion in the input data. This is
     * used do determine the number of required iterations.
     */
    public static final double DEFAULT_MAX_OUTLIERS_PROPORTION = 0.8;

    /**
     * Minimum allowed value for maximum allowed outliers proportion in the
     * input data.
     */
    public static final double MIN_MAX_OUTLIERS_PROPORTION = 0.0;

    /**
     * Maximum allowed value for maximum allowed outliers proportion in the
     * input data.
     */
    public static final double MAX_MAX_OUTLIERS_PROPORTION = 1.0;

    /**
     * Defines the default value for the maximum probability that a solution
     * with more than inliersNStar in U_nStar exist and was not found after k
     * samples.
     */
    public static final double DEFAULT_ETA0 = 0.05;

    /**
     * Minimum allowed value for eta0.
     */
    public static final double MIN_ETA0 = 0.0;

    /**
     * Maximum allowed value for eta0.
     */
    public static final double MAX_ETA0 = 1.0;

    /**
     * Defines the default value for beta, which is the probability that a
     * match is declared inlier by mistake, i.e. the ratio of the "inlier"
     * surface by the total surface. The inlier surface is a disc with radius
     * 1.96s for homography/displacement computation, or a band with width
     * 1.96*s*2 for epipolar geometry (s is the detection noise), and the total
     * surface is the surface of the image.
     */
    public static final double DEFAULT_BETA = 0.01;

    /**
     * Minimum allowed value for beta.
     */
    public static final double MIN_BETA = 0.0;

    /**
     * Maximum allowed value for beta.
     */
    public static final double MAX_BETA = 1.0;

    /**
     * Default factor to normalize threshold to determine inliers. This factor
     * can be used to increase or lower the dynamically computed threshold so
     * that the algorithm becomes more or less accurate. The stricter the
     * threshold (lower factor), the more time the algorithm will need to
     * converge, if it can converge. By default, the factor is 1.0, which makes
     * the threshold to be computed as the median of residuals.
     */
    public static final double DEFAULT_INLIER_FACTOR = 1.0; //1.5 would also be reasonable

    /**
     * Minimum allowed value for inlier factor.
     */
    public static final double MIN_INLER_FACTOR = 0.0;

    /**
     * Indicates whether the inlier threshold will be used to find inliers along
     * with their median of residuals.
     */
    public static final boolean DEFAULT_USE_INLIER_THRESHOLD = true;

    /**
     * Constant to estimate standard deviation of residuals based on their
     * median.
     */
    public static final double STD_CONSTANT = 1.4826;

    /**
     * Chi squared.
     */
    public static final double CHI_SQUARED = 2.706;

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
     * In this implementation, PROSAC won't stop before having reached the
     * corresponding inliers rate on the complete data set.
     * Maximum allowed outliers proportion in the input data: used to compute
     * nIters (can be as high as 0.95).
     */
    private double maxOutliersProportion;

    /**
     * eta0 is the maximum probability that a solution with more than
     * inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     */
    private double eta0;

    /**
     * beta is the probability that a match is declared inlier by mistake,
     * i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!.
     */
    private double beta;

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
    private PROMedSInliersData bestInliersData;

    /**
     * Indicates whether the algorithm must stop prematurely when dynamically
     * computed threshold using median of residuals has a value lower than
     * provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     */
    private boolean stopThresholdEnabled;

    /**
     * Factor to normalize threshold to determine inliers. This factor can be
     * used to increase or lower the dynamically computed threshold so that the
     * algorithm becomes more or less accurate. The stricter the threshold
     * (lower factor), the more time the algorithm will need to converge, if
     * it can converge. By default, the factor is 1.0, which makes the threshold
     * to be computed as the median of residuals.
     */
    private double inlierFactor;

    /**
     * Flag indicating whether thresholds to determine inliers are used, or if
     * only median of residuals is used. When true, the algorithm will try
     * to fit the best model, otherwise only median of residuals will be used.
     */
    private boolean useInlierThresholds;

    /**
     * Constructor.
     */
    public PROMedSRobustEstimator() {
        super();
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        maxOutliersProportion = DEFAULT_MAX_OUTLIERS_PROPORTION;
        eta0 = DEFAULT_ETA0;
        beta = DEFAULT_BETA;
        iters = maxIterations;
        bestResult = null;
        bestInliersData = null;
        stopThresholdEnabled = DEFAULT_STOP_THRESHOLD_ENABLED;
        inlierFactor = DEFAULT_INLIER_FACTOR;
        useInlierThresholds = DEFAULT_USE_INLIER_THRESHOLD;
    }

    /**
     * Constructor with listener.
     *
     * @param listener listener to be notified of events such as when estimation
     *                 starts, ends or its progress significantly changes, as well as in charge
     *                 of picking samples and doing per-iteration estimations.
     */
    public PROMedSRobustEstimator(final PROMedSRobustEstimatorListener<T> listener) {
        super(listener);
        confidence = DEFAULT_CONFIDENCE;
        maxIterations = DEFAULT_MAX_ITERATIONS;
        maxOutliersProportion = DEFAULT_MAX_OUTLIERS_PROPORTION;
        eta0 = DEFAULT_ETA0;
        beta = DEFAULT_BETA;
        iters = maxIterations;
        bestResult = null;
        bestInliersData = null;
        stopThresholdEnabled = DEFAULT_STOP_THRESHOLD_ENABLED;
        inlierFactor = DEFAULT_INLIER_FACTOR;
        useInlierThresholds = DEFAULT_USE_INLIER_THRESHOLD;
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
     * Returns boolean indicating whether the algorithm must stop prematurely
     * when dynamically computed threshold using median of residuals has a value
     * lower than provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     *
     * @return true if stop threshold is enabled, false otherwise.
     */
    public boolean isStopThresholdEnabled() {
        return stopThresholdEnabled;
    }

    /**
     * Sets boolean indicating whether the algorithm must stop prematurely when
     * dynamically computed threshold using median of residuals has a value
     * lower than provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     *
     * @param stopThresholdEnabled true if stop threshold is enabled, false
     *                             otherwise.
     * @throws LockedException if this estimator is locked because an estimation
     *                         is being computed.
     */
    public void setStopThresholdEnabled(final boolean stopThresholdEnabled) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.stopThresholdEnabled = stopThresholdEnabled;
    }

    /**
     * Returns factor to normalize or adjust threshold to determine inliers.
     * This factor can be used to increase or lower the dynamically computed
     * threshold so that the algorithm becomes more or less accurate. The
     * stricter the threshold (lower factor), the more time the algorithm will
     * need to converge, if it can converge. By default, the factor is 1.0, which
     * makes the threshold to be computed as the median of residuals.
     *
     * @return factor to normalize threshold to determine inliers.
     */
    public double getInlierFactor() {
        return inlierFactor;
    }

    /**
     * Sets factor to normalize or adjust threshold to determine inliers.
     * This factor can be used to increase or lower the dynamically computed
     * threshold so that the algorithm becomes more or less accurate. The
     * stricter the threshold (lower factor), the more time the algorithm will
     * need to converge, if it can converge. By default, the factor is 1.0, which
     * makes the threshold to be computed as the median of residuals.
     *
     * @param inlierFactor inlier factor to be set.
     * @throws IllegalArgumentException if provided value is less or equal than
     *                                  0.0.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setInlierFactor(final double inlierFactor) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (inlierFactor <= MIN_INLER_FACTOR) {
            throw new IllegalArgumentException();
        }
        this.inlierFactor = inlierFactor;
    }

    /**
     * Returns flag indicating whether thresholds to determine inliers are used,
     * or if only median of residuals is used. When true, the algorithm will try
     * to fit the best model, otherwise only median of residuals will be used.
     *
     * @return true if best model is used (threshold or median), otherwise only
     * median of residuals will be used.
     */
    public boolean isUseInlierThresholds() {
        return useInlierThresholds;
    }

    /**
     * Sets flag indicating whether thresholds to determine inliers are used, or
     * if only median of residuals is used. When true, the algorithm will try to
     * fit the best model, otherwise only median of residuals will be used.
     *
     * @param useInlierThresholds true if best model is used (threshold or
     *                            median), oitherwise only median of residuals will be used.
     * @throws LockedException if this estimator is locked because an estimation
     *                         is being computed.
     */
    public void setUseInlierThresholds(final boolean useInlierThresholds) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }

        this.useInlierThresholds = useInlierThresholds;
    }

    /**
     * Returns maximum allowed outliers proportion in the input data. This is
     * used to compute number of iterations to be done (nIters). It typically
     * can be as high as 0.95. Higher values, up to 1 are possible but not
     * recommended.
     * In this implementation, PROSAC won't stop before having reached the
     * corresponding inliers rate on the complete data set.
     *
     * @return maximum allowed outliers proportion in the input data.
     */
    public double getMaxOutliersProportion() {
        return maxOutliersProportion;
    }

    /**
     * Sets maximum allowed outliers proportion in the input data. This is used
     * to compute number of iterations to be done (nIters). It typically can be
     * as high as 0.95. Higher values, up to 1 are possible but not recommended.
     * In this implementation, PROSAC won't stop before having reached the
     * corresponding inliers rate on the complete data set.
     *
     * @param maxOutliersProportion maximum allowed outliers proportion in the
     *                              input data.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     *                                  greater than 1.0.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setMaxOutliersProportion(final double maxOutliersProportion) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxOutliersProportion < MIN_MAX_OUTLIERS_PROPORTION
                || maxOutliersProportion > MAX_MAX_OUTLIERS_PROPORTION) {
            throw new IllegalArgumentException();
        }

        this.maxOutliersProportion = maxOutliersProportion;
    }

    /**
     * Return eta0, which is the maximum probability that a solution with more
     * than inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     *
     * @return eta0 value.
     */
    public double getEta0() {
        return eta0;
    }

    /**
     * Sets eta0, which is the maximum probability that a solution with more
     * than inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     *
     * @param eta0 eta0 value to be set.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     *                                  greater than 1.0.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setEta0(final double eta0) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (eta0 < MIN_ETA0 || eta0 > MAX_ETA0) {
            throw new IllegalArgumentException();
        }

        this.eta0 = eta0;
    }

    /**
     * Returns beta, which is the probability that a match is declared inlier by
     * mistake, i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!
     *
     * @return beta value.
     */
    public double getBeta() {
        return beta;
    }

    /**
     * Sets beta, which is the probability that a match is declared inlier by
     * mistake, i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!
     *
     * @param beta beta value to be set.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     *                                  greater than 1.0.
     * @throws LockedException          if this estimator is locked because an estimation
     *                                  is being computed.
     */
    public void setBeta(final double beta) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (beta < MIN_BETA || beta > MAX_BETA) {
            throw new IllegalArgumentException();
        }

        this.beta = beta;
    }

    /**
     * Returns number of iterations to be done to obtain required confidence.
     * This does not need to be equal to the actual number of iterations the
     * algorithm finally required to obtain a solution.
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
     * Returns data related to inliers found for best result.
     *
     * @return data related to inliers found for best result.
     */
    protected PROMedSInliersData getBestInliersData() {
        return bestInliersData;
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
        return (listener instanceof PROMedSRobustEstimatorListener);
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
            var listener = (PROMedSRobustEstimatorListener<T>) this.listener;

            locked = true;

            listener.onEstimateStart(this);

            // N = CORRESPONDENCES
            final var totalSamples = listener.getTotalSamples();
            final var subsetSize = listener.getSubsetSize();

            final var qualityScores = listener.getQualityScores();
            // check for invalid quality scores length
            if (qualityScores.length != totalSamples) {
                throw new RobustEstimatorException();
            }

            var inlierThreshold = 0.0;
            if (useInlierThresholds) {
                inlierThreshold = listener.getThreshold();
            }
            // obtain indices referring to original samples position after sorting
            // quality scores in descending order
            final var sortedIndices = computeSortedQualityIndices(listener.getQualityScores());

            // reusable list that will contain preliminary solutions on each
            // iteration
            final var iterResults = new ArrayList<T>();
            bestResult = null;
            var previousProgress = 0.0f;
            float progress;
            // subset indices obtained from a subset selector
            final var subsetIndices = new int[subsetSize];
            final var residualsTemp = new double[totalSamples];
            // subset indices referred to the real samples positions after taking
            // into account the sorted indices obtained from quality scores
            final var transformedSubsetIndices = new int[subsetSize];
            // array containing inliers efficiently
            final var inliers = new BitSet(totalSamples);

            // T_N
            iters = Math.min(computeIterations(1.0 - maxOutliersProportion, subsetSize, confidence),
                    maxIterations);

            // termination length
            var sampleSizeStar = totalSamples;
            // number of inliers found within the first
            // nStar data points
            var inliersNStar = 0;
            // best number of inliers found so far
            // (store the model that goes with it)
            var inliersBest = -1;
            // threshold to stop algorithm
            var threshold = Double.MAX_VALUE;
            final var inliersMin = (int) ((1.0 - maxOutliersProportion) * totalSamples);
            // iteration number (t)
            var currentIter = 0;
            // (n) we draw samples from the set U_n
            // of the top n (sampleSize) data points
            var sampleSize = subsetSize;
            // average number of samples "{M_i}_{i=1}^{Tn}"
            // that contains samples from U_n only
            double tn = iters;
            // integer version of Tn
            var tnPrime = 1;
            // number of samples to draw to reach the
            // maximality constraint
            var kNStar = iters;

            // initialize Tn
            for (var i = 0; i < subsetSize; i++) {
                tn *= (double) (sampleSize - i) / (double) (totalSamples - i);
            }

            if (subsetSelector == null) {
                // create new subset selector
                subsetSelector = SubsetSelector.create(totalSamples);
            } else {
                // set number of samples to current subset selector
                subsetSelector.setNumSamples(totalSamples);
            }

            // data related to inliers
            var inliersData = new PROMedSInliersData(totalSamples);
            // sorter to compute medians
            final var sorter = Sorter.<Double>create();

            // indicates if result improved
            boolean improved;
            var continueIteration = true;

            // iterate until the expected number of inliers or the estimated
            // number of iterations is reached
            while (continueIteration) {
                if (kNStar > 0) {
                    progress = Math.min((float) currentIter / (float) kNStar, 1.0f);
                } else {
                    progress = 1.0f;
                }
                if (progress - previousProgress > progressDelta) {
                    previousProgress = progress;
                    listener.onEstimateProgressChange(this, progress);
                }
                currentIter++;

                // choice of the hypothesis generation set

                // The growth function is defined as
                // g(t) = min{n : TnPrime > t} where n is sampleSize
                // Thus sampleSize should be incremented if currentIter > TnPrime
                if ((currentIter > tnPrime) && (sampleSize < sampleSizeStar)) {
                    final var TnPlus1 = (tn * (sampleSize + 1)) / (sampleSize + 1 - subsetSize);
                    sampleSize++;
                    tnPrime += (int) Math.ceil(TnPlus1 - tn);
                    tn = TnPlus1;
                }

                // Draw semi-random sample
                if (currentIter > tnPrime) {
                    // during the finishing stage (sampleSize == sampleSizeStar &&
                    // currentIter > TnPrime), draw a standard RANSAC sample
                    // The sample contains subsetSize points selected from U_n at
                    // random
                    subsetSelector.computeRandomSubsets(subsetSize, subsetIndices);
                } else {
                    // The sample contains subsetSize-1 points selected from
                    // U_sampleSize_1 at random and u_sampleSize

                    subsetSelector.computeRandomSubsetsInRange(0, sampleSize, subsetSize, true,
                            subsetIndices);
                }

                transformIndices(subsetIndices, sortedIndices, transformedSubsetIndices);

                // clear list of preliminary solutions before calling listener
                iterResults.clear();
                // compute solution for current iteration
                listener.estimatePreliminarSolutions(transformedSubsetIndices, iterResults);

                // total number of inliers for a
                // given result
                int inliersCurrent;

                // iterate over all solutions that have been found
                improved = false;
                for (final var iterResult : iterResults) {
                    // compute inliers
                    computeInliers(iterResult, subsetSize, inlierFactor, useInlierThresholds, inlierThreshold,
                            residualsTemp, listener, sorter, inliersData);
                    inliersCurrent = inliersData.getNumInliers();

                    if (inliersData.isMedianResidualImproved()) {
                        improved = true;

                        // keep current solution
                        bestResult = iterResult;

                        // update estimated thresholds to be used as stop
                        // criterion
                        threshold = inliersData.getEstimatedThreshold();
                    }

                    if (inliersCurrent > inliersBest) {
                        // update best number of inliers
                        inliersBest = inliersCurrent;
                        // keep current solution
                        bestResult = iterResult;

                        keepInliersData(inliersData, totalSamples);
                        // create new instance after keeping inlier data
                        inliersData = inliersData.createCopy();

                        // select new termination length sampleSizeStar if possible
                        // only when a new sample is better than the others found
                        // so far

                        // best value found so far in terms of inliers ratio
                        var sampleSizeBest = totalSamples;
                        var inliersSampleSizeBest = inliersCurrent;
                        // test value for the termination length
                        int sampleSizeTest;
                        // number of inliers for that test value
                        int inliersSampleSizeTest;
                        var epsilonSampleSizeBest = (double) inliersSampleSizeBest / (double) sampleSizeBest;

                        for (sampleSizeTest = totalSamples, inliersSampleSizeTest = inliersCurrent;
                             sampleSizeTest > subsetSize; sampleSizeTest--) {
                            // Loop invariants:
                            // - inliersSampleSizeTest is the number of inliers
                            //   for the sampleSizeTest first correspondences
                            // - sampleSizeBest is the value between
                            //   sampleSizeTest+1 and totalSamples that maximizes
                            //   the ratio inliersSampleSizeBest/sampleSizeBest

                            // - Non-randomness: In >= imin(n*)
                            // - Maximality: the number of samples that were drawn
                            //   so far must be enough so that the probability of
                            //   having missed a set of inliers is below eta=0.01.
                            //   This is the classical RANSAC termination criterion,
                            //   except that it takes into account only the
                            //   sampleSize first samples (not the total number
                            //   of samples)
                            //   kNStar = log(eta0) / log(1 - (inliersNStar/
                            //   sampleSizeStar)^subsetSize
                            //   We have to minimize kNStar, e.g. maximize
                            //   inliersNStar/sampleSizeStar, a straightforward
                            //   implementation would use the following test:
                            //   if(inliersSampleSizeTest > epsilonSampleSizeBest *
                            //   sampleSizeTest){ ... blah blah blah
                            //   However, since In is binomial, and in the case of
                            //   evenly distributed inliers, a better test would be
                            //   to reduce sampleSizeStar only if there's a
                            //   significant improvement in epsilon. Thus, we use a
                            //   Chi-squared test (P=0.10), together with the normal
                            //   approximation to the binomial (mu =
                            //   epsilonSampleSizeStart * sampleSizeTest, sigma =
                            //   sqrt(sampleSizeTest * epsilonSampleSizeStar * (1 -
                            //   epsilonSampleSizeStar))).
                            //   There is a significant difference between the two
                            //   tests (e.g. with the computeInliers function
                            //   provided)
                            //   We do the cheap test first, and the expensive test
                            //   only if the cheap one passes
                            if ((inliersSampleSizeTest * sampleSizeBest > inliersSampleSizeBest * sampleSizeTest)
                                    && (inliersSampleSizeTest > epsilonSampleSizeBest * sampleSizeTest
                                    + Math.sqrt(sampleSizeTest * epsilonSampleSizeBest
                                    * (1.0 - epsilonSampleSizeBest) * CHI_SQUARED))) {

                                if (inliersSampleSizeTest < imin(subsetSize, sampleSizeTest, beta)) {
                                    // equation not satisfied, no need to test for
                                    // smaller sampleSizeTest values anyway

                                    // jump out of the for sampleSizeTest loop
                                    break;
                                }
                                sampleSizeBest = sampleSizeTest;
                                inliersSampleSizeBest = inliersSampleSizeTest;
                                epsilonSampleSizeBest = (double) inliersSampleSizeBest / (double) sampleSizeBest;
                            }

                            // prepare for next loop iteration
                            inliersSampleSizeTest -= inliers.get(sortedIndices[sampleSizeTest - 1]) ? 1 : 0;
                        }

                        // is the best one we found even better than sampleSizeStar?
                        if (inliersSampleSizeBest * sampleSizeStar > inliersNStar * sampleSizeBest) {

                            // update all values
                            sampleSizeStar = sampleSizeBest;
                            inliersNStar = inliersSampleSizeBest;
                            kNStar = computeIterations((double) inliersNStar
                                            / (double) sampleSizeStar, subsetSize, 1.0 - eta0);
                        }
                    }
                }

                continueIteration = (currentIter < maxIterations);
                if (useInlierThresholds && stopThresholdEnabled) {
                    // if inlier threshold is being used, and stop threshold is
                    // enabled, then stop the algorithm if threshold determined
                    // by median of residuals has a value lower than inlier
                    // threshold
                    continueIteration &= (threshold > inlierThreshold);
                }

                if (!improved) {
                    continueIteration &= ((inliersBest < inliersMin) || (currentIter < kNStar))
                            && (currentIter < iters);
                }

                listener.onEstimateNextIteration(this, currentIter);
            }

            // no solution could be found after completing all iterations
            if (bestResult == null) {
                throw new RobustEstimatorException();
            }

            listener.onEstimateEnd(this);

            return bestResult;
        } catch (final SubsetSelectorException | SortingException e) {
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
        return RobustEstimatorMethod.PROMEDS;
    }

    /**
     * Transforms indices picked by the subset selector into the indices where
     * samples are actually localed by taking into account their original
     * position before sorting quality scores.
     *
     * @param subsetIndices            indices picked by the subset selector. These are
     *                                 positions after sorting. Must have the subset length.
     * @param sortedIndices            indices relating sorted positions to their original
     *                                 positions. Each position i-th in the array refers to the original
     *                                 position before sorting. Must have the number of samples length.
     * @param transformedSubsetIndices array where result is stored. Must have
     *                                 the subset length.
     */
    private static void transformIndices(
            final int[] subsetIndices, final int[] sortedIndices, final int[] transformedSubsetIndices) {
        final var length = transformedSubsetIndices.length;
        for (var i = 0; i < length; i++) {
            transformedSubsetIndices[i] = sortedIndices[subsetIndices[i]];
        }
    }

    /**
     * Computes inliers data for current iteration.
     *
     * @param <T>                 type of result to be estimated.
     * @param iterResult          result to be tested on current iteration.
     * @param subsetSize          subset sample size to be picked on each iteration.
     * @param inlierFactor        factor to adjust threshold to determine whether
     *                            samples are inliers or not.
     * @param useInlierThresholds true to use thresholds to determine inliers,
     *                            false otherwise.
     * @param inlierThreshold     threshold to determine which samples are inliers.
     * @param residualsTemp       temporal array to store residuals, since median
     *                            computation requires modifying the original array.
     * @param listener            listener to obtain residuals for samples.
     * @param sorter              sorter instance to compute median of residuals.
     * @param inliersData         inliers data to be reused on each iteration
     */
    private static <T> void computeInliers(
            final T iterResult, final int subsetSize, final double inlierFactor, final boolean useInlierThresholds,
            final double inlierThreshold, final double[] residualsTemp, final LMedSRobustEstimatorListener<T> listener,
            final Sorter<Double> sorter, final PROMedSInliersData inliersData) {

        final var residuals = inliersData.getResiduals();
        final var lmedsInliers = inliersData.getInliersLMedS();
        final var msacInliers = inliersData.getInliersMSAC();
        var bestMedianResidual = inliersData.getBestMedianResidual();
        var medianResidualImproved = false;

        final var totalSamples = residuals.length;

        for (var i = 0; i < totalSamples; i++) {
            residuals[i] = Math.abs(listener.computeResidual(iterResult, i));
        }

        System.arraycopy(residuals, 0, residualsTemp, 0, residuals.length);
        final var medianResidual = sorter.median(residualsTemp);
        if (medianResidual < bestMedianResidual) {
            bestMedianResidual = medianResidual;
            medianResidualImproved = true;
        }

        final var standardDeviation = STD_CONSTANT * (1.0 + 5.0 / (totalSamples - subsetSize))
                * Math.sqrt(medianResidual);
        final var normEstimatedThreshold = inlierFactor * medianResidual;

        // determine which points are inliers
        var numInliers = 0;
        // by default if thresholds are not used
        var lmedsInlierModelEnabled = true;
        if (useInlierThresholds) {
            var numInliersMsac = 0;
            var numInliersLmedS = 0;
            for (var i = 0; i < totalSamples; i++) {
                if (residuals[i] <= normEstimatedThreshold) {
                    numInliersLmedS++;
                    lmedsInliers.set(i);
                } else {
                    lmedsInliers.clear(i);
                }
                if (residuals[i] <= inlierThreshold) {
                    numInliersMsac++;
                    msacInliers.set(i);
                } else {
                    msacInliers.clear(i);
                }
            }

            // keep model with smaller number of inliers (to be more restrictive)
            lmedsInlierModelEnabled = numInliersLmedS < numInliersMsac;
            numInliers = lmedsInlierModelEnabled ? numInliersLmedS : numInliersMsac;
        } else {

            for (var i = 0; i < totalSamples; i++) {
                if (residuals[i] <= normEstimatedThreshold) {
                    numInliers++;
                    lmedsInliers.set(i);
                } else {
                    lmedsInliers.clear(i);
                }
            }
        }

        // store values in inliers data, only if residuals improve
        if (medianResidualImproved) {
            inliersData.update(bestMedianResidual, standardDeviation, lmedsInlierModelEnabled, residuals, numInliers,
                    medianResidual, normEstimatedThreshold, true);
        } else {
            inliersData.medianResidualImproved = false;
        }
    }

    /**
     * Obtains indices of samples corresponding to samples ordered in descending
     * quality scores.
     *
     * @param qualityScores quality scores associated to each sample to be used
     *                      to obtain indices to sort samples in descending order of quality values.
     * @return indices to sort samples in descending order of quality values.
     * @throws SortingException if sorting fails.
     */
    private static int[] computeSortedQualityIndices(double[] qualityScores) throws SortingException {
        final var sorter = Sorter.<Double>create();
        final var qualityScoresCopy = Arrays.copyOf(qualityScores, qualityScores.length);
        // this method modifies quality scores copy array because it gets sorted
        // in ascending order. Indices contains indices of samples corresponding
        // to quality scores ordered in ascending order
        final var indices = sorter.sortWithIndices(qualityScoresCopy);

        // reverse indices so we have indices of samples ordered in descending
        // order of quality
        reverse(indices);

        return indices;
    }

    /**
     * Reverses provided array.
     *
     * @param array array to be reversed.
     */
    private static void reverse(final int[] array) {
        var length = array.length;
        for (var i = 0; i < length / 2; i++) {
            var temp = array[i];
            var pos = length - 1 - i;
            array[i] = array[pos];
            array[pos] = temp;
        }
    }

    /**
     * Computes number of required iterations to achieve required confidence
     * with current probability of inlier and sample subset size.
     *
     * @param probInlier probability of inlier.
     * @param subsetSize sample subset size.
     * @param confidence required confidence of result.
     * @return number of required iterations.
     */
    private static int computeIterations(final double probInlier, final int subsetSize, final double confidence) {

        // compute number of times the algorithm needs to be executed depending
        // on number of inliers respect total points to achieve with probability
        // confidence that we have all inliers and probability 1 - confidence
        // that we have some outliers
        final var probSubsetAllInliers = Math.pow(probInlier, subsetSize);
        if (Math.abs(probSubsetAllInliers) < Double.MIN_VALUE || Double.isNaN(probSubsetAllInliers)) {
            return Integer.MAX_VALUE;
        } else {
            final var logProbSomeOutliers = Math.log(1.0 - probSubsetAllInliers);
            if (Math.abs(logProbSomeOutliers) < Double.MIN_VALUE || Double.isNaN(logProbSomeOutliers)) {
                return Integer.MAX_VALUE;
            } else {
                return (int) Math.ceil(Math.abs(Math.log(1.0 - confidence) / logProbSomeOutliers));
            }
        }
    }

    /**
     * Non randomness states that i-m (where i is the cardinal of the set of
     * inliers for a wrong model) follows the binomial distribution B(n,beta).
     * For n big enough, B(n,beta) approximates to normal distribution N(mu,
     * sigma^2) by the central limit theorem, with mu = n*beta and sigma =
     * sqrt(n*beta*(1 - beta)).
     * Psi, the probability that In_star out of n_star data points are by chance
     * inliers to an arbitrary incorrect model, is set to 0.05 (5%, as in the
     * original paper), and you must change the Chi2 value if you chose a
     * different value for psi.
     *
     * @param subsetSize sample subset size.
     * @param sampleSize total number of samples.
     * @param beta       beta value.
     * @return i-m.
     */
    private static int imin(final int subsetSize, final int sampleSize, final double beta) {
        final var mu = sampleSize * beta;
        final var sigma = Math.sqrt(sampleSize * beta * (1.0 - beta));

        return (int) Math.ceil(subsetSize + mu + sigma * Math.sqrt(CHI_SQUARED));
    }

    /**
     * Keeps inliers data stored and initializes a new one with proper
     * initial values.
     *
     * @param inliersData  inliers data to be stored.
     * @param totalSamples total number of samples.
     */
    private void keepInliersData(PROMedSInliersData inliersData, final int totalSamples) {
        // keep the best inliers data corresponding to best solution,
        // in case it can be useful along with the result
        bestInliersData = inliersData;

        // create new inliers data instance until a new best solution
        // is found
        final var bestMedianResidual = inliersData.getBestMedianResidual();
        inliersData = new PROMedSInliersData(totalSamples);
        // update the best median residual on new instance so that
        // only better solutions that are found later can update
        // inliers data
        inliersData.update(bestMedianResidual, inliersData.getStandardDeviation(),
                inliersData.isLMedSInlierModelEnabled(), inliersData.getResiduals(), inliersData.getNumInliers(),
                bestMedianResidual, inliersData.getEstimatedThreshold(), false);
    }

    /**
     * Contains data related to inliers estimated in one iteration.
     */
    public static class PROMedSInliersData extends InliersData {
        /**
         * Best median of error found so far taking into account all provided
         * samples.
         */
        private double bestMedianResidual;

        /**
         * Standard deviation of error among all provided samples respect to
         * currently estimated result.
         */
        private double standardDeviation;

        /**
         * Inliers considering LMedS model.
         */
        private BitSet inliersLmeds;

        /**
         * Inliers considering MSAC model.
         */
        private BitSet inliersMsac;

        /**
         * Indicates whether LMedS or MSAC inlier model is enabled.
         */
        private boolean lmedsInlierModelEnabled;

        /**
         * Median of error found on current iteration among all provided
         * samples.
         */
        private double medianResidual;

        /**
         * Estimated threshold to determine whether samples are inliers or not.
         */
        private double estimatedThreshold;

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
        PROMedSInliersData(final int totalSamples) {
            bestMedianResidual = standardDeviation = medianResidual = estimatedThreshold = Double.MAX_VALUE;
            inliersLmeds = new BitSet(totalSamples);
            inliersMsac = new BitSet(totalSamples);
            lmedsInlierModelEnabled = true;
            residuals = new double[totalSamples];
            numInliers = 0;
            medianResidualImproved = false;
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
            return lmedsInlierModelEnabled ? inliersLmeds : inliersMsac;
        }

        /**
         * Creates a copy of inlier data.
         *
         * @return copy of inlier data.
         */
        PROMedSInliersData createCopy() {
            final var result = new PROMedSInliersData(residuals.length);
            result.bestMedianResidual = bestMedianResidual;
            result.standardDeviation = standardDeviation;
            result.medianResidual = medianResidual;
            result.estimatedThreshold = estimatedThreshold;
            result.inliersLmeds = (BitSet) inliersLmeds.clone();
            result.inliersMsac = (BitSet) inliersMsac.clone();
            result.lmedsInlierModelEnabled = lmedsInlierModelEnabled;
            result.residuals = Arrays.copyOf(residuals, residuals.length);
            result.numInliers = numInliers;
            result.medianResidualImproved = medianResidualImproved;

            return result;
        }

        /**
         * Returns best median of error found so far taking into account all
         * provided samples.
         *
         * @return best median of error found so far taking into account all
         * provided samples.
         */
        double getBestMedianResidual() {
            return bestMedianResidual;
        }

        /**
         * Returns standard deviation of error among all provided samples
         * respect to currently estimated result.
         *
         * @return standard deviation of error among all provided samples
         * respect to currently estimated result.
         */
        double getStandardDeviation() {
            return standardDeviation;
        }

        /**
         * Returns inliers considering LMedS model.
         *
         * @return inliers considering LMedS model.
         */
        BitSet getInliersLMedS() {
            return inliersLmeds;
        }

        /**
         * Returns inliers considering MSAC model.
         *
         * @return inliers considering MSAC model.
         */
        BitSet getInliersMSAC() {
            return inliersMsac;
        }

        /**
         * Returns boolean indicating whether LMedS or MSAC inlier model is
         * enabled. If true, estimated threshold was used to determine inliers,
         * if false only median of residuals was used.
         *
         * @return true if LMedS model is used, false if MSAC model is used.
         */
        boolean isLMedSInlierModelEnabled() {
            return lmedsInlierModelEnabled;
        }

        /**
         * Returns estimated threshold to determine whether samples are inliers
         * or not.
         *
         * @return estimated threshold to determine whether samples are inliers
         * or not.
         */
        public double getEstimatedThreshold() {
            return estimatedThreshold;
        }

        /**
         * Returns boolean indicating whether median residual computed in
         * current iteration has improved respect to previous iterations.
         *
         * @return true if median residual improved, false otherwise.
         */
        boolean isMedianResidualImproved() {
            return medianResidualImproved;
        }

        /**
         * Updates data contained in this instance.
         *
         * @param bestMedianResidual      best median of error found so far taking
         *                                into account all provided samples.
         * @param standardDeviation       standard deviation of error among all
         *                                provided samples respect to currently estimated result.
         * @param lmedsInlierModelEnabled indicates whether the LMedS or MSAC
         *                                inlier model is used.
         * @param residuals               residuals obtained for each sample of data.
         * @param numInliers              number of inliers found on current iteration.
         * @param medianResidual          median of error found on current iteration
         *                                among all provided samples.
         * @param estimatedThreshold      estimated threshold to determine whether
         *                                samples are inliers or not.
         * @param medianResidualImproved  indicates whether median residual
         *                                computed in current iteration has improved respect to previous
         *                                iteration.
         */
        protected void update(
                final double bestMedianResidual, final double standardDeviation, final boolean lmedsInlierModelEnabled,
                final double[] residuals, final int numInliers, final double medianResidual,
                final double estimatedThreshold, final boolean medianResidualImproved) {
            this.bestMedianResidual = bestMedianResidual;
            this.standardDeviation = standardDeviation;
            this.lmedsInlierModelEnabled = lmedsInlierModelEnabled;
            this.residuals = residuals;
            this.numInliers = numInliers;
            this.medianResidual = medianResidual;
            this.estimatedThreshold = estimatedThreshold;
            this.medianResidualImproved = medianResidualImproved;
        }
    }

}
