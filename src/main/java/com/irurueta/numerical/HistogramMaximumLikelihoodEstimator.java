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

import com.irurueta.algebra.ArrayUtils;

import java.util.Arrays;

/**
 * Class to estimate the most likely value from a series of samples assumed to
 * be normally distributed.
 * This implementation will compute a histogram of the probability distribution
 * function of all the provided samples.
 * The probability distribution function is computed by aggregating all the
 * samples as a series of Gaussian functions with a small sigma and centered at
 * the exact value of the sample.
 */
public class HistogramMaximumLikelihoodEstimator extends MaximumLikelihoodEstimator {

    /**
     * Default number of bins to be used on the histogram.
     */
    public static final int DEFAULT_NUMBER_OF_BINS = 100;

    /**
     * Minimum number of bins allowed on the histogram.
     */
    public static final int MIN_NUMBER_OF_BINS = 2;

    /**
     * Value to be considered as the machine precision.
     */
    public static final double EPS = 1e-9;

    /**
     * Number of bins to be used on the histogram. The larger the value the
     * more precise results will be but more data will be needed to obtain
     * statistically meaningful data so that enough samples are present on each
     * bin.
     */
    private int numberOfBins;

    /**
     * Number of samples contained on each bin. The number of samples is not
     * an integer because the Gaussians of each sample will only have a value
     * equal to one on their respective centers, but at other locations, only
     * the fractional part is added to each bin. Using Gaussians on each sample
     * behave as an interpolation of the location of each sample around a
     * certain value.
     */
    private double[] bins;

    /**
     * Array containing the Gaussian values corresponding to each bin for a
     * single sample. The values of each array are updated for each sample and
     * aggregated to the values of the bins array.
     */
    private double[] gaussian;

    /**
     * Empty constructor.
     */
    public HistogramMaximumLikelihoodEstimator() {
        super();
        numberOfBins = DEFAULT_NUMBER_OF_BINS;
        bins = null;
        gaussian = null;
    }

    /**
     * Constructor.
     *
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @param numberOfBins  Number of bins to be used on the histogram.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero, or if provided number of bins is smaller than the
     *                                  minimum allowed number of bins.
     */
    public HistogramMaximumLikelihoodEstimator(final double gaussianSigma, final int numberOfBins) {
        super(gaussianSigma);
        internalSetNumberOfBins(numberOfBins);
        bins = null;
        gaussian = null;
    }

    /**
     * Constructor.
     *
     * @param inputData     Array containing input data where most likely value must
     *                      be estimated from.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @param numberOfBins  Number of bins to be used on the histogram.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero, or if provided number of bins is smaller than the
     *                                  minimum allowed number of bins.
     */
    public HistogramMaximumLikelihoodEstimator(final double[] inputData, final double gaussianSigma,
                                               final int numberOfBins) {
        super(inputData, gaussianSigma);
        internalSetNumberOfBins(numberOfBins);
        bins = null;
        gaussian = null;
    }

    /**
     * Constructor.
     *
     * @param minValue      Minimum value assumed to be contained within input data
     *                      array.
     * @param maxValue      Maximum value assumed to be contained within input data
     *                      array.
     * @param inputData     Array containing input data where most likely value must
     *                      be estimated from.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @param numberOfBins  Number of bins to be used on the histogram.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     *                                  negative or zero, or if provided number of bins is smaller than the
     *                                  minimum allowed number of bins or if minValue &lt; maxValue.
     */
    public HistogramMaximumLikelihoodEstimator(
            final double minValue, final double maxValue, final double[] inputData, final double gaussianSigma,
            final int numberOfBins) {
        super(minValue, maxValue, inputData, gaussianSigma);
        internalSetNumberOfBins(numberOfBins);
        bins = null;
        gaussian = null;
    }


    /**
     * Returns method to be used for maximum likelihood estimation, which for
     * this class is MaximumLikelihoodEstimatorMethod.
     * HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR.
     *
     * @return Method for maximum likelihood estimation.
     */
    @Override
    public MaximumLikelihoodEstimatorMethod getMethod() {
        return MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR;
    }

    /**
     * Returns number of bins to be used on the histogram. The larger the value
     * the more precise results will be but more data will be needed to obtain
     * statistically meaningful data so that enough samples are present on each
     * bin.
     *
     * @return Number of bins to be used on the histogram.
     */
    public int getNumberOfBins() {
        return numberOfBins;
    }

    /**
     * Sets number of bins to be used on the histogram. The larger the provided
     * value being set the more precise results will be but more data will be
     * needed to obtain statistically meaningful data so that enough samples are
     * present on each bin.
     *
     * @param numberOfBins Number of bins to be used on the histogram.
     * @throws LockedException          Exception raised if this instance is locked.
     *                                  This method can only be executed when computations finish and this
     *                                  instance becomes unlocked.
     * @throws IllegalArgumentException Raised if provided value is lower than
     *                                  the allowed minimum.
     */
    public void setNumberOfBins(final int numberOfBins) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetNumberOfBins(numberOfBins);
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

        if (!areMinMaxAvailable) {
            computeMinMaxValues();
        }

        if ((maxValue - minValue) < EPS) {
            // min-max limits are almost equal, so we return it as the solution
            locked = false;
            return (minValue + maxValue) * 0.5;
        }

        // Create histogram (initialized to 0.0)
        if (bins == null || bins.length != numberOfBins) {
            bins = new double[numberOfBins];
        }
        Arrays.fill(bins, 0.0);


        // Vector containing gaussians being added to histogram
        // Data is added to histogram as if it was convoluted with a gaussian
        // to add smoothness to results
        if (gaussian == null || gaussian.length != numberOfBins) {
            gaussian = new double[numberOfBins];
        }

        // set of data values between bins
        final var delta = (maxValue - minValue) / (numberOfBins - 1);

        // iterate over input data to add data to histogram
        double gaussianCenterPos;
        for (var data : inputData) {
            gaussianCenterPos = (data - minValue) / delta;
            computeGaussian(gaussian, gaussianCenterPos);

            // add gaussian centered at input data value to histogram bins
            ArrayUtils.sum(bins, gaussian, bins);
        }

        // find location of maximum in histogram
        var maxBin = 0.0;
        final double maxFuncValue;
        int maxPos = 0;
        for (var i = 0; i < numberOfBins; i++) {
            if (bins[i] > maxBin) {
                maxBin = bins[i];
                maxPos = i;
            }
        }

        maxFuncValue = minValue + maxPos * delta;

        locked = false;
        return maxFuncValue;
    }

    /**
     * Internal method to compute values of Gaussian vector assumed to be
     * centered at provided value. Gaussian values are stored in provided array
     * so that array can be reused.
     *
     * @param gaussian  Array containing Gaussian values
     * @param centerPos Value where Gaussian is centered
     */
    protected void computeGaussian(final double[] gaussian, final double centerPos) {
        final var length = gaussian.length;

        double x;
        for (var i = 0; i < length; i++) {
            x = i - centerPos;
            gaussian[i] = Math.exp(-x * x / (2.0 * gaussianSigma * gaussianSigma))
                    / (Math.sqrt(2.0 * Math.PI) * gaussianSigma);
        }
    }

    /**
     * Internal method to set number of bins to be used on the histogram. The
     * larger the value being set the more precise results will be but more data
     * will be needed to obtain statistically meaningful data so that enough
     * samples are present on each bin.
     *
     * @param numberOfBins Number of bins to be used on the histogram.
     * @throws IllegalArgumentException Raised if provided value is lower than
     *                                  the allowed minimum.
     */
    private void internalSetNumberOfBins(final int numberOfBins) {
        if (numberOfBins < MIN_NUMBER_OF_BINS) {
            throw new IllegalArgumentException();
        }

        this.numberOfBins = numberOfBins;
    }
}
