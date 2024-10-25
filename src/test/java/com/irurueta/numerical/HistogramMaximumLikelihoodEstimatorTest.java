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

import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HistogramMaximumLikelihoodEstimatorTest {

    private static final int MIN_BINS = 10;
    private static final int MAX_BINS = 100;

    private static final int NUMBER_OF_SAMPLES = 100000;

    private static final double MIN_MEAN = 1.0;
    private static final double MAX_MEAN = 10.0;

    private static final double MIN_STD = 1.0;
    private static final double MAX_STD = 5.0;

    private static final double MIN_GAUSSIAN_SIGMA = 0.5;
    private static final double MAX_GAUSSIAN_SIGMA = 2.0;

    @Test
    void testConstructor() throws NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);
        final var gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, MAX_GAUSSIAN_SIGMA);

        final var inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final var gaussianRandomizer = new GaussianRandomizer();
        var minValue = Double.MAX_VALUE;
        var maxValue = -Double.MAX_VALUE;
        for (var i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        // test 1st constructor
        var estimator = new HistogramMaximumLikelihoodEstimator();
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(HistogramMaximumLikelihoodEstimator.DEFAULT_NUMBER_OF_BINS, estimator.getNumberOfBins());
        assertEquals(HistogramMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(),
                0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertThrows(NotAvailableException.class, estimator::getInputData);
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::estimate);

        // Test 2nd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(gaussianSigma, numberOfBins);
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(numberOfBins, estimator.getNumberOfBins());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertThrows(NotAvailableException.class, estimator::getInputData);
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::estimate);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(gaussianSigma,
                1));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(0.0,
                numberOfBins));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(0.0,
                1));

        // Test 3rd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(inputData, gaussianSigma, numberOfBins);
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(numberOfBins, estimator.getNumberOfBins());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(inputData,
                gaussianSigma, 1));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(inputData,
                0.0, numberOfBins));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(inputData,
                0.0, 1));

        // Test 4th constructor
        estimator = new HistogramMaximumLikelihoodEstimator(minValue, maxValue, inputData, gaussianSigma, numberOfBins);
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(numberOfBins, estimator.getNumberOfBins());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertEquals(minValue, estimator.getMinValue(), 0.0);
        assertEquals(maxValue, estimator.getMaxValue(), 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        final var minValue2 = minValue;
        final var maxValue2 = maxValue;
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(maxValue2,
                minValue2, inputData, gaussianSigma, numberOfBins));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(minValue2,
                maxValue2, inputData, gaussianSigma, 1));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(maxValue2,
                minValue2, inputData, gaussianSigma, 1));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(minValue2,
                maxValue2, inputData, 0.0, numberOfBins));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(maxValue2,
                minValue2, inputData, 0.0, numberOfBins));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(minValue2,
                maxValue2, inputData, 0.0, 1));
        assertThrows(IllegalArgumentException.class, () -> new HistogramMaximumLikelihoodEstimator(maxValue2,
                minValue2, inputData, 0.0, 1));
    }

    @Test
    void testGetMethod() {
        final var estimator = new HistogramMaximumLikelihoodEstimator();
        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
    }

    @Test
    void testGetSetNumberOfBins() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);

        final var estimator = new HistogramMaximumLikelihoodEstimator();

        assertEquals(HistogramMaximumLikelihoodEstimator.DEFAULT_NUMBER_OF_BINS, estimator.getNumberOfBins());

        // set new number of bins
        estimator.setNumberOfBins(numberOfBins);

        // check correctness
        assertEquals(numberOfBins, estimator.getNumberOfBins());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setNumberOfBins(1));
    }

    @Test
    void testGetSetGaussianSigma() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, MAX_GAUSSIAN_SIGMA);

        final var estimator = new HistogramMaximumLikelihoodEstimator();

        assertEquals(HistogramMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(),
                0.0);

        // set new gaussian sigma
        estimator.setGaussianSigma(gaussianSigma);

        // check correctness
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setGaussianSigma(0.0));
    }

    @Test
    void testGetSetMinMaxValuesAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final var standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);

        final var inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final var gaussianRandomizer = new GaussianRandomizer(mean, standardDeviation);
        var minValue = Double.MAX_VALUE;
        var maxValue = -Double.MAX_VALUE;
        for (var i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final var estimator = new HistogramMaximumLikelihoodEstimator();

        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());

        // set min max values
        estimator.setMinMaxValues(minValue, maxValue);

        // check correctness
        assertEquals(minValue, estimator.getMinValue(), 0.0);
        assertEquals(maxValue, estimator.getMaxValue(), 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());

        // Force IllegalArgumentException
        final var minValue2 = minValue;
        final var maxValue2 = maxValue;
        assertThrows(IllegalArgumentException.class, () -> estimator.setMinMaxValues(maxValue2, minValue2));
    }

    @Test
    void testIsLocked() {
        final var estimator = new HistogramMaximumLikelihoodEstimator();
        assertFalse(estimator.isLocked());
    }

    @Test
    void testGetSetInputDataAndAvailability() throws LockedException, NotAvailableException {

        final var inputData = new double[NUMBER_OF_SAMPLES];

        final var gaussianRandomizer = new GaussianRandomizer();
        var minValue = Double.MAX_VALUE;
        var maxValue = -Double.MAX_VALUE;
        for (var i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final var estimator = new HistogramMaximumLikelihoodEstimator();

        assertThrows(NotAvailableException.class, estimator::getInputData);
        assertFalse(estimator.isInputDataAvailable());

        // set input data
        estimator.setInputData(inputData);

        // check correctness
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
    }

    @Test
    void testIsReady() throws LockedException {

        final var inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final var gaussianRandomizer = new GaussianRandomizer();
        var minValue = Double.MAX_VALUE;
        var maxValue = -Double.MAX_VALUE;
        for (var i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final var estimator = new HistogramMaximumLikelihoodEstimator();

        assertFalse(estimator.isReady());

        // set input data
        estimator.setInputData(inputData);

        // check correctness
        assertTrue(estimator.isReady());
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException {

        final var randomizer = new UniformRandomizer();
        final var mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final var standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);

        final var inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final var gaussianRandomizer = new GaussianRandomizer(mean, standardDeviation);
        var minValue = Double.MAX_VALUE;
        var maxValue = -Double.MAX_VALUE;
        for (var i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        var estimator = new HistogramMaximumLikelihoodEstimator();
        estimator.setInputData(inputData);

        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());

        final var estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());

        var binSize = (maxValue - minValue) * estimator.getNumberOfBins();
        assertEquals(mean, estimatedMean, 2.0 * binSize);

        // Force NotReadyException
        estimator = new HistogramMaximumLikelihoodEstimator();
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::estimate);
    }
}
