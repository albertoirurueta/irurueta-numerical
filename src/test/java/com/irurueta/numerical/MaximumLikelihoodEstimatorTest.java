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

class MaximumLikelihoodEstimatorTest {

    private static final int NUMBER_OF_SAMPLES = 1000000;

    private static final double MIN_MEAN = 1.0;
    private static final double MAX_MEAN = 10.0;

    private static final double MIN_STD = 1.0;
    private static final double MAX_STD = 5.0;

    private static final double MIN_GAUSSIAN_SIGMA = 0.5;
    private static final double MAX_GAUSSIAN_SIGMA = 2.0;

    private static final double RELATIVE_ERROR = 0.2;

    @Test
    void testCreate() throws NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final var standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        final var gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, MAX_GAUSSIAN_SIGMA);

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

        // instantiate with no parameters
        var estimator = MaximumLikelihoodEstimator.create();
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertThrows(NotAvailableException.class, estimator::getInputData);
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::estimate);

        // Instantiate with gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma);
        assertNotNull(estimator);

        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
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
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(0.0));

        // Instantiate with gaussian sigma and Histogram method
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma, MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
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
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(0.0,
                MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR));

        // Instantiate with gaussian sigma and Accurate method
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
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
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(0.0,
                MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR));

        // Instantiate with input data
        estimator = MaximumLikelihoodEstimator.create(inputData);
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Instantiate with input data and gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma);
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class,
                () -> MaximumLikelihoodEstimator.create(inputData, 0.0));

        // Instantiate with input data, gaussian sigma and HISTOGRAM method
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma,
                MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(inputData,
                0.0, MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR));

        // Instantiate with input data, gaussian sigma and ACCURATE method
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma);
        assertEquals(MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(inputData,
                0.0, MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR));

        // instantiate with min value, max value and input data
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue, inputData);
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(), 0.0);
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
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(maxValue2, minValue2,
                inputData));

        // instantiate with min value, max value, input data and gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue, inputData, gaussianSigma);
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertEquals(minValue, estimator.getMinValue(), 0.0);
        assertEquals(maxValue, estimator.getMaxValue(), 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(maxValue2, minValue2,
                inputData, gaussianSigma));
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(minValue2, maxValue2,
                inputData, 0.0));

        // instantiate with min value, max value, input data, gaussian sigma and
        // HISTOGRAM method
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue, inputData, gaussianSigma,
                MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertEquals(minValue, estimator.getMinValue(), 0.0);
        assertEquals(maxValue, estimator.getMaxValue(), 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(maxValue2, minValue2,
                inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR));
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(minValue2, maxValue2,
                inputData, 0.0, MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR));

        // instantiate with min value, max value, input data, gaussian sigma and
        // ACCURATE method
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue, inputData, gaussianSigma,
                MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
        assertEquals(gaussianSigma, estimator.getGaussianSigma(), 0.0);
        assertEquals(minValue, estimator.getMinValue(), 0.0);
        assertEquals(maxValue, estimator.getMaxValue(), 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(inputData, estimator.getInputData());
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(maxValue2, minValue2,
                inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR));
        assertThrows(IllegalArgumentException.class, () -> MaximumLikelihoodEstimator.create(minValue2, maxValue2,
                inputData, 0.0, MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR));
    }

    @Test
    void testGetMethod() {

        final var randomizer = new UniformRandomizer();
        final var gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, MAX_GAUSSIAN_SIGMA);

        var estimator = MaximumLikelihoodEstimator.create(gaussianSigma);
        assertEquals(MaximumLikelihoodEstimator.DEFAULT_METHOD, estimator.getMethod());

        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());

        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR, estimator.getMethod());
    }

    @Test
    void testGetSetGaussianSigma() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, MAX_GAUSSIAN_SIGMA);

        final var estimator = MaximumLikelihoodEstimator.create();

        assertEquals(MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, estimator.getGaussianSigma(), 0.0);

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

        final var estimator = MaximumLikelihoodEstimator.create();

        assertThrows(NotAvailableException.class, estimator::getMinValue);
        assertThrows(NotAvailableException.class, estimator::getMaxValue);
        assertFalse(estimator.areMinMaxValuesAvailable());

        // set min max values
        estimator.setMinMaxValues(minValue, maxValue);

        // check correctness
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());

        // Force IllegalArgumentException
        final var minValue2 = minValue;
        final var maxValue2 = maxValue;
        assertThrows(IllegalArgumentException.class, () -> estimator.setMinMaxValues(maxValue2, minValue2));
    }

    @Test
    void testIsLocked() {
        final var estimator = MaximumLikelihoodEstimator.create();
        assertFalse(estimator.isLocked());
    }

    @Test
    void testGetSetInputDataAndAvailability() throws LockedException, NotAvailableException {

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

        final var estimator = MaximumLikelihoodEstimator.create();

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

        final var estimator = MaximumLikelihoodEstimator.create();

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

        var estimator = MaximumLikelihoodEstimator.create(inputData);

        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());

        final var estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());

        assertEquals(mean, estimatedMean, RELATIVE_ERROR * estimatedMean);

        // Force NotReadyException
        estimator = MaximumLikelihoodEstimator.create();
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::estimate);
    }
}
