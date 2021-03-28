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
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class HistogramMaximumLikelihoodEstimatorTest {

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
    public void testConstructor() throws LockedException, NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);
        final double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);

        final double[] inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        HistogramMaximumLikelihoodEstimator estimator;

        // test 1st constructor
        estimator = new HistogramMaximumLikelihoodEstimator();
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_NUMBER_OF_BINS);
        assertEquals(estimator.getGaussianSigma(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA,
                0.0);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // Test 2nd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(gaussianSigma,
                numberOfBins);
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(gaussianSigma,
                    1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(0.0,
                    numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // Test 3rd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(inputData,
                gaussianSigma, numberOfBins);
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(inputData,
                    gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(inputData, 0.0,
                    numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(inputData, 0.0,
                    1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // Test 4th constructor
        estimator = new HistogramMaximumLikelihoodEstimator(minValue, maxValue,
                inputData, gaussianSigma, numberOfBins);
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, gaussianSigma, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, 0.0, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, 0.0, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, 0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, 0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        assertNull(estimator);
    }

    @Test
    public void testGetMethod() {
        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
    }

    @Test
    public void testGetSetNumberOfBins() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);

        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        assertEquals(estimator.getNumberOfBins(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_NUMBER_OF_BINS);

        // set new number of bins
        estimator.setNumberOfBins(numberOfBins);

        // check correctness
        assertEquals(estimator.getNumberOfBins(), numberOfBins);

        // Force IllegalArgumentException
        try {
            estimator.setNumberOfBins(1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetGaussianSigma() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);

        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        assertEquals(estimator.getGaussianSigma(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA,
                0.0);

        // set new gaussian sigma
        estimator.setGaussianSigma(gaussianSigma);

        // check correctness
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setGaussianSigma(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetMinMaxValuesAndAvailability() throws LockedException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);

        final double[] inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.areMinMaxValuesAvailable());

        // set min max values
        estimator.setMinMaxValues(minValue, maxValue);

        // check correctness
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());

        // Force IllegalArgumentException
        try {
            estimator.setMinMaxValues(maxValue, minValue);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsLocked() {
        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        assertFalse(estimator.isLocked());
    }

    @Test
    public void testGetSetInputDataAndAvailability() throws LockedException,
            NotAvailableException {

        final double[] inputData = new double[NUMBER_OF_SAMPLES];

        final GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isInputDataAvailable());

        // set input data
        estimator.setInputData(inputData);

        // check correctness
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
    }

    @Test
    public void testIsReady() throws LockedException {

        final double[] inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        final HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();

        assertFalse(estimator.isReady());

        // set input data
        estimator.setInputData(inputData);

        // check correctness
        assertTrue(estimator.isReady());
    }

    @Test
    public void testEstimate() throws LockedException, NotReadyException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);

        final double[] inputData = new double[NUMBER_OF_SAMPLES];

        // initialize input data with gaussian data
        final GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if (inputData[i] < minValue) {
                minValue = inputData[i];
            }
            if (inputData[i] > maxValue) {
                maxValue = inputData[i];
            }
        }

        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        estimator.setInputData(inputData);

        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());

        final double estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());

        double binSize = (maxValue - minValue) * estimator.getNumberOfBins();
        assertEquals(mean, estimatedMean, 2.0 * binSize);

        // Force NotReadyException
        estimator = new HistogramMaximumLikelihoodEstimator();
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }
    }
}
