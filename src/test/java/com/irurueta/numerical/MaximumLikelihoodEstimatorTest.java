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

public class MaximumLikelihoodEstimatorTest {

    private static final int NUMBER_OF_SAMPLES = 1000000;

    private static final double MIN_MEAN = 1.0;
    private static final double MAX_MEAN = 10.0;

    private static final double MIN_STD = 1.0;
    private static final double MAX_STD = 5.0;

    private static final double MIN_GAUSSIAN_SIGMA = 0.5;
    private static final double MAX_GAUSSIAN_SIGMA = 2.0;

    private static final double RELATIVE_ERROR = 0.2;

    @Test
    public void testCreate() throws LockedException, NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        final double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        final double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);

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

        MaximumLikelihoodEstimator estimator;

        // instantiate with no parameters
        estimator = MaximumLikelihoodEstimator.create();
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
        assertEquals(estimator.getGaussianSigma(),
                MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 0.0);
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


        // Instantiate with gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma);
        assertNotNull(estimator);

        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
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
        try {
            MaximumLikelihoodEstimator.create(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }


        // Instantiate with gaussian sigma and Histogram method
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
        try {
            MaximumLikelihoodEstimator.create(0.0,
                    MaximumLikelihoodEstimatorMethod.
                            HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }


        // Instantiate with gaussian sigma and Accurate method
        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.
                        ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
        try {
            MaximumLikelihoodEstimator.create(0.0,
                    MaximumLikelihoodEstimatorMethod.
                            ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }


        // Instantiate with input data
        estimator = MaximumLikelihoodEstimator.create(inputData);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
        assertEquals(estimator.getGaussianSigma(),
                MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 0.0);
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


        // Instantiate with input data and gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
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
            estimator = MaximumLikelihoodEstimator.create(inputData, 0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // Instantiate with input data, gaussian sigma and HISTOGRAM method
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma,
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
            estimator = MaximumLikelihoodEstimator.create(inputData, 0.0,
                    MaximumLikelihoodEstimatorMethod.
                            HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // Instantiate with input data, gaussian sigma and ACCURATE method
        estimator = MaximumLikelihoodEstimator.create(inputData, gaussianSigma);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
            estimator = MaximumLikelihoodEstimator.create(inputData, 0.0,
                    MaximumLikelihoodEstimatorMethod.
                            ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // instantiate with min value, max value and input data
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                inputData);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
        assertEquals(estimator.getGaussianSigma(),
                MaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 0.0);
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
            estimator = MaximumLikelihoodEstimator.create(maxValue, minValue,
                    inputData);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // instantiate with min value, max value, input data and gaussian sigma
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                inputData, gaussianSigma);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);
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
            estimator = MaximumLikelihoodEstimator.create(maxValue, minValue,
                    inputData, gaussianSigma);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                    inputData, 0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // instantiate with min value, max value, input data, gaussian sigma and
        // HISTOGRAM method
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
            estimator = MaximumLikelihoodEstimator.create(maxValue, minValue,
                    inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.
                            HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                    inputData, 0.0, MaximumLikelihoodEstimatorMethod.
                            HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // instantiate with min value, max value, input data, gaussian sigma and
        // ACCURATE method
        estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.
                        ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                        ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
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
            estimator = MaximumLikelihoodEstimator.create(maxValue, minValue,
                    inputData, gaussianSigma, MaximumLikelihoodEstimatorMethod.
                            ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = MaximumLikelihoodEstimator.create(minValue, maxValue,
                    inputData, 0.0, MaximumLikelihoodEstimatorMethod.
                            ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testGetMethod() {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);

        MaximumLikelihoodEstimator estimator;

        estimator = MaximumLikelihoodEstimator.create(gaussianSigma);
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimator.DEFAULT_METHOD);

        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.
                        HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);

        estimator = MaximumLikelihoodEstimator.create(gaussianSigma,
                MaximumLikelihoodEstimatorMethod.
                        ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
    }

    @Test
    public void testGetSetGaussianSigma() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);

        final MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create();

        assertEquals(estimator.getGaussianSigma(), MaximumLikelihoodEstimator.
                DEFAULT_GAUSSIAN_SIGMA, 0.0);

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

        final MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create();

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

        final MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create();

        assertFalse(estimator.isLocked());
    }

    @Test
    public void testGetSetInputDataAndAvailability()
            throws LockedException, NotAvailableException {

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

        final MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create();

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

        final MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create();

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

        MaximumLikelihoodEstimator estimator =
                MaximumLikelihoodEstimator.create(inputData);

        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());

        final double estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());

        assertEquals(mean, estimatedMean, RELATIVE_ERROR * estimatedMean);


        // Force NotReadyException
        estimator = MaximumLikelihoodEstimator.create();
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }
    }
}
