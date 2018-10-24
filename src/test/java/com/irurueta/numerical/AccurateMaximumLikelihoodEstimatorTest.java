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
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.*;

@SuppressWarnings("Duplicates")
public class AccurateMaximumLikelihoodEstimatorTest {
    
    private static final int NUMBER_OF_SAMPLES = 100000;
    
    private static final double MIN_MEAN = 1.0;
    private static final double MAX_MEAN = 10.0;
    
    private static final double MIN_STD = 1.0;
    private static final double MAX_STD = 5.0;
    
    private static final double MIN_GAUSSIAN_SIGMA = 0.5;
    private static final double MAX_GAUSSIAN_SIGMA = 2.0;
    
    private static final double RELATIVE_ERROR = 0.2;
    
    public AccurateMaximumLikelihoodEstimatorTest() { }

    @BeforeClass
    public static void setUpClass() { }

    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() throws NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);        
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        boolean useHistogramInitialSolution = randomizer.nextBoolean();
        
        AccurateMaximumLikelihoodEstimator estimator;
        
        //instantiate with empty constructor
        estimator = new AccurateMaximumLikelihoodEstimator();
        assertNotNull(estimator);
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getGaussianSigma(), 
                AccurateMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 0.0);
        assertEquals(estimator.isHistogramInitialSolutionUsed(),
                AccurateMaximumLikelihoodEstimator.
                DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (LockedException e) {
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        
        //Instantiate with gaussian sigma and use histogram initial solution
        estimator = new AccurateMaximumLikelihoodEstimator(gaussianSigma, 
                useHistogramInitialSolution);
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        assertEquals(estimator.isHistogramInitialSolutionUsed(), 
                useHistogramInitialSolution);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (LockedException ignore) {
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }

        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new AccurateMaximumLikelihoodEstimator(0.0, 
                    useHistogramInitialSolution);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
                
        
        //Instantiate with input data, gaussian sigma and use histogram initial
        //solution
        estimator = new AccurateMaximumLikelihoodEstimator(inputData, 
                gaussianSigma, useHistogramInitialSolution);
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        assertEquals(estimator.isHistogramInitialSolutionUsed(), 
                useHistogramInitialSolution);
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());        
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new AccurateMaximumLikelihoodEstimator(inputData, 0.0, 
                    useHistogramInitialSolution);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        

        //Instantiate with min, max values, input data, gaussian sigma and 
        //use histogram initial solution
        estimator = new AccurateMaximumLikelihoodEstimator(minValue, maxValue,
                inputData, gaussianSigma, useHistogramInitialSolution);
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        assertEquals(estimator.isHistogramInitialSolutionUsed(), 
                useHistogramInitialSolution);
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());        
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());

        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new AccurateMaximumLikelihoodEstimator(maxValue, 
                    minValue, inputData, gaussianSigma, 
                    useHistogramInitialSolution);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
        estimator = null;
        try {
            estimator = new AccurateMaximumLikelihoodEstimator(minValue, 
                    maxValue, inputData, 0.0, useHistogramInitialSolution);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
    }
    
    @Test
    public void testGetMethod() {
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR);
    }
    
    @Test
    public void testGetSetGaussianSigma() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);
        
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        
        assertEquals(estimator.getGaussianSigma(),
                AccurateMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 0.0);
        
        //set new gaussian sigma
        estimator.setGaussianSigma(gaussianSigma);
        
        //check correctness
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setGaussianSigma(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetMinMaxValuesAndAvailability() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);        
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        
        try {
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.areMinMaxValuesAvailable());
        
        //set min max values
        estimator.setMinMaxValues(minValue, maxValue);
        
        //check correctness
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        
        //Force IllegalArgumentException
        try {
            estimator.setMinMaxValues(maxValue, minValue);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testIsLocked() {
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        
        assertFalse(estimator.isLocked());
    }
    
    @Test
    public void testGetSetInputDataAndAvailability() throws LockedException, NotAvailableException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        
        try {
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isInputDataAvailable());
        
        //set input data
        estimator.setInputData(inputData);
        
        //check correctness
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
    }
    
    @Test
    public void testIsReady() throws LockedException {
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        AccurateMaximumLikelihoodEstimator estimator =
                new AccurateMaximumLikelihoodEstimator();
        
        assertFalse(estimator.isReady());
        
        estimator.setInputData(inputData);
        
        //check correctness
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testGetSetHistogramInitialSolutionUsed() throws LockedException {
        AccurateMaximumLikelihoodEstimator estimator =
                new AccurateMaximumLikelihoodEstimator();
        
        assertEquals(estimator.isHistogramInitialSolutionUsed(),
                AccurateMaximumLikelihoodEstimator.
                DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION);
        
        //disable
        estimator.setHistogramInitialSolutionUsed(false);
        
        //check correctness
        assertFalse(estimator.isHistogramInitialSolutionUsed());
        
        //enable
        estimator.setHistogramInitialSolutionUsed(true);
        
        //check correctness
        assertTrue(estimator.isHistogramInitialSolutionUsed());
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        AccurateMaximumLikelihoodEstimator estimator = 
                new AccurateMaximumLikelihoodEstimator();
        estimator.setInputData(inputData);
        
        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());
        
        double estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());
        
        assertEquals(mean, estimatedMean, RELATIVE_ERROR * estimatedMean);
        
        //Force NotReadyException
        estimator = new AccurateMaximumLikelihoodEstimator();
        assertFalse(estimator.isReady());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
    }
}
