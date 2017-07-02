/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.HistogramMaximumLikelihoodEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 22, 2012
 */
package com.irurueta.numerical;

import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.*;
import org.junit.*;

public class HistogramMaximumLikelihoodEstimatorTest {
    
    public static final int MIN_BINS = 10;
    public static final int MAX_BINS = 100;
    
    public static final int NUMBER_OF_SAMPLES = 100000;
    
    public static final double MIN_MEAN = 1.0;
    public static final double MAX_MEAN = 10.0;
    
    public static final double MIN_STD = 1.0;
    public static final double MAX_STD = 5.0;
    
    public static final double MIN_GAUSSIAN_SIGMA = 0.5;
    public static final double MAX_GAUSSIAN_SIGMA = 2.0;
    
    public static final double RELATIVE_ERROR = 5.0;
    
    public HistogramMaximumLikelihoodEstimatorTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testConstructor() throws LockedException, NotAvailableException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA, 
                MAX_GAUSSIAN_SIGMA);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for(int i = 0; i < NUMBER_OF_SAMPLES; i++){
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        HistogramMaximumLikelihoodEstimator estimator;
        
        //test 1st constructor
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
        try{
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try{
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try{
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        
        
        //Test 2nd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(gaussianSigma,
                numberOfBins);
        assertNotNull(estimator);
        
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        try{
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        try{
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isInputDataAvailable());
        assertFalse(estimator.isReady());
        try{
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}

        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(gaussianSigma, 
                    1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(0.0, 
                    numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(estimator);
        
        
        //Test 3rd constructor
        estimator = new HistogramMaximumLikelihoodEstimator(inputData,
                gaussianSigma, numberOfBins);
        assertNotNull(estimator);
        
        assertEquals(estimator.getMethod(),
                MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        try{
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.areMinMaxValuesAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
        assertTrue(estimator.isReady());
        
        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(inputData,
                    gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(inputData, 0.0,
                    numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(inputData, 0.0, 
                    1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(estimator);
        
        
        
        //Test 4th constructor
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
        
        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, gaussianSigma, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, gaussianSigma, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, 0.0, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, 0.0, numberOfBins);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(minValue,
                    maxValue, inputData, 0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator = new HistogramMaximumLikelihoodEstimator(maxValue,
                    minValue, inputData, 0.0, 1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        
        assertNull(estimator);        
    }
    
    @Test
    public void testGetMethod(){
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        assertEquals(estimator.getMethod(), MaximumLikelihoodEstimatorMethod.
                HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR);
    }
    
    @Test
    public void testGetSetNumberOfBins() throws LockedException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numberOfBins = randomizer.nextInt(MIN_BINS, MAX_BINS);
        
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        assertEquals(estimator.getNumberOfBins(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_NUMBER_OF_BINS);
        
        //set new number of bins
        estimator.setNumberOfBins(numberOfBins);
        
        //check correctness
        assertEquals(estimator.getNumberOfBins(), numberOfBins);
        
        //Force IllegalArgumentException
        try{
            estimator.setNumberOfBins(1);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetGaussianSigma() throws LockedException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double gaussianSigma = randomizer.nextDouble(MIN_GAUSSIAN_SIGMA,
                MAX_GAUSSIAN_SIGMA);
        
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        assertEquals(estimator.getGaussianSigma(),
                HistogramMaximumLikelihoodEstimator.DEFAULT_GAUSSIAN_SIGMA, 
                0.0);
        
        //set new gaussian sigma
        estimator.setGaussianSigma(gaussianSigma);
        
        //check correctness
        assertEquals(estimator.getGaussianSigma(), gaussianSigma, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setGaussianSigma(0.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetMinMaxValuesAndAvailability() throws LockedException, NotAvailableException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for(int i = 0; i < NUMBER_OF_SAMPLES; i++){
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        HistogramMaximumLikelihoodEstimator estimator = 
                new HistogramMaximumLikelihoodEstimator();
        
        try{
            estimator.getMinValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getMaxValue();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.areMinMaxValuesAvailable());
        
        //set min max values
        estimator.setMinMaxValues(minValue, maxValue);
        
        //check correctness
        assertEquals(estimator.getMinValue(), minValue, 0.0);
        assertEquals(estimator.getMaxValue(), maxValue, 0.0);
        assertTrue(estimator.areMinMaxValuesAvailable());
        
        //Force IllegalArgumentException
        try{
            estimator.setMinMaxValues(maxValue, minValue);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testIsLocked(){
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        assertFalse(estimator.isLocked());
    }
    
    @Test
    public void testGetSetInputDataAndAvailability() throws LockedException, 
        NotAvailableException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for(int i = 0; i < NUMBER_OF_SAMPLES; i++){
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        try{
            estimator.getInputData();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isInputDataAvailable());
        
        //set input data
        estimator.setInputData(inputData);
        
        //check correctness
        assertEquals(estimator.getInputData(), inputData);
        assertTrue(estimator.isInputDataAvailable());
    }
    
    @Test
    public void testIsReady() throws LockedException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random());
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for(int i = 0; i < NUMBER_OF_SAMPLES; i++){
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        
        assertFalse(estimator.isReady());
        
        //set input data
        estimator.setInputData(inputData);
        
        //check correctness
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double mean = randomizer.nextDouble(MIN_MEAN, MAX_MEAN);
        double standardDeviation = randomizer.nextDouble(MIN_STD, MAX_STD);
        
        double[] inputData = new double[NUMBER_OF_SAMPLES];
        
        //initialize input data with gaussian data
        GaussianRandomizer gaussianRandomizer = new GaussianRandomizer(
                new Random(), mean, standardDeviation);
        double minValue = Double.MAX_VALUE;
        double maxValue = -Double.MAX_VALUE;
        for(int i = 0; i < NUMBER_OF_SAMPLES; i++){
            inputData[i] = gaussianRandomizer.nextDouble();
            if(inputData[i] < minValue) minValue = inputData[i];
            if(inputData[i] > maxValue) maxValue = inputData[i];
        }
        
        HistogramMaximumLikelihoodEstimator estimator =
                new HistogramMaximumLikelihoodEstimator();
        estimator.setInputData(inputData);
        
        assertTrue(estimator.isReady());
        assertFalse(estimator.isLocked());
        
        double estimatedMean = estimator.estimate();
        assertFalse(estimator.isLocked());
        
        double binSize = (maxValue - minValue) * estimator.getNumberOfBins();
        assertEquals(mean, estimatedMean, 2.0 * binSize);
        
        //Force NotReadyException
        estimator = new HistogramMaximumLikelihoodEstimator();
        assertFalse(estimator.isReady());
        try{
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
    }
}
