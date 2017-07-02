/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.LMedSRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 6, 2015
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.robust.LMedSRobustEstimator.LMedSInliersData;
import com.irurueta.statistics.UniformRandomizer;
import java.util.List;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class LMedSRobustEstimatorTest {
    
    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double MIN_ERROR = 1e-5;
    public static final double MAX_ERROR = 1.0;
    
    public static final double MIN_CONFIDENCE = 0.95;
    public static final double MAX_CONFIDENCE = 0.99;
    
    public static final int MIN_MAX_ITERATIONS = 500;
    public static final int MAX_MAX_ITERATIONS = 5000;
    
    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public static final double ABSOLUTE_ERROR = 1e-6;
    
    public static final int PERCENTAGE_OUTLIER = 20;
    
    public static final int NUM_PARAMS = 2;
    
    public static final int TIMES = 100;
    
    public LMedSRobustEstimatorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}

    @Test
    public void testConstructor(){
        //test empty constructor
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getConfidence(),
                LMedSRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(), 
                LMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.getNIters(),
                LMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
        assertEquals(estimator.getStopThreshold(),
                LMedSRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getInlierFactor(),
                LMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
        
        //test constructor with listener
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        TestLMedSRobustEstimatorListener listener =
                new TestLMedSRobustEstimatorListener(numSamples, 
                PERCENTAGE_OUTLIER);
        estimator = new LMedSRobustEstimator<double[]>(listener);
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertEquals(estimator.isReady(), listener.isReady());
        assertTrue(estimator.isReady());
        assertEquals(estimator.getConfidence(),
                LMedSRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(), 
                LMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.getNIters(),
                LMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
        assertEquals(estimator.getStopThreshold(),
                LMedSRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getInlierFactor(),
                LMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
    }
    
    @Test
    public void testGetSetListenerAvailabilityAndIsReady()
            throws LockedException{
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());
        
        //set listener
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        TestLMedSRobustEstimatorListener listener =
                new TestLMedSRobustEstimatorListener(numSamples, 
                PERCENTAGE_OUTLIER);
        
        estimator.setListener(listener);
        
        //check correctness
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testGetSetProgressDelta() throws IllegalArgumentException,
            LockedException{
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        float progressDelta = randomizer.nextFloat(0.0f, 1.0f);
        estimator.setProgressDelta(progressDelta);
        
        //check correctness
        assertEquals(estimator.getProgressDelta(), progressDelta, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setProgressDelta(-1.0f);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator.setProgressDelta(2.0f);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetConfidence() throws IllegalArgumentException, 
            LockedException{
        LMedSRobustEstimator<double[]> estimator = 
                new LMedSRobustEstimator<double[]>();
        assertEquals(estimator.getConfidence(), 
                LMedSRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double confidence = randomizer.nextDouble(0.0, 1.0);
        estimator.setConfidence(confidence);
        
        //check correctness
        assertEquals(estimator.getConfidence(), confidence, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setConfidence(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator.setConfidence(2.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetMaxIterations() throws IllegalArgumentException, 
            LockedException{
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertEquals(estimator.getMaxIterations(), 
                LMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int maxIterations = randomizer.nextInt(MIN_MAX_ITERATIONS, 
                MAX_MAX_ITERATIONS);
        estimator.setMaxIterations(maxIterations);
        
        //check correctness
        assertEquals(estimator.getMaxIterations(), maxIterations);
        
        //Force IllegalArgumentException
        try{
            estimator.setMaxIterations(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }    
    
    @Test
    public void testGetSetStopThreshold() throws IllegalArgumentException,
            LockedException{
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertEquals(estimator.getStopThreshold(),
                LMedSRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        
        //set new value
        estimator.setStopThreshold(2.0);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(), 2.0, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setStopThreshold(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetInlierFactor() throws IllegalArgumentException,
            LockedException{
        LMedSRobustEstimator<double[]> estimator =
                new LMedSRobustEstimator<double[]>();
        assertEquals(estimator.getInlierFactor(),
                LMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
        
        //set new value
        estimator.setInlierFactor(2.0);
        
        //check correctness
        assertEquals(estimator.getInlierFactor(), 2.0, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setInlierFactor(0.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException,
            RobustEstimatorException{
        for(int i = 0; i < TIMES; i++){
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            TestLMedSRobustEstimatorListener listener =
                    new TestLMedSRobustEstimatorListener(numSamples, 
                    PERCENTAGE_OUTLIER);
            LMedSRobustEstimator<double[]> estimator =
                    new LMedSRobustEstimator<double[]>();
            
            //Force NotReadyException
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}

            //set listener
            estimator.setListener(listener);
            listener.reset();
            assertEquals(listener.getStartCounter(), 0);
            assertEquals(listener.getEndCounter(), 0);
            assertFalse(estimator.isLocked());

            //estimate
            double[] params = estimator.estimate();

            assertNotNull(estimator.getBestResult());
            assertNotNull(estimator.getBestInliersData());

            //check status after estimation
            assertFalse(estimator.isLocked());
            assertEquals(listener.getStartCounter(), 1);
            assertEquals(listener.getEndCounter(), 1);

            //check correctness of estimation
            assertEquals(params.length, listener.getParams().length);
            assertEquals(params.length, NUM_PARAMS);

            assertArrayEquals(params, listener.getParams(), ABSOLUTE_ERROR);
            
            LMedSInliersData inliersData = estimator.getBestInliersData();
            assertNotNull(inliersData);
            assertSame(inliersData, estimator.getInliersData());
            assertTrue(inliersData.getNumInliers() > 0);
            assertNotNull(inliersData.getInliers());
            assertNotNull(inliersData.getResiduals());
        }
    }
    
    private double[] computeParams(){
        //we will estimate parameters a and b for equation y = a*x + b
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] params = new double[NUM_PARAMS];
        //a parameter
        params[0] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        //b parameter
        params[1] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        return params;
    }
    
    private void computeSamples(double[] params, int numSamples, 
            int percentageOutliers, double[] ys, double[] xs){
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for(int i = 0; i < numSamples; i++){
            //compute x values
            xs[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            //compute exact y values
            ys[i] = params[0] * xs[i] + params[1];
            if(randomizer.nextInt(0, 100) < percentageOutliers){
                //is outlier, so we add a certain amount of error
                double error = randomizer.nextDouble(MIN_ERROR, MAX_ERROR);
                ys[i] += error;
            }
        }        
    }    
    
    public class TestLMedSRobustEstimatorListener implements
            LMedSRobustEstimatorListener<double[]>{
        
        private double[] params;
        private double[] xs;
        private double[] ys;
        private int numSamples;
        
        private int startCounter;
        private int endCounter;
        private float previousProgress;
        
        public TestLMedSRobustEstimatorListener(int numSamples,
                int percentageOutliers){
            this.numSamples = numSamples;
            params = computeParams();
            xs = new double[numSamples];
            ys = new double[numSamples];
            computeSamples(params, numSamples, percentageOutliers, ys, xs);
            reset();
        }
        
        public double[] getParams(){
            return params;
        }
        
        public double[] getXs(){
            return xs;
        }
        
        public double[] getYs(){
            return ys;
        }
        
        public int getStartCounter(){
            return startCounter;
        }
        
        public int getEndCounter(){
            return endCounter;
        }

        @Override
        public int getTotalSamples() {
            return numSamples;
        }

        @Override
        public int getSubsetSize() {
            //only two matches x and y are required to determine a and b as 
            //follows:
            //y1 = a* x1 + b --> b = y1 - a * x1
            //y2 = a * x2 + b --> y2 = a * x2 + y1 - a * x1 -->
            //y2 - y1 = (x2 - x1)*a --> a = (y2 - y1) / (x2 - x1)
            
            //Hence:
            //a = (y2 - y1) / (x2 - x1)
            //b = y1 - (y2 - y1) / (x2 - x1) * x1
            return NUM_PARAMS; 
        }

        @Override
        public void estimatePreliminarSolutions(int[] samplesIndices, 
                List<double[]> solutions){
            
            if(samplesIndices.length != NUM_PARAMS) 
                throw new IllegalArgumentException();
            int index1 = samplesIndices[0];
            int index2 = samplesIndices[1];
            
            double y1 = ys[index1];
            double y2 = ys[index2];
            double x1 = xs[index1];
            double x2 = xs[index2];
            
            double a = (y2 - y1) / (x2 - x1);
            double b = y1 - a * x1;
            
            double[] solution = new double[NUM_PARAMS];
            solution[0] = a;
            solution[1] = b;
            
            solutions.add(solution);
        }

        @Override
        public double computeResidual(double[] currentEstimation, int i) {
            double a = currentEstimation[0];
            double b = currentEstimation[1];

            double estimatedY = a * xs[i] + b;
            double y = ys[i];
            
            return  Math.abs(estimatedY - y);
        }

        @Override
        public boolean isReady() {
            return params != null && xs != null && ys != null;
        }

        @Override
        public void onEstimateStart(RobustEstimator<double[]> estimator) {
            testIsLocked((LMedSRobustEstimator<double[]>)estimator);
            startCounter++;
        }

        @Override
        public void onEstimateEnd(RobustEstimator<double[]> estimator) {
            testIsLocked((LMedSRobustEstimator<double[]>)estimator);
            endCounter++;
        }

        @Override
        public void onEstimateNextIteration(RobustEstimator<double[]> estimator, int iteration) {
            LMedSRobustEstimator<double[]> lmedsEstimator = 
                    (LMedSRobustEstimator<double[]>)estimator;
            testIsLocked(lmedsEstimator);
            assertTrue(iteration > 0);
            assertTrue(lmedsEstimator.getNIters() >= 0);
            assertNotNull(lmedsEstimator.getBestResult());
            assertNotNull(lmedsEstimator.getBestInliersData());
        }

        @Override
        public void onEstimateProgressChange(RobustEstimator<double[]> estimator, float progress) {
            testIsLocked((LMedSRobustEstimator<double[]>)estimator);
            assertTrue(progress >= 0.0f);
            assertTrue(progress <= 1.0f);
            assertTrue(progress >= previousProgress);
            previousProgress = progress;
        }
        
        private void testIsLocked(LMedSRobustEstimator<double[]> estimator){
            assertTrue(estimator.isLocked());
            //test that estimator cannot be modified while locked
            try{
                estimator.setConfidence(0.5);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setListener(this);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setMaxIterations(1);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setProgressDelta(0.5f);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setStopThreshold(1.0);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setInlierFactor(1.0);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
        }
        
        public final void reset(){
            startCounter = endCounter = 0;
            previousProgress = 0.0f;            
        }        
    }
}
