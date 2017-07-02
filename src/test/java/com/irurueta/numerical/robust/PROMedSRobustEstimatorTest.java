/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.PROMedSRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 7, 2015
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.robust.PROMedSRobustEstimator.PROMedSInliersData;
import com.irurueta.statistics.UniformRandomizer;
import java.util.List;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class PROMedSRobustEstimatorTest {

    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double THRESHOLD = 1e-6;
    
    //error added to samples and related to quality scores
    public static final double MIN_ERROR = 1e-5;
    public static final double MAX_ERROR = 1.0;
    
    //error added to quality scores so they are not totally related to sample
    //error
    public static final double MIN_SCORE_ERROR = -0.30;
    public static final double MAX_SCORE_ERROR = 0.30;
    
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
    
    public PROMedSRobustEstimatorTest() {}
    
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
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getConfidence(), 
                PROMedSRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(), 
                PROMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isStopThresholdEnabled(),
                PROMedSRobustEstimator.DEFAULT_STOP_THRESHOLD_ENABLED);
        assertEquals(estimator.getInlierFactor(),
                PROMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
        assertEquals(estimator.isUseInlierThresholds(),
                PROMedSRobustEstimator.DEFAULT_USE_INLIER_THRESHOLD);
        assertEquals(estimator.getNIters(), 
                PROMedSRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertEquals(estimator.getMaxOutliersProportion(),
                PROMedSRobustEstimator.DEFAULT_MAX_OUTLIERS_PROPORTION, 0.0);
        assertEquals(estimator.getEta0(),
                PROSACRobustEstimator.DEFAULT_ETA0, 0.0);
        assertEquals(estimator.getBeta(),
                PROSACRobustEstimator.DEFAULT_BETA, 0.0);
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
        
        //test constructor with listener
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);        
        TestPROMedSCRobustEstimatorListener listener =
                new TestPROMedSCRobustEstimatorListener(numSamples, 
                PERCENTAGE_OUTLIER, THRESHOLD);
        estimator = new PROMedSRobustEstimator<double[]>(listener);
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        assertTrue(estimator.isReady());
        assertEquals(estimator.getConfidence(), 
                PROSACRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(), 
                PROSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isStopThresholdEnabled(),
                PROMedSRobustEstimator.DEFAULT_STOP_THRESHOLD_ENABLED);        
        assertEquals(estimator.getInlierFactor(),
                PROMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
        assertEquals(estimator.isUseInlierThresholds(),
                PROMedSRobustEstimator.DEFAULT_USE_INLIER_THRESHOLD);        
        assertEquals(estimator.getNIters(), 
                PROSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertEquals(estimator.getMaxOutliersProportion(),
                PROSACRobustEstimator.DEFAULT_MAX_OUTLIERS_PROPORTION, 0.0);
        assertEquals(estimator.getEta0(),
                PROSACRobustEstimator.DEFAULT_ETA0, 0.0);
        assertEquals(estimator.getBeta(),
                PROSACRobustEstimator.DEFAULT_BETA, 0.0);   
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
    }
    
    @Test
    public void testGetSetListenerAvailabilityAndIsReady() 
            throws LockedException{
        PROMedSRobustEstimator<double[]> estimator = 
                new PROMedSRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());
        
        //set listener
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);        
        TestPROMedSCRobustEstimatorListener listener = 
                new TestPROMedSCRobustEstimatorListener(numSamples, 
                PERCENTAGE_OUTLIER, THRESHOLD);

        estimator.setListener(listener);
        
        //check correctness
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());        
    }
    
    @Test
    public void testGetSetProgressDelta() throws IllegalArgumentException, 
            LockedException{
        PROMedSRobustEstimator<double[]> estimator = 
                new PROMedSRobustEstimator<double[]>();
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
        PROMedSRobustEstimator<double[]> estimator = 
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.getConfidence(), 
                PROSACRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        
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
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.getMaxIterations(), 
                PROSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        
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
    public void testIsSetStopThresholdEnabled() throws LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.isStopThresholdEnabled(),
                PROMedSRobustEstimator.DEFAULT_STOP_THRESHOLD_ENABLED);
        
        //set new value
        estimator.setStopThresholdEnabled(
                !PROMedSRobustEstimator.DEFAULT_STOP_THRESHOLD_ENABLED);
        
        //check correctness
        assertEquals(estimator.isStopThresholdEnabled(),
                !PROMedSRobustEstimator.DEFAULT_STOP_THRESHOLD_ENABLED);
    }
        
    @Test
    public void testGetSetInlierFactor() throws IllegalArgumentException,
            LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.getInlierFactor(),
                PROMedSRobustEstimator.DEFAULT_INLIER_FACTOR, 0.0);
        
        //set new value 
        estimator.setInlierFactor(0.5);
        
        //check correctness
        assertEquals(estimator.getInlierFactor(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setInlierFactor(0.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testIsSetUseInlierThresholds() throws LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.isUseInlierThresholds(),
                PROMedSRobustEstimator.DEFAULT_USE_INLIER_THRESHOLD);
        
        //set new value
        estimator.setUseInlierThresholds(
                !PROMedSRobustEstimator.DEFAULT_USE_INLIER_THRESHOLD);
        
        //check correctness
        assertEquals(estimator.isUseInlierThresholds(),
                !PROMedSRobustEstimator.DEFAULT_USE_INLIER_THRESHOLD);
    }
    
    @Test
    public void testGetSetMaxOutliersProportion() 
            throws IllegalArgumentException, LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.getMaxOutliersProportion(),
                PROMedSRobustEstimator.DEFAULT_MAX_OUTLIERS_PROPORTION, 0.0);
        
        //set new value
        estimator.setMaxOutliersProportion(0.5);
        
        //check correctness
        assertEquals(estimator.getMaxOutliersProportion(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setMaxOutliersProportion(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator.setMaxOutliersProportion(2.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetEta0() throws IllegalArgumentException, 
            LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        assertEquals(estimator.getEta0(), PROMedSRobustEstimator.DEFAULT_ETA0, 
                0.0);
        
        //set new value
        estimator.setEta0(0.5);
        
        //check correctness
        assertEquals(estimator.getEta0(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setEta0(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator.setEta0(2.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetBeta() throws IllegalArgumentException,
            LockedException{
        PROMedSRobustEstimator<double[]> estimator =
                new PROMedSRobustEstimator<double[]>();
        
        assertEquals(estimator.getBeta(), PROMedSRobustEstimator.DEFAULT_BETA, 
                0.0);
        
        //set new value
        estimator.setBeta(0.5);
        
        //check correctness
        assertEquals(estimator.getBeta(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try{
            estimator.setBeta(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            estimator.setBeta(2.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException,
            RobustEstimatorException{
        for(int i = 0; i < TIMES; i++){
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            TestPROMedSCRobustEstimatorListener listener =
                    new TestPROMedSCRobustEstimatorListener(numSamples,
                    PERCENTAGE_OUTLIER, THRESHOLD);
            PROMedSRobustEstimator<double[]> estimator =
                    new PROMedSRobustEstimator<double[]>();            

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

            //check status after estimation
            assertFalse(estimator.isLocked());
            assertEquals(listener.getStartCounter(), 1);
            assertEquals(listener.getEndCounter(), 1);

            //check correctness of estimation
            assertEquals(params.length, listener.getParams().length);
            assertEquals(params.length, NUM_PARAMS);

            assertArrayEquals(params, listener.getParams(), ABSOLUTE_ERROR);
            
            PROMedSInliersData inliersData = estimator.getBestInliersData();
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
    
    private void computeSamplesAndQualityScores(double[] params, int numSamples, 
            int percentageOutliers, double[] ys, double[] xs, 
            double[] qualityScores){
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for(int i = 0; i < numSamples; i++){
            //compute x values
            xs[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            //compute exact y values
            ys[i] = params[0] * xs[i] + params[1];
            double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, 
                    MAX_SCORE_ERROR);          
            //inliers score can also have error
            qualityScores[i] = 1.0 + scoreError;
            if(randomizer.nextInt(0, 100) < percentageOutliers){
                //is outlier, so we add a certain amount of error
                double error = randomizer.nextDouble(MIN_ERROR, MAX_ERROR);
                ys[i] += error; //add sample error
                //quality score is (1 / (1 + error)) + scoreError
                qualityScores[i] = 1.0 / (1.0 + error) + scoreError; 
            }
        }        
    } 
    
    public class TestPROMedSCRobustEstimatorListener implements
            PROMedSRobustEstimatorListener<double[]>{
        
        private double[] params;
        private double[] xs;
        private double[] ys;
        private double[] qualityScores;
        private int numSamples;
        private double threshold;
        
        private int startCounter;
        private int endCounter;
        private float previousProgress;
        
        public TestPROMedSCRobustEstimatorListener(int numSamples,
                int percentageOutliers, double threshold){
            this.numSamples = numSamples;
            params = computeParams();
            xs = new double[numSamples];
            ys = new double[numSamples];
            qualityScores = new double[numSamples];
            computeSamplesAndQualityScores(params, numSamples, 
                    percentageOutliers, ys, xs, qualityScores);
            this.threshold = threshold;
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
        public double[] getQualityScores() {
            return qualityScores;
        }

        @Override
        public double getThreshold() {
            return threshold;
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
            testIsLocked((PROMedSRobustEstimator<double[]>)estimator);
            startCounter++;
        }

        @Override
        public void onEstimateEnd(RobustEstimator<double[]> estimator) {
            testIsLocked((PROMedSRobustEstimator<double[]>)estimator);
            endCounter++;
        }

        @Override
        public void onEstimateNextIteration(RobustEstimator<double[]> estimator, int iteration) {
            PROMedSRobustEstimator<double[]> ransacEstimator = 
                    (PROMedSRobustEstimator<double[]>)estimator;
            testIsLocked(ransacEstimator);
            assertTrue(iteration > 0);
            assertTrue(ransacEstimator.getNIters() >= 0);
            assertNotNull(ransacEstimator.getBestResult());
        }

        @Override
        public void onEstimateProgressChange(
                RobustEstimator<double[]> estimator, float progress) {
            testIsLocked((PROMedSRobustEstimator<double[]>)estimator);
            assertTrue(progress >= 0.0f);
            assertTrue(progress <= 1.0f);
            assertTrue(progress >= previousProgress);
            previousProgress = progress;
        } 
        
        private void testIsLocked(PROMedSRobustEstimator<double[]> estimator){
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
                estimator.setMaxOutliersProportion(0.5);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setEta0(0.5);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setBeta(0.5);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
            try{
                estimator.setStopThresholdEnabled(false);
                fail("LockedException expected but not thrown");
            }catch(LockedException e){}
        }
        
        public final void reset(){
            startCounter = endCounter = 0;
            previousProgress = 0.0f;            
        }        
    }
}
