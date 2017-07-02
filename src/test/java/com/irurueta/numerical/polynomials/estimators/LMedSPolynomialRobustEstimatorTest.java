/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.LMedSPolynomialRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 14, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.numerical.robust.RobustEstimatorException;
import com.irurueta.numerical.robust.RobustEstimatorMethod;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class LMedSPolynomialRobustEstimatorTest implements 
        PolynomialRobustEstimatorListener{
    
    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public static final int MIN_DEGREE = 1;
    public static final int MAX_DEGREE = 5;
    
    public static final double ABSOLUTE_ERROR = 1e-8;
    
    public static final double PERCENTAGE_OUTLIER = 20;
    
    public static final int MIN_EVALUATIONS = 500;
    public static final int MAX_EVALUATIONS = 1000;
    
    public static final double STD_ERROR = 100.0;
    
    public static final int TIMES = 10;
    
    private int estimateStart;
    private int estimateEnd;
    private int estimateNextIteration;
    private int estimateProgressChange;
    
    public LMedSPolynomialRobustEstimatorTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() {
        //test empty constructor
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        
        //test constructor with degree
        estimator = new LMedSPolynomialRobustEstimator(2);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());    
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>(); 
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new LMedSPolynomialRobustEstimator(evaluations);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());     
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        //test constructor with listener
        estimator = new LMedSPolynomialRobustEstimator(this);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        
        
        //test constructor with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new LMedSPolynomialRobustEstimator(2, evaluations);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores()); 
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new LMedSPolynomialRobustEstimator(2, wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with degree and listener
        estimator = new LMedSPolynomialRobustEstimator(2, this);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());    
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with evaluations and listener
        estimator = new LMedSPolynomialRobustEstimator(evaluations, this);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());     
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(wrongEvals, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with degree, evaluations and listener
        estimator = new LMedSPolynomialRobustEstimator(2, evaluations, this);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(), 
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());    
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMedSPolynomialRobustEstimator(0, evaluations, 
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new LMedSPolynomialRobustEstimator(2, wrongEvals, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
    }
    
    @Test
    public void testGetSetStopThreshold() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.getStopThreshold(),
                LMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD, 0.0);
        
        //set new value
        estimator.setStopthreshold(1.0);
        
        //check correctness
        assertEquals(estimator.getStopThreshold(), 1.0, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setStopthreshold(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetEvaluations() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertNull(estimator.getEvaluations());
        
        //set new value
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>(); 
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator.setEvaluations(evaluations);
        
        //check correctness
        assertSame(estimator.getEvaluations(), evaluations);
        
        //Force IllegalArgumentException
        try {
            estimator.setEvaluations(null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setEvaluations(new ArrayList<PolynomialEvaluation>());
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }        
    }
    
    @Test
    public void testGetSetListener() {
        LMedSPolynomialRobustEstimator estimator = 
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertNull(estimator.getListener());
        
        //set new value
        estimator.setListener(this);
        
        //check correctness
        assertSame(estimator.getListener(), this);
    }    
    
    @Test
    public void testGetSetProgressDelta() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.getProgressDelta(), 
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        
        //set new value
        estimator.setProgressDelta(0.5f);
        
        //check correctness
        assertEquals(estimator.getProgressDelta(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setProgressDelta(-1.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setProgressDelta(2.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }   
    
    @Test
    public void testGetSetConfidence() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        
        //set new value
        estimator.setConfidence(0.5);
        
        //check correctness
        assertEquals(estimator.getConfidence(), 0.5, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setConfidence(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setConfidence(2.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }    
    
    @Test
    public void testGetSetMaxIterations() throws LockedException {
        LMedSPolynomialRobustEstimator estimator = 
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        
        //set new value
        estimator.setMaxIterations(10);
        
        //check correctness
        assertEquals(estimator.getMaxIterations(), 10);
        
        //Force IllegalArgumentException
        try {
            estimator.setMaxIterations(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }    
    
    @Test
    public void testIsSetGeometricDistanceUsed() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        
        //set new value
        estimator.setGeometricDistanceUsed(
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        
        //check correctness
        assertEquals(estimator.isGeometricDistanceUsed(),
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
    }   
    
    @Test
    public void testGetSetDegree() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        
        //set new value
        estimator.setDegree(2);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        
        //Force IllegalArgumentException
        try {
            estimator.setDegree(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }    
    
    @Test
    public void testGetSetQualityScores() throws LockedException {
        LMedSPolynomialRobustEstimator estimator =
                new LMedSPolynomialRobustEstimator();
        
        //check default value
        assertNull(estimator.getQualityScores());
        
        //set new value
        estimator.setQualityScores(null);
        
        //check correctness
        assertNull(estimator.getQualityScores());
    }    
    
    @Test
    public void testEstimateDirectEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException, 
            RobustEstimatorException {                
        
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluate(x);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                        valueWithError);
                evaluations.add(eval);            
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);            
        }
    }
    
    @Test
    public void testEstimateDirectAndDerivativeEvaluationsAlgebraicDistance() 
            throws LockedException, NotReadyException, 
            RobustEstimatorException {   

        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations / 2; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluate(x);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                        valueWithError);
                evaluations.add(eval);            
            }
            for (int i = 0; i < numEvaluations / 2; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluateDerivative(x);
                
                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            
                
                DerivativePolynomialEvaluation eval =
                        new DerivativePolynomialEvaluation(x, valueWithError, 
                        1);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            boolean failed = false;
            for (int i = 0; i < polynomial2.getPolyParams().length; i++) {
                if (Math.abs(polynomial2.getPolyParams()[i] - polyParams[i]) > ABSOLUTE_ERROR) {
                    failed = true;
                    break;
                }
            }
            if (failed) {
                continue;
            }

            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }
    
    @Test
    public void testEstimateIntegralEvaluationsAlgebraicDistance() 
            throws LockedException, NotReadyException, 
            RobustEstimatorException {   
        
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double constant = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                Polynomial integral = polynomial.integrationAndReturnNew(
                        constant);
                double value = integral.evaluate(x);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                IntegralPolynomialEvaluation eval =
                        new IntegralPolynomialEvaluation(x, valueWithError, 
                                new double[]{constant}, 1);
                evaluations.add(eval);            
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);            
        }
    }
    
    @Test
    public void testEstimateIntegralIntervalEvaluationsAlgebraicDistance() 
            throws LockedException, NotReadyException, 
            RobustEstimatorException {   
        
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations; i++) {
                double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);            
                double value = polynomial.integrateInterval(startX, endX);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                IntegralIntervalPolynomialEvaluation eval = 
                        new IntegralIntervalPolynomialEvaluation(startX, endX, 
                        valueWithError, 1);
                evaluations.add(eval);            
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);            
        }
    }
    
    @Test
    public void testEstimateDirectEvaluationsGeometricDistance()
            throws LockedException, NotReadyException, RobustEstimatorException {                
        
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluate(x);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                        valueWithError);
                evaluations.add(eval);            
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);            
        }
    }
    
    @Test
    public void testEstimateDirectAndDerivativeEvaluationsGeometricDistance() 
            throws LockedException, NotReadyException, 
            RobustEstimatorException {   
        
        for (int t = 0; t < TIMES; t++) {
            LMedSPolynomialRobustEstimator estimator =
                    new LMedSPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            //check default values
            assertEquals(estimator.getDegree(), 1);
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

            //Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException e) { }

            //create random 1st degree polynomial
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            Polynomial polynomial = new Polynomial(polyParams);

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, 
                    MAX_EVALUATIONS);
            GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            List<PolynomialEvaluation> evaluations =
                    new ArrayList<PolynomialEvaluation>();
            for (int i = 0; i < numEvaluations / 2; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluate(x);

                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            

                DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                        valueWithError);
                evaluations.add(eval);            
            }
            for (int i = 0; i < numEvaluations / 2; i++) {
                double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                        MAX_RANDOM_VALUE);
                double value = polynomial.evaluateDerivative(x);
                
                double valueWithError;
                if(randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    //evaluation is outlier
                    double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                } else {
                    valueWithError = value;
                }            
                
                DerivativePolynomialEvaluation eval =
                        new DerivativePolynomialEvaluation(x, valueWithError, 
                        1);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);

            estimator.setListener(this);
            reset();

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);        
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());        

            //estimate
            Polynomial polynomial2 = estimator.estimate();

            //check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);            
        }
    }
    
    
    private void reset() {
        estimateStart = estimateEnd = estimateNextIteration =
                estimateProgressChange = 0;
    }

    @Override
    public void onEstimateStart(PolynomialRobustEstimator estimator) {
        estimateStart++;
    }

    @Override
    public void onEstimateEnd(PolynomialRobustEstimator estimator) {
        estimateEnd++;
    }

    @Override
    public void onEstimateNextIteration(PolynomialRobustEstimator estimator, 
            int iteration) {
        estimateNextIteration++;
    }

    @Override
    public void onEstimateProgressChange(PolynomialRobustEstimator estimator, 
            float progress) {
        estimateProgressChange++;
    }
}
