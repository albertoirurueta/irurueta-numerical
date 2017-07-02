/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 16, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.robust.RobustEstimatorMethod;
import java.util.ArrayList;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class PolynommialRobustEstimatorTest implements 
        PolynomialRobustEstimatorListener{
    
    public PolynommialRobustEstimatorTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testCreate() {
        //test empty creator
        PolynomialRobustEstimator estimator = 
                PolynomialRobustEstimator.create();
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        
        //test creator with degree
        estimator = PolynomialRobustEstimator.create(2);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener
        estimator = PolynomialRobustEstimator.create(this);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        
        //test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
    }
    
    @Test
    public void testCreteRANSAC() {
        //test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        
        //test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, 
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);        
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, 
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and method
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, 
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, 
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        
        //test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, 
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);   
        
        
        //test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.RANSAC);
        
        //check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.RANSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);        
    }
    
    @Test
    public void testCreateLMedS() {
        //test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);        
        
        
        //test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, 
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);        
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, 
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and method
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, 
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, 
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        
        //test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, 
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);   
        
        
        //test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.LMedS);
        
        //check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.LMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.LMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);        
    }
    
    @Test
    public void testCreateMSAC() {
        //test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);        
        
        
        //test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, 
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);        
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, 
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and method
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, 
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, 
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        
        //test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, 
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);   
        
        
        //test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.MSAC);
        
        //check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);        
    }
    
    @Test
    public void testCreatePROSAC() {
        //test empty creator
        PolynomialRobustEstimator estimator = 
                PolynomialRobustEstimator.create(RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        
        //test creator with degree
        estimator = PolynomialRobustEstimator.create(2, 
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, 
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        
        //test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this, 
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.PROSAC);
        
        //check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
    }

    @Test
    public void testCreatePROMedS() {
        //test empty creator
        PolynomialRobustEstimator estimator = 
                PolynomialRobustEstimator.create(RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        
        //test creator with degree
        estimator = PolynomialRobustEstimator.create(2, 
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, 
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with listener
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        
        //test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this, 
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
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
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.PROMedS);
        
        //check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
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
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROMedS);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.PROMedS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
    }    

    @Override
    public void onEstimateStart(PolynomialRobustEstimator estimator) { }

    @Override
    public void onEstimateEnd(PolynomialRobustEstimator estimator) { }

    @Override
    public void onEstimateNextIteration(PolynomialRobustEstimator estimator, 
            int iteration) { }

    @Override
    public void onEstimateProgressChange(PolynomialRobustEstimator estimator, 
            float progress) { }
}
