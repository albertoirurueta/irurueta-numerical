/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.PolynomialEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 10, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import java.util.ArrayList;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class PolynomialEstimatorTest implements PolynomialEstimatorListener {
    
    public PolynomialEstimatorTest() { }
    
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
        //default type
        PolynomialEstimator estimator = PolynomialEstimator.create();
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        //default type with degree
        estimator = PolynomialEstimator.create(2);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        estimator = PolynomialEstimator.create(evaluations);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with listener
        estimator = PolynomialEstimator.create(this);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with degree and listener
        estimator = PolynomialEstimator.create(2, this);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
    }
    
    @Test
    public void testCreateLMSE() {
        //default type
        PolynomialEstimator estimator = PolynomialEstimator.create(
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        //default type with degree
        estimator = PolynomialEstimator.create(2, 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        estimator = PolynomialEstimator.create(evaluations, 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with listener
        estimator = PolynomialEstimator.create(this, 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with degree and listener
        estimator = PolynomialEstimator.create(2, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);        
        
        
        //default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);                
    }

    @Test
    public void testCreateWeighted() {
        //default type
        PolynomialEstimator estimator = PolynomialEstimator.create(
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);        
        
        //default type with degree
        estimator = PolynomialEstimator.create(2,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        estimator = PolynomialEstimator.create(evaluations,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);        
        
        
        //default type with listener
        estimator = PolynomialEstimator.create(this,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);        
        
        
        //default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        
        //default type with degree and listener
        estimator = PolynomialEstimator.create(2, this,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        
        //default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);        
        
        
        //default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this,
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);             
    }
    
    @Override
    public void onEstimateStart(PolynomialEstimator estimator) { }

    @Override
    public void onEstimateEnd(PolynomialEstimator estimator) { }
}
