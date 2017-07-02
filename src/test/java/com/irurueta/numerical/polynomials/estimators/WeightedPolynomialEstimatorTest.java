/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.WeightedPolynomialEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 12, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
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

public class WeightedPolynomialEstimatorTest implements 
        PolynomialEstimatorListener {
    
    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public static final int MIN_DEGREE = 1;
    public static final int MAX_DEGREE = 5;
    
    public static final double ABSOLUTE_ERROR = 1e-8;
    
    private int estimateStart;
    private int estimateEnd; 
    
    public WeightedPolynomialEstimatorTest() { }
    
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
        //empty constructor
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(), 
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(), 
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        
        //test constructor with degre
        estimator = new WeightedPolynomialEstimator(2);
        
        //check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with evaluations and weights
        List<PolynomialEvaluation> evaluations =
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        double[] weights = new double[1];
        estimator = new WeightedPolynomialEstimator(evaluations, weights);
        
        //check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);

        
        //Force IllegalArgumentException
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(
                    (List<PolynomialEvaluation>)null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations, 
                    (double[])null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(wrongEvals, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations, 
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with listener
        estimator = new WeightedPolynomialEstimator(this);
        
        //check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(), 
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(), 
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        
        //test constructor with degree, evaluations and weights
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights);
        
        //check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);

        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(2,
                    (List<PolynomialEvaluation>)null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations, 
                    (double[])null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, wrongEvals, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations, 
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator); 
        
        
        //test degree and listener
        estimator = new WeightedPolynomialEstimator(2, this);
        
        //check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(evaluations, weights, this);
        
        //check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);

        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(
                    (List<PolynomialEvaluation>)null, weights, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations, 
                    (double[])null, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(wrongEvals, weights, 
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations, 
                    new double[2], this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);
        
        
        //test constructor with degree, evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights, 
                this);
        
        //check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(estimator.getMaxEvaluations(),
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        assertEquals(estimator.isSortWeightsEnabled(),
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR);

        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(2,
                    (List<PolynomialEvaluation>)null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations, 
                    (double[])null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, wrongEvals, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations, 
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(estimator);         
    }
    
    @Test
    public void testGetSetMaxEvaluations() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.getMaxEvaluations(), 
                WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS);
        
        //set new value
        estimator.setMaxEvaluations(100);
        
        //check correctness
        assertEquals(estimator.getMaxEvaluations(), 100);
        
        //Force IllegalArgumentException
        try {
            estimator.setMaxEvaluations(1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testIsSetSortWeightsEnabled() throws LockedException {
        WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.isSortWeightsEnabled(), 
                WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        
        //set new value
        estimator.setSortWeightsEnabled(
                !WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
        
        //check correctness
        assertEquals(estimator.isSortWeightsEnabled(),
                !WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
    }
    
    @Test
    public void testGetSetEvaluationsAndWeights() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());
        
        //set new values
        List<PolynomialEvaluation> evaluations =
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        double[] weights = new double[1];        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        //check correctness
        assertSame(estimator.getEvaluations(), evaluations);
        assertSame(estimator.getWeights(), weights);
        
        //Force IllegalArgumentException
        try {
            estimator.setEvaluations(evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeAndEvaluations(2, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        try {
            estimator.setEvaluationsAndWeights((List<PolynomialEvaluation>)null, 
                    weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setEvaluationsAndWeights(evaluations, (double[])null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setEvaluationsAndWeights(wrongEvals, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setEvaluationsAndWeights(evaluations, new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }                
    }
    
    @Test
    public void testGetSetDegree() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.getDegree(), 1);
        
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
    public void testGetSetDegreeEvaluationsAndWeights() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        assertEquals(estimator.getDegree(), 
                WeightedPolynomialEstimator.MIN_DEGREE);
        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());
        
        //set new values
        List<PolynomialEvaluation> evaluations =
                new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        double[] weights = new double[1];        
        estimator.setDegreeEvaluationsAndWeights(2, evaluations, weights);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertSame(estimator.getWeights(), weights);
        
        //Force IllegalArgumentException
        try {
            estimator.setEvaluations(evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeAndEvaluations(2, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeEvaluationsAndWeights(0, evaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }

        
        List<PolynomialEvaluation> wrongEvals = 
                new ArrayList<PolynomialEvaluation>();
        try {
            estimator.setDegreeEvaluationsAndWeights(2, 
                    (List<PolynomialEvaluation>)null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, evaluations, 
                    (double[])null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, wrongEvals, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, evaluations, 
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }                
    }

    @Test
    public void testIsReady() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = new double[2];
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
            weights[i] = 1.0;
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testGetMinNumberOfEvaluations() {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int degree = randomizer.nextInt(MIN_DEGREE, MAX_DEGREE);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(degree), 
                degree + 1);
        
        //Force IllegalArgumentException
        try {
            PolynomialEstimator.getMinNumberOfEvaluations(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator(degree);
        assertEquals(estimator.getMinNumberOfEvaluations(), degree + 1);
    }
    
    @Test
    public void testGetSetListener() throws LockedException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default value
        assertNull(estimator.getListener());
        
        //set new value
        estimator.setListener(this);
        
        //check correctness
        assertSame(estimator.getListener(), this);
    }
    
    @Test
    public void testEstimateWithDirectEvaluations() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
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
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);    
    }
    
    @Test
    public void testEstimateWithDirectAndDerivativeEvaluations() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
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
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[4 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluateDerivative(x);
            DerivativePolynomialEvaluation eval = 
                    new DerivativePolynomialEvaluation(x, value, 1);
            evaluations.add(eval);
            weights[2 * estimator.getMinNumberOfEvaluations() + i] = 
                    randomizer.nextDouble(0.5, 1.0);
        }
        
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);     
    }
    
    @Test
    public void testEstimateWithIntegralEvaluations() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
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
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double constant = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.integrationAndReturnNew(constant);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, 
                            new double[]{constant}, 1);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);     
    }    

    @Test
    public void testEstimateWithIntegralIntervalEvaluations() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
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
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);            
            double value = polynomial.integrateInterval(startX, endX);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);     
    }    
    
    @Test
    public void testEstimateWithDirectEvaluationsSecondDegree() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator(2);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException e) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);      
    }
    
    @Test
    public void testEstimateWithSecondOrderIntegralEvaluationsSecondDegree() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator(2);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException e) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.nthIntegrationAndReturnNew(2, 
                    constants);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);       
    }

    @Test
    public void testEstimateWithSecondOrderIntegralIntervalEvaluationsSecondDegree() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        WeightedPolynomialEstimator estimator = 
                new WeightedPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMaxEvaluations(), 50);
        
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
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<PolynomialEvaluation>();
        double[] weights = 
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);  
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            
            double value = polynomial.nthOrderIntegrateInterval(startX, endX, 2,
                    constants);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);
            
            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        
        estimator.setEvaluationsAndWeights(evaluations, weights);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);      
    }    
    
    private void reset() {
        estimateStart = estimateEnd = 0;
    }

    @Override
    public void onEstimateStart(PolynomialEstimator estimator) {
        estimateStart++;
        checkIsLocked(estimator);
    }

    @Override
    public void onEstimateEnd(PolynomialEstimator estimator) {
        estimateEnd++;
        checkIsLocked(estimator);
    }
    
    private void checkIsLocked(PolynomialEstimator estimator) {
        assertTrue(estimator.isLocked());
        
        //Force LockedException
        try {
            estimator.setDegree(2);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }
        try {
            ((WeightedPolynomialEstimator)estimator).setEvaluationsAndWeights(
                    null, null);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }
        try {
            ((WeightedPolynomialEstimator)estimator).
                    setDegreeEvaluationsAndWeights(2, null, null);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }     
        try {
            ((WeightedPolynomialEstimator)estimator).setMaxEvaluations(1);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }             
        try {
            ((WeightedPolynomialEstimator)estimator).setSortWeightsEnabled(
                    false);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }                     
        try {
            estimator.setListener(null);
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { }
        try {
            estimator.estimate();
            fail("LockedException expected but not thrown");
        } catch (LockedException e) { 
        } catch (Exception e) { 
            fail("LockedException expected but not thrown");
        }
    }    
}
