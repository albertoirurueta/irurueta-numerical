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
package com.irurueta.numerical.roots;

import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.*;

@SuppressWarnings("Duplicates")
public class SecondDegreePolynomialRootsEstimatorTest {
    
    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;
    
    private static final double TOLERANCE = 3e-8;
    
    private static final int TIMES = 100;
    
    public SecondDegreePolynomialRootsEstimatorTest() { }

    @BeforeClass
    public static void setUpClass() { }

    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }

    @Test
    public void testConstructor() throws LockedException, 
            RootEstimationException, NotAvailableException,
            NotReadyException {
        
        double[] polyParams = new double[3];
        polyParams[2] = 1.0; //to ensure it's second degree
        
        SecondDegreePolynomialRootsEstimator estimator;
        
        //test 1st constructor
        estimator = new SecondDegreePolynomialRootsEstimator();
        assertNotNull(estimator);
        
        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        try {
            //noinspection deprecation
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        try {
            estimator.isSecondDegree();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        try {
            estimator.hasTwoDistinctRealRoots();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        try {
            estimator.hasDoubleRoot();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        try {
            estimator.hasTwoComplexConjugateRoots();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //test 2nd constructor
        estimator = new SecondDegreePolynomialRootsEstimator(polyParams);
        assertNotNull(estimator);
        
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams, 
                TOLERANCE);
        try {
            //noinspection deprecation
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertTrue(estimator.isSecondDegree());
        estimator.hasTwoDistinctRealRoots();
        estimator.hasDoubleRoot();
        estimator.hasTwoComplexConjugateRoots();
        
        //Force IllegalArgumentException
        double[] badPolyParams = new double[1];
        
        estimator = null;
        try {
            estimator = new SecondDegreePolynomialRootsEstimator(badPolyParams);
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
    }
    
    @Test
    public void testGetSetPolynomialParameters() throws LockedException, 
        NotAvailableException{
        
        double[] polyParams = new double[3];
        polyParams[2] = 1.0;
        double[] polyParams2 = new double[4];
        polyParams2[2] = 1.0;
        polyParams2[3] = 0.0;
        double[] badPolyParams = new double[2];
        double[] badPolyParams2 = new double[3];
        badPolyParams2[2] = 0.0;
        
        SecondDegreePolynomialRootsEstimator estimator = 
                new SecondDegreePolynomialRootsEstimator();
        
        //check default values
        try {
            //noinspection deprecation
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getRealPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.arePolynomialParametersAvailable());
        
        //set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        //check correctness
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams, 
                TOLERANCE);
        try {
            //noinspection deprecation
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        estimator.setPolynomialParameters(polyParams2);
        //check correctness
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams2, 
                TOLERANCE);
        try {
            //noinspection deprecation
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        
        //Force IllegalArgumentException
        try{
            estimator.setPolynomialParameters(badPolyParams);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        
        try {
            estimator.setPolynomialParameters(badPolyParams2);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        
        //attempting to use complex parameters will also raise an 
        //IllegalArgumentException
        try {
            estimator.setPolynomialParameters(new Complex[3]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException, 
        RootEstimationException, NotAvailableException {
        
        for (int t = 0; t < TIMES; t++) {
        
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, 
                    0.5 * MAX_EVAL_POINT);
            double realRoot2 = randomizer.nextDouble(0.5 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT);
        
            Complex root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, 
                    0.5 * MAX_EVAL_POINT));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, 
                    0.5 * MAX_EVAL_POINT));
        
            Complex root2 = new Complex();
            root2.setReal(randomizer.nextDouble(0.5 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT));
            root2.setImaginary(randomizer.nextDouble(0.5 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT));
        
        
            //also compute conjugates
            Complex conjRoot1 = root1.conjugateAndReturnNew();
        
            SecondDegreePolynomialRootsEstimator estimator = 
                    new SecondDegreePolynomialRootsEstimator();
        
            double[] polyParams;
            Complex[] roots;
        
            //attempt set parameters for constant
            polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(SecondDegreePolynomialRootsEstimator.isSecondDegree(
                    polyParams));
            try {
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            } catch (IllegalArgumentException ignore) { }
            try {
                estimator.isSecondDegree();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException ignore) { }
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException ignore) { }
            //check correctness
            assertFalse(estimator.areRootsAvailable());
            
        
            //set parameters for first degree polynomial with real root
            polyParams = generateFirstDegreePolynomialParams(realRoot1);
            assertFalse(SecondDegreePolynomialRootsEstimator.isSecondDegree(
                    polyParams));
            try {
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            } catch (IllegalArgumentException ignore) { }
            try {
                estimator.isSecondDegree();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException ignore) { }
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException ignore) { }
            //check correctness
            assertFalse(estimator.areRootsAvailable());
        
            //set parameters for second degree polynomial with real roots
            polyParams = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2));
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(
                    polyParams));
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isSecondDegree());
            assertTrue(estimator.hasTwoDistinctRealRoots());
            assertFalse(estimator.hasDoubleRoot());
            assertFalse(estimator.hasTwoComplexConjugateRoots());
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 2);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));

            
            //set parameters for second degree polynomial with complex conjugate
            //roots (and real coeficients)
            polyParams = generateSecondDegreePolynomialParams(root1, conjRoot1);
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(
                    polyParams));
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isSecondDegree());
            assertFalse(estimator.hasTwoDistinctRealRoots());
            assertFalse(estimator.hasDoubleRoot());
            assertTrue(estimator.hasTwoComplexConjugateRoots());
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 2);
            //because root[0] and root[1] might be exchanged, we check for their
            //real parts and absolute value of their imaginary parts (which are
            //the same but with opposite sign because they are complex 
            //conjugates)
            assertEquals(roots[0].getReal(), root1.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[0].getImaginary()), 
                    Math.abs(root1.getImaginary()), TOLERANCE);
            assertEquals(roots[1].getReal(), conjRoot1.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[1].getImaginary()), 
                    Math.abs(conjRoot1.getImaginary()), TOLERANCE);
        
            //set parameters for second degree polynomial with double real roots
            polyParams = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1));
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(
                    polyParams));
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isSecondDegree());
            assertFalse(estimator.hasTwoDistinctRealRoots());
            assertTrue(estimator.hasDoubleRoot());
            assertFalse(estimator.hasTwoComplexConjugateRoots());
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 2);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));

        }
    }

    private double vectorNorm(double[] v) {
        double normValue = 0.0;
        for (double value : v) {
            normValue += value * value; //square norm
        }

        return Math.sqrt(normValue);
    }

    private double[] generateConstantPolynomialParams(double param) {

        double[] out = new double[1];
        out[0] = param;
        return out;
    }

    private double[] generateFirstDegreePolynomialParams(double root1) {

        double[] out = new double[2];
        //p(x) = x - root1
        out[1] = 1.0;
        out[0] = -root1;

        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);
        return out;
    }

    private double[] generateSecondDegreePolynomialParams(Complex root1,
                                                          Complex root2) {

        double[] out = new double[3];
        //p(x) = (x - root1) * (x - root2) = x * x - (root1 + root2) * x +
        //root1 * root2
        out[2] = 1.0;
        out[1] = root1.addAndReturnNew(root2).multiplyByScalarAndReturnNew(
                -1.0).getReal();
        out[0] = root1.multiplyAndReturnNew(root2).getReal();

        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);
        return out;
    }
}
