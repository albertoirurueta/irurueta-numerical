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

public class LaguerrePolynomialRootsEstimatorTest {
    
    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;
    
    private static final double TOLERANCE = 3e-8;
    
    private static final int TIMES = 100;
    
    public LaguerrePolynomialRootsEstimatorTest() { }

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
        RootEstimationException, NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        boolean polishRoots = randomizer.nextBoolean();
        
        Complex[] polyParams = new Complex[2];
        Complex[] badPolyParams = new Complex[1];
        LaguerrePolynomialRootsEstimator estimator;
        
        
        //test 1st constructor
        estimator = new LaguerrePolynomialRootsEstimator();
        assertNotNull(estimator);
        
        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.areRootsPolished(),
                LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS);
        try{
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        
        
        //test 2nd constructor
        estimator = new LaguerrePolynomialRootsEstimator(polishRoots);
        assertNotNull(estimator);
        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.areRootsPolished(), polishRoots);
        try{
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());        
        
        //test 3rd constructor
        estimator = new LaguerrePolynomialRootsEstimator(polyParams);
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.areRootsPolished(),
                LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS);
        assertSame(estimator.getPolynomialParameters(), polyParams);
        try{
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        
        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new LaguerrePolynomialRootsEstimator(badPolyParams);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(estimator);
        
        
        //test 4th constructor
        estimator = new LaguerrePolynomialRootsEstimator(polyParams, 
                polishRoots);
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.areRootsPolished(), polishRoots);
        assertSame(estimator.getPolynomialParameters(), polyParams);
        try{
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        
        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new LaguerrePolynomialRootsEstimator(badPolyParams,
                    polishRoots);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(estimator);        
    }
    
    @Test
    public void testGetSetPolishRoots() throws LockedException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        boolean polishRoots = randomizer.nextBoolean();
        
        LaguerrePolynomialRootsEstimator estimator =
                new LaguerrePolynomialRootsEstimator();
        
        //check default value
        assertEquals(estimator.areRootsPolished(),
                LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS);
        
        //set new value
        estimator.setPolishRootsEnabled(polishRoots);
        //check correctness
        assertEquals(estimator.areRootsPolished(), polishRoots);
    }
    
    @Test
    public void testGetSetPolynomialParameters() throws LockedException, 
        NotAvailableException{
        
        Complex[] polyParams = new Complex[2];
        Complex[] badPolyParams = new Complex[1];
        
        LaguerrePolynomialRootsEstimator estimator = 
                new LaguerrePolynomialRootsEstimator();
        
        //check default values
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");           
        }catch(NotAvailableException e){}
        assertFalse(estimator.arePolynomialParametersAvailable());
        
        //set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        //check correctness
        assertSame(estimator.getPolynomialParameters(), polyParams);
        
        //Force IllegalArgumentException
        try{
            estimator.setPolynomialParameters(badPolyParams);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException, 
        RootEstimationException, NotAvailableException{
        
        for(int t = 0; t < TIMES; t++){
        
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT / 3.0);
            double realRoot2 = randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 
                    2.0 / 3.0 * MAX_EVAL_POINT);
            double realRoot3 = randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT);
        
            Complex root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT / 3.0));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT / 3.0));
        
            Complex root2 = new Complex();
            root2.setReal(randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 
                    2.0 / 3.0 * MAX_EVAL_POINT));
            root2.setImaginary(randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 
                    2.0 / 3.0 * MAX_EVAL_POINT));
        
            Complex root3 = new Complex();
            root3.setReal(randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT));
            root3.setImaginary(randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, 
                    MAX_EVAL_POINT));
        
            //also compute conjugates
            Complex conjRoot1 = root1.conjugateAndReturnNew();
            Complex conjRoot2 = root2.conjugateAndReturnNew();
        
        
            LaguerrePolynomialRootsEstimator estimator = 
                    new LaguerrePolynomialRootsEstimator();
        
            Complex[] polyParams;
            Complex[] roots;
        
            //attempt set parameters for constant
            polyParams = generateConstantPolynomialParams(new Complex(
                    realRoot1));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
        
            //set parameters for first degree polynomial with real root
            polyParams = generateFirstDegreePolynomialParams(new Complex(
                    realRoot1));
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 1);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
        
            //set parameters for first degree polynomial with complex root
            polyParams = generateFirstDegreePolynomialParams(root1);
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 1);
            assertTrue(roots[0].equals(root1, TOLERANCE));        

            //set parameters for second degree polynomial with real roots
            polyParams = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2));
            estimator.setPolynomialParameters(polyParams);
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
            estimator.setPolynomialParameters(polyParams);
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
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 2);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));

            //set parameters for third degree polynomial with real roots
            polyParams = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2), 
                    new Complex(realRoot3));
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 3);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot3), TOLERANCE));

            //set parameters for third degree polynomial with real root and two
            //complex conjugate roots
            if(realRoot1 < root2.getReal()){
                polyParams = generateThirdDegreePolynomialParams(
                        new Complex(realRoot1), root2, conjRoot2);
            }else{
                polyParams = generateThirdDegreePolynomialParams(root2, 
                        conjRoot2, new Complex(realRoot1));
            }
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 3);
            //because roots1 might be exchanged, we check for their
            //real parts and absolute value of their imaginary parts (which are
            //the same but with opposite sign because they are complex 
            //conjugates)
            assertEquals(roots[0].getReal(), realRoot1, TOLERANCE);
            assertEquals(Math.abs(roots[0].getImaginary()), 0.0, TOLERANCE);
            assertEquals(roots[1].getReal(), root2.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[1].getImaginary()), 
                    Math.abs(root2.getImaginary()), TOLERANCE);
            assertEquals(roots[2].getReal(), conjRoot2.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[2].getImaginary()), 
                    Math.abs(conjRoot2.getImaginary()), TOLERANCE);                    

            //set parameters for third degree polynomial with two double real 
            //roots
            polyParams = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2), 
                    new Complex(realRoot2));
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
            
            assertEquals(roots.length, 3);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot2), TOLERANCE));

            //set parameters for third degree polynomial with one triple real 
            //roots
            polyParams = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1), 
                    new Complex(realRoot1));
            estimator.setPolynomialParameters(polyParams);
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
        
            assertEquals(roots.length, 3);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot1), TOLERANCE));
        }
    }

    private double vectorNorm(Complex[] v) {
        double normValue = 0.0, real, imag;
        for (Complex value : v) {
            real = value.getReal();
            imag = value.getImaginary();
            normValue += real * real + imag * imag; //square norm
        }

        return Math.sqrt(normValue);
    }

    private Complex[] generateConstantPolynomialParams(Complex param) {

        Complex[] out = new Complex[1];
        out[0] = param;
        return out;
    }

    private Complex[] generateFirstDegreePolynomialParams(Complex root1) {

        Complex[] out = new Complex[2];
        //p(x) = x - root1
        out[1] = new Complex(1.0, 0.0);
        out[0] = new Complex(-root1.getReal(), -root1.getImaginary());

        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }

    private Complex[] generateSecondDegreePolynomialParams(Complex root1,
            Complex root2) {

        Complex[] out = new Complex[3];
        //p(x) = (x - root1) * (x - root2) = x * x - (root1 * root2) * x +
        //root1 + root2
        out[2] = new Complex(1.0, 0.0);
        out[1] = root1.addAndReturnNew(root2).multiplyByScalarAndReturnNew(
                -1.0);
        out[0] = root1.multiplyAndReturnNew(root2);

        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }

    private Complex[] generateThirdDegreePolynomialParams(Complex root1,
            Complex root2, Complex root3) {

        Complex[] out = new Complex[4];
        //p(x) = (x - root1) * (x - root2) * (x - root3) =
        //(x * x - (root1 + root2) * x + root1 * root2) * (x - root3) =
        //(x * x * x - (root1 + root2) * x * x + (root1 + root2) * x
        //- root3 * x * x + (root1 + root2) * root3 * x
        //- (root1 + root2) * root3 =

        //x * x * x - (root1 + root2 + root3) * x * x +
        //((root1 * root2) + (root1 + root2) * root3) * x
        //- root1 * root2 * root3

        out[3] = new Complex(1.0, 0.0);
        out[2] = root1.addAndReturnNew(root2).addAndReturnNew(root3).
                multiplyByScalarAndReturnNew(-1.0);
        out[1] = root1.multiplyAndReturnNew(root2).addAndReturnNew(
                root1.addAndReturnNew(root2).multiplyAndReturnNew(root3));
        out[0] = root1.multiplyAndReturnNew(root2).multiplyAndReturnNew(root3).
                multiplyByScalarAndReturnNew(-1.0);

        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }

}
