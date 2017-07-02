/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.roots.ThirdDegreePolynomialRootsEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date June 2, 2012
 */
package com.irurueta.numerical.roots;

import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.*;
import org.junit.*;

public class ThirdDegreePolynomialRootsEstimatorTest {
    
    public static final double MIN_EVAL_POINT = 0.0;
    public static final double MAX_EVAL_POINT = 1.0;
    
    public static final double TOLERANCE = 3e-7;
    
    public static final int TIMES = 100;
    
    public ThirdDegreePolynomialRootsEstimatorTest() {
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
    
    public double vectorNorm(double[] v){
        double normValue = 0.0, value;
        for(int i = 0; i < v.length; i++){
            value = v[i];
            normValue += value * value; //square norm
        }
        
        return Math.sqrt(normValue);
    }
    
    public double[] generateConstantPolynomialParams(double param){
        
        double[] out = new double[1];
        out[0] = param;
        return out;
    }
    
    public double[] generateFirstDegreePolynomialParams(double root1){
        
        double[] out = new double[2];
        //p(x) = x - root1
        out[1] = 1.0;
        out[0] = -root1;
        
        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);        
        return out;
    }
    
    public double[] generateSecondDegreePolynomialParams(Complex root1, 
            Complex root2){
        
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
    
    public double[] generateThirdDegreePolynomialParams(Complex root1,
            Complex root2, Complex root3){
        
        double[] out = new double[4];
        //p(x) = (x - root1) * (x - root2) * (x - root3) =
        //(x * x - (root1 + root2) * x + root1 * root2) * (x - root3) =
        //(x * x * x - (root1 + root2) * x * x + (root1 + root2) * x
        //- root3 * x * x + (root1 + root2) * root3 * x
        //- (root1 + root2) * root3 =
        
        //x * x * x - (root1 + root2 + root3) * x * x +
        //((root1 * root2) + (root1 + root2) * root3) * x 
        //- root1 * root2 * root3
        
        out[3] = 1.0;
        out[2] = root1.addAndReturnNew(root2).addAndReturnNew(root3).
                multiplyByScalarAndReturnNew(-1.0).getReal();
        out[1] = root1.multiplyAndReturnNew(root2).addAndReturnNew(
                root1.addAndReturnNew(root2).multiplyAndReturnNew(root3)).
                getReal();
        out[0] = root1.multiplyAndReturnNew(root2).multiplyAndReturnNew(root3).
                multiplyByScalarAndReturnNew(-1.0).getReal();
        
        double normValue = vectorNorm(out);
        //normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);        
        return out;
    }
    
    @Test
    public void testConstructor() throws LockedException, 
        RootEstimationException, 
        NotAvailableException, 
        NotReadyException{
        
        double[] polyParams = new double[4];
        polyParams[3] = 1.0; //to ensure it's third degree
        
        ThirdDegreePolynomialRootsEstimator estimator;
        
        //test 1st constructor
        estimator = new ThirdDegreePolynomialRootsEstimator();
        assertNotNull(estimator);
        
        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
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
        try{
            estimator.isThirdDegree();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        try{
            estimator.hasThreeDistinctRealRoots();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        try{
            estimator.hasMultipleRealRoot();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        try{
            estimator.hasOneRealRootAndTwoComplexConjugateRoots();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        
        
        
        //test 2nd constructor
        estimator = new ThirdDegreePolynomialRootsEstimator(polyParams);
        assertNotNull(estimator);
        
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.getRealPolynomialParameters(), polyParams);
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getRoots();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertTrue(estimator.isThirdDegree());
        estimator.hasThreeDistinctRealRoots();
        estimator.hasMultipleRealRoot();
        estimator.hasOneRealRootAndTwoComplexConjugateRoots();
        
        //Force IllegalArgumentException
        double[] badPolyParams = new double[1];
        
        estimator = null;
        try{
            estimator = new ThirdDegreePolynomialRootsEstimator(badPolyParams);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(estimator);
    }
    
    @Test
    public void testGetSetPolynomialParameters() throws LockedException, 
        NotAvailableException{
        
        double[] polyParams = new double[4];
        polyParams[3] = 1.0;
        double[] polyParams2 = new double[5];
        polyParams2[3] = 1.0;
        polyParams2[4] = 0.0;
        double[] badPolyParams = new double[3];
        double[] badPolyParams2 = new double[4];
        badPolyParams2[3] = 0.0;
        
        ThirdDegreePolynomialRootsEstimator estimator =
                new ThirdDegreePolynomialRootsEstimator();
        
        //check default values
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        try{
            estimator.getRealPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}        
        assertFalse(estimator.arePolynomialParametersAvailable());
        
        //set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        //check correctness
        assertEquals(estimator.getRealPolynomialParameters(), polyParams);
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        
        estimator.setPolynomialParameters(polyParams2);
        //check correctness
        assertEquals(estimator.getRealPolynomialParameters(), polyParams2);
        try{
            estimator.getPolynomialParameters();
            fail("NotAvailableException expected but not thrown");
        }catch(NotAvailableException e){}
        
        
        //Force IllegalArgumentException
        try{
            estimator.setPolynomialParameters(badPolyParams);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        
        try{
            estimator.setPolynomialParameters(badPolyParams2);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        
        //attempting to use complex parameters will also raise an 
        //IllegalArgumentException
        try{
            estimator.setPolynomialParameters(new Complex[4]);
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

            ThirdDegreePolynomialRootsEstimator estimator = 
                    new ThirdDegreePolynomialRootsEstimator();
            
            double[] polyParams;
            Complex[] roots;
        
            //attempt set parameters for constant
            polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));        
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isThirdDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());
        
        
            //attempt estimate root for first degree polynomial
            polyParams = generateFirstDegreePolynomialParams(realRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isThirdDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());
        
            //set parameters for second degree polynomial with real roots
            polyParams = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2));
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isThirdDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());
            
            //set parameters for second degree polynomial with complex conjugate
            //roots (and real coeficients)
            polyParams = generateSecondDegreePolynomialParams(root1, conjRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isThirdDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());
        
            //set parameters for second degree polynomial with double real roots
            polyParams = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1));
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isThirdDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());

            //set parameters for third degree polynomial with real roots
            polyParams = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2), 
                    new Complex(realRoot3));
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isThirdDegree());
            assertTrue(estimator.hasThreeDistinctRealRoots());
            assertFalse(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());        
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
            polyParams = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), root2, conjRoot2);
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));        
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertFalse(estimator.hasMultipleRealRoot());
            assertTrue(estimator.hasOneRealRootAndTwoComplexConjugateRoots());                
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
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));                
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertTrue(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());        
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
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(
                    polyParams));                        
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertTrue(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());        
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
}
