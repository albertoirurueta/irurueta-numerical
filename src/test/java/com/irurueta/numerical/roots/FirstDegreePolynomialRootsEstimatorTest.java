/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.roots.FirstDegreePolynomialRootsEstimator
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

public class FirstDegreePolynomialRootsEstimatorTest {

    public static final double MIN_EVAL_POINT = 0.0;
    public static final double MAX_EVAL_POINT = 1.0;
    
    public static final double TOLERANCE = 3e-8;
    
    public static final int TIMES = 100;
    
    public FirstDegreePolynomialRootsEstimatorTest() {
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
        
    @Test
    public void testConstructor() throws LockedException, 
        RootEstimationException, 
        NotAvailableException, 
        NotReadyException{
        
        double[] polyParams = new double[2];
        polyParams[1] = 1.0; //to ensure it's second degree
        
        FirstDegreePolynomialRootsEstimator estimator;
        
        //test 1st constructor
        estimator = new FirstDegreePolynomialRootsEstimator();
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
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        try{
            estimator.isFirstDegree();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
        assertTrue(estimator.isRealSolution());
        
        
        //test 2nd constructor
        estimator = new FirstDegreePolynomialRootsEstimator(polyParams);
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
        assertTrue(estimator.isFirstDegree());
        assertTrue(estimator.isReady());
        
        //Force IllegalArgumentException
        double[] badPolyParams = new double[1];
        
        estimator = null;
        try{
            estimator = new FirstDegreePolynomialRootsEstimator(badPolyParams);
        }catch(IllegalArgumentException e){}
        assertNull(estimator);
    }
    
    @Test
    public void testGetSetPolynomialParameters() throws LockedException, 
        NotAvailableException{
        
        double[] polyParams = new double[2];
        polyParams[1] = 1.0;
        double[] polyParams2 = new double[3];
        polyParams2[1] = 1.0;
        polyParams2[2] = 0.0;
        double[] badPolyParams = new double[1];
        double[] badPolyParams2 = new double[2];
        badPolyParams2[1] = 0.0;
        
        FirstDegreePolynomialRootsEstimator estimator = 
                new FirstDegreePolynomialRootsEstimator();
        
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
            estimator.setPolynomialParameters(new Complex[3]);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException, 
        RootEstimationException, NotAvailableException{
        
        for(int t = 0; t < TIMES; t++){
        
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            double realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT);
        
            Complex root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, 
                    MAX_EVAL_POINT));
                
        
            FirstDegreePolynomialRootsEstimator estimator = 
                    new FirstDegreePolynomialRootsEstimator();
        
            double[] polyParams;
            Complex[] roots;
        
            //attempt set parameters for constant
            polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(FirstDegreePolynomialRootsEstimator.isFirstDegree(
                    polyParams));
            try{
                estimator.setPolynomialParameters(polyParams);
                fail("IllegalArgumentException expected but not thrown");
            }catch(IllegalArgumentException e){}
            try{
                estimator.isFirstDegree();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}            
            try{
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            }catch(NotReadyException e){}
            //check correctness
            assertFalse(estimator.areRootsAvailable());
            
        
            //set parameters for first degree polynomial with real root
            polyParams = generateFirstDegreePolynomialParams(realRoot1);
            assertTrue(FirstDegreePolynomialRootsEstimator.isFirstDegree(
                    polyParams));
            estimator.setPolynomialParameters(polyParams);
            assertTrue(estimator.isFirstDegree());
            assertTrue(estimator.isRealSolution());
            estimator.estimate();
            //check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();
            
            assertEquals(roots.length, 1);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
        }
    }    
}
