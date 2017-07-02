/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.SvdMultiDimensionLinearFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 27, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class SvdMultiDimensionLinearFitterTest {
    
    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    
    public static final double MIN_SIGMA_VALUE = 1e-4;
    public static final double MAX_SIGMA_VALUE = 1.0;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    //For functions like: a*1 + b*x + c*y + d*x*y + e*x^2 + f*y^2
    public static final int NUM_QUADRATIC_PARAMS = 6;
    
    public SvdMultiDimensionLinearFitterTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}
    
    @Test
    public void testConstructor() throws FittingException, WrongSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        //test empty constructor
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdMultiDimensionLinearFitter.TOL, 0.0);
        
        //test constructor with input data
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        fitter = new SvdMultiDimensionLinearFitter(x, y, sig);
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdMultiDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        double[] shortY = new double[nPoints - 1];
        double[] shortSig = new double[nPoints - 1];
        
        fitter = null;
        try{
            fitter = new SvdMultiDimensionLinearFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(x, y, 1.0);
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++){
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdMultiDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdMultiDimensionLinearFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with evaluator
        LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }
                    
            @Override
            public double[] createResultArray() {
                return new double[nPoints];
            }

            @Override
            public void evaluate(double[] point, double[] result) 
                    throws Throwable {
            }
        };
        
        fitter = new SvdMultiDimensionLinearFitter(evaluator);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //test constructor with evaluator and input data
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, sig);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdMultiDimensionLinearFitter(evaluator, shortX, y, 
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, shortY, 
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, 
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with evaluator and input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertTrue(fitter.isReady());
        assertEquals(fitter.getTol(), SvdMultiDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdMultiDimensionLinearFitter(evaluator, shortX, y, 
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, shortY, 
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);        
    }

    @Test
    public void testGetSetFunctionEvaluator() throws FittingException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        
        //set new value
        LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }
                    
            @Override
            public double[] createResultArray() {
                return new double[nPoints];
            }

            @Override
            public void evaluate(double[] point, double[] result) 
                    throws Throwable {
            }
        };
        
        //set new value
        fitter.setFunctionEvaluator(evaluator);
        
        //check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
    }
    
    @Test
    public void testGetSetInputData() throws WrongSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        //set input data
        fitter.setInputData(x, y, sig);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        
        //Force IllegalArgumentException
        Matrix wrongX = new Matrix(nPoints - 1, 2);
        double[] wrongY = new double[nPoints - 1];
        double[] wrongSig = new double[nPoints - 1];
        
        try{
            fitter.setInputData(wrongX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, wrongY, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, y, wrongSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }   
    
    @Test
    public void testGetSetInputDataWithConstantSigma() 
            throws WrongSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        
        //set input data
        fitter.setInputData(x, y, 1.0);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++){
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        
        //Force IllegalArgumentException
        Matrix wrongX = new Matrix(nPoints - 1, 2);
        double[] wrongY = new double[nPoints - 1];
        
        try{
            fitter.setInputData(wrongX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, wrongY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testIsReady() throws FittingException, WrongSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default value
        assertFalse(fitter.isReady());
        
        //set new values
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        fitter.setInputData(x, y, sig);
        
        assertFalse(fitter.isReady());
        
        LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }
                    
            @Override
            public double[] createResultArray() {
                return new double[nPoints];
            }

            @Override
            public void evaluate(double[] point, double[] result) 
                    throws Throwable {
            }
        };
        
        fitter.setFunctionEvaluator(evaluator);
        
        assertTrue(fitter.isReady());
        
        //test bad evaluator
        LinearFitterMultiDimensionFunctionEvaluator badEvaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 3;
            }
                    
            @Override
            public double[] createResultArray() {
                return new double[nPoints];
            }

            @Override
            public void evaluate(double[] point, double[] result) 
                    throws Throwable {
            }
        };

        fitter.setFunctionEvaluator(badEvaluator);
        
        assertFalse(fitter.isReady());        
    }
    
    @Test
    public void testGetSetTol(){
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter();
        
        //check default values
        assertEquals(fitter.getTol(), SvdMultiDimensionLinearFitter.TOL, 0.0);
        
        //set new value
        fitter.setTol(1e-3);
        
        //check correctness
        assertEquals(fitter.getTol(), 1e-3, 0.0);
    }
    
    @Test
    public void testFitQuadratic() throws FittingException, NotReadyException, 
            WrongSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        double[] params = new double[NUM_QUADRATIC_PARAMS];
        for(int i = 0; i < NUM_QUADRATIC_PARAMS; i++){
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
        
        Matrix x = Matrix.createWithUniformRandomValues(npoints, 2, 
                MIN_RANDOM_VALUE, MAX_RANDOM_VALUE, new Random());
        double[] y = new double[npoints];
        GaussianRandomizer errorRandomizer = new GaussianRandomizer(
            new Random(), 0.0, sigma);
        double error, x1, x2;
        for(int i = 0; i < npoints; i++){
            x1 = x.getElementAt(i, 0);
            x2 = x.getElementAt(i, 1);
            //function is: a*1 + b*x + c*y + d*x*y + e*x^2 + f*y^2
            y[i] = params[0] + params[1]*x1 + params[2]*x2 + params[3]*x1*x2 +
                    params[4]*x1*x1 + params[5]*x2*x2;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2; //x, y coordinates
            }
                    
            @Override
            public double[] createResultArray() {
                //parameters a, b, c, d, e, f for function:
                //a + b*x + c*y + d*x*y + e*x*x + f*y*y
                return new double[NUM_QUADRATIC_PARAMS];
            }

            @Override
            public void evaluate(double[] point, double[] result) 
                    throws Throwable {
                double x = point[0];
                double y = point[1];
                
                result[0] = 1.0;
                result[1] = x;
                result[2] = y;
                result[3] = x*y;
                result[4] = x*x;
                result[5] = y*y;
            }
        };
        
        SvdMultiDimensionLinearFitter fitter = 
                new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        //fit
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, NUM_QUADRATIC_PARAMS);
        for(int i = 0; i < NUM_QUADRATIC_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);                  
    }
}
