/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.StraightLineFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 24, 2015
 */
package com.irurueta.numerical.fitting;

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

public class StraightLineFitterTest {
    public static final double MIN_PARAM_VALUE = -100.0;
    public static final double MAX_PARAM_VALUE = 100.0;
    
    public static final double MIN_SIGMA_VALUE = 1e-3;
    public static final double MAX_SIGMA_VALUE = 3.0;
    
    public static final int MIN_POINTS = 1000;
    public static final int MAX_POINTS = 10000;
    
    public static final double MIN_DATA_VALUE = -100.0;
    public static final double MAX_DATA_VALUE = 100.0;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    public StraightLineFitterTest() {}
    
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
        //test constructor without arguments
        StraightLineFitter fitter = new StraightLineFitter();
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(fitter.getA(), 0.0, 0.0);
        assertEquals(fitter.getB(), 0.0, 0.0);
        assertEquals(fitter.getSigA(), 0.0, 0.0);
        assertEquals(fitter.getSigB(), 0.0, 0.0);
        assertEquals(fitter.getChi2(), 0.0, 0.0);
        assertEquals(fitter.getQ(), 1.0, 0.0);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        
        //test constructor with input data
        double[] x = new double[2];
        double[] y = new double[2];
        
        fitter = new StraightLineFitter(x, y);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(fitter.getA(), 0.0, 0.0);
        assertEquals(fitter.getB(), 0.0, 0.0);
        assertEquals(fitter.getSigA(), 0.0, 0.0);
        assertEquals(fitter.getSigB(), 0.0, 0.0);
        assertEquals(fitter.getChi2(), 0.0, 0.0);
        assertEquals(fitter.getQ(), 1.0, 0.0);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        
        //Force IllegalArgumentException
        double[] shortX = new double[1];
        
        fitter = null;
        try{
            fitter = new StraightLineFitter(shortX, y);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with input data and sigmas
        double[] sig = new double[2];
        fitter = new StraightLineFitter(x, y, sig);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(fitter.getA(), 0.0, 0.0);
        assertEquals(fitter.getB(), 0.0, 0.0);
        assertEquals(fitter.getSigA(), 0.0, 0.0);
        assertEquals(fitter.getSigB(), 0.0, 0.0);
        assertEquals(fitter.getChi2(), 0.0, 0.0);
        assertEquals(fitter.getQ(), 1.0, 0.0);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        
        //Force IllegalArgumentException
        double[] shortSig = new double[1];
        
        fitter = null;
        try{
            fitter = new StraightLineFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
    }
    
    @Test
    public void testGetSetInputDataAndIsReady(){
        StraightLineFitter fitter = new StraightLineFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        
        //set new values
        double[] x = new double[2];
        double[] y = new double[2];
        
        fitter.setInputData(x, y);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        
        //Force IllegalArgumentException
        double[] shortX = new double[1];
        
        try{
            fitter.setInputData(shortX, y);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetInputDataAndStandardDeviations(){
        StraightLineFitter fitter = new StraightLineFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        
        //set new values
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];
        
        fitter.setInputDataAndStandardDeviations(x, y, sig);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        
        fitter.setInputDataAndStandardDeviations(x, y, null);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);        
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        
        //Force IllegalArgumentException
        double[] shortX = new double[1];
        double[] shortSig = new double[1];
        try{
            fitter.setInputDataAndStandardDeviations(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputDataAndStandardDeviations(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }

    @Test
    public void testFitNoSig() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        double b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);                
        
        GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        
        
        //generate data
        int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double[] x = new double[nPoints];
        double[] y = new double[nPoints];
        for(int i = 0; i < nPoints; i++){
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
        }
        
        StraightLineFitter fitter = new StraightLineFitter(x, y);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(fitter.getA(), 0.0, 0.0);
        assertEquals(fitter.getB(), 0.0, 0.0);
        assertEquals(fitter.getSigA(), 0.0, 0.0);
        assertEquals(fitter.getSigB(), 0.0, 0.0);
        assertEquals(fitter.getChi2(), 0.0, 0.0);
        assertEquals(fitter.getQ(), 1.0, 0.0);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        
        //fit data
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertEquals(fitter.getA(), a, ABSOLUTE_ERROR);
        assertEquals(fitter.getB(), b, ABSOLUTE_ERROR);
        assertEquals(fitter.getSigdat(), sigma, ABSOLUTE_ERROR);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(), 
                fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);
        
        //Force NotReadyException
        fitter = new StraightLineFitter();
        
        assertFalse(fitter.isReady());
        try{
            fitter.fit();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
    }
    
    @Test
    public void testFitSig() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        double b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);                
        
        GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        
        
        //generate data
        int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double[] x = new double[nPoints];
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        for(int i = 0; i < nPoints; i++){
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
            sig[i] = sigma;
        }
        
        StraightLineFitter fitter = new StraightLineFitter(x, y, sig);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(fitter.getA(), 0.0, 0.0);
        assertEquals(fitter.getB(), 0.0, 0.0);
        assertEquals(fitter.getSigA(), 0.0, 0.0);
        assertEquals(fitter.getSigB(), 0.0, 0.0);
        assertEquals(fitter.getChi2(), 0.0, 0.0);
        assertEquals(fitter.getQ(), 1.0, 0.0);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        
        //fit data
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertEquals(fitter.getA(), a, ABSOLUTE_ERROR);
        assertEquals(fitter.getB(), b, ABSOLUTE_ERROR);
        assertEquals(fitter.getSigdat(), 0.0, 0.0);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(), 
                fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);
        
        //Force NotReadyException
        fitter = new StraightLineFitter();
        
        assertFalse(fitter.isReady());
        try{
            fitter.fit();
            fail("NotReadyException expected but not thrown");
        }catch(NotReadyException e){}
    }    
}
