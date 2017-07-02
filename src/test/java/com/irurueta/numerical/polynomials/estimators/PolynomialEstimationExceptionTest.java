/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.PolynomialEstimationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 8, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class PolynomialEstimationExceptionTest {
    
    public PolynomialEstimationExceptionTest() { }
    
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
        PolynomialEstimationException ex;
        assertNotNull(ex = new PolynomialEstimationException());
        
        ex = null;
        assertNotNull(ex = new PolynomialEstimationException("message"));
        
        ex = null;
        assertNotNull(ex = new PolynomialEstimationException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new PolynomialEstimationException("message", 
                new Exception()));
    }
}
