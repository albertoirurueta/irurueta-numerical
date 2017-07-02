/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.PolynomialsException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 8, 2016.
 */
package com.irurueta.numerical.polynomials;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class PolynomialsExceptionTest {
    
    public PolynomialsExceptionTest() { }
    
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
        PolynomialsException ex;
        assertNotNull(ex = new PolynomialsException());
        
        ex = null;
        assertNotNull(ex = new PolynomialsException("message"));
        
        ex = null;
        assertNotNull(ex = new PolynomialsException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new PolynomialsException("message", 
                new Exception()));
    }
}
