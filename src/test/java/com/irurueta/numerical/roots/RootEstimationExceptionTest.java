/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.roots.RootEstimationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.roots;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class RootEstimationExceptionTest {
    
    public RootEstimationExceptionTest() {}
    
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
        RootEstimationException ex;
        assertNotNull(ex = new RootEstimationException());
        
        ex = null;
        assertNotNull(ex = new RootEstimationException("message"));
        
        ex = null;
        assertNotNull(ex = new RootEstimationException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new RootEstimationException("message", 
                new Exception()));        
    }                  
}