/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.RobustEstimatorException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class RobustEstimatorExceptionTest {
    
    public RobustEstimatorExceptionTest() {}
    
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
        RobustEstimatorException ex;
        assertNotNull(ex = new RobustEstimatorException());
        
        ex = null;
        assertNotNull(ex = new RobustEstimatorException("message"));
        
        ex = null;
        assertNotNull(ex = new RobustEstimatorException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new RobustEstimatorException("message", 
                new Exception()));        
    }              
}