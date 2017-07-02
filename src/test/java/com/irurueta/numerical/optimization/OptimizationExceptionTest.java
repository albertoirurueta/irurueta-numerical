/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.optimization.OptimizationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 25, 2012
 */
package com.irurueta.numerical.optimization;

import static org.junit.Assert.assertNotNull;
import org.junit.*;

public class OptimizationExceptionTest {
    
    public OptimizationExceptionTest() { }

    @BeforeClass
    public static void setUpClass() throws Exception { }

    @AfterClass
    public static void tearDownClass() throws Exception { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() {
        OptimizationException ex;
        assertNotNull(ex = new OptimizationException());
        
        ex = null;
        assertNotNull(ex = new OptimizationException("message"));
        
        ex = null;
        assertNotNull(ex = new OptimizationException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new OptimizationException("message", 
                new Exception()));        
    }          
}
