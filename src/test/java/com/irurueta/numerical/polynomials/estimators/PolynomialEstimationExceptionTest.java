/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.polynomials.estimators;

import org.junit.*;

import static org.junit.Assert.assertNotNull;

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
        PolynomialEstimationException ex = new PolynomialEstimationException();
        assertNotNull(ex);

        ex = new PolynomialEstimationException("message");
        assertNotNull(ex);

        ex = new PolynomialEstimationException(new Exception());
        assertNotNull(ex);

        ex = new PolynomialEstimationException("message",
                new Exception());
        assertNotNull(ex);
    }
}
