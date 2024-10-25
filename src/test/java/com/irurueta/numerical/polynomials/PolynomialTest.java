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
package com.irurueta.numerical.polynomials;

import com.irurueta.algebra.ArrayUtils;
import com.irurueta.numerical.NumericalException;
import com.irurueta.numerical.SerializationHelper;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class PolynomialTest {

    private static final int MIN_LENGTH = 1;
    private static final int MAX_LENGTH = 10;

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final double ABSOLUTE_ERROR = 1e-8;
    private static final double LARGE_ABSOLUTE_ERROR = 1e-5;
    private static final double VERY_LARGE_ABSOLUTE_ERROR = 1e-3;

    private static final int TIMES = 50;

    @Test
    void testConstructor() throws NumericalException {
        // test empty constructor
        var p = new Polynomial();

        // check correctness
        assertArrayEquals(new double[1], p.getPolyParams(), 0.0);
        assertEquals(0, p.getDegree());
        assertNull(p.getRoots());
        assertNull(p.getMaxima());
        assertNull(p.getMaxima(1.0));
        assertNull(p.getMinima());
        assertNull(p.getMinima(1.0));
        assertNull(p.getExtrema());
        assertNull(p.getExtrema(1.0));

        // test constructor with number of parameters
        p = new Polynomial(2);

        // check correctness
        assertArrayEquals(new double[2], p.getPolyParams(), 0.0);
        assertEquals(0, p.getDegree());
        assertNull(p.getRoots());
        assertNull(p.getMaxima());
        assertNull(p.getMaxima(1.0));
        assertNull(p.getMinima());
        assertNull(p.getMinima(1.0));
        assertNull(p.getExtrema());
        assertNull(p.getExtrema(1.0));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new Polynomial(0));

        // test constructor with array of polynomial parameters
        final var polyParams = new double[1];
        p = new Polynomial(polyParams);

        // check correctness
        assertSame(p.getPolyParams(), polyParams);
        assertEquals(0, p.getDegree());
        assertNull(p.getRoots());
        assertNull(p.getMaxima());
        assertNull(p.getMaxima(1.0));
        assertNull(p.getMinima());
        assertNull(p.getMinima(1.0));
        assertNull(p.getExtrema());
        assertNull(p.getExtrema(1.0));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new Polynomial(new double[0]));
    }

    @Test
    void testGetSetPolyParams() {
        final var polyParams = new double[2];

        final var p = new Polynomial(polyParams);

        // check correctness
        assertSame(p.getPolyParams(), polyParams);

        // set new value
        final var polyParams2 = new double[3];
        p.setPolyParams(polyParams2);

        // check correctness
        assertSame(polyParams2, p.getPolyParams());

        // Force IllegalArgumentException
        final var wrong = new double[0];
        assertThrows(IllegalArgumentException.class, () -> p.setPolyParams(wrong));
    }

    @Test
    void testGetDegree() {
        final var p = new Polynomial(1.0);
        assertEquals(0, p.getDegree());

        p.setPolyParams(1.0, 0.0);
        assertEquals(0, p.getDegree());

        p.setPolyParams(1.0, 1.0);
        assertEquals(1, p.getDegree());

        p.setPolyParams(1.0, 1.0, 0.0);
        assertEquals(1, p.getDegree());

        p.setPolyParams(1.0, 1.0, 1.0);
        assertEquals(2, p.getDegree());

        p.setPolyParams(1.0, 1.0, 1.0, 0.0);
        assertEquals(2, p.getDegree());
    }

    @Test
    void testAddWithResult() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        var p3 = new Polynomial(2);
        p1.add(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0}, p3.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p3 = new Polynomial();
        p1.add(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{5.0, 7.0, 3.0}, p3.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p3 = new Polynomial(4);
        p1.add(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0, 5.0}, p3.getPolyParams(), 0.0);
    }

    @Test
    void testAdd() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        p1.add(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0}, p1.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p1.add(p2);

        // check correctness
        assertArrayEquals(new double[]{5.0, 7.0, 3.0}, p1.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p1.add(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0, 5.0}, p1.getPolyParams(), 0.0);
    }

    @Test
    void testAddAndReturnNew() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        var p3 = p1.addAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0}, p3.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p3 = p1.addAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{5.0, 7.0, 3.0}, p3.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p3 = p1.addAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 6.0, 5.0}, p3.getPolyParams(), 0.0);
    }

    @Test
    void testSubtractWithResult() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        var p3 = new Polynomial(2);
        p1.subtract(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0}, p3.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p3 = new Polynomial();
        p1.subtract(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-3.0, -3.0, 3.0}, p3.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p3 = new Polynomial(4);
        p1.subtract(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0, -5.0}, p3.getPolyParams(), 0.0);
    }

    @Test
    void testSubtract() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        p1.subtract(p2);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0}, p1.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p1.subtract(p2);

        // check correctness
        assertArrayEquals(new double[]{-3.0, -3.0, 3.0}, p1.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p1.subtract(p2);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0, -5.0}, p1.getPolyParams(), 0.0);
    }

    @Test
    void testSubtractAndReturnNew() {
        // test equal length
        var p1 = new Polynomial(1.0, 2.0);
        var p2 = new Polynomial(3.0, 4.0);

        var p3 = p1.subtractAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0}, p3.getPolyParams(), 0.0);

        // test p2 length < p1 length
        p1 = new Polynomial(1.0, 2.0, 3.0);
        p2 = new Polynomial(4.0, 5.0);

        p3 = p1.subtractAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-3.0, -3.0, 3.0}, p3.getPolyParams(), 0.0);

        // test p2 length > p1 length
        p1 = new Polynomial(1.0, 2.0);
        p2 = new Polynomial(3.0, 4.0, 5.0);

        p3 = p1.subtractAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-2.0, -2.0, -5.0}, p3.getPolyParams(), 0.0);
    }

    @Test
    void testMultiplyWithResult() {
        // (x + 1) * (x - 1) = x^2 - 1 --> (1,1)*(-1,1) = (-1, 0, 1)
        var p1 = new Polynomial(1.0, 1.0);
        var p2 = new Polynomial(-1.0, 1.0);

        var p3 = new Polynomial(2);
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-1.0, 0.0, 1.0}, p3.getPolyParams(), 0.0);

        // (x - 1) * (x - 1) = x^2 - 2*x + 1 --> (-1,1)*(-1,1) = (1,-2,1)
        p1 = new Polynomial(-1.0, 1.0);
        p2 = new Polynomial(-1.0, 1.0);

        p3 = new Polynomial(2);
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{1.0, -2.0, 1.0}, p3.getPolyParams(), 0.0);

        // (x + 1) * 2 = 2*x + 2 --> (1,1)*(2) = (2,2)
        p1 = new Polynomial(1.0, 1.0);
        p2 = new Polynomial(2.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{2.0, 2.0}, p3.getPolyParams(), 0.0);


        // (x + 3)* (x + 2) = x^2 + 5*x + 6 --> (3,1)*(2,1) = (6,5,1)
        p1 = new Polynomial(3.0, 1.0);
        p2 = new Polynomial(2.0, 1.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{6.0, 5.0, 1.0}, p3.getPolyParams(), 0.0);


        // (4*x^2 - 4*x -7)*(x + 3) = 4*x^3 + 8*x^2 - 19*x - 21 -->
        // (-7,-4,4)*(3,1)=(-21,-19,8,4)
        p1 = new Polynomial(-7.0, -4.0, 4.0);
        p2 = new Polynomial(3.0, 1.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-21.0, -19.0, 8.0, 4.0}, p3.getPolyParams(), 0.0);


        // (x+2)(x^3 + 3*x^2 + 4*x - 17) = x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        // (2, 1)*(-17, 4, 3, 1) = (-34, -9, 10, 5, 1)
        p1 = new Polynomial(2.0, 1.0);
        p2 = new Polynomial(-17.0, 4.0, 3.0, 1.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, p3.getPolyParams(), 0.0);

        // (3*x^2 - 9*x + 5)(2*x^ + 4*x -7) = 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        // (5, -9, 3)*(-7, 4, 2) = (-35, 83, -47, -6, 6)
        p1 = new Polynomial(5.0, -9.0, 3.0);
        p2 = new Polynomial(-7.0, 4.0, 2.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, p3.getPolyParams(), 0.0);

        // (x^3 + 2*x^2 + 4)(2*x^3 + x + 1) = 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        // (4, 0, 2, 1)*(1, 1, 0, 2) = (4, 4, 2, 11, 1, 4, 2)
        p1 = new Polynomial(4.0, 0.0, 2.0, 1.0);
        p2 = new Polynomial(1.0, 1.0, 0.0, 2.0);

        p3 = new Polynomial();
        p1.multiply(p2, p3);

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, p3.getPolyParams(), 0.0);
    }

    @Test
    void testMultiply() {
        // (x + 1) * (x - 1) = x^2 - 1 --> (1,1)*(-1,1) = (-1, 0, 1)
        var p1 = new Polynomial(1.0, 1.0);
        var p2 = new Polynomial(-1.0, 1.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{-1.0, 0.0, 1.0}, p1.getPolyParams(), 0.0);

        // (x - 1) * (x - 1) = x^2 - 2*X + 1 --> (-1,1)*(-1,1) = (1,-2,1)
        p1 = new Polynomial(-1.0, 1.0);
        p2 = new Polynomial(-1.0, 1.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{1.0, -2.0, 1.0}, p1.getPolyParams(), 0.0);

        // (x + 1) * 2 = 2*x + 2 --> (1,1)*(2) = (2,2)
        p1 = new Polynomial(1.0, 1.0);
        p2 = new Polynomial(2.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{2.0, 2.0}, p1.getPolyParams(), 0.0);


        // (x + 3)* (x + 2) = x^2 + 5*X + 6 --> (3,1)*(2,1) = (6,5,1)
        p1 = new Polynomial(3.0, 1.0);
        p2 = new Polynomial(2.0, 1.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{6.0, 5.0, 1.0}, p1.getPolyParams(), 0.0);


        // (4*x^2 - 4*x -7)*(x + 3) = 4*x^3 + 8*x2 - 19*x -21 -->
        // (-7,-4,4)*(3,1)=(-21,-19,8,4)
        p1 = new Polynomial(-7.0, -4.0, 4.0);
        p2 = new Polynomial(3.0, 1.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{-21.0, -19.0, 8.0, 4.0}, p1.getPolyParams(), 0.0);


        // (x+2)(x^3 + 3*x^2 + 4*x - 17) = x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        // (2, 1)*(-17, 4, 3, 1) = (-34, -9, 10, 5, 1)
        p1 = new Polynomial(2.0, 1.0);
        p2 = new Polynomial(-17.0, 4.0, 3.0, 1.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, p1.getPolyParams(), 0.0);

        // (3*x^2 - 9*x + 5)(2*x^ + 4*x -7) = 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        // (5, -9, 3)*(-7, 4, 2) = (-35, 83, -47, -6, 6)
        p1 = new Polynomial(5.0, -9.0, 3.0);
        p2 = new Polynomial(-7.0, 4.0, 2.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, p1.getPolyParams(), 0.0);

        // (x^3 + 2*x^2 + 4)(2*x^3 + x + 1) = 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        // (4, 0, 2, 1)*(1, 1, 0, 2) = (4, 4, 2, 11, 1, 4, 2)
        p1 = new Polynomial(4.0, 0.0, 2.0, 1.0);
        p2 = new Polynomial(1.0, 1.0, 0.0, 2.0);

        p1.multiply(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, p1.getPolyParams(), 0.0);
    }

    @Test
    void testMultiplyAndReturnNew() {
        // (x + 1) * (x - 1) = x^2 - 1 --> (1,1)*(-1,1) = (-1, 0, 1)
        var p1 = new Polynomial(1.0, 1.0);
        var p2 = new Polynomial(-1.0, 1.0);

        var p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-1.0, 0.0, 1.0}, p3.getPolyParams(), 0.0);

        // (x - 1) * (x - 1) = x^2 - 2*X + 1 --> (-1,1)*(-1,1) = (1,-2,1)
        p1 = new Polynomial(-1.0, 1.0);
        p2 = new Polynomial(-1.0, 1.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{1.0, -2.0, 1.0}, p3.getPolyParams(), 0.0);

        // (x + 1) * 2 = 2*x + 2 --> (1,1)*(2) = (2,2)
        p1 = new Polynomial(1.0, 1.0);
        p2 = new Polynomial(2.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{2.0, 2.0}, p3.getPolyParams(), 0.0);


        // (x + 3)* (x + 2) = x^2 + 5*X + 6 --> (3,1)*(2,1) = (6,5,1)
        p1 = new Polynomial(3.0, 1.0);
        p2 = new Polynomial(2.0, 1.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{6.0, 5.0, 1.0}, p3.getPolyParams(), 0.0);


        // (4*x^2 - 4*x -7)*(x + 3) = 4*x^3 + 8*x2 - 19*x -21 -->
        // (-7,-4,4)*(3,1)=(-21,-19,8,4)
        p1 = new Polynomial(-7.0, -4.0, 4.0);
        p2 = new Polynomial(3.0, 1.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-21.0, -19.0, 8.0, 4.0}, p3.getPolyParams(), 0.0);


        // (x+2)(x^3 + 3*x^2 + 4*x - 17) = x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        // (2, 1)*(-17, 4, 3, 1) = (-34, -9, 10, 5, 1)
        p1 = new Polynomial(2.0, 1.0);
        p2 = new Polynomial(-17.0, 4.0, 3.0, 1.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, p3.getPolyParams(), 0.0);

        // (3*x^2 - 9*x + 5)(2*x^ + 4*x -7) = 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        // (5, -9, 3)*(-7, 4, 2) = (-35, 83, -47, -6, 6)
        p1 = new Polynomial(5.0, -9.0, 3.0);
        p2 = new Polynomial(-7.0, 4.0, 2.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, p3.getPolyParams(), 0.0);

        // (x^3 + 2*x^2 + 4)(2*x^3 + x + 1) = 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        // (4, 0, 2, 1)*(1, 1, 0, 2) = (4, 4, 2, 11, 1, 4, 2)
        p1 = new Polynomial(4.0, 0.0, 2.0, 1.0);
        p2 = new Polynomial(1.0, 1.0, 0.0, 2.0);

        p3 = p1.multiplyAndReturnNew(p2);

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0},
                p3.getPolyParams(), 0.0);
    }

    @Test
    void testMultiplyByScalarWithResult() {
        final var randomizer = new UniformRandomizer();

        final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final var polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var p1 = new Polynomial(polyParams);

        final var scalar = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var result = ArrayUtils.multiplyByScalarAndReturnNew(polyParams, scalar);

        final var p2 = new Polynomial();

        p1.multiplyByScalar(scalar, p2);

        // check correctness
        assertArrayEquals(result, p2.getPolyParams(), 0.0);
    }

    @Test
    void testMultiplyByScalar() {
        final var randomizer = new UniformRandomizer();

        final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final var polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var p1 = new Polynomial(polyParams);

        final var scalar = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var result = ArrayUtils.multiplyByScalarAndReturnNew(polyParams, scalar);

        p1.multiplyByScalar(scalar);

        // check correctness
        assertArrayEquals(result, p1.getPolyParams(), 0.0);
    }

    @Test
    void testMultiplyByScalarAndReturnNew() {
        final var randomizer = new UniformRandomizer();

        final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final var polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var p1 = new Polynomial(polyParams);

        final var scalar = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var result = ArrayUtils.multiplyByScalarAndReturnNew(polyParams, scalar);

        final var p2 = p1.multiplyByScalarAndReturnNew(scalar);

        // check correctness
        assertArrayEquals(result, p2.getPolyParams(), 0.0);
    }

    @Test
    void testGetRoots() throws NumericalException {
        var numValid = 0;
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // test degree 0
            var p = new Polynomial(2.0);

            assertNull(p.getRoots());


            // test degree 1
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            p = new Polynomial(-root1, 1.0);

            var roots = p.getRoots();

            // check correctness
            assertEquals(1, roots.length);
            assertEquals(root1, roots[0].getReal(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[0].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[0].getReal()), ABSOLUTE_ERROR);


            // test degree 2
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            roots = p.getRoots();

            // check correctness
            assertEquals(2, roots.length);
            assertTrue(Math.abs(roots[0].getReal() - root1) <= ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root2) <= ABSOLUTE_ERROR);
            assertTrue(Math.abs(roots[1].getReal() - root1) <= ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root2) <= ABSOLUTE_ERROR);
            assertEquals(0.0, roots[0].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[1].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[0].getReal()), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[1].getReal()), ABSOLUTE_ERROR);


            // test degree 3
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);
            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3);

            roots = p.getRoots();

            // check correctness
            assertEquals(3, roots.length);
            assertTrue(Math.abs(roots[0].getReal() - root1) <= ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root2) <= ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root3) <= ABSOLUTE_ERROR);
            assertTrue(Math.abs(roots[1].getReal() - root1) <= ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root2) <= ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root3) <= ABSOLUTE_ERROR);
            assertTrue(Math.abs(roots[2].getReal() - root1) <= ABSOLUTE_ERROR
                    || Math.abs(roots[2].getReal() - root2) <= ABSOLUTE_ERROR
                    || Math.abs(roots[2].getReal() - root3) <= ABSOLUTE_ERROR);
            assertEquals(0.0, roots[0].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[1].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[2].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[0].getReal()), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[1].getReal()), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[2].getReal()), ABSOLUTE_ERROR);

            // test degree 4
            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);
            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            roots = p.getRoots();

            // check correctness
            assertEquals(4, roots.length);
            final var cond1 = Math.abs(roots[0].getReal() - root1) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root2) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root3) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[0].getReal() - root4) <= LARGE_ABSOLUTE_ERROR;
            if (!cond1) {
                continue;
            }

            final var cond2 = Math.abs(roots[1].getReal() - root1) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root2) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root3) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[1].getReal() - root4) <= LARGE_ABSOLUTE_ERROR;
            if (!cond2) {
                continue;
            }

            final var cond3 = Math.abs(roots[2].getReal() - root1) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[2].getReal() - root2) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[2].getReal() - root3) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[2].getReal() - root4) <= LARGE_ABSOLUTE_ERROR;
            if (!cond3) {
                continue;
            }

            final var cond4 = Math.abs(roots[3].getReal() - root1) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[3].getReal() - root2) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[3].getReal() - root3) <= LARGE_ABSOLUTE_ERROR
                    || Math.abs(roots[3].getReal() - root4) <= LARGE_ABSOLUTE_ERROR;
            if (!cond4) {
                continue;
            }

            assertEquals(0.0, roots[0].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[1].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[2].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, roots[3].getImaginary(), ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[0].getReal()), VERY_LARGE_ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[1].getReal()), VERY_LARGE_ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[2].getReal()), VERY_LARGE_ABSOLUTE_ERROR);
            assertEquals(0.0, p.evaluate(roots[3].getReal()), VERY_LARGE_ABSOLUTE_ERROR);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testEvaluate() {
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(x * x - 1, p.evaluate(x), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(x * x - 2.0 * x + 1.0, p.evaluate(x), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(2.0 * x + 2.0, p.evaluate(x), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(x * x + 5 * x + 6, p.evaluate(x), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(4 * x * x * x + 8 * x * x - 19 * x - 21, p.evaluate(x), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(x * x * x * x + 5 * x * x * x + 10 * x * x - 9 * x - 34, p.evaluate(x), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(6 * x * x * x * x - 6 * x * x * x - 47 * x * x + 83 * x - 35, p.evaluate(x),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(2 * x * x * x * x * x * x + 4 * x * x * x * x * x + x * x * x * x + 11 * x * x * x
                + 2 * x * x + 4 * x + 4, p.evaluate(x), ABSOLUTE_ERROR);
    }

    @Test
    void testDerivativeWithResult() {
        // 5 --> 0
        var p = new Polynomial(5.0);
        var d = new Polynomial();
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{0.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1 --> 2*x (0, 2)
        p = new Polynomial(-1.0, 0.0, 1.0);
        d = new Polynomial();
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{0.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1 --> 2*x - 2 (-2, 2)
        p = new Polynomial(1.0, -2.0, 1.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{-2.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2 --> 2
        p = new Polynomial(2.0, 2.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{2.0}, d.getPolyParams(), ABSOLUTE_ERROR);


        // x^2 + 5*x + 6 --> 2*x + 5 (5, 2)
        p = new Polynomial(6.0, 5.0, 1.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{5.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21 --> 12*x^2 + 16*x - 19 (-19, 16, 12)
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{-19.0, 16.0, 12.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34 --> 4*x^3 + 15*x^2 + 20*x - 9
        // (-9, 20, 15, 4)
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{-9.0, 20.0, 15.0, 4.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35 --> 24*x^3 - 18*x^2 - 94*x + 83
        // (83, -94, -18, 24)
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{83.0, -94.0, -18.0, 24.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4 -->
        // 12*x^5 + 20*x^4 + 4*x^3 + 33*x^2 + 4*x + 4 (4, 4, 33, 4, 20, 12)
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.derivative(d);

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 33.0, 4.0, 20.0, 12.0}, d.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testDerivative() {
        // 5 --> 0
        var p = new Polynomial(5.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{0.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1 --> 2*x (0, 2)
        p = new Polynomial(-1.0, 0.0, 1.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{0.0, 2.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1 --> 2*x - 2 (-2, 2)
        p = new Polynomial(1.0, -2.0, 1.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{-2.0, 2.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2 --> 2
        p = new Polynomial(2.0, 2.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{2.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6 --> 2*x + 5 (5, 2)
        p = new Polynomial(6.0, 5.0, 1.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{5.0, 2.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21 --> 12*x^2 + 16*x - 19 (-19, 16, 12)
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{-19.0, 16.0, 12.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34 --> 4*x^3 + 15*x^2 + 20*x - 9
        // (-9, 20, 15, 4)
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{-9.0, 20.0, 15.0, 4.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35 --> 24*x^3 - 18*x^2 - 94*x + 83
        // (83, -94, -18, 24)
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{83.0, -94.0, -18.0, 24.0}, p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4 -->
        // 12*x^5 + 20*x^4 + 4*x^3 + 33*x^2 + 4*x + 4 (4, 4, 33, 4, 20, 12)
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.derivative();

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 33.0, 4.0, 20.0, 12.0}, p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testDerivativeAndReturnNew() {
        // 5 --> 0
        var p = new Polynomial(5.0);
        var d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{0.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1 --> 2*x (0, 2)
        p = new Polynomial(-1.0, 0.0, 1.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{0.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1 --> 2*x - 2 (-2, 2)
        p = new Polynomial(1.0, -2.0, 1.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{-2.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2 --> 2
        p = new Polynomial(2.0, 2.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6 --> 2*x + 5 (5, 2)
        p = new Polynomial(6.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{5.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21 --> 12*x^2 + 16*x - 19 (-19, 16, 12)
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{-19.0, 16.0, 12.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34 --> 4*x^3 + 15*x^2 + 20*x - 9
        // (-9, 20, 15, 4)
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{-9.0, 20.0, 15.0, 4.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35 --> 24*x^3 - 18*x^2 - 94*x + 83
        // (83, -94, -18, 24)
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{83.0, -94.0, -18.0, 24.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4 -->
        // 12*x^5 + 20*x^4 + 4*x^3 + 33*x^2 + 4*x + 4 (4, 4, 33, 4, 20, 12)
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d = p.derivativeAndReturnNew();

        // check correctness
        assertArrayEquals(new double[]{4.0, 4.0, 33.0, 4.0, 20.0, 12.0}, d.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testEvaluateDerivative() {
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(p.derivativeAndReturnNew().evaluate(x), p.evaluateDerivative(x), ABSOLUTE_ERROR);
    }

    @Test
    void testSecondDerivativeWithResult() {
        final var dd = new Polynomial();

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative(dd);
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testSecondDerivative() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d = p.derivativeAndReturnNew();
        p.secondDerivative();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testSecondDerivativeAndReturnNew() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var d = p.derivativeAndReturnNew();
        var dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d = p.derivativeAndReturnNew();
        dd = p.secondDerivativeAndReturnNew();
        assertArrayEquals(d.derivativeAndReturnNew().getPolyParams(), dd.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testEvaluateSecondDerivative() {
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(p.secondDerivativeAndReturnNew().evaluate(x), p.evaluateSecondDerivative(x), ABSOLUTE_ERROR);
    }

    @Test
    void testNthDerivativeWithResult() {
        // test first order derivatives

        // 5
        var p = new Polynomial(5.0);
        var d1 = p.derivativeAndReturnNew();
        final var d2 = new Polynomial();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1, d2);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // test second order derivatives

        // 5
        p = new Polynomial(5.0);
        var dd1 = p.secondDerivativeAndReturnNew();
        final var dd2 = new Polynomial();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2, dd2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order derivatives

        // 5
        p = new Polynomial(5.0);
        var ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        final var ddd2 = new Polynomial();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3, ddd2);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        final var p2 = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertThrows(IllegalArgumentException.class, () -> p2.nthDerivative(0, ddd2));
    }

    @Test
    void testNthDerivative() {
        // test first order derivatives

        // 5
        var p = new Polynomial(5.0);
        var d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        p.nthDerivative(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test second order derivatives

        // 5
        p = new Polynomial(5.0);
        var dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        p.nthDerivative(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order derivatives

        // 5
        p = new Polynomial(5.0);
        var ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        p.nthDerivative(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        final var p2 = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertThrows(IllegalArgumentException.class, () -> p2.nthDerivative(0));
    }

    @Test
    void testNthDerivativeAndReturnNew() {
        // test first order derivatives

        // 5
        var p = new Polynomial(5.0);
        var d1 = p.derivativeAndReturnNew();
        var d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        d1 = p.derivativeAndReturnNew();
        d2 = p.nthDerivativeAndReturnNew(1);

        // check correctness
        assertArrayEquals(d1.getPolyParams(), d2.getPolyParams(), ABSOLUTE_ERROR);

        // test second order derivatives

        // 5
        p = new Polynomial(5.0);
        Polynomial dd1 = p.secondDerivativeAndReturnNew();
        Polynomial dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        dd1 = p.secondDerivativeAndReturnNew();
        dd2 = p.nthDerivativeAndReturnNew(2);

        // check correctness
        assertArrayEquals(dd1.getPolyParams(), dd2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order derivatives

        // 5
        p = new Polynomial(5.0);
        var ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        var ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35 
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        ddd1 = p.secondDerivativeAndReturnNew().derivativeAndReturnNew();
        ddd2 = p.nthDerivativeAndReturnNew(3);

        // check correctness
        assertArrayEquals(ddd1.getPolyParams(), ddd2.getPolyParams(), ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        final var p2 = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertThrows(IllegalArgumentException.class, () -> p2.nthDerivativeAndReturnNew(0));
    }

    @Test
    void testEvaluateNthDerivative() {
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // test first order derivatives

        // 5
        var p = new Polynomial(5.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(p.evaluateDerivative(x), p.evaluateNthDerivative(x, 1), ABSOLUTE_ERROR);

        // test second order derivatives

        // 5
        p = new Polynomial(5.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(p.evaluateSecondDerivative(x), p.evaluateNthDerivative(x, 2), ABSOLUTE_ERROR);

        // test 3rd order derivatives

        // 5
        p = new Polynomial(5.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        assertEquals(p.nthDerivativeAndReturnNew(3).evaluate(x), p.evaluateNthDerivative(x, 3),
                ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrationWithResultAndConstant() {
        final var randomizer = new UniformRandomizer();
        final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var i = new Polynomial();

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.integration(i, constant);
        var id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.integration(i, constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrationWithResult() {

        final var i = new Polynomial();

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.integration(i);
        var id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.integration(i);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrationWithConstant() {
        final var randomizer = new UniformRandomizer();
        final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.integration(constant);
        var d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-1.0, 0.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(d.getPolyParams(), new double[]{1.0, -2.0, 1.0}, ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{2.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{6.0, 5.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-21.0, -19.0, 8.0, 4.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.integration(constant);
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegration() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.integration();
        var d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-1.0, 0.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{1.0, -2.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{2.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{6.0, 5.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-21.0, -19.0, 8.0, 4.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, d.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.integration();
        d = p.derivativeAndReturnNew();
        assertArrayEquals(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, d.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrationAndReturnNewWithConstant() {
        final var randomizer = new UniformRandomizer();
        final double constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i = p.integrationAndReturnNew(constant);
        var id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.integrationAndReturnNew(constant);
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrationAndReturnNew() {

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i = p.integrationAndReturnNew();
        var id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.integrationAndReturnNew();
        id = i.derivativeAndReturnNew();
        assertArrayEquals(p.getPolyParams(), id.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testIntegrateInterval() {
        final var randomizer = new UniformRandomizer();
        final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var endX = randomizer.nextDouble(startX, startX + MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.integrationAndReturnNew();
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.integrateInterval(startX, endX), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegrationWithResultAndConstant() {
        final var randomizer = new UniformRandomizer();

        // test 1st order
        var order = 1;
        var constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        var i2 = new Polynomial();

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);


        // test 2nd order
        order = 2;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, i2, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegrationWithResult() {
        // test 1st order
        var order = 1;

        var i2 = new Polynomial();

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;

        i2 = new Polynomial();

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, i2);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegrationWithConstant() {
        final var randomizer = new UniformRandomizer();

        // test 1st order
        var order = 1;
        var constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        p.nthIntegration(order, constants);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegration() {
        // test 1st order
        var order = 1;

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order, p);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        p.nthIntegration(order);
        assertArrayEquals(i1.getPolyParams(), p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegrationAndReturnNewWithConstant() {
        final var randomizer = new UniformRandomizer();

        // test 1st order
        var order = 1;
        var constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew(constants[0]);
        var i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[2]).integrationAndReturnNew(constants[1]).integrationAndReturnNew(
                constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew(constants[3]).integrationAndReturnNew(constants[2]).integrationAndReturnNew(
                constants[1]).integrationAndReturnNew(constants[0]);
        i2 = p.nthIntegrationAndReturnNew(order, constants);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthIntegrationAndReturnNew() {
        // test 1st order
        var order = 1;

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i1 = p.integrationAndReturnNew();
        var i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i1 = p.integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew().integrationAndReturnNew();
        i2 = p.nthIntegrationAndReturnNew(order);
        assertArrayEquals(i1.getPolyParams(), i2.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNthOrderIntegrateIntervalWithConstants() {
        final var randomizer = new UniformRandomizer();
        final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var endX = randomizer.nextDouble(startX, startX + MAX_RANDOM_VALUE);

        // 1st order integration
        var order = 1;
        var constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // test 4th order
        order = 4;
        constants = new double[order];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order, constants);
        assertEquals(i.evaluate(endX) - i.evaluate(startX),
                p.nthOrderIntegrateInterval(startX, endX, order, constants), ABSOLUTE_ERROR);
    }

    @Test
    void testNthOrderIntegrateInterval() {
        final var randomizer = new UniformRandomizer();
        final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var endX = randomizer.nextDouble(startX, startX + MAX_RANDOM_VALUE);

        // 1st order integration
        var order = 1;

        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // test 2nd order
        order = 2;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // test 3rd order
        order = 3;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // test 4th order
        order = 4;

        // x^2 - 1
        p = new Polynomial(-1.0, 0.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        i = p.nthIntegrationAndReturnNew(order);
        assertEquals(i.evaluate(endX) - i.evaluate(startX), p.nthOrderIntegrateInterval(startX, endX, order),
                ABSOLUTE_ERROR);
    }

    @Test
    void testTrim() {
        final var p = new Polynomial(1.0, 2.0, 0.0);

        assertArrayEquals(new double[]{1.0, 2.0, 0.0}, p.getPolyParams(), 0.0);

        // trim
        var result = new Polynomial();
        p.trim(result);

        // check
        assertArrayEquals(new double[]{1.0, 2.0}, result.getPolyParams(), 0.0);

        // trim again
        p.trim(result);

        // check
        assertArrayEquals(new double[]{1.0, 2.0}, result.getPolyParams(), 0.0);

        // trim and return new
        result = p.trimAndReturnNew();

        // check
        assertArrayEquals(new double[]{1.0, 2.0}, result.getPolyParams(), 0.0);

        // trim this instance
        p.trim();

        // check
        assertArrayEquals(new double[]{1.0, 2.0}, p.getPolyParams(), 0.0);
    }

    @Test
    void testNormalizeWithResult() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        final var result = new Polynomial();
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-1.0, 0.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{1.0, -2.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{2.0, 2.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{6.0, 5.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-21.0, -19.0, 8.0, 4.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.normalize(result);

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNormalize() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-1.0, 0.0, 1.0}), p.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{1.0, -2.0, 1.0}), p.getPolyParams(),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{2.0, 2.0}), p.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{6.0, 5.0, 1.0}), p.getPolyParams(),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-21.0, -19.0, 8.0, 4.0}), p.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.normalize();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}),
                p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNormalizeAndReturnNew() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-1.0, 0.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{1.0, -2.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{2.0, 2.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{6.0, 5.0, 1.0}), result.getPolyParams(),
                ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-21.0, -19.0, 8.0, 4.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        result = p.normalizeAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.normalizeAndReturnNew(new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}),
                result.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNormalizeHighestDegreeTermWithResult() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var result = new Polynomial();
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-1.0, 0.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{1.0, -2.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{2.0, 2.0}, 1.0 / 2.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{6.0, 5.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-21.0, -19.0, 8.0, 4.0}, 1.0 / 4.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(
                ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, 1.0 / 6.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.normalizeHighestDegreeTerm(result);

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(
                new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, 1.0 / 2.0), result.getPolyParams(),
                ABSOLUTE_ERROR);
    }

    @Test
    void testNormalizeHighestDegreeTerm() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-1.0, 0.0, 1.0}, 1.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{1.0, -2.0, 1.0}, 1.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{2.0, 2.0}, 1.0 / 2.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{6.0, 5.0, 1.0}, 1.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-21.0, -19.0, 8.0, 4.0}, 1.0 / 4.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, 1.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(
                ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, 1.0 / 6.0),
                p.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        p.normalizeHighestDegreeTerm();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(
                new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, 1.0 / 2.0), p.getPolyParams(), ABSOLUTE_ERROR);
    }

    @Test
    void testNormalizeHighestDegreeTermAndReturnNew() {
        // x^2 - 1
        var p = new Polynomial(-1.0, 0.0, 1.0);
        var result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-1.0, 0.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 - 2*x + 1
        p = new Polynomial(1.0, -2.0, 1.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{1.0, -2.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x + 2
        p = new Polynomial(2.0, 2.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{2.0, 2.0}, 1.0 / 2.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // x^2 + 5*x + 6
        p = new Polynomial(6.0, 5.0, 1.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(new double[]{6.0, 5.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 4*x^3 + 8*x^2 - 19*x - 21
        p = new Polynomial(-21.0, -19.0, 8.0, 4.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(
                new double[]{-21.0, -19.0, 8.0, 4.0}, 1.0 / 4.0), result.getPolyParams(), ABSOLUTE_ERROR);

        // x^4 + 5*x^3 + 10*x^2 - 9*x - 34
        p = new Polynomial(-34.0, -9.0, 10.0, 5.0, 1.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(
                ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-34.0, -9.0, 10.0, 5.0, 1.0}, 1.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 6*x^4 - 6*x^3 - 47*x^2 + 83*x - 35
        p = new Polynomial(-35.0, 83.0, -47.0, -6.0, 6.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(
                ArrayUtils.multiplyByScalarAndReturnNew(new double[]{-35.0, 83.0, -47.0, -6.0, 6.0}, 1.0 / 6.0),
                result.getPolyParams(), ABSOLUTE_ERROR);

        // 2*x^6 + 4*x^5 + x^4 + 11*x^3 + 2*x^2 + 4*x + 4
        p = new Polynomial(4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0);
        result = p.normalizeHighestDegreeTermAndReturnNew();

        // check correctness
        assertArrayEquals(ArrayUtils.multiplyByScalarAndReturnNew(
                new double[]{4.0, 4.0, 2.0, 11.0, 1.0, 4.0, 2.0}, 1.0 / 2.0), result.getPolyParams(),
                ABSOLUTE_ERROR);
    }

    @Test
    void testGetMaxima() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getMaxima());

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            assertNull(p.getMaxima());

            // reverse sign
            p.multiplyByScalar(-1.0);

            var maxima = p.getMaxima();
            assertEquals(1, maxima.length);
            assertEquals(maxima[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            var p3 = new Polynomial(-root3, 1.0);

            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            maxima = p.getMaxima();
            assertEquals(1, maxima.length);
            assertTrue(maxima[0] >= root2 && maxima[0] <= root3);

            // reverse sign
            p.multiplyByScalar(-1.0);

            maxima = p.getMaxima();
            assertEquals(2, maxima.length);

            assertTrue(maxima[0] >= root1 && maxima[0] <= root2);
            assertTrue(maxima[1] >= root3 && maxima[1] <= root4);
        }
    }

    @Test
    void testGetMaximaWithThreshold() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getMaxima(ABSOLUTE_ERROR));

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            assertNull(p.getMaxima(ABSOLUTE_ERROR));

            // reverse sign
            p.multiplyByScalar(-1.0);

            var maxima = p.getMaxima(ABSOLUTE_ERROR);
            assertEquals(1, maxima.length);
            assertEquals(maxima[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);

            var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            maxima = p.getMaxima(ABSOLUTE_ERROR);
            assertEquals(1, maxima.length);
            assertTrue(maxima[0] >= root2 && maxima[0] <= root3);

            // reverse sign
            p.multiplyByScalar(-1.0);

            maxima = p.getMaxima();
            assertEquals(2, maxima.length);

            assertTrue(maxima[0] >= root1 && maxima[0] <= root2);
            assertTrue(maxima[1] >= root3 && maxima[1] <= root4);

            // Force IllegalArgumentException
            final var p5 = p;
            assertThrows(IllegalArgumentException.class, () -> p5.getMaxima(-1.0));
        }
    }

    @Test
    void testGetMinima() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getMinima());

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            var minima = p.getMinima();
            assertEquals(1, minima.length);
            assertEquals(minima[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // reverse sign
            p.multiplyByScalar(-1.0);

            assertNull(p.getMinima());

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);

            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            minima = p.getMinima();
            assertEquals(2, minima.length);

            assertTrue(minima[0] >= root1 && minima[0] <= root2);
            assertTrue(minima[1] >= root3 && minima[1] <= root4);

            // reverse sign
            p.multiplyByScalar(-1.0);

            minima = p.getMinima();
            assertEquals(1, minima.length);
            assertTrue(minima[0] >= root2 && minima[0] <= root3);
        }
    }

    @Test
    void testGetMinimaWithThreshold() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getMinima(ABSOLUTE_ERROR));

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            var minima = p.getMinima(ABSOLUTE_ERROR);
            assertEquals(1, minima.length);
            assertEquals(minima[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // reverse sign
            p.multiplyByScalar(-1.0);

            assertNull(p.getMinima(ABSOLUTE_ERROR));

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);

            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            minima = p.getMinima(ABSOLUTE_ERROR);
            assertEquals(2, minima.length);

            assertTrue(minima[0] >= root1 && minima[0] <= root2);
            assertTrue(minima[1] >= root3 && minima[1] <= root4);

            // reverse sign
            p.multiplyByScalar(-1.0);

            minima = p.getMinima(ABSOLUTE_ERROR);
            assertEquals(1, minima.length);
            assertTrue(minima[0] >= root2 && minima[0] <= root3);
        }
    }

    @Test
    void testGetExtrema() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getExtrema());

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            var extrema = p.getExtrema();
            assertEquals(1, extrema.length);
            assertEquals(extrema[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);

            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            extrema = p.getExtrema();
            assertEquals(3, extrema.length);
            assertTrue(extrema[0] >= root1 && extrema[0] <= root2);
            assertTrue(extrema[1] >= root2 && extrema[1] <= root3);
            assertTrue(extrema[2] >= root3 && extrema[2] <= root4);
        }
    }

    @Test
    void testGetExtremaWithThreshold() throws NumericalException {
        final var randomizer = new UniformRandomizer();

        for (var t = 0; t < TIMES; t++) {
            // 1st degree
            final var root1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var p1 = new Polynomial(-root1, 1.0);
            var p = p1;

            assertNull(p.getExtrema(ABSOLUTE_ERROR));

            // 2nd degree
            final var root2 = randomizer.nextDouble(root1, root1 + MAX_RANDOM_VALUE);
            final var p2 = new Polynomial(-root2, 1.0);
            p = p1.multiplyAndReturnNew(p2);

            var extrema = p.getExtrema(ABSOLUTE_ERROR);
            assertEquals(1, extrema.length);
            assertEquals(extrema[0], 0.5 * (root1 + root2), ABSOLUTE_ERROR);

            // 4th degree
            final var root3 = randomizer.nextDouble(root2, root2 + MAX_RANDOM_VALUE);
            final var p3 = new Polynomial(-root3, 1.0);

            final var root4 = randomizer.nextDouble(root3, root3 + MAX_RANDOM_VALUE);
            final var p4 = new Polynomial(-root4, 1.0);

            p = p1.multiplyAndReturnNew(p2).multiplyAndReturnNew(p3).multiplyAndReturnNew(p4);

            extrema = p.getExtrema(ABSOLUTE_ERROR);
            assertEquals(3, extrema.length);
            assertTrue(extrema[0] >= root1 && extrema[0] <= root2);
            assertTrue(extrema[1] >= root2 && extrema[1] <= root3);
            assertTrue(extrema[2] >= root3 && extrema[2] <= root4);
        }
    }

    @Test
    void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final var randomizer = new UniformRandomizer();
        final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        final var polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var p1 = new Polynomial(polyParams);

        final var bytes = SerializationHelper.serialize(p1);
        final var p2 = SerializationHelper.<Polynomial>deserialize(bytes);

        assertArrayEquals(p1.getPolyParams(), p2.getPolyParams(), 0.0);
        assertEquals(p1.getDegree(), p2.getDegree());
    }
}
