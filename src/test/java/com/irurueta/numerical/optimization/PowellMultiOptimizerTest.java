/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.optimization;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;

public class PowellMultiOptimizerTest implements OnIterationCompletedListener {

    public static final int MIN_DIMS = 2;
    public static final int MAX_DIMS = 4;

    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;

    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;

    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;

    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 3e-6;

    public static final double ABSOLUTE_ERROR = 1e-6;

    public static final int TIMES = 100;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    private int iterations = 0;

    private final MultiDimensionFunctionEvaluatorListener listener;

    public PowellMultiOptimizerTest() {
        listener = new MultiDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double[] point) {
                final int dims = Math.min(Math.min(point.length, minimum.length),
                        width.length);

                double value = 1.0;
                for (int i = 0; i < dims; i++) {
                    value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
                }

                value += offset;

                return value;
            }
        };
    }

    @Test
    public void testConstructor() throws WrongSizeException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        minimum = new double[ndims];
        final double[] point = new double[ndims];
        final Matrix setsOfDirections = new Matrix(ndims, ndims);

        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        PowellMultiOptimizer optimizer;

        // Test 1st constructor
        optimizer = new PowellMultiOptimizer();
        assertNotNull(optimizer);

        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isDirectionAvailable());
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(),
                PowellMultiOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());


        // Test 2nd constructor
        optimizer = new PowellMultiOptimizer(listener, tolerance);
        assertNotNull(optimizer);

        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new PowellMultiOptimizer(listener, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);


        // Test 3rd constructor
        optimizer = new PowellMultiOptimizer(listener, point, setsOfDirections,
                tolerance);
        assertNotNull(optimizer);

        assertEquals(optimizer.getDirections(), setsOfDirections);
        assertTrue(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isDirectionAvailable());
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), point, 0.0);
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());
    }

    @Test
    public void testGetSetDirectionsAndAvailability()
            throws WrongSizeException, LockedException, NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final Matrix directions = new Matrix(nDims, nDims);
        final double[] point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areDirectionsAvailable());

        // set directions
        optimizer.setPointAndDirections(point, directions);
        // check correctness
        assertEquals(optimizer.getDirections(), directions);
        assertTrue(optimizer.areDirectionsAvailable());
    }

    @Test
    public void testGetSetDirection() throws LockedException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final double[] direction = new double[nDims];


        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isDirectionAvailable());

        // set direction
        optimizer.setStartPointAndDirection(point, direction);

        // Check correctness
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);
        assertTrue(optimizer.isDirectionAvailable());
    }

    @Test
    public void testIsReady() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] point = new double[dims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final double[] direction = new double[dims];

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(listener);
        assertFalse(optimizer.isReady());

        // set start point
        optimizer.setStartPointAndDirection(point, direction);
        assertTrue(optimizer.isReady());
    }

    @Test
    public void testGetSetOnIterationCompletedListener() throws LockedException {
        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    public void testGetSetTolerance() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        // get tolerance
        optimizer.setTolerance(tolerance);
        // check correctness
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);

        // Force IllegalArgumentException
        try {
            optimizer.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetStartPointAndAvailability() throws LockedException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[dims];
        final double[] direction = new double[dims];

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        // get start point
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isStartPointAvailable());

        // set start point
        optimizer.setStartPointAndDirection(startPoint, direction);
        // check correctness
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertTrue(optimizer.isStartPointAvailable());
    }

    @Test
    public void testGetSetListenerAndAvailability() throws LockedException {

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        // get listener
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isListenerAvailable());

        // set listener
        optimizer.setListener(listener);
        assertTrue(optimizer.isListenerAvailable());
    }

    @Test
    public void testIsLocked() {

        final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    public void testMinimize() throws Throwable {

        for (int i = 0; i < TIMES; i++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            final double[] point = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

            final double[] startPoint = Arrays.copyOf(point, ndims);

            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);


            final PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
            optimizer.setListener(listener);
            optimizer.setStartPoint(startPoint);
            optimizer.setOnIterationCompletedListener(this);

            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());

            reset();
            assertEquals(0, iterations);

            // minimize
            optimizer.minimize();

            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());

            // get result
            final double[] result = optimizer.getResult();

            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());
            assertTrue(iterations > 0);

            // check correctness of estimated result by comparing the difference of
            // function value at true minimum and the estimated minimum. If both are
            // very similar, then the algorithm will have converged

            // value at estimated result
            final double value1 = listener.evaluate(result);
            // returned value by optimizer (must be equal to value1)
            final double value2 = optimizer.getEvaluationAtResult();
            // value at true minimum
            final double value3 = listener.evaluate(minimum);
            // value at start point
            final double value4 = listener.evaluate(startPoint);

            assertEquals(value1, value2, 0.0);
            assertEquals(value1, value3, ABSOLUTE_ERROR);

            // check that indeed function has been minimized
            assertTrue(value1 <= value4);
        }
    }

    @Override
    public void onIterationCompleted(
            final Optimizer optimizer, final int iteration, final Integer maxIterations) {
        assertNotNull(optimizer);
        assertTrue(optimizer.isLocked());
        assertTrue(iteration >= 0);
        assertNotNull(maxIterations);
        assertTrue(iteration <= maxIterations);
        iterations++;
    }

    private void reset() {
        iterations = 0;
    }
}
