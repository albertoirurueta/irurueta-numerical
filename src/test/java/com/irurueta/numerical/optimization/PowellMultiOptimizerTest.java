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
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

class PowellMultiOptimizerTest implements OnIterationCompletedListener {

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
        listener = point -> {
            final var dims = Math.min(Math.min(point.length, minimum.length), width.length);

            var value = 1.0;
            for (var i = 0; i < dims; i++) {
                value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
            }

            value += offset;

            return value;
        };
    }

    @Test
    void testConstructor() throws WrongSizeException, NotAvailableException {

        final var randomizer = new UniformRandomizer(new Random());

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        minimum = new double[ndims];
        final var point = new double[ndims];
        final var setsOfDirections = new Matrix(ndims, ndims);

        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        PowellMultiOptimizer optimizer;

        // Test 1st constructor
        optimizer = new PowellMultiOptimizer();
        assertNotNull(optimizer);

        assertThrows(NotAvailableException.class, optimizer::getDirections);
        assertFalse(optimizer.areDirectionsAvailable());
        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertFalse(optimizer.isDirectionAvailable());
        assertFalse(optimizer.isReady());
        assertEquals(PowellMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Test 2nd constructor
        optimizer = new PowellMultiOptimizer(listener, tolerance);
        assertNotNull(optimizer);

        assertThrows(NotAvailableException.class, optimizer::getDirections);
        assertFalse(optimizer.areDirectionsAvailable());
        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertFalse(optimizer.isDirectionAvailable());
        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertFalse(optimizer.isReady());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertEquals(listener, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PowellMultiOptimizer(listener, -tolerance));

        // Test 3rd constructor
        optimizer = new PowellMultiOptimizer(listener, point, setsOfDirections, tolerance);
        assertNotNull(optimizer);

        assertEquals(setsOfDirections, optimizer.getDirections());
        assertTrue(optimizer.areDirectionsAvailable());
        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertFalse(optimizer.isDirectionAvailable());
        assertTrue(optimizer.isReady());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(point, optimizer.getStartPoint(), 0.0);
        assertEquals(listener, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testGetSetDirectionsAndAvailability() throws WrongSizeException, LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer(new Random());
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var directions = new Matrix(nDims, nDims);
        final var point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final var optimizer = new PowellMultiOptimizer();

        assertThrows(NotAvailableException.class, optimizer::getDirections);
        assertFalse(optimizer.areDirectionsAvailable());

        // set directions
        optimizer.setPointAndDirections(point, directions);
        // check correctness
        assertEquals(directions, optimizer.getDirections());
        assertTrue(optimizer.areDirectionsAvailable());
    }

    @Test
    void testGetSetDirection() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer(new Random());
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final var direction = new double[nDims];

        final var optimizer = new PowellMultiOptimizer();

        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertFalse(optimizer.isDirectionAvailable());

        // set direction
        optimizer.setStartPointAndDirection(point, direction);

        // Check correctness
        assertArrayEquals(direction, optimizer.getDirection(), 0.0);
        assertTrue(optimizer.isDirectionAvailable());
    }

    @Test
    void testIsReady() throws LockedException {

        final var randomizer = new UniformRandomizer(new Random());
        final var dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var point = new double[dims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var direction = new double[dims];

        final var optimizer = new PowellMultiOptimizer();
        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(listener);
        assertFalse(optimizer.isReady());

        // set start point
        optimizer.setStartPointAndDirection(point, direction);
        assertTrue(optimizer.isReady());
    }

    @Test
    void testGetSetOnIterationCompletedListener() throws LockedException {
        final var optimizer = new PowellMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testGetSetTolerance() throws LockedException {

        final var randomizer = new UniformRandomizer(new Random());
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var optimizer = new PowellMultiOptimizer();

        // get tolerance
        optimizer.setTolerance(tolerance);
        // check correctness
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> optimizer.setTolerance(-tolerance));
    }

    @Test
    void testGetSetStartPointAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer(new Random());
        final var dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[dims];
        final var direction = new double[dims];

        final var optimizer = new PowellMultiOptimizer();

        // get start point
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertFalse(optimizer.isStartPointAvailable());

        // set start point
        optimizer.setStartPointAndDirection(startPoint, direction);
        // check correctness
        assertArrayEquals(startPoint, optimizer.getStartPoint(), 0.0);
        assertTrue(optimizer.isStartPointAvailable());
    }

    @Test
    void testGetSetListenerAndAvailability() throws LockedException {

        final var optimizer = new PowellMultiOptimizer();

        // get listener
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());

        // set listener
        optimizer.setListener(listener);
        assertTrue(optimizer.isListenerAvailable());
    }

    @Test
    void testIsLocked() {

        final var optimizer = new PowellMultiOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    void testMinimize() throws Throwable {

        for (var i = 0; i < TIMES; i++) {
            final var randomizer = new UniformRandomizer(new Random());

            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            final var point = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);

            final var startPoint = Arrays.copyOf(point, ndims);

            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            final var optimizer = new PowellMultiOptimizer();
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
            final var result = optimizer.getResult();

            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());
            assertTrue(iterations > 0);

            // check correctness of estimated result by comparing the difference of
            // function value at true minimum and the estimated minimum. If both are
            // very similar, then the algorithm will have converged

            // value at estimated result
            final var value1 = listener.evaluate(result);
            // returned value by optimizer (must be equal to value1)
            final var value2 = optimizer.getEvaluationAtResult();
            // value at true minimum
            final var value3 = listener.evaluate(minimum);
            // value at start point
            final var value4 = listener.evaluate(startPoint);

            assertEquals(value1, value2, 0.0);
            assertEquals(value1, value3, ABSOLUTE_ERROR);

            // check that indeed function has been minimized
            assertTrue(value1 <= value4);
        }
    }

    @Override
    public void onIterationCompleted(final Optimizer optimizer, final int iteration, final Integer maxIterations) {
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
