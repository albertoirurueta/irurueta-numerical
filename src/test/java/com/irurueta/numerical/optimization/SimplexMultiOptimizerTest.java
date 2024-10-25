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

import static org.junit.jupiter.api.Assertions.*;

class SimplexMultiOptimizerTest implements MultiDimensionFunctionEvaluatorListener, OnIterationCompletedListener {

    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 4;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 1e-5;

    private static final double MIN_OFFSET = -10.0;
    private static final double MAX_OFFSET = 10.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;

    private static final double MIN_DELTA = -1.0;
    private static final double MAX_DELTA = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-4;

    private static final int TIMES = 100;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    private int iterations = 0;

    @Test
    void testConstructor() throws NotAvailableException, WrongSizeException {

        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var delta = randomizer.nextDouble(MIN_DELTA, MAX_DELTA);

        final var deltas = new double[ndims];
        randomizer.fill(deltas, MIN_DELTA, MAX_DELTA);

        final var badDeltas = new double[ndims + 1];

        final var simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims, MIN_EVAL_POINT,
                MAX_EVAL_POINT);

        SimplexMultiOptimizer optimizer;

        // Test 1st constructor
        optimizer = new SimplexMultiOptimizer();
        assertNotNull(optimizer);

        assertThrows(NotAvailableException.class, optimizer::getSimplex);
        assertFalse(optimizer.isSimplexAvailable());
        assertThrows(NotAvailableException.class, optimizer::getEvaluationsAtSimplex);
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(SimplexMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isReady());
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Test 2nd constructor
        optimizer = new SimplexMultiOptimizer(this, startPoint, delta, tolerance);
        assertNotNull(optimizer);

        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());
        assertThrows(NotAvailableException.class, optimizer::getEvaluationsAtSimplex);
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(this, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SimplexMultiOptimizer(this, startPoint, delta,
                -tolerance));

        // Test 3rd constructor
        optimizer = new SimplexMultiOptimizer(this, startPoint, deltas, tolerance);
        assertNotNull(optimizer);

        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());
        assertThrows(NotAvailableException.class, optimizer::getEvaluationsAtSimplex);
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(this, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SimplexMultiOptimizer(this, startPoint, deltas,
                -tolerance));
        assertThrows(IllegalArgumentException.class, () ->new SimplexMultiOptimizer(this, startPoint, badDeltas,
                tolerance) );

        // Test 4th constructor
        optimizer = new SimplexMultiOptimizer(this, simplex, tolerance);
        assertNotNull(optimizer);

        assertEquals(optimizer.getSimplex(), simplex);
        assertTrue(optimizer.isSimplexAvailable());
        assertThrows(NotAvailableException.class, optimizer::getEvaluationsAtSimplex);
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(this, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SimplexMultiOptimizer(this, simplex,
                -tolerance));
    }

    @Test
    void testGetSetSimplexAndAvailability() throws LockedException, NotAvailableException, WrongSizeException {

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var delta = randomizer.nextDouble(MIN_DELTA, MAX_DELTA);

        final var deltas = new double[ndims];
        randomizer.fill(deltas, MIN_DELTA, MAX_DELTA);

        final var badDeltas = new double[ndims + 1];

        final var simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims, MIN_EVAL_POINT,
                MAX_EVAL_POINT);

        var optimizer = new SimplexMultiOptimizer();
        assertThrows(NotAvailableException.class, optimizer::getSimplex);
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex using start point and delta
        optimizer.setSimplex(startPoint, delta);

        // check correctness
        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());

        optimizer = new SimplexMultiOptimizer();
        assertThrows(NotAvailableException.class, optimizer::getSimplex);
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex using start point and deltas
        optimizer.setSimplex(startPoint, deltas);

        // check correctness
        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());

        // Force IllegalArgumentException
        final var finalOptimizer = optimizer;
        assertThrows(IllegalArgumentException.class, () -> finalOptimizer.setSimplex(startPoint, badDeltas));

        optimizer = new SimplexMultiOptimizer();
        assertThrows(NotAvailableException.class, optimizer::getSimplex);
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex
        optimizer.setSimplex(simplex);

        // check correctness
        assertEquals(simplex, optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());
    }

    @Test
    void testGetSetTolerance() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var optimizer = new SimplexMultiOptimizer();

        assertEquals(SimplexMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);

        // set tolerance
        optimizer.setTolerance(tolerance);

        // check correctness
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> optimizer.setTolerance(-tolerance));
    }

    @Test
    void testIsReady() throws LockedException, WrongSizeException {

        final var optimizer = new SimplexMultiOptimizer();

        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(this);
        assertFalse(optimizer.isReady());

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims, MIN_EVAL_POINT,
                MAX_EVAL_POINT);

        // set simplex
        optimizer.setSimplex(simplex);

        // now optimizer is ready
        assertTrue(optimizer.isReady());
    }

    @Test
    void testGetSetOnIterationCompletedListener() throws LockedException {
        final var optimizer = new SimplexMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testGetSetListenerAndAvailability() throws LockedException {

        final var optimizer = new SimplexMultiOptimizer();

        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());

        // set listener
        optimizer.setListener(this);
        assertTrue(optimizer.isListenerAvailable());
    }

    @Test
    void testMinimize() throws Throwable {

        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            final var startPoint = Arrays.copyOf(minimum, ndims);
            // add some noise to start point
            for (var i = 0; i < ndims; i++) {
                startPoint[i] += randomizer.nextDouble(MIN_DELTA, MAX_DELTA);
            }

            final var delta = randomizer.nextDouble(2.0 * MIN_DELTA, 2.0 * MAX_DELTA);

            final var optimizer = new SimplexMultiOptimizer(this, startPoint, delta,
                    SimplexMultiOptimizer.DEFAULT_TOLERANCE);
            optimizer.setOnIterationCompletedListener(this);

            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            assertTrue(optimizer.isListenerAvailable());
            assertTrue(optimizer.isSimplexAvailable());
            assertFalse(optimizer.areFunctionEvaluationsAvailable());
            assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
            assertThrows(NotAvailableException.class, optimizer::getResult);

            reset();
            assertEquals(0, iterations);

            // optimize
            optimizer.minimize();

            // check correctness
            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertTrue(optimizer.isListenerAvailable());
            assertTrue(optimizer.isSimplexAvailable());
            assertTrue(optimizer.areFunctionEvaluationsAvailable());
            assertTrue(iterations >= 0);

            final var evaluationsAtSimplex = optimizer.getEvaluationsAtSimplex();
            final var evaluationAtResult = optimizer.getEvaluationAtResult();

            final var result = optimizer.getResult();

            // check correctness

            // check that obtained function value is close to that on the true
            // minimum
            final var evaluationAtMinimum = evaluate(minimum);

            assertEquals(evaluationAtResult, evaluationAtMinimum, ABSOLUTE_ERROR);

            assertNotNull(evaluationsAtSimplex);
            assertNotNull(result);
        }
    }

    @Override
    public double evaluate(final double[] point) {
        final var dims = Math.min(Math.min(point.length, minimum.length), width.length);

        var value = 1.0;
        for (var i = 0; i < dims; i++) {
            value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
        }

        value += offset;
        return value;
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
