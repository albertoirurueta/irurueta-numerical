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
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;

public class SimplexMultiOptimizerTest implements
        MultiDimensionFunctionEvaluatorListener,
        OnIterationCompletedListener {

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
    public void testConstructor() throws NotAvailableException,
            WrongSizeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final double delta = randomizer.nextDouble(MIN_DELTA, MAX_DELTA);

        final double[] deltas = new double[ndims];
        randomizer.fill(deltas, MIN_DELTA, MAX_DELTA);

        final double[] badDeltas = new double[ndims + 1];

        final Matrix simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims,
                MIN_EVAL_POINT, MAX_EVAL_POINT);

        SimplexMultiOptimizer optimizer;

        // Test 1st constructor
        optimizer = new SimplexMultiOptimizer();
        assertNotNull(optimizer);

        try {
            optimizer.getSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isSimplexAvailable());
        try {
            optimizer.getEvaluationsAtSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(optimizer.getTolerance(),
                SimplexMultiOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isReady());
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
        optimizer = new SimplexMultiOptimizer(this, startPoint, delta,
                tolerance);
        assertNotNull(optimizer);

        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());
        try {
            optimizer.getEvaluationsAtSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getListener(), this);
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
            optimizer = new SimplexMultiOptimizer(this, startPoint, delta,
                    -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);

        // Test 3rd constructor
        optimizer = new SimplexMultiOptimizer(this, startPoint, deltas,
                tolerance);
        assertNotNull(optimizer);

        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());
        try {
            optimizer.getEvaluationsAtSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getListener(), this);
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
            optimizer = new SimplexMultiOptimizer(this, startPoint, deltas,
                    -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            optimizer = new SimplexMultiOptimizer(this, startPoint, badDeltas,
                    tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);

        // Test 4th constructor
        optimizer = new SimplexMultiOptimizer(this, simplex, tolerance);
        assertNotNull(optimizer);

        assertEquals(optimizer.getSimplex(), simplex);
        assertTrue(optimizer.isSimplexAvailable());
        try {
            optimizer.getEvaluationsAtSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areFunctionEvaluationsAvailable());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getListener(), this);
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
            optimizer = new SimplexMultiOptimizer(this, simplex, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);
    }

    @Test
    public void testGetSetSimplexAndAvailability() throws LockedException,
            NotAvailableException, WrongSizeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final double delta = randomizer.nextDouble(MIN_DELTA, MAX_DELTA);

        final double[] deltas = new double[ndims];
        randomizer.fill(deltas, MIN_DELTA, MAX_DELTA);

        final double[] badDeltas = new double[ndims + 1];

        final Matrix simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims,
                MIN_EVAL_POINT, MAX_EVAL_POINT);

        SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer();

        try {
            optimizer.getSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex using start point and delta
        optimizer.setSimplex(startPoint, delta);

        // check correctness
        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());


        optimizer = new SimplexMultiOptimizer();
        try {
            optimizer.getSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex using start point and deltas
        optimizer.setSimplex(startPoint, deltas);

        // check correctness
        assertNotNull(optimizer.getSimplex());
        assertTrue(optimizer.isSimplexAvailable());

        // Force IllegalArgumentException
        try {
            optimizer.setSimplex(startPoint, badDeltas);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }


        optimizer = new SimplexMultiOptimizer();
        try {
            optimizer.getSimplex();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isSimplexAvailable());

        // set simplex
        optimizer.setSimplex(simplex);

        // check correctness
        assertEquals(optimizer.getSimplex(), simplex);
        assertTrue(optimizer.isSimplexAvailable());
    }

    @Test
    public void testGetSetTolerance() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer();

        assertEquals(optimizer.getTolerance(),
                SimplexMultiOptimizer.DEFAULT_TOLERANCE, 0.0);

        // set tolerance
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
    public void testIsReady() throws LockedException, WrongSizeException {

        final SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer();

        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(this);
        assertFalse(optimizer.isReady());

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final Matrix simplex = Matrix.createWithUniformRandomValues(ndims + 1, ndims,
                MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set simplex
        optimizer.setSimplex(simplex);

        // now optimizer is ready
        assertTrue(optimizer.isReady());
    }

    @Test
    public void testGetSetOnIterationCompletedListener() throws LockedException {
        final SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    public void testGetSetListenerAndAvailability() throws LockedException {

        final SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer();

        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isListenerAvailable());

        // set listener
        optimizer.setListener(this);
        assertTrue(optimizer.isListenerAvailable());
    }

    @Test
    public void testMinimize() throws Throwable {

        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            final double[] startPoint = Arrays.copyOf(minimum, ndims);
            // add some noise to start point
            for (int i = 0; i < ndims; i++) {
                startPoint[i] += randomizer.nextDouble(MIN_DELTA, MAX_DELTA);
            }

            final double delta = randomizer.nextDouble(2.0 * MIN_DELTA,
                    2.0 * MAX_DELTA);


            final SimplexMultiOptimizer optimizer = new SimplexMultiOptimizer(this,
                    startPoint, delta, SimplexMultiOptimizer.DEFAULT_TOLERANCE);
            optimizer.setOnIterationCompletedListener(this);

            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            assertTrue(optimizer.isListenerAvailable());
            assertTrue(optimizer.isSimplexAvailable());
            assertFalse(optimizer.areFunctionEvaluationsAvailable());
            try {
                optimizer.getEvaluationAtResult();
                fail("NotAvailableException expected but not thrown");
            } catch (final NotAvailableException ignore) {
            }
            try {
                optimizer.getResult();
                fail("NotAvailableException expected but not thrown");
            } catch (final NotAvailableException ignore) {
            }

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

            final double[] evaluationsAtSimplex = optimizer.getEvaluationsAtSimplex();
            final double evaluationAtResult = optimizer.getEvaluationAtResult();

            final double[] result = optimizer.getResult();

            // check correctness

            // check that obtained function value is close to that on the true
            // minimum
            final double evaluationAtMinimum = evaluate(minimum);

            assertEquals(evaluationAtResult, evaluationAtMinimum,
                    ABSOLUTE_ERROR);

            assertNotNull(evaluationsAtSimplex);
            assertNotNull(result);
        }
    }

    @Override
    public double evaluate(final double[] point) throws EvaluationException {
        final int dims = Math.min(Math.min(point.length, minimum.length),
                width.length);

        double value = 1.0;
        for (int i = 0; i < dims; i++) {
            value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
        }

        value += offset;
        return value;
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
