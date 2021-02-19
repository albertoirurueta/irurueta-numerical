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

import com.irurueta.numerical.GradientFunctionEvaluatorListener;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class DerivativeConjugateGradientMultiOptimizerTest implements
        OnIterationCompletedListener {

    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 4;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_OFFSET = 0.0;
    private static final double MAX_OFFSET = 1.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;

    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 3e-6;

    private static final double ABSOLUTE_ERROR = 1e-3;

    private static final double MIN_DIRECTION = -1.0;
    private static final double MAX_DIRECTION = 1.0;

    private static final int TIMES = 100;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    private int iterations = 0;

    private final MultiDimensionFunctionEvaluatorListener listener;
    private final GradientFunctionEvaluatorListener gradientListener;

    public DerivativeConjugateGradientMultiOptimizerTest() {
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

        gradientListener = new GradientFunctionEvaluatorListener() {

            @Override
            public void evaluateGradient(final double[] params, final double[] result) {
                final int dims = Math.min(Math.min(params.length, minimum.length),
                        width.length);

                double value;
                for (int j = 0; j < dims; j++) {
                    value = 1.0;
                    for (int i = 0; i < dims; i++) {
                        if (i != j) {
                            value *= Math.pow(params[i] - minimum[i], 2.0) /
                                    width[i];
                        } else {
                            value *= 2.0 * (params[i] - minimum[i]) / width[i];
                        }
                    }

                    result[j] = value;
                }
            }
        };
    }

    @Test
    public void testConstructor() throws NotAvailableException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final boolean usePolakRibiere = randomizer.nextBoolean();

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        final double[] direction = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);

        DerivativeConjugateGradientMultiOptimizer optimizer;

        // test 1st constructor
        optimizer = new DerivativeConjugateGradientMultiOptimizer();
        assertNotNull(optimizer);
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(),
                DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE,
                0.0);
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.isPolakRibiereEnabled(),
                DerivativeConjugateGradientMultiOptimizer.
                        DEFAULT_USE_POLAK_RIBIERE);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getDirection();
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
        assertNull(optimizer.getOnIterationCompletedListener());


        //test 2nd constructor
        optimizer = new DerivativeConjugateGradientMultiOptimizer(listener,
                gradientListener, startPoint, direction, tolerance,
                usePolakRibiere);
        assertNotNull(optimizer);
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.isPolakRibiereEnabled(), usePolakRibiere);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);
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
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new DerivativeConjugateGradientMultiOptimizer(listener,
                    gradientListener, startPoint, direction, -tolerance,
                    usePolakRibiere);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);
    }

    @Test
    public void testIsReady() throws LockedException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(listener);

        assertFalse(optimizer.isReady());

        // set gradient listener
        optimizer.setGradientListener(gradientListener);

        assertFalse(optimizer.isReady());

        // set start point
        optimizer.setStartPoint(startPoint);

        // once both listeners and start point are available, optimizer becomes
        // ready
        assertTrue(optimizer.isReady());
    }

    @Test
    public void testGetSetOnIterationCompletedListener() throws LockedException {
        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

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

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        // get tolerance
        assertEquals(optimizer.getTolerance(),
                DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE,
                0.0);

        // set new tolerance
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
    public void testGetSetGradientListenerAndAvailability()
            throws LockedException, NotAvailableException {

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        // Get gradient listener
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isGradientListenerAvailable());

        // set gradient listener
        optimizer.setGradientListener(gradientListener);

        // check correctness
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
    }

    @Test
    public void testGetSetUsePolakRibiere() throws LockedException {
        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        // get initial status
        assertEquals(optimizer.isPolakRibiereEnabled(),
                DerivativeConjugateGradientMultiOptimizer.
                        DEFAULT_USE_POLAK_RIBIERE);

        // disable
        optimizer.setUsePolakRibiere(false);

        // check correctness
        assertFalse(optimizer.isPolakRibiereEnabled());

        // enable
        optimizer.setUsePolakRibiere(true);

        // check correctness
        assertTrue(optimizer.isPolakRibiereEnabled());
    }

    @Test
    public void testGetSetStartPointAndAvailability() throws LockedException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        // get start point
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }

        // set start point
        optimizer.setStartPoint(startPoint);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
    }

    @Test
    public void testGetSetStartPointAndDirectionAndAvailability()
            throws LockedException, NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final double[] startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final double[] direction = new double[nDims];
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        assertFalse(optimizer.isStartPointAvailable());
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }

        // set start point and direction
        optimizer.setStartPointAndDirection(startPoint, direction);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);

        // Force IllegalArgumentException
        final double[] wrongStartPoint = new double[nDims + 1];
        try {
            optimizer.setStartPointAndDirection(wrongStartPoint, direction);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetListenerAndAvailability() throws LockedException,
            NotAvailableException {

        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        // check initial status
        assertFalse(optimizer.isListenerAvailable());
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }

        // set listener
        optimizer.setListener(listener);

        // check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
    }

    @Test
    public void testIsLocked() {
        final DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    public void testMinimize() throws Throwable {

        // we repeat the process because depending on the start point this
        // algorithm is not very accurate
        for (int i = 0; i < TIMES; i++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            final double[] startPoint = new double[ndims];
            minimum = new double[ndims];
            width = new double[ndims];
            randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            final DerivativeConjugateGradientMultiOptimizer optimizer =
                    new DerivativeConjugateGradientMultiOptimizer();
            optimizer.setListener(listener);
            optimizer.setGradientListener(gradientListener);
            optimizer.setStartPoint(startPoint);
            optimizer.setOnIterationCompletedListener(this);

            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            try {
                optimizer.getResult();
                fail("NotAvailableException expected but not thrown");
            } catch (final NotAvailableException ignore) {
            }

            reset();
            assertEquals(0, iterations);

            // minimize
            optimizer.minimize();

            // check correctness of result
            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertTrue(iterations > 0);

            final double[] result = optimizer.getResult();

            // check that function at estimated result and at tue minimum have
            // almost the same value
            final double valueMin = listener.evaluate(minimum);
            final double valueResult = listener.evaluate(result);

            assertEquals(valueMin, valueResult, ABSOLUTE_ERROR);
        }
    }

    @Override
    public void onIterationCompleted(Optimizer optimizer, int iteration, Integer maxIterations) {
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
