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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class DerivativeConjugateGradientMultiOptimizerTest implements OnIterationCompletedListener {

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

    private final MultiDimensionFunctionEvaluatorListener listener = point -> {
        final var dims = Math.min(Math.min(point.length, minimum.length),
                width.length);

        var value = 1.0;
        for (var i = 0; i < dims; i++) {
            value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
        }

        value += offset;
        return value;
    };

    private final GradientFunctionEvaluatorListener gradientListener = (params, result) -> {
        final var dims = Math.min(Math.min(params.length, minimum.length), width.length);

        for (var j = 0; j < dims; j++) {
            var value = 1.0;
            for (var i = 0; i < dims; i++) {
                if (i != j) {
                    value *= Math.pow(params[i] - minimum[i], 2.0) / width[i];
                } else {
                    value *= 2.0 * (params[i] - minimum[i]) / width[i];
                }
            }

            result[j] = value;
        }
    };

    @Test
    void testConstructor() throws NotAvailableException {
        final var randomizer = new UniformRandomizer();

        final var usePolakRibiere = randomizer.nextBoolean();

        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        final var direction = new double[ndims];
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
        assertEquals(DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);
        assertThrows(NotAvailableException.class, optimizer::getGradientListener);
        assertFalse(optimizer.isGradientListenerAvailable());
        assertEquals(DerivativeConjugateGradientMultiOptimizer.DEFAULT_USE_POLAK_RIBIERE,
                optimizer.isPolakRibiereEnabled());
        assertFalse(optimizer.isStartPointAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertFalse(optimizer.isDirectionAvailable());
        assertThrows(NotAvailableException.class, optimizer::getDirection);
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertNull(optimizer.getOnIterationCompletedListener());

        //test 2nd constructor
        optimizer = new DerivativeConjugateGradientMultiOptimizer(listener, gradientListener, startPoint, direction,
                tolerance, usePolakRibiere);
        assertNotNull(optimizer);
        assertTrue(optimizer.isReady());
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertEquals(gradientListener, optimizer.getGradientListener());
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(usePolakRibiere, optimizer.isPolakRibiereEnabled());
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(startPoint, optimizer.getStartPoint(), 0.0);
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(direction, optimizer.getDirection(), 0.0);
        assertEquals(listener, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new DerivativeConjugateGradientMultiOptimizer(listener,
                gradientListener, startPoint, direction, -tolerance, usePolakRibiere));
    }

    @Test
    void testIsReady() throws LockedException {
        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

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
    void testGetSetOnIterationCompletedListener() throws LockedException {
        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testGetSetTolerance() throws LockedException {
        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        // get tolerance
        assertEquals(DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);

        // set new tolerance
        optimizer.setTolerance(tolerance);

        // check correctness
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> optimizer.setTolerance(-tolerance));
    }

    @Test
    void testGetSetGradientListenerAndAvailability() throws LockedException, NotAvailableException {

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        // Get gradient listener
        assertThrows(NotAvailableException.class, optimizer::getGradientListener);
        assertFalse(optimizer.isGradientListenerAvailable());

        // set gradient listener
        optimizer.setGradientListener(gradientListener);

        // check correctness
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
    }

    @Test
    void testGetSetUsePolakRibiere() throws LockedException {
        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        // get initial status
        assertEquals(DerivativeConjugateGradientMultiOptimizer.DEFAULT_USE_POLAK_RIBIERE,
                optimizer.isPolakRibiereEnabled());

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
    void testGetSetStartPointAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        // get start point
        assertFalse(optimizer.isStartPointAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);

        // set start point
        optimizer.setStartPoint(startPoint);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(startPoint, optimizer.getStartPoint(), 0.0);
    }

    @Test
    void testGetSetStartPointAndDirectionAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var direction = new double[nDims];
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        assertFalse(optimizer.isStartPointAvailable());
        assertFalse(optimizer.isDirectionAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertThrows(NotAvailableException.class, optimizer::getDirection);

        // set start point and direction
        optimizer.setStartPointAndDirection(startPoint, direction);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(startPoint, optimizer.getStartPoint(), 0.0);
        assertArrayEquals(direction, optimizer.getDirection(), 0.0);

        // Force IllegalArgumentException
        final var wrongStartPoint = new double[nDims + 1];
        assertThrows(IllegalArgumentException.class,
                () -> optimizer.setStartPointAndDirection(wrongStartPoint, direction));
    }

    @Test
    void testGetSetListenerAndAvailability() throws LockedException, NotAvailableException {

        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        // check initial status
        assertFalse(optimizer.isListenerAvailable());
        assertThrows(NotAvailableException.class, optimizer::getListener);

        // set listener
        optimizer.setListener(listener);

        // check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
    }

    @Test
    void testIsLocked() {
        final var optimizer = new DerivativeConjugateGradientMultiOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    void testMinimize() throws Throwable {

        // we repeat the process because depending on the start point this
        // algorithm is not very accurate
        for (var i = 0; i < TIMES; i++) {
            final var randomizer = new UniformRandomizer();
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            final var startPoint = new double[ndims];
            minimum = new double[ndims];
            width = new double[ndims];
            randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            final var optimizer = new DerivativeConjugateGradientMultiOptimizer();
            optimizer.setListener(listener);
            optimizer.setGradientListener(gradientListener);
            optimizer.setStartPoint(startPoint);
            optimizer.setOnIterationCompletedListener(this);

            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            assertThrows(NotAvailableException.class, optimizer::getResult);

            reset();
            assertEquals(0, iterations);

            // minimize
            optimizer.minimize();

            // check correctness of result
            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertTrue(iterations > 0);

            final var result = optimizer.getResult();

            // check that function at estimated result and at tue minimum have
            // almost the same value
            final var valueMin = listener.evaluate(minimum);
            final var valueResult = listener.evaluate(result);

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
