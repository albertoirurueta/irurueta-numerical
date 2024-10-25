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

class QuasiNewtonMultiOptimizerTest implements OnIterationCompletedListener {

    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 4;

    private static final double MIN_EVAL_POINT = -1e3;
    private static final double MAX_EVAL_POINT = 1e3;

    private static final double MIN_OFFSET = -1e3;
    private static final double MAX_OFFSET = 1e3;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;

    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 3e-6;

    private static final double ABSOLUTE_ERROR = 1e-6;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    private int iterations = 0;

    private final MultiDimensionFunctionEvaluatorListener listener = point -> {
        final var dims = Math.min(Math.min(point.length, minimum.length), width.length);

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
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var startPoint = new double[ndims];
        minimum = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        width = new double[ndims];
        randomizer.fill(width, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        QuasiNewtonMultiOptimizer optimizer;

        // Test 1st constructor
        optimizer = new QuasiNewtonMultiOptimizer();
        assertNotNull(optimizer);

        assertEquals(QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertEquals(0, optimizer.getIterations());

        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertThrows(NotAvailableException.class, optimizer::getGradientListener);
        assertFalse(optimizer.isGradientListenerAvailable());
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Test 2nd constructor
        optimizer = new QuasiNewtonMultiOptimizer(listener, gradientListener, tolerance);
        assertNotNull(optimizer);
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertEquals(0, optimizer.getIterations());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);
        assertEquals(gradientListener, optimizer.getGradientListener());
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(listener, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new QuasiNewtonMultiOptimizer(listener, gradientListener,
                -tolerance));

        // Test 3rd constructor
        optimizer = new QuasiNewtonMultiOptimizer(listener, gradientListener, startPoint, tolerance);
        assertNotNull(optimizer);
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(startPoint, optimizer.getStartPoint(), ABSOLUTE_ERROR);
        assertEquals(gradientListener, optimizer.getGradientListener());
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(listener, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        assertEquals(0, optimizer.getIterations());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new QuasiNewtonMultiOptimizer(listener, gradientListener,
                startPoint, -tolerance));
    }

    @Test
    void testIsReady() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[nDims];

        final var optimizer = new QuasiNewtonMultiOptimizer();

        // optimizer will be ready when listener, gradient listener and start
        // point are provided
        assertFalse(optimizer.isReady());

        // set start point
        optimizer.setStartPoint(startPoint);

        // optimizer is not ready yet
        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(listener);

        // optimizer is not ready yet
        assertFalse(optimizer.isReady());

        // set gradient listener
        optimizer.setGradientListener(gradientListener);

        // optimizer is now ready
        assertTrue(optimizer.isReady());
    }

    @Test
    void testGetSetTolerance() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var optimizer = new QuasiNewtonMultiOptimizer();

        // get default tolerance
        assertEquals(QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);

        // set new tolerance
        optimizer.setTolerance(tolerance);

        // check correctness
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
    }

    @Test
    void testGetSetStartPointAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[nDims];

        final var optimizer = new QuasiNewtonMultiOptimizer();

        // initially start point is not available
        assertFalse(optimizer.isStartPointAvailable());
        assertThrows(NotAvailableException.class, optimizer::getStartPoint);

        // set start point
        optimizer.setStartPoint(startPoint);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, ABSOLUTE_ERROR);

        // set start point
        optimizer.setStartPoint(startPoint);

        // check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, ABSOLUTE_ERROR);
    }

    @Test
    void testGetSetGradientListenerAndAvailability() throws LockedException, NotAvailableException {

        final var optimizer = new QuasiNewtonMultiOptimizer();

        // initially gradient listener is not available
        assertFalse(optimizer.isGradientListenerAvailable());
        assertThrows(NotAvailableException.class, optimizer::getGradientListener);

        // set new listener
        optimizer.setGradientListener(gradientListener);

        // check correctness
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.getGradientListener(), gradientListener);
    }

    @Test
    void testGetSetListenerAndAvailability() throws LockedException, NotAvailableException {

        final var optimizer = new QuasiNewtonMultiOptimizer();

        // initially listener is not available
        assertFalse(optimizer.isListenerAvailable());
        assertThrows(NotAvailableException.class, optimizer::getListener);

        // set new listener
        optimizer.setListener(listener);

        // check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
    }

    @Test
    void testIsLocked() {

        final var optimizer = new QuasiNewtonMultiOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    void testGetSetOnIterationCompletedListener() throws LockedException {
        final var optimizer = new QuasiNewtonMultiOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testMinimize() throws Throwable {

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        final var startPoint = new double[ndims];
        minimum = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);

        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        final var optimizer = new QuasiNewtonMultiOptimizer(listener, gradientListener, startPoint,
                QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE);
        optimizer.setOnIterationCompletedListener(this);

        // check status
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isLocked());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);

        reset();
        assertEquals(0, iterations);

        // minimize
        optimizer.minimize();

        // check status
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isLocked());
        assertTrue(optimizer.isResultAvailable());
        assertTrue(optimizer.getIterations() >= 0);
        assertTrue(iterations >= 0);

        // pick result
        final var result = optimizer.getResult();
        final var valueAtResult = optimizer.getEvaluationAtResult();

        // check correctness of result

        // compute function value at true minimum
        final var valueAtMinimum = listener.evaluate(minimum);

        // compare values
        assertNotNull(result);
        assertEquals(valueAtResult, valueAtMinimum, ABSOLUTE_ERROR);
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
