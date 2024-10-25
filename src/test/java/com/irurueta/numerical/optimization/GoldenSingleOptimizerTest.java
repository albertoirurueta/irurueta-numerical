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

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class GoldenSingleOptimizerTest implements SingleDimensionFunctionEvaluatorListener,
        OnIterationCompletedListener {

    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;

    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 1e-5;

    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;

    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;

    public static final double ABSOLUTE_ERROR = 1e-5;
    public static final int TIMES = 100;

    private double minimum;
    private double offset;
    private double width;

    private int iterations = 0;

    @Test
    void testConstructor() throws NotAvailableException, InvalidBracketRangeException {

        final var randomizer = new UniformRandomizer();
        final var minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var middleEvalPoint = randomizer.nextDouble(minEvalPoint, MAX_EVAL_POINT);
        final var maxEvalPoint = randomizer.nextDouble(middleEvalPoint, MAX_EVAL_POINT);
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        // test 1st constructor
        var optimizer = new GoldenSingleOptimizer();
        assertNotNull(optimizer);
        assertEquals(GoldenSingleOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);
        assertFalse(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(GoldenSingleOptimizer.DEFAULT_MIN_EVAL_POINT, optimizer.getMinEvaluationPoint(), 0.0);
        assertEquals(GoldenSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, optimizer.getMiddleEvaluationPoint(), 0.0);
        assertEquals(GoldenSingleOptimizer.DEFAULT_MAX_EVAL_POINT, optimizer.getMaxEvaluationPoint(), 0.0);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMin);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMiddle);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMax);
        assertFalse(optimizer.areBracketEvaluationsAvailable());
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // test 2nd constructor
        optimizer = new GoldenSingleOptimizer(this, minEvalPoint, middleEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(optimizer);
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);
        assertTrue(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(minEvalPoint, optimizer.getMinEvaluationPoint(), 0.0);
        assertEquals(middleEvalPoint, optimizer.getMiddleEvaluationPoint(), 0.0);
        assertEquals(maxEvalPoint, optimizer.getMaxEvaluationPoint(), 0.0);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMin);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMiddle);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMax);
        assertFalse(optimizer.areBracketEvaluationsAvailable());
        assertEquals(this, optimizer.getListener());
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        assertThrows(NotAvailableException.class, optimizer::getResult);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);
        assertFalse(optimizer.isLocked());
        assertNull(optimizer.getOnIterationCompletedListener());

        // Force InvalidBracketRangeException
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, maxEvalPoint,
                middleEvalPoint, minEvalPoint, tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, minEvalPoint,
                maxEvalPoint, middleEvalPoint, tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, maxEvalPoint,
                minEvalPoint, middleEvalPoint, tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, middleEvalPoint,
                minEvalPoint, maxEvalPoint, tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, middleEvalPoint,
                maxEvalPoint, minEvalPoint, tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, maxEvalPoint,
                middleEvalPoint, minEvalPoint, -tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, minEvalPoint,
                maxEvalPoint, middleEvalPoint, -tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, maxEvalPoint,
                minEvalPoint, middleEvalPoint, -tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, middleEvalPoint,
                minEvalPoint, maxEvalPoint, -tolerance));
        assertThrows(InvalidBracketRangeException.class, () -> new GoldenSingleOptimizer(this, middleEvalPoint,
                maxEvalPoint, minEvalPoint, -tolerance));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new GoldenSingleOptimizer(this, minEvalPoint,
                middleEvalPoint, maxEvalPoint, -tolerance));
    }

    @Test
    void testGetSetTolerance() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var optimizer = new GoldenSingleOptimizer();

        assertEquals(GoldenSingleOptimizer.DEFAULT_TOLERANCE, optimizer.getTolerance(), 0.0);

        // set new tolerance
        optimizer.setTolerance(tolerance);

        // check correctness
        assertEquals(tolerance, optimizer.getTolerance(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> optimizer.setTolerance(-tolerance));
    }

    @Test
    void testGetSetBracketAndAvailability() throws NotAvailableException, LockedException,
            InvalidBracketRangeException {

        final var randomizer = new UniformRandomizer();
        final var minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var middleEvalPoint = randomizer.nextDouble(minEvalPoint, MAX_EVAL_POINT);
        final var maxEvalPoint = randomizer.nextDouble(middleEvalPoint, MAX_EVAL_POINT);

        final var optimizer = new GoldenSingleOptimizer();

        assertTrue(optimizer.isBracketAvailable());

        assertEquals(GoldenSingleOptimizer.DEFAULT_MIN_EVAL_POINT, optimizer.getMinEvaluationPoint(), 0.0);
        assertEquals(GoldenSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, optimizer.getMiddleEvaluationPoint(), 0.0);
        assertEquals(GoldenSingleOptimizer.DEFAULT_MAX_EVAL_POINT, optimizer.getMaxEvaluationPoint(), 0.0);

        // set new bracket
        optimizer.setBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);

        // check correctness
        assertEquals(optimizer.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(), middleEvalPoint, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(), maxEvalPoint, 0.0);

        // Force InvalidBracketRangeException
        assertThrows(InvalidBracketRangeException.class, () -> optimizer.setBracket(maxEvalPoint, middleEvalPoint,
                minEvalPoint));
        assertThrows(InvalidBracketRangeException.class, () -> optimizer.setBracket(minEvalPoint, maxEvalPoint,
                middleEvalPoint));
        assertThrows(InvalidBracketRangeException.class, () -> optimizer.setBracket(maxEvalPoint, minEvalPoint,
                middleEvalPoint));
        assertThrows(InvalidBracketRangeException.class, () -> optimizer.setBracket(middleEvalPoint, minEvalPoint,
                maxEvalPoint));
        assertThrows(InvalidBracketRangeException.class, () -> optimizer.setBracket(middleEvalPoint, maxEvalPoint,
                minEvalPoint));
    }

    @Test
    void testGetEvaluationsAndEvaluateBracket() throws InvalidBracketRangeException, LockedException, NotReadyException,
            OptimizationException {

        final var randomizer = new UniformRandomizer();

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set bracket
        final var minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var middleEvalPoint = randomizer.nextDouble(minEvalPoint, MAX_EVAL_POINT);
        final var maxEvalPoint = randomizer.nextDouble(middleEvalPoint, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final var optimizer = new GoldenSingleOptimizer(this, minEvalPoint, middleEvalPoint, maxEvalPoint,
                GoldenSingleOptimizer.DEFAULT_TOLERANCE);

        // attempting to retrieve evaluations fails because although bracket is
        // available, evaluations have not yet been computed
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMin);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMiddle);
        assertThrows(NotAvailableException.class, optimizer::getEvaluationAtMax);

        assertFalse(optimizer.areBracketEvaluationsAvailable());

        // we compute evaluations
        optimizer.evaluateBracket();

        assertTrue(optimizer.areBracketEvaluationsAvailable());
    }

    @Test
    void testComputeBracket() throws LockedException, NotReadyException, OptimizationException, NotAvailableException {

        final var randomizer = new UniformRandomizer();

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final var optimizer = new GoldenSingleOptimizer();
        optimizer.setListener(this);

        assertTrue(optimizer.isBracketAvailable());
        assertFalse(optimizer.areBracketEvaluationsAvailable());

        optimizer.computeBracket();

        assertTrue(optimizer.isBracketAvailable());
        assertTrue(optimizer.areBracketEvaluationsAvailable());

        // after computing bracket we can only ensure that ax < bx < cx and also
        // that fa > fb and fc > fb
        assertTrue(optimizer.getMinEvaluationPoint() <= optimizer.getMiddleEvaluationPoint());
        assertTrue(optimizer.getMiddleEvaluationPoint() <= optimizer.getMaxEvaluationPoint());
        assertTrue(optimizer.getEvaluationAtMin() >= optimizer.getEvaluationAtMiddle());
        assertTrue(optimizer.getEvaluationAtMax() >= optimizer.getEvaluationAtMiddle());

        // also bracket limits must surround the real minimum location, that is
        // ax < minimum and cx > minimum
        assertTrue(optimizer.getMinEvaluationPoint() <= minimum);
        assertTrue(optimizer.getMaxEvaluationPoint() >= minimum);
    }

    @Test
    void testGetSetListenerAndAvailability() throws LockedException, NotAvailableException {

        final var randomizer = new UniformRandomizer();

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final var optimizer = new GoldenSingleOptimizer();

        assertFalse(optimizer.isListenerAvailable());
        assertThrows(NotAvailableException.class, optimizer::getListener);
        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(this);

        // check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), this);
        assertTrue(optimizer.isReady());
    }

    @Test
    void testIsLocked() {
        final var optimizer = new GoldenSingleOptimizer();
        assertFalse(optimizer.isLocked());
    }

    @Test
    void testIsReady() throws LockedException {

        final var randomizer = new UniformRandomizer();

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);


        final var optimizer = new GoldenSingleOptimizer();

        assertFalse(optimizer.isReady());

        // Set listener
        optimizer.setListener(this);

        // check correctness
        assertTrue(optimizer.isReady());
    }

    @Test
    void testGetSetOnIterationCompletedListener() throws LockedException {
        final var optimizer = new GoldenSingleOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    void testMinimizeGetResultAndAvailability() throws Throwable {

        for (var i = 0; i < TIMES; i++) {
            final var randomizer = new UniformRandomizer();

            // set minimum to be estimated
            minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

            // set offset
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            // set width
            width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

            final var optimizer = new GoldenSingleOptimizer();
            optimizer.setListener(this);
            optimizer.setOnIterationCompletedListener(this);

            assertFalse(optimizer.isResultAvailable());
            assertThrows(NotAvailableException.class, optimizer::getResult);
            assertThrows(NotAvailableException.class, optimizer::getEvaluationAtResult);

            // minimize
            if (minimum > 0.0) {
                optimizer.setBracket(0.8 * minimum, 1.1 * minimum, 1.5 * minimum);
            } else {
                optimizer.setBracket(1.5 * minimum, 1.1 * minimum, 0.8 * minimum);
            }

            reset();
            assertEquals(0, iterations);

            optimizer.minimize();

            // test correctness
            assertTrue(optimizer.isResultAvailable());
            assertEquals(optimizer.getResult(), minimum, ABSOLUTE_ERROR);
            assertEquals(optimizer.getEvaluationAtResult(), evaluate(optimizer.getResult()), 0.0);
            assertTrue(iterations > 0);
        }
    }

    @Override
    public double evaluate(final double point) {
        return (point - minimum) * (point - minimum) / width + offset;
    }

    @Override
    public void onIterationCompleted(Optimizer optimizer, int iteration, Integer maxIterations) {
        assertNotNull(optimizer);
        assertTrue(optimizer.isLocked());
        assertTrue(iteration >= 0);
        assertNull(maxIterations);
        iterations++;
    }

    private void reset() {
        iterations = 0;
    }
}
