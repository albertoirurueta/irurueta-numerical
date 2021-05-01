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

import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class GoldenSingleOptimizerTest implements
        SingleDimensionFunctionEvaluatorListener,
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
    public void testConstructor() throws NotAvailableException,
            InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double middleEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(middleEvalPoint,
                MAX_EVAL_POINT);
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        // test 1st constructor
        GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(),
                GoldenSingleOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(optimizer.getMinEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MAX_EVAL_POINT, 0.0);
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areBracketEvaluationsAvailable());
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

        // test 2nd constructor
        optimizer = new GoldenSingleOptimizer(this, minEvalPoint,
                middleEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(optimizer.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(), middleEvalPoint,
                0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.areBracketEvaluationsAvailable());
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

        // Force InvalidBracketRangeException
        optimizer = null;
        try {
            optimizer = new GoldenSingleOptimizer(this, maxEvalPoint,
                    middleEvalPoint, minEvalPoint, tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, minEvalPoint,
                    maxEvalPoint, middleEvalPoint, tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, maxEvalPoint,
                    minEvalPoint, middleEvalPoint, tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, middleEvalPoint,
                    minEvalPoint, maxEvalPoint, tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, middleEvalPoint,
                    maxEvalPoint, minEvalPoint, tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }

        try {
            optimizer = new GoldenSingleOptimizer(this, maxEvalPoint,
                    middleEvalPoint, minEvalPoint, -tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, minEvalPoint,
                    maxEvalPoint, middleEvalPoint, -tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, maxEvalPoint,
                    minEvalPoint, middleEvalPoint, -tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, middleEvalPoint,
                    minEvalPoint, maxEvalPoint, -tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer = new GoldenSingleOptimizer(this, middleEvalPoint,
                    maxEvalPoint, minEvalPoint, -tolerance);
        } catch (final InvalidBracketRangeException ignore) {
        }

        // Force IllegalArgumentException
        try {
            optimizer = new GoldenSingleOptimizer(this, minEvalPoint,
                    middleEvalPoint, maxEvalPoint, -tolerance);
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(optimizer);
    }

    @Test
    public void testGetSetTolerance() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertEquals(optimizer.getTolerance(),
                GoldenSingleOptimizer.DEFAULT_TOLERANCE, 0.0);

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
    public void testGetSetBracketAndAvailability() throws NotAvailableException,
            LockedException, InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double middleEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(middleEvalPoint,
                MAX_EVAL_POINT);

        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertTrue(optimizer.isBracketAvailable());

        assertEquals(optimizer.getMinEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(),
                GoldenSingleOptimizer.DEFAULT_MAX_EVAL_POINT, 0.0);

        // set new bracket
        optimizer.setBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);

        // check correctness
        assertEquals(optimizer.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(), middleEvalPoint,
                0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(), maxEvalPoint, 0.0);

        // Force InvalidBracketRangeException
        try {
            optimizer.setBracket(maxEvalPoint, middleEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer.setBracket(minEvalPoint, maxEvalPoint, middleEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer.setBracket(maxEvalPoint, minEvalPoint, middleEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer.setBracket(middleEvalPoint, minEvalPoint, maxEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        try {
            optimizer.setBracket(middleEvalPoint, maxEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
    }

    @Test
    public void testGetEvaluationsAndEvaluateBracket()
            throws InvalidBracketRangeException, LockedException,
            NotReadyException, OptimizationException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set bracket
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double middleEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(middleEvalPoint,
                MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer(this,
                minEvalPoint, middleEvalPoint, maxEvalPoint,
                GoldenSingleOptimizer.DEFAULT_TOLERANCE);

        // attempting to retrieve evaluations fails because although bracket is
        // available, evaluations have not yet been computed
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }

        assertFalse(optimizer.areBracketEvaluationsAvailable());

        // we compute evaluations
        optimizer.evaluateBracket();

        assertTrue(optimizer.areBracketEvaluationsAvailable());
    }

    @Test
    public void testComputeBracket() throws LockedException, NotReadyException,
            OptimizationException, NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();
        optimizer.setListener(this);

        assertTrue(optimizer.isBracketAvailable());
        assertFalse(optimizer.areBracketEvaluationsAvailable());

        optimizer.computeBracket();

        assertTrue(optimizer.isBracketAvailable());
        assertTrue(optimizer.areBracketEvaluationsAvailable());

        // after computing bracket we can only ensure that ax < bx < cx and also
        // that fa > fb and fc > fb
        assertTrue(optimizer.getMinEvaluationPoint() <=
                optimizer.getMiddleEvaluationPoint());
        assertTrue(optimizer.getMiddleEvaluationPoint() <=
                optimizer.getMaxEvaluationPoint());
        assertTrue(optimizer.getEvaluationAtMin() >=
                optimizer.getEvaluationAtMiddle());
        assertTrue(optimizer.getEvaluationAtMax() >=
                optimizer.getEvaluationAtMiddle());

        // also bracket limits must surround the real minimum location, that is
        // ax < minimum and cx > minimum
        assertTrue(optimizer.getMinEvaluationPoint() <= minimum);
        assertTrue(optimizer.getMaxEvaluationPoint() >= minimum);
    }

    @Test
    public void testGetSetListenerAndAvailability() throws LockedException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertFalse(optimizer.isListenerAvailable());
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(optimizer.isReady());

        // set listener
        optimizer.setListener(this);

        // check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), this);
        assertTrue(optimizer.isReady());
    }

    @Test
    public void testIsLocked() {
        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertFalse(optimizer.isLocked());
    }

    @Test
    public void testIsReady() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        // set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        // set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);


        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertFalse(optimizer.isReady());

        // Set listener
        optimizer.setListener(this);

        // check correctness
        assertTrue(optimizer.isReady());
    }

    @Test
    public void testGetSetOnIterationCompletedListener() throws LockedException {
        final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();

        assertNull(optimizer.getOnIterationCompletedListener());

        // set new value
        optimizer.setOnIterationCompletedListener(this);

        // check
        assertSame(this, optimizer.getOnIterationCompletedListener());
    }

    @Test
    public void testMinimizeGetResultAndAvailability() throws Throwable {

        for (int i = 0; i < TIMES; i++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            // set minimum to be estimated
            minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

            // set offset
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            // set width
            width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

            final GoldenSingleOptimizer optimizer = new GoldenSingleOptimizer();
            optimizer.setListener(this);
            optimizer.setOnIterationCompletedListener(this);

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

            // minimize
            if (minimum > 0.0) {
                optimizer.setBracket(0.8 * minimum, 1.1 * minimum,
                        1.5 * minimum);
            } else {
                optimizer.setBracket(1.5 * minimum, 1.1 * minimum,
                        0.8 * minimum);
            }

            reset();
            assertEquals(0, iterations);

            optimizer.minimize();

            // test correctness
            assertTrue(optimizer.isResultAvailable());
            assertEquals(optimizer.getResult(), minimum, ABSOLUTE_ERROR);
            assertEquals(optimizer.getEvaluationAtResult(), evaluate(
                    optimizer.getResult()), 0.0);
            assertTrue(iterations > 0);
        }
    }

    @Override
    public double evaluate(final double point) throws EvaluationException {
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
