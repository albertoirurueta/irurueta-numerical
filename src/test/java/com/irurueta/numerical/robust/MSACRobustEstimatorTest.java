/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.robust.MSACRobustEstimator.MSACInliersData;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

public class MSACRobustEstimatorTest {

    private static final int MIN_POINTS = 500;
    private static final int MAX_POINTS = 1000;

    private static final double THRESHOLD = 1e-3;

    private static final double MIN_ERROR = 1e-5;
    private static final double MAX_ERROR = 1.0;

    private static final int MIN_MAX_ITERATIONS = 500;
    private static final int MAX_MAX_ITERATIONS = 5000;

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final double ABSOLUTE_ERROR = 1e-6;

    private static final int PERCENTAGE_OUTLIER = 20;

    private static final int NUM_PARAMS = 2;

    private static final int TIMES = 100;

    @Test
    public void testConstructor() {
        // test empty constructor
        MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getConfidence(),
                MSACRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                MSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.getNIters(),
                MSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestResultInliersData());
        assertNull(estimator.getBestNumberInliersData());

        // test constructor with listener
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final TestMSACRobustEstimatorListener listener =
                new TestMSACRobustEstimatorListener(numSamples,
                        PERCENTAGE_OUTLIER, THRESHOLD);
        estimator = new MSACRobustEstimator<>(listener);
        assertSame(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.MSAC);
        assertTrue(estimator.isReady());
        assertEquals(estimator.getConfidence(),
                MSACRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                MSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.getNIters(),
                MSACRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestResultInliersData());
        assertNull(estimator.getBestNumberInliersData());
    }

    @Test
    public void testGetSetListenerAvailabilityAndIsReady()
            throws LockedException {
        final MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());

        // set listener
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final TestMSACRobustEstimatorListener listener =
                new TestMSACRobustEstimatorListener(numSamples,
                        PERCENTAGE_OUTLIER, THRESHOLD);

        estimator.setListener(listener);

        // check correctness
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());
    }

    @Test
    public void testGetSetProgressDelta() throws IllegalArgumentException,
            LockedException {
        final MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();
        assertEquals(estimator.getProgressDelta(),
                RobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final float progressDelta = randomizer.nextFloat(0.0f, 1.0f);
        estimator.setProgressDelta(progressDelta);

        // check correctness
        assertEquals(estimator.getProgressDelta(), progressDelta, 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setProgressDelta(-1.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setProgressDelta(2.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetConfidence() throws IllegalArgumentException,
            LockedException {
        final MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();
        assertEquals(estimator.getConfidence(),
                MSACRobustEstimator.DEFAULT_CONFIDENCE, 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double confidence = randomizer.nextDouble(0.0, 1.0);
        estimator.setConfidence(confidence);

        // check correctness
        assertEquals(estimator.getConfidence(), confidence, 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setConfidence(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setConfidence(2.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetMaxIterations() throws IllegalArgumentException,
            LockedException {
        final MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();
        assertEquals(estimator.getMaxIterations(),
                MSACRobustEstimator.DEFAULT_MAX_ITERATIONS);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int maxIterations = randomizer.nextInt(MIN_MAX_ITERATIONS,
                MAX_MAX_ITERATIONS);
        estimator.setMaxIterations(maxIterations);

        // check correctness
        assertEquals(estimator.getMaxIterations(), maxIterations);

        // Force IllegalArgumentException
        try {
            estimator.setMaxIterations(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testEstimate() throws LockedException, NotReadyException,
            RobustEstimatorException {
        for (int i = 0; i < TIMES; i++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            final int numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final TestMSACRobustEstimatorListener listener =
                    new TestMSACRobustEstimatorListener(numSamples,
                            PERCENTAGE_OUTLIER, THRESHOLD);
            final MSACRobustEstimator<double[]> estimator = new MSACRobustEstimator<>();

            // Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (final NotReadyException ignore) {
            }

            // set listener
            estimator.setListener(listener);
            listener.reset();
            assertEquals(listener.getStartCounter(), 0);
            assertEquals(listener.getEndCounter(), 0);
            assertFalse(estimator.isLocked());

            // estimate
            final double[] params = estimator.estimate();

            assertNotNull(estimator.getBestResult());
            assertNotNull(estimator.getBestResultInliersData());
            assertNotNull(estimator.getBestNumberInliersData());

            // check status after estimation
            assertFalse(estimator.isLocked());
            assertEquals(listener.getStartCounter(), 1);
            assertEquals(listener.getEndCounter(), 1);

            // check correctness of estimation
            assertEquals(params.length, listener.getParams().length);
            assertEquals(params.length, NUM_PARAMS);

            assertArrayEquals(params, listener.getParams(), ABSOLUTE_ERROR);

            assertNotNull(estimator.getBestResultInliersData());

            final MSACInliersData inliersData = estimator.getBestNumberInliersData();
            assertNotNull(inliersData);
            assertSame(inliersData, estimator.getInliersData());
            assertTrue(inliersData.getNumInliers() > 0);
            assertNotNull(inliersData.getInliers());
            assertNotNull(inliersData.getResiduals());
        }
    }

    private double[] computeParams() {
        // we will estimate parameters a and b for equation y = a*x + b
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] params = new double[NUM_PARAMS];
        // a parameter
        params[0] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        // b parameter
        params[1] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        return params;
    }

    private void computeSamples(final double[] params, final int numSamples,
                                final int percentageOutliers, final double[] ys, final double[] xs) {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for (int i = 0; i < numSamples; i++) {
            // compute x values
            xs[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            // compute exact y values
            ys[i] = params[0] * xs[i] + params[1];
            if (randomizer.nextInt(0, 100) < percentageOutliers) {
                // is outlier, so we add a certain amount of error
                final double error = randomizer.nextDouble(MIN_ERROR, MAX_ERROR);
                ys[i] += error;
            }
        }
    }

    public class TestMSACRobustEstimatorListener implements
            MSACRobustEstimatorListener<double[]> {

        private final double[] params;
        private final double[] xs;
        private final double[] ys;
        private final int numSamples;
        private final double threshold;

        private int startCounter;
        private int endCounter;
        private float previousProgress;

        TestMSACRobustEstimatorListener(final int numSamples,
                                        final int percentageOutliers, final double threshold) {
            this.numSamples = numSamples;
            params = computeParams();
            xs = new double[numSamples];
            ys = new double[numSamples];
            computeSamples(params, numSamples, percentageOutliers, ys, xs);
            this.threshold = threshold;
            reset();
        }

        public double[] getParams() {
            return params;
        }

        @Override
        public int getTotalSamples() {
            return numSamples;
        }

        @Override
        public int getSubsetSize() {
            // only two matches x and y are required to determine a and b as
            // follows:
            // y1 = a* x1 + b --> b = y1 - a * x1
            // y2 = a * x2 + b --> y2 = a * x2 + y1 - a * x1 -->
            // y2 - y1 = (x2 - x1)*a --> a = (y2 - y1) / (x2 - x1)

            // Hence:
            // a = (y2 - y1) / (x2 - x1)
            // b = y1 - (y2 - y1) / (x2 - x1) * x1
            return NUM_PARAMS;
        }

        @Override
        public double getThreshold() {
            return threshold;
        }

        @Override
        public void estimatePreliminarSolutions(final int[] samplesIndices,
                                                final List<double[]> solutions) {

            if (samplesIndices.length != NUM_PARAMS) {
                throw new IllegalArgumentException();
            }
            final int index1 = samplesIndices[0];
            final int index2 = samplesIndices[1];

            final double y1 = ys[index1];
            final double y2 = ys[index2];
            final double x1 = xs[index1];
            final double x2 = xs[index2];

            final double a = (y2 - y1) / (x2 - x1);
            final double b = y1 - a * x1;

            final double[] solution = new double[NUM_PARAMS];
            solution[0] = a;
            solution[1] = b;

            solutions.add(solution);
        }

        @Override
        public double computeResidual(final double[] currentEstimation, final int i) {
            final double a = currentEstimation[0];
            final double b = currentEstimation[1];

            final double estimatedY = a * xs[i] + b;
            final double y = ys[i];

            return Math.abs(estimatedY - y);
        }

        @Override
        public boolean isReady() {
            return params != null && xs != null && ys != null;
        }

        @Override
        public void onEstimateStart(final RobustEstimator<double[]> estimator) {
            testIsLocked((MSACRobustEstimator<double[]>) estimator);
            startCounter++;
        }

        @Override
        public void onEstimateEnd(final RobustEstimator<double[]> estimator) {
            testIsLocked((MSACRobustEstimator<double[]>) estimator);
            endCounter++;
        }

        @Override
        public void onEstimateNextIteration(final RobustEstimator<double[]> estimator,
                                            final int iteration) {
            final MSACRobustEstimator<double[]> ransacEstimator =
                    (MSACRobustEstimator<double[]>) estimator;
            testIsLocked(ransacEstimator);
            assertTrue(iteration > 0);
            assertTrue(ransacEstimator.getNIters() >= 0);
            assertNotNull(ransacEstimator.getBestResult());
        }

        @Override
        public void onEstimateProgressChange(
                final RobustEstimator<double[]> estimator, final float progress) {
            testIsLocked((MSACRobustEstimator<double[]>) estimator);
            assertTrue(progress >= 0.0f);
            assertTrue(progress <= 1.0f);
            assertTrue(progress >= previousProgress);
            previousProgress = progress;
        }

        int getStartCounter() {
            return startCounter;
        }

        int getEndCounter() {
            return endCounter;
        }

        private void testIsLocked(final MSACRobustEstimator<double[]> estimator) {
            assertTrue(estimator.isLocked());
            // test that estimator cannot be modified while locked
            try {
                estimator.setConfidence(0.5);
                fail("LockedException expected but not thrown");
            } catch (final LockedException ignore) {
            }
            try {
                estimator.setListener(this);
                fail("LockedException expected but not thrown");
            } catch (final LockedException ignore) {
            }
            try {
                estimator.setMaxIterations(1);
                fail("LockedException expected but not thrown");
            } catch (final LockedException ignore) {
            }
            try {
                estimator.setProgressDelta(0.5f);
                fail("LockedException expected but not thrown");
            } catch (final LockedException ignore) {
            }
        }

        public final void reset() {
            startCounter = endCounter = 0;
            previousProgress = 0.0f;
        }
    }
}
