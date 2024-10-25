/*
 * Copyright (C) 2013 Alberto Irurueta Carro (alberto@irurueta.com)
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
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class RANSACRobustEstimatorTest {

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

    private static final int PERCENTAGE_OUTLIER = 15;

    private static final int NUM_PARAMS = 2;

    private static final int TIMES = 100;

    @Test
    void testConstructor() {
        // test empty constructor
        var estimator = new RANSACRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(RobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());
        assertFalse(estimator.isReady());
        assertEquals(RANSACRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(RANSACRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(RANSACRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getNIters());
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
        assertEquals(RANSACRobustEstimator.DEFAULT_COMPUTE_AND_KEEP_INLIERS,
                estimator.isComputeAndKeepInliersEnabled());
        assertEquals(RANSACRobustEstimator.DEFAULT_COMPUTE_AND_KEEP_RESIDUALS,
                estimator.isComputeAndKeepResidualsEnabled());

        // test constructor with listener
        final var randomizer = new UniformRandomizer();
        final var numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var listener = new TestRANSACRobustEstimatorListener(numSamples, PERCENTAGE_OUTLIER, THRESHOLD);
        estimator = new RANSACRobustEstimator<>(listener);
        assertEquals(estimator.getListener(), listener);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertEquals(RobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());
        assertEquals(estimator.isReady(), listener.isReady());
        assertTrue(estimator.isReady());
        assertEquals(RANSACRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(RANSACRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(RANSACRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getNIters());
        assertNull(estimator.getBestResult());
        assertNull(estimator.getInliersData());
        assertNull(estimator.getBestInliersData());
        assertEquals(RANSACRobustEstimator.DEFAULT_COMPUTE_AND_KEEP_INLIERS,
                estimator.isComputeAndKeepInliersEnabled());
        assertEquals(RANSACRobustEstimator.DEFAULT_COMPUTE_AND_KEEP_RESIDUALS,
                estimator.isComputeAndKeepResidualsEnabled());
    }

    @Test
    void testGetSetListenerAvailabilityAndIsReady() throws LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();
        assertNull(estimator.getListener());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());

        // set listener
        final var randomizer = new UniformRandomizer();
        final var numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var listener = new TestRANSACRobustEstimatorListener(numSamples, PERCENTAGE_OUTLIER, THRESHOLD);

        estimator.setListener(listener);

        // check correctness
        assertEquals(listener, estimator.getListener());
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());
    }

    @Test
    void testGetSetProgressDelta() throws IllegalArgumentException, LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();
        assertEquals(RobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();
        final var progressDelta = randomizer.nextFloat(0.0f, 1.0f);
        estimator.setProgressDelta(progressDelta);

        // check correctness
        assertEquals(estimator.getProgressDelta(), progressDelta, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setProgressDelta(-1.0f));
        assertThrows(IllegalArgumentException.class, () -> estimator.setProgressDelta(2.0f));
    }

    @Test
    void testGetSetConfidence() throws IllegalArgumentException, LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();
        assertEquals(RANSACRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();
        final var confidence = randomizer.nextDouble(0.0, 1.0);
        estimator.setConfidence(confidence);

        // check correctness
        assertEquals(confidence, estimator.getConfidence(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setConfidence(-1.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setConfidence(2.0));
    }

    @Test
    void testGetSetMaxIterations() throws IllegalArgumentException, LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();
        assertEquals(RANSACRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());

        // set new value
        final var randomizer = new UniformRandomizer();
        final var maxIterations = randomizer.nextInt(MIN_MAX_ITERATIONS, MAX_MAX_ITERATIONS);
        estimator.setMaxIterations(maxIterations);

        // check correctness
        assertEquals(maxIterations, estimator.getMaxIterations());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setMaxIterations(0));
    }

    @Test
    void testIsSetComputeAndKeepInliersEnabled() throws LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();

        // check default value
        assertFalse(estimator.isComputeAndKeepInliersEnabled());

        // set new value
        estimator.setComputeAndKeepInliersEnabled(true);

        // check correctness
        assertTrue(estimator.isComputeAndKeepInliersEnabled());
    }

    @Test
    void testIsSetComputeAndKeepResidualsEnabled() throws LockedException {
        final var estimator = new RANSACRobustEstimator<double[]>();

        // check default value
        assertFalse(estimator.isComputeAndKeepResidualsEnabled());

        // set new value
        estimator.setComputeAndKeepResidualsEnabled(true);

        // check correctness
        assertTrue(estimator.isComputeAndKeepResidualsEnabled());
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException, RobustEstimatorException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            final var numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final var listener = new TestRANSACRobustEstimatorListener(numSamples, PERCENTAGE_OUTLIER, THRESHOLD);
            final var estimator = new RANSACRobustEstimator<double[]>();

            estimator.setComputeAndKeepInliersEnabled(false);
            estimator.setComputeAndKeepResidualsEnabled(false);

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // set listener
            estimator.setListener(listener);
            listener.reset();
            assertEquals(0, listener.getStartCounter());
            assertEquals(0, listener.getEndCounter());
            assertFalse(estimator.isLocked());

            // estimate
            final var params = estimator.estimate();

            // check status after estimation
            assertFalse(estimator.isLocked());
            assertEquals(1, listener.getStartCounter());
            assertEquals(1, listener.getEndCounter());

            // check correctness of estimation
            assertEquals(params.length, listener.getParams().length);
            assertEquals(NUM_PARAMS, params.length);

            var failed = false;
            for (var i = 0; i < params.length; i++) {
                if (Math.abs(params[i] - listener.getParams()[i]) > 10.0 * ABSOLUTE_ERROR) {
                    failed = true;
                    break;
                }
            }

            if (failed) {
                continue;
            }
            assertArrayEquals(params, listener.getParams(), 10.0 * ABSOLUTE_ERROR);

            assertNull(estimator.getBestInliersData());
            assertNull(estimator.getInliersData());

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testEstimateWithInliersData() throws LockedException, NotReadyException, RobustEstimatorException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            final var numSamples = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final var listener = new TestRANSACRobustEstimatorListener(numSamples, PERCENTAGE_OUTLIER, THRESHOLD);
            final var estimator = new RANSACRobustEstimator<double[]>();

            estimator.setComputeAndKeepInliersEnabled(true);
            estimator.setComputeAndKeepResidualsEnabled(true);

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // set listener
            estimator.setListener(listener);
            listener.reset();
            assertEquals(0, listener.getStartCounter());
            assertEquals(0, listener.getEndCounter());
            assertFalse(estimator.isLocked());

            // estimate
            final var params = estimator.estimate();

            // check status after estimation
            assertFalse(estimator.isLocked());
            assertEquals(1, listener.getStartCounter());
            assertEquals(1, listener.getEndCounter());

            // check correctness of estimation
            assertEquals(params.length, listener.getParams().length);
            assertEquals(NUM_PARAMS, params.length);

            var failed = false;
            for (var i = 0; i < params.length; i++) {
                if (Math.abs(params[i] - listener.getParams()[i]) > ABSOLUTE_ERROR) {
                    failed = true;
                    break;
                }
            }
            if (failed) {
                continue;
            }
            assertArrayEquals(params, listener.getParams(), ABSOLUTE_ERROR);

            final var inliersData = estimator.getBestInliersData();
            assertNotNull(inliersData);
            assertSame(inliersData, estimator.getInliersData());
            assertTrue(inliersData.getNumInliers() > 0);
            assertNotNull(inliersData.getInliers());
            assertNotNull(inliersData.getResiduals());

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    private static double[] computeParams() {
        // we will estimate parameters a and b for equation y = a*x + b
        final var randomizer = new UniformRandomizer();
        final var params = new double[NUM_PARAMS];
        // a parameter
        params[0] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        // b parameter
        params[1] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        return params;
    }

    private static void computeSamples(final double[] params, final int numSamples, final int percentageOutliers,
                                       final double[] ys, final double[] xs) {

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < numSamples; i++) {
            // compute x values
            xs[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            // compute exact y values
            ys[i] = params[0] * xs[i] + params[1];
            if (randomizer.nextInt(0, 100) < percentageOutliers) {
                // is outlier, so we add a certain amount of error
                final var error = randomizer.nextDouble(MIN_ERROR, MAX_ERROR);
                ys[i] += error;
            }
        }
    }

    private static class TestRANSACRobustEstimatorListener implements RANSACRobustEstimatorListener<double[]> {

        private final double[] params;
        private final double[] xs;
        private final double[] ys;
        private final int numSamples;
        private final double threshold;

        private int startCounter;
        private int endCounter;
        private float previousProgress;

        TestRANSACRobustEstimatorListener(final int numSamples, final int percentageOutliers, final double threshold) {
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
        public void estimatePreliminarSolutions(final int[] samplesIndices, final List<double[]> solutions) {

            if (samplesIndices.length != NUM_PARAMS) {
                throw new IllegalArgumentException();
            }
            final var index1 = samplesIndices[0];
            final var index2 = samplesIndices[1];

            final var y1 = ys[index1];
            final var y2 = ys[index2];
            final var x1 = xs[index1];
            final var x2 = xs[index2];

            final var a = (y2 - y1) / (x2 - x1);
            final var b = y1 - a * x1;

            final var solution = new double[NUM_PARAMS];
            solution[0] = a;
            solution[1] = b;

            solutions.add(solution);
        }

        @Override
        public double computeResidual(final double[] currentEstimation, final int i) {
            final var a = currentEstimation[0];
            final var b = currentEstimation[1];

            final var estimatedY = a * xs[i] + b;
            final var y = ys[i];

            return Math.abs(estimatedY - y);
        }

        @Override
        public boolean isReady() {
            return params != null && xs != null && ys != null;
        }

        @Override
        public void onEstimateStart(final RobustEstimator<double[]> estimator) {
            testIsLocked((RANSACRobustEstimator<double[]>) estimator);
            startCounter++;
        }

        @Override
        public void onEstimateEnd(final RobustEstimator<double[]> estimator) {
            testIsLocked((RANSACRobustEstimator<double[]>) estimator);
            endCounter++;
        }

        @Override
        public void onEstimateNextIteration(final RobustEstimator<double[]> estimator, final int iteration) {
            final var ransacEstimator = (RANSACRobustEstimator<double[]>) estimator;
            testIsLocked(ransacEstimator);
            assertTrue(iteration > 0);
            assertTrue(ransacEstimator.getNIters() >= 0);
            assertNotNull(ransacEstimator.getBestResult());
        }

        @Override
        public void onEstimateProgressChange(final RobustEstimator<double[]> estimator, final float progress) {
            testIsLocked((RANSACRobustEstimator<double[]>) estimator);
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

        private void testIsLocked(final RANSACRobustEstimator<double[]> estimator) {
            assertTrue(estimator.isLocked());
            // test that estimator cannot be modified while locked
            assertThrows(LockedException.class, () -> estimator.setConfidence(0.5));
            assertThrows(LockedException.class, () -> estimator.setListener(this));
            assertThrows(LockedException.class, () -> estimator.setMaxIterations(1));
            assertThrows(LockedException.class, () -> estimator.setProgressDelta(0.5f));
        }

        private void reset() {
            startCounter = endCounter = 0;
            previousProgress = 0.0f;
        }
    }
}