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

import com.irurueta.statistics.UniformRandomizer;

import java.util.HashSet;
import java.util.Set;

/**
 * This class computes indices of subsets of samples using a uniform randomizer
 * to pick random samples as fast as possible.
 * This class ensures that samples are not repeated within a single subset.
 */
public class FastRandomSubsetSelector extends SubsetSelector {

    /**
     * Constant defining whether randomizer needs to be initialized with system
     * timer. By disabling this option it is ensured that the same result is
     * obtained on each execution.
     */
    public static final boolean DEFAULT_SEED_RANDOMIZER_WITH_TIME = false;

    /**
     * Randomizer to pick random indexes.
     */
    private UniformRandomizer randomizer;

    /**
     * Set containing selected indices on a given run.
     * This is kept around between executions to avoid excessive calls to
     * garbage collector, as random subset generation is likely to be called
     * a large number of times during robust estimations. Because of that,
     * if a robust estimator is capable of spanning multiple executions at once
     * on different threads, then each thread needs to have a different subset
     * selector instance.
     */
    private final Set<Integer> selectedIndices;

    /**
     * Constructor.
     *
     * @param numSamples number of samples to select subsets from
     * @throws IllegalArgumentException if provided number of samples is zero
     *                                  or negative
     */
    public FastRandomSubsetSelector(final int numSamples) {
        this(numSamples, DEFAULT_SEED_RANDOMIZER_WITH_TIME);
    }

    /**
     * Constructor.
     *
     * @param numSamples             number of samples to select subsets from.
     * @param seedRandomizerWithTime true if randomizer seed must be initialized
     *                               to system timer to obtain more random results.
     * @throws IllegalArgumentException if provided number of samples is zero
     *                                  or negative.
     */
    public FastRandomSubsetSelector(final int numSamples, final boolean seedRandomizerWithTime) {
        super(numSamples);
        createRandomizer(seedRandomizerWithTime);
        selectedIndices = new HashSet<>();
    }

    /**
     * Returns type of this subset selector.
     *
     * @return type of this subset selector.
     */
    @Override
    public SubsetSelectorType getType() {
        return SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR;
    }

    /**
     * Computes a random subset of indices within range of number of samples to
     * be used on robust estimators.
     *
     * @param subsetSize subset size to be computed. This value must be smaller
     *                   than total number of samples.
     * @param result     array containing indices to be picked. Provided array must
     *                   be at least of length subsetSize. The former subsetSize entries of the
     *                   array will be modified by this method.
     * @throws NotEnoughSamplesException  if subset size is greater than the
     *                                    total number of samples.
     * @throws InvalidSubsetSizeException if subset size is zero or if result
     *                                    array does not have at least a length of subsetSize.
     */
    @Override
    public void computeRandomSubsets(final int subsetSize, final int[] result) throws NotEnoughSamplesException,
            InvalidSubsetSizeException {
        if (subsetSize == 0) {
            throw new InvalidSubsetSizeException();
        }
        if (result.length < subsetSize) {
            throw new InvalidSubsetSizeException();
        }
        if (numSamples < subsetSize) {
            throw new NotEnoughSamplesException();
        }

        // On start set of selected indices is empty
        selectedIndices.clear();

        var counter = 0;
        int index;
        do {
            index = randomizer.nextInt(0, numSamples);

            // check whether this index has already been selected
            if (selectedIndices.contains(index)) {
                continue;
            }

            // if not selected, pick it now
            selectedIndices.add(index);
            result[counter] = index;
            counter++;
        } while (counter < subsetSize);

        selectedIndices.clear();
    }

    /**
     * Computes a random subset of indices within provided range of positions to
     * be used on robust estimators.
     *
     * @param minPos     minimum position to be picked. This value must be greater
     *                   or equal than zero and smaller than the total number of samples and less
     *                   than maxPos.
     * @param maxPos     maximum position to be picked. This value must be greater
     *                   or equal than zero and smaller than the total number of samples and
     *                   greater than minPos.
     * @param subsetSize subset size to be computed. This value must be smaller
     *                   than total number of samples.
     * @param pickLast   true indicates that last sample in range must always be
     *                   picked within subset. This is done to obtain faster execution times and
     *                   greater stability on some algorithms.
     * @param result     array containing indices to be picked. Provided array must
     *                   be at least of length subsetSize. The former subsetSize entries of the
     *                   array will be modified by this method.
     * @throws NotEnoughSamplesException   if subset size is greater than the
     *                                     total number of samples or if maxPos is greater than the total number of
     *                                     samples.
     * @throws InvalidSubsetSizeException  if subset size is zero or if result
     *                                     array does not have at least a length of subsetSize, or ig subset size
     *                                     is greater than the allowed range of positions to be picked.
     * @throws InvalidSubsetRangeException if maximum position is smaller than
     *                                     minimum position or maximum or minimum position are negative.
     */
    @Override
    public void computeRandomSubsetsInRange(
            final int minPos, final int maxPos, final int subsetSize, final boolean pickLast, final int[] result)
            throws NotEnoughSamplesException, InvalidSubsetSizeException, InvalidSubsetRangeException {
        if (subsetSize == 0) {
            throw new InvalidSubsetSizeException();
        }
        if (result.length < subsetSize) {
            throw new InvalidSubsetSizeException();
        }
        if (minPos >= maxPos || maxPos < 0 || minPos < 0) {
            throw new InvalidSubsetRangeException();
        }
        if ((maxPos - minPos) < subsetSize) {
            throw new InvalidSubsetSizeException();
        }
        if (numSamples < subsetSize) {
            throw new NotEnoughSamplesException();
        }
        if (maxPos > numSamples) {
            throw new NotEnoughSamplesException();
        }

        // On start set of selected indices is empty
        selectedIndices.clear();

        var counter = 0;
        int index;

        if (pickLast) {
            // this is done to accelerate computations and obtain more
            // stable results in some cases
            // pick last element in range
            index = maxPos - 1;
            selectedIndices.add(index);
            result[counter] = index;
            counter++;
        }

        if (!pickLast || counter < result.length) {
            // keep selecting only if not all elements have already been selected
            do {
                index = randomizer.nextInt(minPos, maxPos);

                // check whether this index has already been selected
                if (selectedIndices.contains(index)) {
                    continue;
                }

                // if not selected, pick it now
                selectedIndices.add(index);
                result[counter] = index;
                counter++;
            } while (counter < subsetSize);
        }

        selectedIndices.clear();
    }

    /**
     * Returns internal randomizer to generate uniformly distributed random
     * values.
     *
     * @return internal randomizer.
     */
    protected UniformRandomizer getRandomizer() {
        return randomizer;
    }

    /**
     * Initializes randomizer for an instance of this class.
     *
     * @param seedWithTime if true randomizer will be initialized using current
     *                     time as seed. If false, randomizer will always generate the same
     *                     pseudo-random sequence on each JVM execution.
     */
    private void createRandomizer(final boolean seedWithTime) {
        randomizer = new UniformRandomizer();
        if (seedWithTime) {
            randomizer.setSeed(System.currentTimeMillis());
        }
    }
}
