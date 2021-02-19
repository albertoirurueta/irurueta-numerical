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

/**
 * Base class to pick subsets of samples.
 */
@SuppressWarnings("WeakerAccess")
public abstract class SubsetSelector {

    /**
     * Constant defining minimum amount of allowed samples.
     */
    public static final int MIN_NUM_SAMPLES = 1;

    /**
     * Defines default subset selector type.
     */
    public static final SubsetSelectorType DEFAULT_SUBSET_SELECTOR_TYPE =
            SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR;

    /**
     * Total number of samples to pick subsets from.
     * Subsets need to be always smaller or equal than total number of samples.
     */
    protected int mNumSamples;

    /**
     * Constructor.
     *
     * @param numSamples number of samples to select subsets from.
     * @throws IllegalArgumentException if provided number of samples is zero
     *                                  or negative.
     */
    public SubsetSelector(final int numSamples) {
        setNumSamples(numSamples);
    }

    /**
     * Returns number of samples to select subsets from.
     *
     * @return number of samples to select subsets from.
     */
    public int getNumSamples() {
        return mNumSamples;
    }

    /**
     * Sets number of samples to select subsets from.
     *
     * @param numSamples number of samples to select subsets from.
     * @throws IllegalArgumentException if provided number of samples is zero or
     *                                  negative.
     */
    public final void setNumSamples(final int numSamples) {
        if (numSamples < MIN_NUM_SAMPLES) {
            throw new IllegalArgumentException();
        }
        mNumSamples = numSamples;
    }

    /**
     * Computes a random subset of indices within range of number of samples to
     * be used on robust estimators.
     * If subsets need to be computed repeatedly in a small span of time then it
     * is suggested to use computeRandomSubsets(int, int[]) for better memory
     * usage.
     *
     * @param subsetSize subset size to be computed. This value must be smaller
     *                   than total number of samples.
     * @return array containing indices to be picked.
     * @throws NotEnoughSamplesException  if subset size is greater than the
     *                                    total number of samples.
     * @throws InvalidSubsetSizeException if subset size is zero or if result
     *                                    array does not have at least a length of subsetSize.
     * @see #computeRandomSubsets(int, int[])
     */
    public int[] computeRandomSubsets(final int subsetSize)
            throws NotEnoughSamplesException, InvalidSubsetSizeException {
        final int[] result = new int[subsetSize];
        computeRandomSubsets(subsetSize, result);
        return result;
    }

    /**
     * Computes a random subset of indices within provided range of positions to
     * be used on robust estimators.
     * If subsets need to be computed repeatedly in a small span of time then it
     * is suggested to use computeRandomSubsets(int, int, int, bool, int[]) for
     * better memory usage.
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
     * @return array containing indices to be picked.
     * @throws NotEnoughSamplesException   if subset size is greater than the
     *                                     total number of samples or if maxPos is greater than the total number of
     *                                     samples.
     * @throws InvalidSubsetSizeException  if subset size is zero, or if subset
     *                                     size is greater than the allowed range of positions to be picked.
     * @throws InvalidSubsetRangeException if maximum position is smaller than
     *                                     minimum position or maximum or minimum position are negative.
     */
    public int[] computeRandomSubsetsInRange(final int minPos, final int maxPos,
                                             final int subsetSize, final boolean pickLast)
            throws NotEnoughSamplesException, InvalidSubsetSizeException,
            InvalidSubsetRangeException {
        final int[] result = new int[subsetSize];
        computeRandomSubsetsInRange(minPos, maxPos, subsetSize, pickLast,
                result);
        return result;
    }

    /**
     * Returns type of this subset selector.
     *
     * @return type of this subset selector.
     */
    public abstract SubsetSelectorType getType();

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
    public abstract void computeRandomSubsets(final int subsetSize, final int[] result)
            throws NotEnoughSamplesException, InvalidSubsetSizeException;

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
     *                                     array does not have at least a length of subsetSize, or if subset size
     *                                     is greater than the allowed range of positions to be picked.
     * @throws InvalidSubsetRangeException if maximum position is smaller than
     *                                     minimum position or maximum or minimum position are negative.
     */
    public abstract void computeRandomSubsetsInRange(
            final int minPos, final int maxPos, final int subsetSize,
            final boolean pickLast, final int[] result)
            throws NotEnoughSamplesException, InvalidSubsetSizeException,
            InvalidSubsetRangeException;

    /**
     * Creates a new subset selector instance using provided total number of
     * samples and subset selector type.
     *
     * @param numSamples number of samples to select subsets from.
     * @param type       subset selector type.
     * @return a subset selector.
     * @throws IllegalArgumentException if provided number of samples is zero
     *                                  or negative.
     */
    public static SubsetSelector create(
            final int numSamples, final SubsetSelectorType type) {
        //noinspection SwitchStatementWithTooFewBranches
        switch (type) {
            case FAST_RANDOM_SUBSET_SELECTOR:
            default:
                return new FastRandomSubsetSelector(numSamples);
        }
    }

    /**
     * Creates a new subset selector instance using provided total number of
     * samples and default subset selector type.
     *
     * @param numSamples number of samples to select subsets from.
     * @return a subset selector.
     * @throws IllegalArgumentException if provided number of samples is zero or
     *                                  negative.
     */
    public static SubsetSelector create(final int numSamples) {
        return create(numSamples, DEFAULT_SUBSET_SELECTOR_TYPE);
    }
}
