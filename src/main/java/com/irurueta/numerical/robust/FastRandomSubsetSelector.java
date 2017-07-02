/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.FastRandomSubsetSelector
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

import com.irurueta.statistics.UniformRandomizer;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * This class computes indices of subsets of samples using a uniform randomizer
 * to pick random samples as fast as possible.
 * This class ensures that samples are not repeated within a single subset.
 */
public class FastRandomSubsetSelector extends SubsetSelector{
    
    /**
     * Constant defining whether randomizer needs to be initialized with system
     * timer. By disabling this option it is ensured that the same result is
     * obtained on each execution
     */
    public static final boolean DEFAULT_SEED_RANDOMIZER_WITH_TIME = false;
    
    private UniformRandomizer mRandomizer;
    
    /**
     * Set containing selected indices on a given run.
     * This is kept around between executions to avoid excessive calls to 
     * garbage collector, as random subset generation is likely to be called
     * a large number of times during robust estimations. Because of that,
     * if a robust estimator is capable of spanning multiple executions at once
     * on different threads, then each thread needs to have a different subset
     * selector instance.
     */
    private Set<Integer> mSelectedIndices;
    
    /**
     * Constructor
     * @param numSamples number of samples to select subsets from
     * @throws IllegalArgumentException if provided number of samples is zero
     * or negative
     */
    public FastRandomSubsetSelector(int numSamples) 
            throws IllegalArgumentException{
        this(numSamples, DEFAULT_SEED_RANDOMIZER_WITH_TIME);
    }
    
    /**
     * Constructor
     * @param numSamples number of samples to select subsets from
     * @param seedRandomizerWithTime true if randomizer seed must be initialized
     * to system timer to obtain more random results
     * @throws IllegalArgumentException if provided number of samples is zero
     * or negative
     */
    public FastRandomSubsetSelector(int numSamples, 
            boolean seedRandomizerWithTime) throws IllegalArgumentException{
        super(numSamples);
        createRandomizer(seedRandomizerWithTime);
        mSelectedIndices = new HashSet<Integer>();
    }
    
    /**
     * Returns type of this subset selector
     * @return type of this subset selector
     */
    @Override
    public SubsetSelectorType getType(){
        return SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR;
    }

    /**
     * Returns internal randomizer to generate uniformly distributed random
     * values
     * @return internal randomizer
     */
    protected UniformRandomizer getRandomizer(){
        return mRandomizer;
    }
    
    /**
     * Computes a random subset of indices within range of number of samples to
     * be used on robust estimators
     * @param subsetSize subset size to be computed. This value must be smaller
     * than total number of samples
     * @param result array containing indices to be picked. Provided array must
     * be at least of length subsetSize. The former subsetSize entries of the
     * array will be modified by this method
     * @throws NotEnoughSamplesException if subset size is greater than the 
     * total number of samples
     * @throws InvalidSubsetSizeException if subset size is zero or if result
     * array does not have at least a length of subsetSize
     */
    @Override
    public void computeRandomSubsets(int subsetSize, int[] result) 
            throws NotEnoughSamplesException, InvalidSubsetSizeException {
        if(subsetSize == 0) throw new InvalidSubsetSizeException();
        if(result.length < subsetSize) throw new InvalidSubsetSizeException();
        if(mNumSamples < subsetSize) throw new NotEnoughSamplesException();
        
        //On start set of selected indices is empty
        mSelectedIndices.clear();        
        
        int counter = 0;
        int index;
        do{
            index = mRandomizer.nextInt(0, mNumSamples);
            
            //check whether this index has already been selected
            Integer intIndex = index;
            if(mSelectedIndices.contains(intIndex)) continue;
            
            //if not selected, pick it now
            mSelectedIndices.add(intIndex);
            result[counter] = index;
            counter++;
        }while(counter < subsetSize);
        
        mSelectedIndices.clear();
    }

    /**
     * Computes a random subset of indices within provided range of positions to
     * be used on robust estimators
     * @param minPos minimum position to be picked. This value must be greater 
     * or equal than zero and smaller than the total number of samples and less
     * than maxPos
     * @param maxPos maximum position to be picked. This value must be greater
     * or equal than zero and smaller than the total number of samples and 
     * greater than minPos
     * @param subsetSize subset size to be computed. This value must be smaller
     * than total number of samples
     * @param pickLast true indicates that last sample in range must always be
     * picked within subset. This is done to obtain faster execution times and
     * greater stability on some algorithms
     * @param result array containing indices to be picked. Provided array must
     * be at least of length subsetSize. The former subsetSize entries of the
     * array will be modified by this method
     * @throws NotEnoughSamplesException if subset size is greater than the 
     * total number of samples or if maxPos is greater than the total number of
     * samples
     * @throws InvalidSubsetSizeException if subset size is zero or if result
     * array does not have at least a length of subsetSize, or ig subset size
     * is greater than the allowed range of positions to be picked
     * @throws InvalidSubsetRangeException if maximum position is smaller than
     * minimum position or maximum or minimum position are negative
     */
    @Override
    public void computeRandomSubsetsInRange(int minPos, int maxPos, 
            int subsetSize, boolean pickLast, int[] result) 
            throws NotEnoughSamplesException, InvalidSubsetSizeException, 
            InvalidSubsetRangeException {
        if(subsetSize == 0) throw new InvalidSubsetSizeException();
        if(result.length < subsetSize) throw new InvalidSubsetSizeException();
        if(minPos >= maxPos || minPos < 0 || maxPos < 0) 
            throw new InvalidSubsetRangeException();        
        if((maxPos - minPos) < subsetSize) 
            throw new InvalidSubsetSizeException();        
        if(mNumSamples < subsetSize) throw new NotEnoughSamplesException();
        if(maxPos > mNumSamples) throw new NotEnoughSamplesException();
        
        //On start set of selected indices is empty
        mSelectedIndices.clear();        
        
        int counter = 0;
        int index;
        
        Integer intIndex;
        if(pickLast){ //this is done to accelerate computations and obtain more
                      //stable results in some cases
            //pick last element in range
            index = maxPos - 1;
            intIndex = index;
            mSelectedIndices.add(intIndex);
            result[counter] = index;
            counter++;
        }
        
        if(!pickLast || (pickLast && counter < result.length)){
            //keep selecting only if not all elements have already been selected
            do{
                index = mRandomizer.nextInt(minPos, maxPos);

                //check whether this index has already been selected
                intIndex = index;
                if(mSelectedIndices.contains(intIndex)) continue;

                //if not selected, pick it now
                mSelectedIndices.add(intIndex);
                result[counter] = index;
                counter++;
            }while(counter < subsetSize);
        }
        
        mSelectedIndices.clear();
    }
    
    /**
     * Initializes randomizer for an instance of this class.
     * @param seedWithTime if true randomizer will be initialized using current
     * time as seed. If false, randomizer will always generate the same 
     * pseudo-random sequence on each JVM execution
     */
    private void createRandomizer(boolean seedWithTime){
        mRandomizer = new UniformRandomizer(new Random());
        if(seedWithTime) mRandomizer.setSeed(System.currentTimeMillis());
    }
}
