/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.LMedSRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 4, 2015
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.sorting.Sorter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * This class implements LMedS (Least Median of Squares) algorithm to robustly
 * estimate a data model.
 * LMedS is based on the idea that a given proportion of outliers exists in the
 * total amount of samples provided. This algorithm tries to iteratively find 
 * the beast subset of samples picking the onea with the least median of error.
 * To determine whether a sample is an outlier or not, and the estimated error 
 * for each sample, provided listener must be used.
 * Contrary to RANSAC, this algorithm does not require a fixed threshold to be
 * set to determine whether samples are inliers or not. Instead, threshold is
 * computed dynamically. Because of that LMedS typically produces results with
 * larger error than RANSAC having a similar computational cost, becaus esamples 
 * usually contain a large error. Hence, if threshold is known in advance for a 
 * given estimation, RANSAC should be preferred rather than LMedS.
 * On the contrary, if it can be ensured that samples are very accurate except
 * for some outliers, then LMedS becomes much more accurate than RANSAC because 
 * it typically converges to a solution with a very small threshold. However,
 * typically inlier samples tend to have certain error, and in practice LMedS
 * produces results with a similar accuracy and computational cost than RANSAC.
 * @param <T> type of object to be estimated.
 */
public class LMedSRobustEstimator<T> extends RobustEstimator<T> {

    /**
     * Constant defining default confidence of the estimated result, which is
     * 99%. This means that with a probability of 99% estimation will be 
     * accurate because chosen subsamples will be inliers.
     */
    public static final double DEFAULT_CONFIDENCE = 0.99;
    
    /**
     * Default maximum allowed number of iterations.
     */
    public static final int DEFAULT_MAX_ITERATIONS = 5000;
    
    /**
     * Minimum allowed confidence value.
     */
    public static final double MIN_CONFIDENCE = 0.0;
    
    /**
     * Maximum allowed confidence value.
     */
    public static final double MAX_CONFIDENCE = 1.0;    
    
    /**
     * Minimum allowed number of iterations.
     */
    public static final int MIN_ITERATIONS = 1;        
    
    /**
     * Default value to be used for stop threshold. Stop threshold can be used 
     * to keep the algorithm iterating in case that best threshold is not small 
     * enough. Once a better solution is found yielding a threshold smaller than 
     * this value, the algorithm will stop.
     */
    public static final double DEFAULT_STOP_THRESHOLD = 0.0;
    
    /**
     * Minimum allowed stop threshold value.
     */
    public static final double MIN_STOP_THRESHOLD = 0.0;
    
    /**
     * Default factor to normalize threshold to determine inliers. This factor 
     * can be used to increase or lower the dynamically computed threshold so 
     * that the algorithm becomes more or less accurate. The stricter the 
     * threshold (lower factor), the more time the algorithm will need to 
     * converge, if it can converge. By default the factor is 1.0, which makes 
     * the threshold to be computed as the median of residuals.
     */
    public static final double DEFAULT_INLIER_FACTOR = 1.0; //1.5 would also be reasonable
    
    /**
     * Minimum allowed value for inlier factor.
     */
    public static final double MIN_INLER_FACTOR = 0.0;
    
    /**
     * Constant to estimate standard deviation of residuals based on their 
     * median.
     */
    public static final double STD_CONSTANT = 1.4826;
    
    /**
     * Amount of confidence expressed as a value between 0 and 1.0 (which is
     * equivalent to 100%). The amount of confidence indicates the probability
     * that the estimated result is correct. Usually this value will be close
     * to 1.0, but not exactly 1.0.
     */
    private double mConfidence;
    
    /**
     * Maximum allowed number of iterations. When the maximum number of 
     * iterations is exceeded, result will not be available, however an 
     * approximate result will be available for retrieval.
     */
    private int mMaxIterations;
    
    /**
     * Instance in charge of picking random subsets of samples.
     */
    private SubsetSelector subsetSelector;   
    
    /**
     * Number of iterations to be done to obtain required confidence.
     */
    private int nIters;
    
    /**
     * Best solution that has been found so far during an estimation.
     */
    private T bestResult;    
    
    /**
     * Data related to inliers found for best result.
     */
    private LMedSInliersData mBestInliersData;
    
    /**
     * Threshold to be used to keep the algorithm iterating in case that
     * best threshold is not small enough. Once a better solution is found
     * yielding a threshold smaller than this value, the algorithm will stop.
     */
    private double mStopThreshold;
    
    /**
     * Factor to normalize threshold to determine inliers. This factor can be 
     * used to increase or lower the dynamically computed threshold so that the
     * algorithm becomes more or less accurate. The stricter the threshold 
     * (lower factor), the more time the algorithm will need to converge, if
     * it can converge. By default the factor is 1.0, which makes the threshold
     * to be computed as the median of residuals.
     */    
    private double mInlierFactor;
    
    
    /**
     * Constructor.
     */
    public LMedSRobustEstimator() {
        super();
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        nIters = mMaxIterations;
        bestResult = null;
        mBestInliersData = null;
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        mInlierFactor = DEFAULT_INLIER_FACTOR;
    }
    
    /**
     * Constructor with listener.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes, as well as in charge
     * of picking samples and doing per-iteration estimations.
     */
    public LMedSRobustEstimator(LMedSRobustEstimatorListener<T> listener) {
        super(listener);
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        nIters = mMaxIterations;
        bestResult = null;      
        mBestInliersData = null;
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        mInlierFactor = DEFAULT_INLIER_FACTOR;
    }
    
    /**
     * Returns amount of confidence expressed as a value between 0 and 1.0 
     * (which is equivalent to 100%). The amount of confidence indicates the
     * probability that the estimated result is correct. Usually this value will
     * be close to 1.0, but not exactly 1.0.
     * @return  amount of confidence as a value between 0.0 and 1.0.
     */
    public double getConfidence() {
        return mConfidence;
    }

    /**
     * Sets amount of confidence expressed as a value between 0 and 1.0 (which 
     * is equivalent to 100%). The amount of confidence indicates the 
     * probability that the estimated result is correct. Usually this value will 
     * be close to 1.0, but not exactly 1.0.
     * @param confidence confidence to be set as a value between 0.0 and 1.0.
     * @throws IllegalArgumentException if provided value is not between 0.0 and
     * 1.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setConfidence(double confidence) 
            throws IllegalArgumentException, LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (confidence < MIN_CONFIDENCE || confidence > MAX_CONFIDENCE) {
            throw new IllegalArgumentException();
        }
        mConfidence = confidence;
    }    
    
    /**
     * Maximum allowed number of iterations. When the maximum number of 
     * iterations is exceeded, result will not be available, however an 
     * approximate result will be available for retrieval.
     * @return maximum allowed number of iterations.
     */
    public int getMaxIterations() {
        return mMaxIterations;
    }
    
    /**
     * Sets maximum allowed number of iterations. When the maximum number of 
     * iterations is exceeded, result will not be available, however an 
     * approximate result will be available for retrieval.
     * @param maxIterations maximum allowed number of iterations to be set.
     * @throws IllegalArgumentException if provided value is less than 1.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setMaxIterations(int maxIterations) 
            throws IllegalArgumentException, LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxIterations < MIN_ITERATIONS) {
            throw new IllegalArgumentException();
        }
        mMaxIterations = maxIterations;
    }    
    
    /**
     * Returns threshold to be used to keep the algorithm iterating in case that
     * best threshold is not small enough. Once a better solution is found 
     * yielding a threshold smaller than this value, the algorithm will stop.
     * @return threshold to be used to keep the algorithm iterating in case that
     * best threshold is not small enough.
     */
    public double getStopThreshold() {
        return mStopThreshold;
    }
    
    /**
     * Sets threshold to be used to keep the algorithm iterating in case that
     * best threshold is not small enough. Once a better solution is found 
     * yielding a threshold smaller than this vlaue, the algorithm will stop.
     * @param stopThreshold threshold to be used to keep the algorithm iterating
     * in case that best threshold is not small enough.
     * @throws IllegalArgumentException if provided value is less or equal than 
     * 0.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setStopThreshold(double stopThreshold) 
            throws IllegalArgumentException, LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (stopThreshold < MIN_STOP_THRESHOLD) {
            throw new IllegalArgumentException();
        }
        
        mStopThreshold = stopThreshold;
    }        
    
    /**
     * Returns factor to normalize or adjust threshold to determine inliers. 
     * This factor can be used to increase or lower the dynamically computed 
     * threshold so that the algorithm becomes more or less accurate. The 
     * stricter the threshold (lower factor), the more time the algorithm will 
     * need to converge, if it can converge. By default the factor is 1.0, which 
     * makes the threshold to be computed as the median of residuals.
     * @return factor to normalize threshold to determine inliers.
     */    
    public double getInlierFactor() {
        return mInlierFactor;
    }
    
    /**
     * Sets factor to normalize or adjust threshold to determine inliers. 
     * This factor can be used to increase or lower the dynamically computed 
     * threshold so that the algorithm becomes more or less accurate. The 
     * stricter the threshold (lower factor), the more time the algorithm will 
     * need to converge, if it can converge. By default the factor is 1.0, which 
     * makes the threshold to be computed as the median of residuals.
     * @param inlierFactor inlier factor to be set.
     * @throws IllegalArgumentException if provided value is less or equal than 
     * 0.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setInlierFactor(double inlierFactor)
            throws IllegalArgumentException, LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (inlierFactor <= MIN_INLER_FACTOR) {
            throw new IllegalArgumentException();
        }
        
        mInlierFactor = inlierFactor;
    }
    
    /**
     * Returns number of iterations to be done to obtain required confidence.
     * @return number of iterations to be done to obtain required confidence.
     */
    public int getNIters() {
        return nIters;
    }
    
    /**
     * Returns best solution that has been found so far during an estimation.
     * @return best solution that has been found so far during an estimation.
     */
    public T getBestResult() {
        return bestResult;
    }
    
    /**
     * Returns data related to inliers found for best result.
     * @return data related to inliers found for best result.
     */
    public LMedSInliersData getBestInliersData() {
        return mBestInliersData;
    }
        
    /**
     * Indicates if estimator is ready to start the estimation process.
     * @return true if ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        if (!super.isReady()) {
            return false;
        }
        return (mListener instanceof LMedSRobustEstimatorListener);
    }        
    
    /**
     * Robustly estimates an instance of T.
     * @return estimated object.
     * @throws LockedException if robust estimator is locked.
     * @throws NotReadyException if provided input data is not enough to start
     * the estimation.
     * @throws RobustEstimatorException if estimation fails for any reason
     * (i.e. numerical instability, no solution available, etc).
     */       
    @Override
    public T estimate() throws LockedException, NotReadyException, 
            RobustEstimatorException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        try {
            LMedSRobustEstimatorListener<T> listener = 
                    (LMedSRobustEstimatorListener<T>)mListener;
            
            mLocked = true;
            
            listener.onEstimateStart(this);
        
            int totalSamples = listener.getTotalSamples();
            int subsetSize = listener.getSubsetSize();    
            int bestNumInliers;
            double threshold = Double.MAX_VALUE;
            nIters = Integer.MAX_VALUE;
            int newNIters;
            int currentIter = 0;                     
            //reusable list that will contain preliminar solutions on each 
            //iteration
            List<T> iterResults = new ArrayList<T>();            
            bestResult = null; //best result found so far
            //progress and previous progress to determine when progress
            //notification must occur
            float previousProgress = 0.0f, progress;
            //indices of subset picked in one iteration
            int[] subsetIndices = new int[subsetSize];
            double[] residualsTemp = new double[totalSamples];
            boolean improved;    //indicates if result improved
            //indicates whether algorithm must continue iterating
            boolean continueIteration = true;  
        
            if (subsetSelector == null) {
                //create new subset selector
                subsetSelector = SubsetSelector.create(totalSamples);
            } else {
                //set number of samples to current subset selector
                subsetSelector.setNumSamples(totalSamples);
            }
            
            //data related to inliers
            LMedSInliersData inliersData = new LMedSInliersData(totalSamples);
            //sorter to compute medians
            Sorter sorter = Sorter.create(); 
            
            while (continueIteration) {
                //generate a random subset of samples
                subsetSelector.computeRandomSubsets(subsetSize, subsetIndices);
                
                //clear list of preliminar solutions before calling listener
                iterResults.clear(); 
                //compute solution for current iteration
                listener.estimatePreliminarSolutions(subsetIndices, 
                        iterResults);    
                
                //iterate over all solutions that have been found
                improved = false;
                if (iterResults != null) {
                    for (T iterResult : iterResults) {
                        //compute inliers
                        computeInliers(iterResult, subsetSize, mInlierFactor, 
                                residualsTemp, listener, sorter, inliersData);

                        //save solution that produces the best residual
                        if (inliersData.isMedianResidualImproved()) {
                            improved = true;

                            //keep current solution
                            bestResult = iterResult;

                            //keep best inliers data corresponding to best solution,
                            //in case it can be useful along with the result
                            mBestInliersData = inliersData;

                            //recompute number of times the algorithm needs to be
                            //executed depending on current number of inliers to 
                            //achieve with probability mConfidence that we have 
                            //inliers and probability 1 - mConfidence that we have
                            //outliers
                            bestNumInliers = inliersData.getNumInliers();
                            double probInlier = ((double)bestNumInliers) /
                                    ((double)totalSamples);

                            double probSubsetAllInliers =
                                    Math.pow(probInlier, (double)subsetSize);

                            if (Math.abs(probSubsetAllInliers) < Double.MIN_VALUE ||
                                    Double.isNaN(probSubsetAllInliers)) {
                                newNIters = Integer.MAX_VALUE;
                            } else {
                                double logProbSomeOutliers = 
                                        Math.log(1.0 - probSubsetAllInliers);
                                if (Math.abs(logProbSomeOutliers) < 
                                        Double.MIN_VALUE ||
                                        Double.isNaN(logProbSomeOutliers)) {
                                    newNIters = Integer.MAX_VALUE;
                                } else {
                                    newNIters = (int)Math.ceil(Math.abs(
                                            Math.log(1.0 - mConfidence) / 
                                            logProbSomeOutliers));
                                }
                            }

                            if (newNIters < nIters) {
                                nIters = newNIters;
                            }

                            threshold = inliersData.getEstimatedThreshold();
                            
                            //create new inliers data instance until a new best 
                            //solution is found
                            double bestMedianResidual = 
                                    inliersData.getBestMedianResidual();
                            inliersData = new LMedSInliersData(totalSamples);
                            //update best median residual on new instance so 
                            //that only better solutions that are found later 
                            //can update inliers data
                            inliersData.update(bestMedianResidual, 
                                    inliersData.getStandardDeviation(),
                                    inliersData.getInliers(), 
                                    inliersData.getResiduals(), 
                                    inliersData.getNumInliers(), 
                                    bestMedianResidual, 
                                    inliersData.getEstimatedThreshold(), false);                            
                        }
                    }
                }
                
                if (nIters > 0) {
                    progress = Math.min((float)currentIter / (float)nIters, 
                            1.0f);
                } else {
                    progress = 1.0f;
                }
                if (progress - previousProgress > mProgressDelta) {
                    previousProgress = progress;
                    listener.onEstimateProgressChange(this, progress);                    
                }                
                currentIter++;
                continueIteration =  (currentIter < mMaxIterations) &&
                        (threshold > mStopThreshold);
                if (!improved) {
                    continueIteration &= (currentIter < nIters);
                } 
            
                listener.onEstimateNextIteration(this, currentIter);
            }
                        
            //no solution could be found after completing all iterations
            if (bestResult == null) {
                throw new RobustEstimatorException();
            }
                        
            listener.onEstimateEnd(this);
            
            return bestResult;
        } catch (SubsetSelectorException e) {
            throw new RobustEstimatorException(e);
        } finally {
            mLocked = false;
        }
    }
    
    /**
     * Returns data about inliers once estimation has been done.
     * @return data about inliers or null if estimation has not been done.
     */
    @Override
    public InliersData getInliersData() {
        return getBestInliersData();
    }    
    
    /**
     * Returns method being used for robust estimation.
     * @return method being used for robust estimation.
     */        
    @Override
    public RobustEstimatorMethod getMethod() {
        return RobustEstimatorMethod.LMedS;
    }    
    
    /**
     * Computes inliers data for current iteration.
     * @param <T> type of result to be estimated.
     * @param iterResult result to be tested on current iteration.
     * @param subsetSize subset sample size to be picked on each iteration.
     * @param inlierFactor factor to adjust threshold to determine whether 
     * samples are inliers or not.
     * @param residualsTemp temporal array to store residuals, since median
     * computation requires modifying the original array.
     * @param listener listener to obtain residuals for samples.
     * @param sorter sorter instance to compute median of residuals.
     * @param inliersData inliers data to be reused on each iteration.
     * @return inliers data.
     */
    private static <T> LMedSInliersData computeInliers(T iterResult, 
            int subsetSize, double inlierFactor, double[] residualsTemp, 
            LMedSRobustEstimatorListener<T> listener, 
            Sorter sorter, LMedSInliersData inliersData) {

        double[] residuals = inliersData.getResiduals();
        BitSet inliers = inliersData.getInliers();
        double bestMedianResidual = inliersData.getBestMedianResidual();
        boolean medianResidualImproved = false;
        
        int totalSamples = residuals.length;
        
        for (int i = 0; i < totalSamples; i++) {
            residuals[i] = Math.abs(listener.computeResidual(iterResult, i));
        }
        System.arraycopy(residuals, 0, residualsTemp, 0, residuals.length);
        double medianResidual = sorter.median(residualsTemp);
        if (medianResidual < bestMedianResidual) {
            bestMedianResidual = medianResidual;
            medianResidualImproved = true;
        }
        
        double standardDeviation = STD_CONSTANT * (1.0 + 5.0 / 
                (double)(totalSamples - subsetSize)) * Math.sqrt(
                medianResidual);
        double normEstimatedThreshold = inlierFactor * medianResidual;
        
        //determine which points are inliers
        int numInliers = 0;
        for (int i = 0; i < totalSamples; i++) {
            if (residuals[i] <= normEstimatedThreshold) {
                numInliers++;
                inliers.set(i);
            } else {
                inliers.clear(i);
            }
        }
        
        //store values in inliers data, only if residuals improve
        if (medianResidualImproved) {
            inliersData.update(bestMedianResidual, standardDeviation, inliers, 
                    residuals, numInliers, medianResidual, 
                    normEstimatedThreshold, medianResidualImproved);
        }
        
        return inliersData;
    }
    
    /**
     * Contains data related to inliers estimated in one iteration.
     */
    public static class LMedSInliersData extends InliersData {
        /**
         * Best median of error found so far taking into account all provided
         * samples.
         */
        private double mBestMedianResidual;
        
        /**
         * Standard deviation of error among all provided samples respect to
         * currently estimated result.
         */
        private double mStandardDeviation;
        
        /**
         * Efficiently stores which samples are considered inliers and which
         * ones aren't.
         */
        private BitSet mInliers;
                
        /**
         * Median of error found on current iteration among all provided 
         * samples.
         */
        private double mMedianResidual;
        
        /**
         * Estimated threshold to determine whether samples are inliers or not.
         */
        private double mEstimatedThreshold;
        
        /**
         * Indicates whether median residual computed in current iteration has
         * improved respect to previous iterations.
         */
        private boolean mMedianResidualImproved;
        
        /**
         * Constructor.
         * @param totalSamples total number of samples.
         */
        protected LMedSInliersData(int totalSamples) {
            mBestMedianResidual = mStandardDeviation = mMedianResidual = 
                    mEstimatedThreshold = Double.MAX_VALUE;
            mInliers = new BitSet(totalSamples);
            mResiduals = new double[totalSamples];
            mNumInliers = 0;     
            mMedianResidualImproved = false;
        }
        
        /**
         * Updates data contained in this instance.
         * @param bestMedianResidual best median of error found so far taking 
         * into account all provided samples.
         * @param standardDeviation standard deviation of error among all 
         * provided samples respect to currently estimated result.
         * @param inliers efficiently stores which samples are considered 
         * inliers and which ones aren't.
         * @param residuals residuals obtained for each sample of data.
         * @param numInliers number of inliers found on current iteration.
         * @param medianResidual median of error found on current iteration 
         * among all provided samples.
         * @param estimatedThreshold estimated threshold to determine whether 
         * samples are inliers or not.
         * @param medianResidualImproved indicates whether median residual 
         * computed in current iteration has improved respect to previous
         * iteration.
         */
        protected void update(double bestMedianResidual, double standardDeviation,
                BitSet inliers, double[] residuals, int numInliers, 
                double medianResidual, double estimatedThreshold, 
                boolean medianResidualImproved) {
            mBestMedianResidual = bestMedianResidual;
            mStandardDeviation = standardDeviation;
            mInliers = inliers;
            mResiduals = residuals;
            mNumInliers = numInliers;
            mMedianResidual = medianResidual;
            mEstimatedThreshold = estimatedThreshold;
            mMedianResidualImproved = medianResidualImproved;
        }
        
        /**
         * Returns best median of error found so far taking into account all 
         * provided samples.
         * @return best median of error found so far taking into account all
         * provided samples.
         */
        public double getBestMedianResidual() {
            return mBestMedianResidual;
        }
        
        /**
         * Returns standard deviation of error among all provided samples 
         * respect to currently estimated result.
         * @return standard deviation of error among all provided samples
         * respect to currently estimated result.
         */
        public double getStandardDeviation() {
            return mStandardDeviation;
        }

        /**
         * Returns efficient array indicating which samples are considered 
         * inliers and which ones aren't.
         * @return array indicating which samples are considered inliers and 
         * which ones aren't.
         */
        @Override
        public BitSet getInliers() {
            return mInliers;
        }

        /**
         * Returns median of error found on current iteration among all provided
         * samples.
         * @return median of error found on current iteration among all provided
         * samples.
         */
        public double getMedianResidual() {
            return mMedianResidual;
        }

        /**
         * Returns estimated threshold to determine whether samples are inliers 
         * or not.
         * @return estimated threshold to determine whether samples are inliers
         * or not.
         */
        public double getEstimatedThreshold() {
            return mEstimatedThreshold;
        }

        /**
         * Returns boolean indicating whether median residual computed in 
         * current iteration has improved respect to previous iterations.
         * @return true if median residual improved, false otherwise.
         */
        public boolean isMedianResidualImproved() {
            return mMedianResidualImproved;
        }
    }
}
