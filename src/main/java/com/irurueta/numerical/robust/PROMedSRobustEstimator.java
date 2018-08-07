/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.PROMedSRobustEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 7, 2015
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * This class implements PROMedS (PROgressive least Median Sample) algorithm 
 * to robustly estimate a data model.
 * This algorithm is a mixture between LMedS and PROSAC, taking the best of
 * both.
 * Firstly, it has the advantage that no threshold is required to be set
 * beforehand, the same as LMedS. Threshold to determine inliers is computed
 * dynamically, which helps for an easier setup that is problem independent and
 * depending on the accuracy of the inliers, results will be more acurate than
 * RANSAC or PROSAC, just the same as LMedS.
 * On the other hand, if a certain information about the quality of the samples
 * is available, as in PROSAC, the algorithm takes advantage of this additional
 * information to prioritize the samples with higher quality in order to find
 * a solution much faster than RANSAC or LMedS.
 * Finally, if by any chance a threshold to determine inliers is also used, the
 * algorithm will try to get the solution that better fits in a pure median of
 * residuals model or in a threshold based one to determine inliers.
 * Hence, PROMedS can be as fast as PROSAC (which is typically about 100x faster
 * than RANSAC or LMedS), can obtain the same accuracy than LMedS (which can be
 * much better than RANSAC or PROSAC in certain scenarios), and has an easier
 * setup, which is problem independent because no threshold is required to be
 * known beforehand although one can be provided as well.
 * @param <T> type of object to be estimated.
 */
public class PROMedSRobustEstimator<T> extends RobustEstimator<T> {
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
     * Indicates whether the algorithm must stop prematurely when dynamically
     * computed threshold using median of residuals has a value lower than
     * provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less 
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     * By default stop threshold is enabled, so that computational cost is 
     * similar to RANSAC and only accuracy gets better if inliers are more
     * accurate.
     */
    public static final boolean DEFAULT_STOP_THRESHOLD_ENABLED = true;
    
    /**
     * Default maximum allowed outliers proportion in the input data. This is
     * used do determine the number of required iterations.
     */
    public static final double DEFAULT_MAX_OUTLIERS_PROPORTION = 0.8;
    
    /**
     * Minimum allowed value for maximum allowed outliers proportion in the
     * input data.
     */
    public static final double MIN_MAX_OUTLIERS_PROPORTION = 0.0;
    
    /**
     * Maximum allowed value for maximum allowed outliers proportion in the
     * input data.
     */
    public static final double MAX_MAX_OUTLIERS_PROPORTION = 1.0;
    
    /**
     * Defines the default value for the maximum probability that a solution
     * with more than inliersNStar in U_nStar exist and was not found after k
     * samples.
     */
    public static final double DEFAULT_ETA0 = 0.05;
    
    /**
     * Minimum allowed value for eta0.
     */
    public static final double MIN_ETA0 = 0.0;
    
    /**
     * Maximum allowed value for eta0.
     */
    public static final double MAX_ETA0 = 1.0;
    
    /**
     * Defines the default value for beta, which is the probability that a
     * match is declared inlier by mistake, i.e. the ratio of the "inlier"
     * surface by the total surface. The inlier surface is a disc with radius
     * 1.96s for homography/displacement computation, or a band with width 
     * 1.96*s*2 for epipolar geometry (s is the detection noise), and the total 
     * surface is the surface of the image.
     */
    public static final double DEFAULT_BETA = 0.01;
    
    /**
     * Minimum allowed value for beta.
     */
    public static final double MIN_BETA = 0.0;
    
    /**
     * Maximum allowed value for beta.
     */
    public static final double MAX_BETA = 1.0;
    
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
     * Minimum allowed value for inlier threshold to determine whether residuals
     * are inliers.
     */
    public static final double MIN_INLER_THRESHOLD = 0.0;
    
    /**
     * Indicates whether the inlier threshold will be used to find inliers along
     * with their median of residuals.
     */
    public static final boolean DEFAULT_USE_INLIER_THRESHOLD = true;
    
    /**
     * Constant to estimate standard deviation of residuals based on their 
     * median.
     */
    public static final double STD_CONSTANT = 1.4826;    
    
    /**
     * Chi squared.
     */
    public static final double CHI_SQUARED = 2.706;
        
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
     * In this implementation, PROSAC won't stop before having reached the 
     * corresponding inliers rate on the complete data set.
     * Maximum allowed outliers proportion in the input data: used to compute
     * nIters (can be as high as 0.95).
     */
    private double mMaxOutliersProportion;
    
    /**
     * eta0 is the maximum probability that a solution with more than 
     * inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     */
    private double mEta0;
    
    /**
     * beta is the probability that a match is declared inlier by mistake,
     * i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!.
     */
    private double mBeta;
    
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
    private PROMedSInliersData mBestInliersData;
    
    /**
     * Indicates whether the algorithm must stop prematurely when dynamically
     * computed threshold using median of residuals has a value lower than
     * provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less 
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     */    
    private boolean mStopThresholdEnabled;
        
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
     * Flag indicating whether thresholds to determine inliers are used, or if
     * only median of residuals is used. When true, the algorithm will try
     * to fit the best model, otherwise only median of residuals will be used.
     */
    private boolean mUseInlierThresholds;
    
    /**
     * Constructor.
     */
    public PROMedSRobustEstimator() {
        super();
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        mMaxOutliersProportion = DEFAULT_MAX_OUTLIERS_PROPORTION;
        mEta0 = DEFAULT_ETA0;
        mBeta = DEFAULT_BETA;
        nIters = mMaxIterations;
        bestResult = null;
        mBestInliersData = null;
        mStopThresholdEnabled = DEFAULT_STOP_THRESHOLD_ENABLED;
        mInlierFactor = DEFAULT_INLIER_FACTOR;  
        mUseInlierThresholds = DEFAULT_USE_INLIER_THRESHOLD;
    }
    
    /**
     * Constructor with listener.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes, as well as in charge
     * of picking samples and doing per-iteration estimations.
     */
    public PROMedSRobustEstimator(PROMedSRobustEstimatorListener<T> listener) {
        super(listener);
        mConfidence = DEFAULT_CONFIDENCE;
        mMaxIterations = DEFAULT_MAX_ITERATIONS;
        mMaxOutliersProportion = DEFAULT_MAX_OUTLIERS_PROPORTION;
        mEta0 = DEFAULT_ETA0;
        mBeta = DEFAULT_BETA;
        nIters = mMaxIterations;
        bestResult = null; 
        mBestInliersData = null;
        mStopThresholdEnabled = DEFAULT_STOP_THRESHOLD_ENABLED;
        mInlierFactor = DEFAULT_INLIER_FACTOR;  
        mUseInlierThresholds = DEFAULT_USE_INLIER_THRESHOLD;
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
        if(confidence < MIN_CONFIDENCE || confidence > MAX_CONFIDENCE) {
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
     * Returns boolean indicating whether the algorithm must stop prematurely 
     * when dynamically computed threshold using median of residuals has a value 
     * lower than provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less 
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     * @return true if stop threshold is enabled, false otherwise.
     */
    public boolean isStopThresholdEnabled() {
        return mStopThresholdEnabled;
    }
    
    /**
     * Sets boolean indicating whether the algorithm must stop prematurely when
     * dynamically computed threshold using median of residuals has a value 
     * lower than provided threshold in listener.
     * When this flag is enabled accuracy of PROMedS worsens to a lever similar
     * to PROSAC but the number of iterations is reduced (i.e. less 
     * computational cost). If more accuracy is desired at the expense of some
     * additional computation cost, then disable this flag.
     * @param stopThresholdEnabled true if stop threshold is enabled, false 
     * otherwise.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setStopThresholdEnabled(boolean stopThresholdEnabled) 
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        mStopThresholdEnabled = stopThresholdEnabled;
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
     * Returns flag indicating whether thresholds to determine inliers are used,
     * or if only median of residuals is used. When true, the algorithm will try
     * to fit the best model, otherwise only median of residuals will be used.
     * @return true if best model is used (threshold or median), otherwise only
     * median of residuals will be used.
     */
    public boolean isUseInlierThresholds() {
        return mUseInlierThresholds;
    }
    
    /**
     * Sets flag indicating whether thresholds to determine inliers are used, or
     * if only median of residuals is used. When true, the algorithm will try to
     * fit the best model, otherwise only median of residuals will be used.
     * @param useInlierThresholds true if best model is used (threshold or 
     * median), oitherwise only median of residuals will be used.
     * @throws LockedException if this estimator is locked because an estimation 
     * is being computed.
     */
    public void setUseInlierThresholds(boolean useInlierThresholds)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        
        mUseInlierThresholds = useInlierThresholds;
    }
        
    /**
     * Returns maximum allowed outliers proportion in the input data. This is
     * used to compute number of iterations to be done (nIters). It typically
     * can be as high as 0.95. Higher values, up to 1 are possible but not
     * recommended.
     * In this implementation, PROSAC won't stop before having reached the 
     * corresponding inliers rate on the complete data set.
     * @return maximum allowed outliers proportion in the input data.
     */
    public double getMaxOutliersProportion() {
        return mMaxOutliersProportion;
    }
    
    /**
     * Sets maximum allowed outliers proportion in the input data. This is used
     * to compute number of iterations to be done (nIters). It typically can be
     * as high as 0.95. Higher values, up to 1 are possible but not recommended.
     * In this implementation, PROSAC won't stop before having reached the
     * corresponding inliers rate on the complete data set.
     * @param maxOutliersProportion maximum allowed outliers proportion in the
     * input data.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     * greater than 1.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setMaxOutliersProportion(double maxOutliersProportion) 
            throws IllegalArgumentException, LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (maxOutliersProportion < MIN_MAX_OUTLIERS_PROPORTION ||
                maxOutliersProportion > MAX_MAX_OUTLIERS_PROPORTION) {
            throw new IllegalArgumentException();
        }
        
        mMaxOutliersProportion = maxOutliersProportion;
    }
    
    /**
     * Return eta0, which is the maximum probability that a solution with more 
     * than inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     * @return eta0 value.
     */
    public double getEta0() {
        return mEta0;
    }
    
    /**
     * Sets eta0, which is the maximum probability that a solution with more 
     * than inliersNStar inliers in U_nStar exists and was not found after k
     * samples (typically set to 5%).
     * @param eta0 eta0 value to be set.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     * greater than 1.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setEta0(double eta0) throws IllegalArgumentException, 
            LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (eta0 < MIN_ETA0 || eta0 > MAX_ETA0) {
            throw new IllegalArgumentException();
        }
        
        mEta0 = eta0;
    }
    
    /**
     * Returns beta, which is the probability that a match is declared inlier by 
     * mistake, i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!
     * @return beta value.
     */
    public double getBeta() {
        return mBeta;
    }
    
    /**
     * Sets beta, which is the probability that a match is declared inlier by 
     * mistake, i.e. the ratio of the "inlier" surface by the total surface. The
     * inlier surface is a disc with radius 1.96s for homography/displacement
     * computation, or a band with width 1.96s*2 for epipolar geometry (s is
     * the detection noise), and the total surface is the surface of the image
     * YOU MUST ADJUST THIS VALUE, DEPENDING ON YOUR PROBLEM!
     * @param beta beta value to be set.
     * @throws IllegalArgumentException if provided value is less than 0.0 or
     * greater than 1.0.
     * @throws LockedException if this estimator is locked because an estimation
     * is being computed.
     */
    public void setBeta(double beta) throws IllegalArgumentException, 
            LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (beta < MIN_BETA || beta > MAX_BETA) {
            throw new IllegalArgumentException();
        }
        
        mBeta = beta;
    }
    
    /**
     * Returns number of iterations to be done to obtain required confidence.
     * This does not need to be equal to the actual number of iterations the
     * algorithm finally required to obtain a solution.
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
    protected PROMedSInliersData getBestInliersData(){
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
        return (mListener instanceof PROMedSRobustEstimatorListener);
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
            PROMedSRobustEstimatorListener<T> listener = 
                    (PROMedSRobustEstimatorListener<T>)mListener;
            
            mLocked = true;
            
            listener.onEstimateStart(this);
        
            int totalSamples = listener.getTotalSamples(); //N = CORRESPONDENCES
            int subsetSize = listener.getSubsetSize();    
            
            double[] qualityScores = listener.getQualityScores();
            //check for invalid quality scores length
            if(qualityScores.length != totalSamples) throw new RobustEstimatorException();
            
            double inlierThreshold = 0.0;
            if (mUseInlierThresholds) {
                inlierThreshold = listener.getThreshold();
            }
            //obtain indices referring to original samples position after sorting
            //quality scores in descending order
            int[] sortedIndices = computeSortedQualityIndices(
                    listener.getQualityScores());
            
            //reusable list that will contain preliminar solutions on each 
            //iteration
            List<T> iterResults = new ArrayList<T>();            
            bestResult = null;
            float previousProgress = 0.0f, progress;
            //subset indices obtained from a subset selector
            int[] subsetIndices = new int[subsetSize];
            double[] residualsTemp = new double[totalSamples];
            //subset indices referred to the real samples positions after taking
            //into account the sorted indices obtained from quality scores
            int[] transformedSubsetIndices = new int[subsetSize];
            //array containing inliers efficiently
            BitSet inliers = new BitSet(totalSamples); 
            
            //T_N
            nIters = Math.min(computeIterations(
                    1.0 - mMaxOutliersProportion, subsetSize, mConfidence), 
                    mMaxIterations);
            
            int sampleSizeStar = totalSamples; //termination length
            int inliersNStar = 0;   //number of inliers found within the first
                                    //nStar data points
            int inliersBest = 0;    //best number of inliers found so far
                                    //(store the model that goes with it)
            double threshold = Double.MAX_VALUE; //threshold to stop algorithm
            int inliersMin = (int)((1.0 - mMaxOutliersProportion) * 
                    (double)totalSamples);
            int currentIter = 0;    //iteration number (t)
            int sampleSize = subsetSize; //(n) we draw samples from the set U_n
                                         //of the top n (sampleSize) data points
            double Tn = (double)nIters; //average number of samples {M_i}_{i=1}^{Tn}
                                        //that contains samples from U_n only
            int TnPrime = 1; //integer version of Tn
            int kNStar = nIters; //number of samples to draw to reach the 
                                //maximality constraint
            
            //initialize Tn
            for (int i = 0; i < subsetSize; i++) {
                Tn *= (double)(sampleSize - i) / (double)(totalSamples - i);
            }
        
            if (subsetSelector == null) {
                //create new subset selector
                subsetSelector = SubsetSelector.create(totalSamples);
            } else {
                //set number of samples to current subset selector
                subsetSelector.setNumSamples(totalSamples);
            }
            
            //data related to inliers
            PROMedSInliersData inliersData = new PROMedSInliersData(
                    totalSamples);
            //sorter to compute medians
            Sorter sorter = Sorter.create(); 
            
            boolean improved;    //indicates if result improved            
            boolean continueIteration = true;

            //iterate until the expected number of inliers or the estimated
            //number of iterations is reached                  
            while (continueIteration) {
                if (kNStar > 0) {
                    progress = Math.min((float)currentIter / (float)kNStar, 
                            1.0f);
                } else {
                    progress = 1.0f;
                }
                if (progress - previousProgress > mProgressDelta) {
                    previousProgress = progress;
                    listener.onEstimateProgressChange(this, progress);                    
                }     
                currentIter++;                
                
                //choice of the hypothesis generation set
                
                //The growth function is defined as
                //g(t) = min{n : TnPrime > t} where n is sampleSize
                //Thus sampleSize should be incremented if currentIter > TnPrime
                if ((currentIter > TnPrime) && (sampleSize < sampleSizeStar)) {
                    double TnPlus1 = (Tn * (double)(sampleSize + 1)) /
                            (double)(sampleSize + 1 - subsetSize);
                    sampleSize++;
                    TnPrime += (int)Math.ceil(TnPlus1 - Tn);
                    Tn = TnPlus1;
                }
                
                //Draw semi-random sample
                if (currentIter > TnPrime) {
                    //during the finishing stage (sampleSize == sampleSizeStar &&
                    //currentIter > TnPrime), draw a standard RANSAC sample
                    //The sample contains subsetSize points seleted from U_n at
                    //random
                    subsetSelector.computeRandomSubsets(subsetSize, 
                            subsetIndices);
                } else {
                    //The sample contains subsetSize-1 points selected from
                    //U_sampleSize_1 at random and u_sampleSize
                    
                    subsetSelector.computeRandomSubsetsInRange(0, sampleSize, 
                            subsetSize, true, subsetIndices);
                }
                
                transformIndices(subsetIndices, sortedIndices, 
                        transformedSubsetIndices);
                                
                //clear list of preliminar solutions before calling listener
                iterResults.clear();                 
                //compute solution for current iteration
                listener.estimatePreliminarSolutions(transformedSubsetIndices, 
                        iterResults);
                
                int inliersCurrent; //total number of inliers for a
                                        //given result                
                //iterate over all solutions that have been found
                improved = false;
                if (iterResults != null) {
                    for (T iterResult : iterResults) {
                        //compute inliers
                        computeInliers(iterResult, subsetSize, mInlierFactor, 
                                mUseInlierThresholds, inlierThreshold, 
                                residualsTemp, listener, sorter, inliersData);
                        inliersCurrent = inliersData.getNumInliers();                                        

                        if (inliersData.isMedianResidualImproved()) {
                            improved = true;

                            //keep current solution
                            bestResult = iterResult;

                            //update estimated thresholds to be used as stop 
                            //criterion
                            threshold = inliersData.getEstimatedThreshold();                        
                        }

                        if (inliersCurrent > inliersBest) {
                            //update best number of inliers
                            inliersBest = inliersCurrent;
                            //keep current solution
                            bestResult = iterResult;                        

                            keepInliersData(inliersData, totalSamples);

                            //select new termination length sampleSizeStar if possible
                            //only when a new sample is better than the others found
                            //so far
                            int sampleSizeBest = totalSamples; //best value found so
                                                    //far in terms of inliers ration
                            int inliersSampleSizeBest = inliersCurrent;

                            int sampleSizeTest; //test value for the termination 
                                                //length
                            int inliersSampleSizeTest;  //number of inliers for that
                                                        //test value
                            double epsilonSampleSizeBest = 
                                    (double)inliersSampleSizeBest /
                                    (double)sampleSizeBest;

                            for (sampleSizeTest = totalSamples, 
                                    inliersSampleSizeTest = inliersCurrent; 
                                    sampleSizeTest > subsetSize; sampleSizeTest--) {
                                //Loop invariants:
                                //- inliersSampleSizeTest is the number of inliers 
                                //  for the sampleSizeTest first correspondences
                                //- sampleSizeBest is the value between 
                                //  sampleSizeTest+1 and totalSamples that maximizes
                                //  the ratio inliersSampleSizeBest/sampleSizeBest

                                //- Non-randomness: In >= Imin(n*)
                                //- Maximality: the number of samples that were drawn 
                                //  so far must be enough so that the probability of
                                //  having missed a set of inliers is below eta=0.01.
                                //  This is the classical RANSAC termination criterion, 
                                //  except that it takes into account only the 
                                //  sampleSize first samples (not the total number 
                                //  of samples)
                                //  kNStar = log(eta0) / log(1 - (inliersNStar/
                                //  sampleSizeStar)^subsetSize
                                //  We have to minimize kNStar, e.g. maximize 
                                //  inliersNStar/sampleSizeStar, a straightforward
                                //  implementation would use the following test:
                                //  if(inliersSampleSizeTest > epsilonSampleSizeBest *
                                //  sampleSizeTest){ ... blah blah blah
                                //  However, since In is binomial, and in the case of
                                //  evenly distributed inliers, a better test would be
                                //  to reduce sampleSizeStar only if there's a 
                                //  significant improvement in epsilon. Thus we use a
                                //  Chi-squared test (P=0.10), together with the normal
                                //  approximation to the binomial (mu = 
                                //  epsilonSampleSizeStart * sampleSizeTest, sigma =
                                //  sqrt(sampleSizeTest * epsilonSampleSizeStar * (1 -
                                //  epsilonSampleSizeStar))).
                                //  There is a significant difference between the two
                                //  tests (e.g. with the computeInliers function 
                                //  provided)
                                //  We do the cheap test first, and the expensive test
                                //  only if the cheap one passes
                                if ((inliersSampleSizeTest * sampleSizeBest > 
                                        inliersSampleSizeBest * sampleSizeTest) &&
                                        (inliersSampleSizeTest > 
                                        epsilonSampleSizeBest * sampleSizeTest +
                                        Math.sqrt(sampleSizeTest * epsilonSampleSizeBest * 
                                        (1.0 - epsilonSampleSizeBest) * CHI_SQUARED))) {

                                    if (inliersSampleSizeTest < 
                                            Imin(subsetSize, sampleSizeTest, mBeta)) {
                                        //equation not satisfied, no need to test for
                                        //smaller sampleSizeTest values anyway
                                        break; //jump out of the for(sampleSizeTest) loop
                                    }
                                    sampleSizeBest = sampleSizeTest;
                                    inliersSampleSizeBest = inliersSampleSizeTest;
                                    epsilonSampleSizeBest = 
                                            (double)inliersSampleSizeBest /
                                            (double)sampleSizeBest;
                                }

                                //prepare for next loop iteration
                                inliersSampleSizeTest -=
                                        inliers.get(sortedIndices[sampleSizeTest - 1]) ? 1 : 0;
                            } //for (sampleSizeTest...

                            //is the best one we found even better than sampleSizeStar?
                            if (inliersSampleSizeBest * sampleSizeStar >
                                    inliersNStar * sampleSizeBest) {

                                //update all values
                                sampleSizeStar = sampleSizeBest;
                                inliersNStar = inliersSampleSizeBest;
                                kNStar = computeIterations((double)inliersNStar / 
                                        (double)sampleSizeStar, subsetSize, 
                                        1.0 - mEta0);
                            }
                        } //if(inliersCurrent > inliersBest)                    
                    }
                }
                
                continueIteration = (currentIter < mMaxIterations); 
                if (mUseInlierThresholds && mStopThresholdEnabled) {
                    //if inlier threshold is being used, and stop threshold is
                    //enabled, then stop the algorithm if threshold determined
                    //by median of residuals has a value lower than inlier 
                    //threshold
                    continueIteration &= (threshold > inlierThreshold);
                }
        
                if (!improved) {
                    continueIteration &= ((inliersBest < inliersMin) || 
                            (currentIter < kNStar)) && (currentIter < nIters);    
                }
                                                
                listener.onEstimateNextIteration(this, currentIter);
            } //while(currentIter < kNStar ...
            
            //no solution could be found after completing all iterations
            if (bestResult == null) {
                throw new RobustEstimatorException();
            }
                        
            listener.onEstimateEnd(this);
            
            return bestResult;
        } catch (SubsetSelectorException | SortingException e) {
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
        return RobustEstimatorMethod.PROMedS;
    }    
    
    /**
     * Transforms indices picked by the subset selector into the indices where
     * samples are actually localed by taking into account their original 
     * position before sorting quality scores.
     * @param subsetIndices indices picked by the subset selector. These are
     * positions after sorting. Must have the subset length.
     * @param sortedIndices indices relating sorted positions to their original
     * positions. Each position i-th in the array refers to the original 
     * position before sorting. Must have the number of samples length.
     * @param transformedSubsetIndices array where result is stored. Must have
     * the subset length.
     */
    private static void transformIndices(int[] subsetIndices, 
            int[] sortedIndices, int[] transformedSubsetIndices) {
        int length = transformedSubsetIndices.length;
        for (int i = 0; i < length; i++) {
            transformedSubsetIndices[i] = sortedIndices[subsetIndices[i]];
        }
    }
    
    /**
     * Computes inliers data for current iteration.
     * @param <T> type of result to be estimated.
     * @param iterResult result to be tested on current iteration.
     * @param subsetSize subset sample size to be picked on each iteration.
     * @param inlierFactor factor to adjust threshold to determine whether 
     * samples are inliers or not.
     * @param useInlierThresholds true to use thresholds to determine inliers,
     *                            false otherwise.
     * @param inlierThreshold threshold to determine which samples are inliers.
     * @param residualsTemp temporal array to store residuals, since median
     * computation requires modifying the original array.
     * @param listener listener to obtain residuals for samples.
     * @param sorter sorter instance to compute median of residuals.
     * @param inliersData inliers data to be reused on each iteration
     * @return inliers data.
     */
    private static <T> PROMedSInliersData computeInliers(T iterResult, 
            int subsetSize, double inlierFactor, boolean useInlierThresholds, 
            double inlierThreshold, double[] residualsTemp,
            LMedSRobustEstimatorListener<T> listener, 
            Sorter sorter, PROMedSInliersData inliersData) {

        double[] residuals = inliersData.getResiduals();
        BitSet lmedsInliers = inliersData.getInliersLMedS();
        BitSet msacInliers = inliersData.getInliersMSAC();
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
        boolean lmedsInlierModelEnabled = true; //by default if thresholds are not used
        if (useInlierThresholds) {
            int numInliersMsac = 0;
            int numInliersLmedS = 0;
            for (int i = 0; i < totalSamples; i++) {
                if (residuals[i] <= normEstimatedThreshold) {
                    numInliersLmedS++;
                    lmedsInliers.set(i);
                } else {
                    lmedsInliers.clear(i);
                }
                if (residuals[i] <= inlierThreshold) {
                    numInliersMsac++;
                    msacInliers.set(i);
                } else {
                    msacInliers.clear(i);
                }
            }
            
            //keep model with smaller number of inliers (to be more restrictive)
            lmedsInlierModelEnabled = numInliersLmedS <numInliersMsac;
            numInliers = lmedsInlierModelEnabled ?  
                    numInliersLmedS : numInliersMsac;
        } else {               

            for (int i = 0; i < totalSamples; i++) {
                if (residuals[i] <= normEstimatedThreshold) {
                    numInliers++;
                    lmedsInliers.set(i);
                } else {
                    lmedsInliers.clear(i);
                }
            }
        }
        
        //store values in inliers data, only if residuals improve
        if (medianResidualImproved) {
            inliersData.update(bestMedianResidual, standardDeviation, 
                    lmedsInliers, msacInliers, lmedsInlierModelEnabled,
                    residuals, numInliers, medianResidual, 
                    normEstimatedThreshold, medianResidualImproved);         
        } else {
            inliersData.mMedianResidualImproved = false;
        }
        
        return inliersData;
    }
    
    /**
     * Obtains indices of samples corresponding to samples ordered in descending
     * quality scores.
     * @param qualityScores quality scores associated to each sample to be used
     * to obtain indices to sort samples in descending order of quality values.
     * @return indices to sort samples in descending order of quality values.
     * @throws SortingException if sorting fails.
     */
    private static int[] computeSortedQualityIndices(double[] qualityScores) 
            throws SortingException {
        Sorter sorter = Sorter.create();
        double[] qualityScoresCopy = Arrays.copyOf(qualityScores, 
                qualityScores.length);
        //this method modifies quality scores copy array because it gets sorted
        //in ascending order. Indices contains indices of samples corresponding
        //to quality scores ordered in ascending order
        int[] indices = sorter.sortWithIndices(qualityScoresCopy);
        
        //reverse indices so we have indices of samples ordered in descening
        //order of quality
        reverse(indices);
        
        return indices;
    }
    
    /**
     * Reverses provided array.
     * @param array array to be reversed.
     */
    private static void reverse(int[] array) {
        int length = array.length;
        for(int i = 0; i < length / 2; i++){
            int temp = array[i];
            int pos = length - 1 - i;
            array[i] = array[pos];
            array[pos] = temp;
        }
    }
    
    /**
     * Computes number of required iterations to achieve required confidence
     * with current probability of inlier and sample subset size.
     * @param probInlier probability of inlier.
     * @param subsetSize sample subset size.
     * @param confidence required confidence of result.
     * @return number of required iterations.
     */
    private static int computeIterations(double probInlier, int subsetSize,
            double confidence) {
        
        //compute number of times the algorithm needs to be executed depending
        //on number of inliers respect total points to achieve with probability
        //confidence that we have all inliers and probability 1 - confidence
        //that we have some outliers
        double probSubsetAllInliers = Math.pow(probInlier, (double)subsetSize);
        if (Math.abs(probSubsetAllInliers) < Double.MIN_VALUE ||
                Double.isNaN(probSubsetAllInliers)) {
            return Integer.MAX_VALUE;
        } else {
            double logProbSomeOutliers = Math.log(1.0 - probSubsetAllInliers);
            if (Math.abs(logProbSomeOutliers) < Double.MIN_VALUE || 
                    Double.isNaN(logProbSomeOutliers)) {
                return Integer.MAX_VALUE;
            } else {
                return (int)Math.ceil(Math.abs(Math.log(1.0 - confidence) / 
                        logProbSomeOutliers));
            }
        }
    }
        
    /**
     * Non randomness states that i-m (where i is the cardinal of the set of
     * inliers for a wrong model) follows the binomial distribution B(n,beta).
     * For n big enough, B(n,beta) aproximates to normal distribution N(mu, 
     * sigma^2) by the central limit theorem, with mu = n*beta and sigma =
     * sqrt(n*beta*(1 - beta)).
     * Psi, the probability that In_star out of n_star data points are by chance
     * inliers to an arbitrary incorrect model, is set to 0.05 (5%, as in the
     * original paper), and you must change the Chi2 value if you chose a
     * different value for psi.
     * @param subsetSize sample subset size.
     * @param sampleSize total number of samples.
     * @param beta beta value.
     * @return i-m.
     */
    private static int Imin(int subsetSize, int sampleSize, double beta) {
        double mu = sampleSize * beta;
        double sigma = Math.sqrt(sampleSize * beta * (1.0 - beta));
        
        return (int)Math.ceil(subsetSize + mu + sigma * Math.sqrt(CHI_SQUARED));
    } 
    
    /**
     * Keeps inliers data stored and initializes a new one with proper
     * initial values.
     * @param inliersData inliers data to be stored.
     * @param totalSamples total number of samples.
     */
    private void keepInliersData(PROMedSInliersData inliersData, 
            int totalSamples) {
        //keep best inliers data corresponding to best solution,
        //in case it can be useful along with the result
        mBestInliersData = inliersData;

        //create new inliers data instance until a new best solution
        //is found
        double bestMedianResidual = inliersData.getBestMedianResidual();
        inliersData = new PROMedSInliersData(totalSamples);
        //update best median residual on new instance so that
        //only better solutions that are found later can update 
        //inliers data
        inliersData.update(bestMedianResidual, 
                inliersData.getStandardDeviation(), 
                inliersData.getInliersLMedS(), 
                inliersData.getInliersMSAC(), 
                inliersData.isLMedSInlierModelEnabled(), 
                inliersData.getResiduals(), 
                inliersData.getNumInliers(), 
                bestMedianResidual, 
                inliersData.getEstimatedThreshold(), false);
    }
    
    /**
     * Contains data related to inliers estimated in one iteration.
     */
    public static class PROMedSInliersData extends InliersData {
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
         * Inliers considering LMedS model.
         */
        private BitSet mInliersLmeds;
        
        /**
         * Inliers considering MSAC model.
         */
        private BitSet mInliersMsac;
        
        /**
         * Indicates whether LMedS or MSAC inlier model is enabled.
         */
        private boolean mLmedsInlierModelEnabled;
                
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
        protected PROMedSInliersData(int totalSamples) {
            mBestMedianResidual = mStandardDeviation = mMedianResidual = 
                    mEstimatedThreshold = Double.MAX_VALUE;
            mInliersLmeds = new BitSet(totalSamples);
            mInliersMsac = new BitSet(totalSamples);
            mLmedsInlierModelEnabled = true;
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
         * @param inliersLmeds stores which samples are considered inliers when
         * LMedS inlier model is used.
         * @param inliersMsac stores which samples are considered inliers when
         * MSAC inlier model is used.
         * @param lmedsInlierModelEnabled indicates whether the LMedS or MSAC
         * inlier model is used.
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
                BitSet inliersLmeds, BitSet inliersMsac, 
                boolean lmedsInlierModelEnabled, double[] residuals, 
                int numInliers, double medianResidual, 
                double estimatedThreshold, boolean medianResidualImproved) {
            mBestMedianResidual = bestMedianResidual;
            mStandardDeviation = standardDeviation;
            mLmedsInlierModelEnabled = lmedsInlierModelEnabled;
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
            return mLmedsInlierModelEnabled ? mInliersLmeds : mInliersMsac;
        }

        /**
         * Returns inliers considering LMedS model.
         * @return inliers considering LMedS model.
         */
        protected BitSet getInliersLMedS() {
            return mInliersLmeds;
        }

        /**
         * Returns inliers considering MSAC model.
         * @return inliers considering MSAC model.
         */
        protected BitSet getInliersMSAC() {
            return mInliersMsac;
        }

        /**
         * Returns boolean indicating whether LMedS or MSAC inlier model is
         * enabled. If true, estimated threshold was used to determine inliers, 
         * if false only median of residuals was used.
         * @return true if LMedS model is used, false if MSAC model is used.
         */
        public boolean isLMedSInlierModelEnabled() {
            return mLmedsInlierModelEnabled;
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
