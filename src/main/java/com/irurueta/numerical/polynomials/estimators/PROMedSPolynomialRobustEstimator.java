/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.numerical.robust.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Finds the best polynomial using PROMedS algorithm.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class PROMedSPolynomialRobustEstimator extends PolynomialRobustEstimator {    
    
    /**
     * Default value to be used for stop threshold. Stop threshold can be used 
     * to keep the algorithm iterating in case that best estimated threshold 
     * using median of residuals is not small enough. Once a solution is found 
     * that generates a threshold below this value, the algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close 
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would 
     * iterate for a long time trying to find the best solution when indeed 
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results
     */
    public static final double DEFAULT_STOP_THRESHOLD = 1e-6;
    
    /**
     * Minimum allowed stop threshold value
     */
    public static final double MIN_STOP_THRESHOLD = 0.0;
    
    /**
     * Threshold to be used to keep the algorithm iterating in case that best 
     * estimated threshold using median of residuals is not small enough. Once 
     * a solution is found that generates a threshold below this value, the 
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close 
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would 
     * iterate for a long time trying to find the best solution when indeed 
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results
     */
    private double mStopThreshold;    
    
    /**
     * Quality scores corresponding to each provided polynomial evaluation.
     * The larger the score value the better the quality of the sample.
     */
    private double[] mQualityScores;
    
    /**
     * Constructor.
     */
    public PROMedSPolynomialRobustEstimator() {
        super();
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public PROMedSPolynomialRobustEstimator(int degree) {
        super(degree);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided number of evaluations is
     * less than the required minimum.
     */
    public PROMedSPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations) {
        super(evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     */
    public PROMedSPolynomialRobustEstimator(
            PolynomialRobustEstimatorListener listener) {
        super(listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided number of evaluations is less than the required minimum for
     * provided degree.
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations) {
        super(degree, evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            PolynomialRobustEstimatorListener listener) {
        super(degree, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided number of evaluations is 
     * less than the required minimum.
     */
    public PROMedSPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations, 
            PolynomialRobustEstimatorListener listener) {
        super(evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided number of evaluations is less than the required minimum for
     * provided degree.
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations,
            PolynomialRobustEstimatorListener listener) {
        super(degree, evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
    }

    /**
     * Constructor.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided quality scores length is
     * smaller than minimum required size for default degree (2 evaluations).
     */
    public PROMedSPolynomialRobustEstimator(double[] qualityScores) {
        super();
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided quality scores length is smaller than minimum required size
     * for provided degree (degree + 1).
     */
    public PROMedSPolynomialRobustEstimator(int degree, double[] qualityScores) {
        super(degree);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided quality scores length or 
     * number of evaluations is smaller than minimum required size for default 
     * degree (2 evaluations).
     */
    public PROMedSPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations, double[] qualityScores) {
        super(evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided quality scores length is
     * smaller than minimum required size for default degree (2 evaluations).
     */
    public PROMedSPolynomialRobustEstimator(
            PolynomialRobustEstimatorListener listener, 
            double[] qualityScores) {
        super(listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided quality scores or evaluations length is smaller than minimum 
     * required size for provided degree (degree + 1).
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations, double[] qualityScores) {
        super(degree, evaluations);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided quality scores length is smaller than minimum required size
     * for provided degree (degree + 1).
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            PolynomialRobustEstimatorListener listener,
            double[] qualityScores) {
        super(degree, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided quality scores or 
     * evaluations length is smaller than minimum required size for default 
     * degree (2 evaluations).
     */
    public PROMedSPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations, 
            PolynomialRobustEstimatorListener listener,
            double[] qualityScores) {
        super(evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @param qualityScores quality scores corresponding to each provided 
     * evaluation.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided quality scores or evaluations length is smaller than minimum 
     * required size for provided degree (degree + 1).
     */
    public PROMedSPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations,
            PolynomialRobustEstimatorListener listener,
            double[] qualityScores) {
        super(degree, evaluations, listener);
        mStopThreshold = DEFAULT_STOP_THRESHOLD;
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Returns threshold to be used to keep the algorithm iterating in case that 
     * best estimated threshold using median of residuals is not small enough. 
     * Once a solution is found that generates a threshold below this value, the 
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close 
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would 
     * iterate for a long time trying to find the best solution when indeed 
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results
     * @return stop threshold to stop the algorithm prematurely when a certain
     * accuracy has been reached
     */
    public double getStopThreshold() {
        return mStopThreshold;
    }
    
    /**
     * Sets threshold to be used to keep the algorithm iterating in case that
     * best estimated threshold using median of residuals is not small enough. 
     * Once a solution is found that generates a threshold below this value, the 
     * algorithm will stop.
     * The stop threshold can be used to prevent the LMedS algorithm iterating
     * too many times in cases where samples have a very similar accuracy.
     * For instance, in cases where proportion of outliers is very small (close 
     * to 0%), and samples are very accurate (i.e. 1e-6), the algorithm would 
     * iterate for a long time trying to find the best solution when indeed 
     * there is no need to do that if a reasonable threshold has already been
     * reached.
     * Because of this behaviour the stop threshold can be set to a value much
     * lower than the one typically used in RANSAC, and yet the algorithm could
     * still produce even smaller thresholds in estimated results
     * @param stopThreshold stop threshold to stop the algorithm prematurely 
     * when a certain accuracy has been reached
     * @throws IllegalArgumentException if provided value is zero or negative
     * @throws LockedException if robust estimator is locked because an 
     * estimation is already in progress
     */
    public void setStopThreshold(double stopThreshold) 
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (stopThreshold <= MIN_STOP_THRESHOLD) {
            throw new IllegalArgumentException();
        }
        
        mStopThreshold = stopThreshold;
    }
            
    /**
     * Returns quality scores corresponding to each provided point.
     * The larger the score value the betther the quality of the sampled point
     * @return quality scores corresponding to each point
     */
    @Override
    public double[] getQualityScores() {
        return mQualityScores;
    }
    
    /**
     * Sets quality scores corresponding to each provided point.
     * The larger the score value the better the quality of the sampled point.
     * @param qualityScores quality scores corresponding to each point.
     * @throws LockedException if robust estimator is locked because an 
     * estimation is already in progress.
     * @throws IllegalArgumentException if provided quality scores length is
     * smaller than required minimum size.
     */
    @Override
    public void setQualityScores(double[] qualityScores) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetQualityScores(qualityScores);
    }
    
    /**
     * Indicates if eatimator is ready to start the polynomial estimation.
     * This is true when input data (i.e. polynomial evaluations and quality 
     * scores) are provided and enough data is available.
     * @return true if estimator is ready, false otherwise
     */
    @Override
    public boolean isReady() {
        return super.isReady() && mQualityScores != null && 
                mQualityScores.length == mEvaluations.size();
    }      
    
    /**
     * Estimates polynomial.
     * @return estimated polynomial.
     * @throws LockedException if robust estimator is locked because an 
     * estimation is already in progress.
     * @throws NotReadyException if provided input data is not enough to start
     * the estimation.
     * @throws RobustEstimatorException if estimation fails for any other reason
     * (i.e. numerical instability, no solution available, etc).
     */    
    @Override
    public Polynomial estimate() throws LockedException, NotReadyException, 
            RobustEstimatorException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        PROMedSRobustEstimator<Polynomial> innerEstimator =
                new PROMedSRobustEstimator<>(
                new PROMedSRobustEstimatorListener<Polynomial>() {
                    
            //subset of evaluations picked on each iteration
            private List<PolynomialEvaluation> mSubsetEvaluations =
                    new ArrayList<>();
                    
            @Override
            public double getThreshold() {
                return mStopThreshold;
            }

            @Override
            public int getTotalSamples() {
                return mEvaluations.size();
            }

            @Override
            public int getSubsetSize() {
                return mPolynomialEstimator.getMinNumberOfEvaluations();
            }

            @Override
            public void estimatePreliminarSolutions(int[] samplesIndices, 
                    List<Polynomial> solutions) {
                mSubsetEvaluations.clear();
                for (int samplesIndex : samplesIndices) {
                    mSubsetEvaluations.add(mEvaluations.get(samplesIndex));
                }
                
                try {
                    mPolynomialEstimator.setLMSESolutionAllowed(false);
                    mPolynomialEstimator.setEvaluations(mSubsetEvaluations);
                    
                    Polynomial polynomial = mPolynomialEstimator.estimate();
                    solutions.add(polynomial);
                } catch (Exception e) {
                    //if anything fails, no solution is added
                }
            }

            @Override
            public double computeResidual(Polynomial currentEstimation, int i) {
                PolynomialEvaluation eval = mEvaluations.get(i);
                return getDistance(eval, currentEstimation);
            }

            @Override
            public boolean isReady() {
                return PROMedSPolynomialRobustEstimator.this.isReady();
            }

            @Override
            public void onEstimateStart(RobustEstimator<Polynomial> estimator) {
                if(mListener != null) {
                    mListener.onEstimateStart(
                            PROMedSPolynomialRobustEstimator.this);
                }
            }

            @Override
            public void onEstimateEnd(RobustEstimator<Polynomial> estimator) {
                if(mListener != null) {
                    mListener.onEstimateEnd(
                            PROMedSPolynomialRobustEstimator.this);
                }
            }

            @Override
            public void onEstimateNextIteration(
                    RobustEstimator<Polynomial> estimator, int iteration) {
                if(mListener != null) {
                    mListener.onEstimateNextIteration(
                            PROMedSPolynomialRobustEstimator.this, iteration);
                }
            }

            @Override
            public void onEstimateProgressChange(
                    RobustEstimator<Polynomial> estimator, float progress) {
                if(mListener != null) {
                    mListener.onEstimateProgressChange(
                            PROMedSPolynomialRobustEstimator.this, progress);
                }
            }
            
            @Override
            public double[] getQualityScores() {
                return mQualityScores;
            }                          
        });
        
        try {
            mLocked = true;
            innerEstimator.setConfidence(mConfidence);
            innerEstimator.setMaxIterations(mMaxIterations);
            innerEstimator.setProgressDelta(mProgressDelta);
            return innerEstimator.estimate();
        } finally {
            mLocked = false;
        }
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
     * Sets quality scores corresponding to each provided polynomial evaluation.
     * This method is used internally and does not check whether instance is
     * locked or not
     * @param qualityScores quality scores to be set
     * @throws IllegalArgumentException if provided quality scores length is
     * smaller than required minimum size.
     */
    private void internalSetQualityScores(double[] qualityScores) {
        if(qualityScores.length < 
                mPolynomialEstimator.getMinNumberOfEvaluations()) {
            throw new IllegalArgumentException();
        }
        
        mQualityScores = qualityScores;        
    }         
}
