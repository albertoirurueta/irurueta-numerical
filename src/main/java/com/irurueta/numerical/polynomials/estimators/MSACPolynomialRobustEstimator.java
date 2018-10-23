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
 * Finds the best polynomial using RANSAC algorithm.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class MSACPolynomialRobustEstimator extends PolynomialRobustEstimator {
    
    /**
     * Constant defining default threshold to determine whether polynomials are
     * inliers or not.
     * Threshold will be used to compare either algebraic or geometric distance
     * of estimated polynomial respect each provided evaluation.
     */
    public static final double DEFAULT_THRESHOLD = 1e-6;
    
    /**
     * Minimum value that can be set as threshold.
     * Threshold must be strictly greater than 0.0.
     */
    public static final double MIN_THRESHOLD = 0.0;
    
    /**
     * Threshold to determine whether polynomial evaluations are inlers or not
     * when testing possible estimation solutions
     */
    private double mThreshold;
    
    /**
     * Constructor.
     */
    public MSACPolynomialRobustEstimator() {
        super();
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public MSACPolynomialRobustEstimator(int degree)
            throws IllegalArgumentException {
        super(degree);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided number of evaluations is
     * less than the required minimum.
     */
    public MSACPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations) 
            throws IllegalArgumentException {
        super(evaluations);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     */
    public MSACPolynomialRobustEstimator(
            PolynomialRobustEstimatorListener listener) {
        super(listener);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1 or if
     * provided number of evaluations is less than the required minimum for
     * provided degree.
     */
    public MSACPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations)
            throws IllegalArgumentException {
        super(degree, evaluations);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public MSACPolynomialRobustEstimator(int degree,
            PolynomialRobustEstimatorListener listener)
            throws IllegalArgumentException {
        super(degree, listener);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events such as when estimation
     * starts, ends or its progress significantly changes.
     * @throws IllegalArgumentException if provided number of evaluations is 
     * less than the required minimum.
     */
    public MSACPolynomialRobustEstimator(
            List<PolynomialEvaluation> evaluations, 
            PolynomialRobustEstimatorListener listener) 
            throws IllegalArgumentException {
        super(evaluations, listener);
        mThreshold = DEFAULT_THRESHOLD;
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
    public MSACPolynomialRobustEstimator(int degree,
            List<PolynomialEvaluation> evaluations,
            PolynomialRobustEstimatorListener listener)
            throws IllegalArgumentException {
        super(degree, evaluations, listener);
        mThreshold = DEFAULT_THRESHOLD;
    }
    
    /**
     * Returns threshold to determine whether polynomials are inliers or not
     * when testing possible estimation solutions.
     * @return threshold to determine whether polynomials are inliers or not
     * when testing possible estimation solutions.
     */
    public double getThreshold() {
        return mThreshold;
    }
    
    /**
     * Sets threshold to determine whether polynomials are inliers or not when
     * testing possible estimation solutions.
     * @param threshold threshold to determine whether polynomials are inliers
     * or not when testing possible estimation solutions.
     * @throws IllegalArgumentException if provided value is equal or less than
     * zero.
     * @throws LockedException if robust estimator is locked.
     */
    public void setThreshold(double threshold) throws IllegalArgumentException,
            LockedException {
        if(isLocked()) throw new LockedException();
        if(threshold <= MIN_THRESHOLD) throw new IllegalArgumentException();
        mThreshold = threshold;
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
        if(isLocked()) throw new LockedException();
        if(!isReady()) throw new NotReadyException();
        
        MSACRobustEstimator<Polynomial> innerEstimator =
                new MSACRobustEstimator<>(
                new MSACRobustEstimatorListener<Polynomial>(){
                    
            //subset of evaluations picked on each iteration
            private List<PolynomialEvaluation> mSubsetEvaluations =
                    new ArrayList<>();
                    
            @Override
            public double getThreshold() {
                return mThreshold;
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
                return MSACPolynomialRobustEstimator.this.isReady();
            }

            @Override
            public void onEstimateStart(RobustEstimator<Polynomial> estimator) {
                if(mListener != null) {
                    mListener.onEstimateStart(
                            MSACPolynomialRobustEstimator.this);
                }
            }

            @Override
            public void onEstimateEnd(RobustEstimator<Polynomial> estimator) {
                if(mListener != null) {
                    mListener.onEstimateEnd(
                            MSACPolynomialRobustEstimator.this);
                }
            }

            @Override
            public void onEstimateNextIteration(
                    RobustEstimator<Polynomial> estimator, int iteration) {
                if(mListener != null) {
                    mListener.onEstimateNextIteration(
                            MSACPolynomialRobustEstimator.this, iteration);
                }
            }

            @Override
            public void onEstimateProgressChange(
                    RobustEstimator<Polynomial> estimator, float progress) {
                if(mListener != null) {
                    mListener.onEstimateProgressChange(
                            MSACPolynomialRobustEstimator.this, progress);
                }
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
        return RobustEstimatorMethod.MSAC;
    }        
}
