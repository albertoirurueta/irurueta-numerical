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

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.Utils;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;

import java.util.List;

/**
 * This class defines an LMSE (Least Mean Square Error) estimator of a 
 * polynomial of a given degree using points where polynomials (or its 
 * derivatives or integrals) are evaluated.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class LMSEPolynomialEstimator extends PolynomialEstimator {
    
    /**
     * Indicates if by default an LMSE (Least Mean Square Error) solution is
     * allowed if more evaluations than the required minimum are provided.
     */
    public static final boolean DEFAULT_ALLOW_LMSE_SOLUTION = false;
    
    /**
     * Indicates if an LMSE (Least Mean Square Error) solution is allowed if
     * more evaluations than the required minimum are provided. If false, the
     * exceeding evaluations are ignored, and only the first minimum required
     * are used.
     */
    private boolean mAllowLMSESolution;
    
    /**
     * Constructor.
     */
    public LMSEPolynomialEstimator() {
        super();
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMSEPolynomialEstimator(int degree) {
        super(degree);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     */
    public LMSEPolynomialEstimator(List<PolynomialEvaluation> evaluations) {
        super(evaluations);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param listener listener to be notified of events.
     */
    public LMSEPolynomialEstimator(PolynomialEstimatorListener listener) {
        super(listener);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMSEPolynomialEstimator(int degree, 
            List<PolynomialEvaluation> evaluations) {
        super(degree, evaluations);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMSEPolynomialEstimator(int degree, 
            PolynomialEstimatorListener listener) {
        super(degree, listener);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events.
     */
    public LMSEPolynomialEstimator(List<PolynomialEvaluation> evaluations,
            PolynomialEstimatorListener listener) {
        super(evaluations, listener);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public LMSEPolynomialEstimator(int degree, 
            List<PolynomialEvaluation> evaluations, 
            PolynomialEstimatorListener listener) {
        super(degree, evaluations, listener);
        mAllowLMSESolution = DEFAULT_ALLOW_LMSE_SOLUTION;
    }
    
    /**
     * Indicates if an LMSE (Least Mean Square Error) solution is allowed if 
     * more evaluations than the required minimum are provided. If false, the
     * exceeding evaluations are ignored, and only the first minimum required 
     * are used.
     * @return true if LMSE solution is allowed, false otherwise.
     */
    public boolean isLMSESolutionAllowed() {
        return mAllowLMSESolution;
    }
    
    /**
     * Specified if an LMSE (Least Mean Square Error) solution is allowed if 
     * more evaluations than the required minimum are provided. If false, the
     * exceeding evaluations are ignored, and only the first minimum required
     * are used.
     * @param allowed true if LMSE solution is allowed, false otherwise.
     * @throws LockedException if estimator is locked.
     */
    public void setLMSESolutionAllowed(boolean allowed) throws LockedException {
        if(isLocked()) throw new LockedException();
        mAllowLMSESolution = allowed;
    }

    /**
     * Estimates a polynomial based on provided evaluations.
     * @return estimated polynomial.
     * @throws LockedException if estimator is locked.
     * @throws NotReadyException if estimator is not ready.
     * @throws PolynomialEstimationException if polynomial estimation fails.
     */
    @Override
    public Polynomial estimate() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        try {
            mLocked = true;
            if(mListener != null) {
                mListener.onEstimateStart(this);
            }
            
            int minNumberOfEvaluations = getMinNumberOfEvaluations();
            
            Matrix a = new Matrix(mEvaluations.size(), mDegree + 1);
            Matrix b = new Matrix(mEvaluations.size(), 1);
            
            int counter = 0;
            for(PolynomialEvaluation evaluation : mEvaluations) {
                switch(evaluation.getType()) {
                    case DIRECT_EVALUATION:
                        fillDirectEvaluation(
                                (DirectPolynomialEvaluation)evaluation, a, b, 
                                counter);
                        break;
                    case DERIVATIVE_EVALUATION:
                        fillDerivativeEvaluation(
                                (DerivativePolynomialEvaluation)evaluation, a, 
                                b, counter);
                        break;
                    case INTEGRAL_EVALUATION:
                        fillIntegralEvaluation(
                                (IntegralPolynomialEvaluation)evaluation, a, b, 
                                counter);
                        break;
                    case INTEGRAL_INTERVAL:
                        fillIntegralIntervalEvaluation(
                                (IntegralIntervalPolynomialEvaluation)evaluation, 
                                a, b, counter);
                        break;
                    default:
                        continue;
                }
                
                normalize(a, b, counter);                
                counter++;
                
                if(!isLMSESolutionAllowed() && 
                        counter >= minNumberOfEvaluations) {
                    break;
                }
            }
            
            Matrix params = Utils.solve(a, b);
            
            Polynomial result = new Polynomial(params.toArray());
            
            if(mListener != null) {
                mListener.onEstimateEnd(this);
            }
            
            return result;
        } catch (AlgebraException e) {
            throw new PolynomialEstimationException(e);
        } finally {
            mLocked = false;
        }
    }

    /**
     * Returns type of polynomial estimator.
     * @return type of polynomial estimator.
     */
    @Override
    public PolynomialEstimatorType getType() {
        return PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR;
    }    
}
