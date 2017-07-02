/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.WeightedPolynomialEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 9, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.robust.WeightSelection;
import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.Utils;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.sorting.SortingException;
import java.util.List;

/**
 * This class implements a polynomial estimator using weigthed evaluations.
 * Weights can be used so that evaluations assumed to have a better quality
 * (i.e. more precisely estimated) are considered to be more relevant.
 * It is discouraged to use a large number of evaluations, even if they are
 * correctly weighted, since as the number of evaluations increase so do the
 * rounding errors.
 */
public class WeightedPolynomialEstimator extends PolynomialEstimator{
    
    /**
     * Default number of evaluations to be weighted and taken into account.
     */
    public static final int DEFAULT_MAX_EVALUATIONS = 50;
    
    /**
     * Indicates if weights are sorted by default so that largest weighted
     * evaluations are used first.
     */
    public static final boolean DEFAULT_SORT_WEIGHTS = true;
    
    /**
     * Maximum number of evaluations to be weighted and taken into account.
     */
    private int mMaxEvaluations = DEFAULT_MAX_EVALUATIONS;
    
    /**
     * Indicates if weights are sorted by default so that largest weighted
     * evaluations are used first.
     */
    private boolean mSortWeights = DEFAULT_SORT_WEIGHTS;
    
    /**
     * Array containing weights for all evaluations.
     */
    private double[] mWeights;
    
    /**
     * Constructor.
     */
    public WeightedPolynomialEstimator() {
        super();
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(int degree) 
            throws IllegalArgumentException {
        super(degree);
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each evaluation. The
     * larger the value of a weight, the most significant the correspondence 
     * will be.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size.
     */
    public WeightedPolynomialEstimator(List<PolynomialEvaluation> evaluations, 
            double[] weights) throws IllegalArgumentException {
        super();
        internalSetEvaluationsAndWeights(evaluations, weights);
    }
    
    /**
     * Constructor.
     * @param listener listener to be notified of events.
     */
    public WeightedPolynomialEstimator(PolynomialEstimatorListener listener) {
        super(listener);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each evaluation. The
     * larger the value of a weight, the most significant the correspondence 
     * will be.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size, or if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(int degree, 
            List<PolynomialEvaluation> evaluations, double[] weights)
            throws IllegalArgumentException {
        super(degree);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(int degree, 
            PolynomialEstimatorListener listener) 
            throws IllegalArgumentException {
        super(degree, listener);
    }
    
    /**
     * Constructor.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each evaluation. The
     * larger the value of a weight, the most significant the correspondence 
     * will be.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size.
     */
    public WeightedPolynomialEstimator(List<PolynomialEvaluation> evaluations, 
            double[] weights, PolynomialEstimatorListener listener)
            throws IllegalArgumentException {
        super(listener);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }
    
    /**
     * Constructor.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each evaluation. The
     * larger the value of a weight, the most significant the correspondence 
     * will be.
     * @param listener listener to be notified of events.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size, or if provided degree is less than 1.
     */
    public WeightedPolynomialEstimator(int degree, 
            List<PolynomialEvaluation> evaluations, double[] weights, 
            PolynomialEstimatorListener listener) 
            throws IllegalArgumentException {
        super(degree, listener);
        internalSetEvaluationsAndWeights(evaluations, weights);
    }
    
    /**
     * Sets collection of polynomial evaluations and their corresponding point
     * of evaluation used to determine a polynomial of required degree.
     * This method override always throws an IllegalArgumentException because it
     * is expected to provide both evaluations and their weights.
     * @param evaluations collection of polynomial evaluations.
     * @throws LockedException if this instance is locked.
     * @throws IllegalArgumentException always thrown.
     */
    @Override
    public void setEvaluations(List<PolynomialEvaluation> evaluations) 
            throws LockedException, IllegalArgumentException {
        throw new IllegalArgumentException(
                "evaluations and weights must be provided at once");
    }    
    
    /**
     * Sets collection of polynomial evaluations along with their corresponding
     * weights.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each polynomial 
     * evaluation. The larger the value of a weight, the most significant the
     * evaluation will be.
     * @throws LockedException if estimator is locked.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size.
     */
    public void setEvaluationsAndWeights(List<PolynomialEvaluation> evaluations,
            double[] weights) throws LockedException, IllegalArgumentException {
        if(isLocked()) throw new LockedException();
        internalSetEvaluationsAndWeights(evaluations, weights);
    }
    
    /**
     * Sets degree of polynomial to be estimated and collection of polynomial
     * evaluations and their corresponding point of evaluation used to determine
     * a polynomial of specified degree.
     * @param degree degree of polynomial to be estimated.
     * @param evaluations collection of polynomial evaluations.
     * @param weights array containing a weight amount for each polynomial 
     * evaluation. The larger the value of a weight, the most significant the
     * evaluation will be.
     * @throws IllegalArgumentException if provided degree is less than 1 or
     * if evaluations or weights are null or don't have the same size.
     * @throws LockedException if this instance is locked.
     */
    public void setDegreeEvaluationsAndWeights(int degree, 
            List<PolynomialEvaluation> evaluations, double[] weights) 
            throws IllegalArgumentException, LockedException {        
        setDegree(degree);
        setEvaluationsAndWeights(evaluations, weights);
    }
    
    
    /**
     * Returns array containing a weight amount for each polynomial evaluation.
     * The larger the value of a weight, the most significant the correspondence
     * will be.
     * @return array containing weights for each correspondence.
     */
    public double[] getWeights() {
        return mWeights;
    }
    
    /**
     * Returns boolean indicating whether weights have been provided and are
     * available for retrieval.
     * @return true if weights are available, false otherwise.
     */
    public boolean areWeightsAvailable() {
        return mWeights != null;
    }

    /**
     * Returns maximum number of evaluations to be weighted and taken into 
     * account.
     * @return maximum number of evaluations to be weighted.
     */
    public int getMaxEvaluations() {
        return mMaxEvaluations;
    }
    
    /**
     * Sets maximum number of evaluations to be weighted and taken into account.
     * This method must be called after setting degree, because the minimum
     * number of required evaluations will be checked based on degree of 
     * polynomial to be estimated.
     * @param maxEvaluations maximum number of evaluations to be weighted.
     * @throws IllegalArgumentException if provided value is less than the
     * minimum number of required evaluations.
     * @throws LockedException if this instance is locked.
     */
    public void setMaxEvaluations(int maxEvaluations) 
            throws IllegalArgumentException, LockedException {
        if(isLocked()) throw new LockedException();
        if(maxEvaluations < getMinNumberOfEvaluations()) {
            throw new IllegalArgumentException();
        }
        mMaxEvaluations = maxEvaluations;
    }
    
    /**
     * Indicates if weights are sorted by so that largest weighted evaluations
     * are used first.
     * @return true if weights are sorted, false otherwise.
     */
    public boolean isSortWeightsEnabled() {
        return mSortWeights;
    }
    
    /**
     * Specifies whether weights are sorted by so that largest weighted 
     * evaluations are used first.
     * @param sortWeights true if weights are sorted, false otherwise.
     * @throws LockedException if this instance is locked.
     */
    public void setSortWeightsEnabled(boolean sortWeights)
            throws LockedException {
        if(isLocked()) throw new LockedException();
        
        mSortWeights = sortWeights;
    }

    /**
     * Indicates if this estimator is ready to start the estimation.
     * Estimator will be ready once enough evaluations and weights are provided.
     * @return true if estimator is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return super.isReady() && areWeightsAvailable() &&
                mEvaluations.size() == mWeights.length;
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
        if(isLocked()) throw new LockedException();
        if(!isReady()) throw new NotReadyException();
        
        try {
            mLocked = true;
            if(mListener != null) {
                mListener.onEstimateStart(this);
            }
            
            WeightSelection selection = WeightSelection.selectWeights(mWeights, 
                    mSortWeights, mMaxEvaluations);
            boolean[] selected = selection.getSelected();
            int nEvaluations = selection.getNumSelected();
            
            
            Matrix a = new Matrix(nEvaluations, mDegree + 1);
            Matrix b = new Matrix(nEvaluations, 1);
            
            int index = 0, counter = 0;
            double weight;
            for(PolynomialEvaluation evaluation : mEvaluations) {
                if(selected[index]) {
                    weight = mWeights[index];
                    
                    switch(evaluation.getType()) {
                        case DIRECT_EVALUATION:
                            fillDirectEvaluation(
                                    (DirectPolynomialEvaluation)evaluation, a, 
                                    b, counter);
                            break;
                        case DERIVATIVE_EVALUATION:
                            fillDerivativeEvaluation(
                                    (DerivativePolynomialEvaluation)evaluation, 
                                    a, b, counter);
                            break;
                        case INTEGRAL_EVALUATION:
                            fillIntegralEvaluation(
                                    (IntegralPolynomialEvaluation)evaluation, a,
                                    b, counter);
                            break;
                        case INTEGRAL_INTERVAL:
                            fillIntegralIntervalEvaluation(
                                    (IntegralIntervalPolynomialEvaluation)evaluation, 
                                    a, b, counter);
                            break;
                        default:
                            continue;
                    }

                    normalize(a, b, counter, weight);
                    counter++;
                }
                       
                index++;
            }
            
            Matrix params = Utils.solve(a, b);
            
            Polynomial result = new Polynomial(params.toArray());
            
            if(mListener != null) {
                mListener.onEstimateEnd(this);
            }
            
            return result;            
        } catch (AlgebraException e) {
            throw new PolynomialEstimationException(e);
        } catch (SortingException e) {
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
        return PolynomialEstimatorType.WEIGHTED_POLyNOMIAL_ESTIMATOR;
    }
    
    /**
     * Normalizes rows of system matrix and values matrix to increase accuracy 
     * of linear system of equations to be solved.
     * @param a system matrix.
     * @param b values matrix.
     * @param row row to normalize.
     * @param weight weight.
     */
    private void normalize(Matrix a, Matrix b, int row, double weight) {
        double sqrNorm = 0.0;
        for(int i = 0; i < a.getColumns(); i++) {
            sqrNorm += Math.pow(a.getElementAt(row, i), 2.0);
        }
        sqrNorm += Math.pow(b.getElementAtIndex(row), 2.0);
        
        double norm = Math.sqrt(sqrNorm);
        double factor = weight / norm;
        
        for(int i = 0; i < a.getColumns(); i++) {
            a.setElementAt(row, i, a.getElementAt(row, i) * factor);
        }
        b.setElementAtIndex(row, b.getElementAtIndex(row) * factor);
    }    
    
    /**
     * Internal method to set evaluations and weights.
     * @param evaluations evaluations.
     * @param weights weights.
     * @throws IllegalArgumentException if evaluations or weights are null or
     * don't have the same size.
     */
    private void internalSetEvaluationsAndWeights(
            List<PolynomialEvaluation> evaluations, double[] weights) 
            throws IllegalArgumentException  {
        if(weights == null || evaluations == null ||
                weights.length != evaluations.size()) {
            throw new IllegalArgumentException();
        }
        mEvaluations = evaluations;
        mWeights = weights;
    }    
}
