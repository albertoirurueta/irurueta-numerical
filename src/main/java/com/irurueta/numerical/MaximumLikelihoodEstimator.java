/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

/**
 * Abstract class to estimate the most likely value from a series of data 
 * assumed to be normally distributed.
 * Assuming such condition, the subclasses of this class should be able to find
 * the maximum of the probability distribution of all provided input data.
 * Such probability distribution function is computed by aggregating all the 
 * samples as a series of Gaussian functions with a small sigma and centered at
 * the exact value of the sample.
 * Subclasses of this class will either build a histogram of such aggregation of
 * Gaussians to find the most likely value, or might even find a more accurate
 * result by refining the maximum using a maximum detector method.
 * 
 * It is suggested to always use the latter (implemented as 
 * AccurateMaximumLikelihoodEstimator), as it will provide much more accurate
 * results at a slightly higher computational cost.
 */
@SuppressWarnings("WeakerAccess")
public abstract class MaximumLikelihoodEstimator {

    /**
     * Default Gaussian sigma assigned to each sample.
     */
    public static final double DEFAULT_GAUSSIAN_SIGMA = 1.0;
    
    /**
     * Minimum allowed Gaussian sigma to be set for each sample. Attempting to
     * set Gaussian sigmas lower than this value will raise an exception.
     */
    public static final double MIN_GAUSSIAN_SIGMA = 0.0;
    
    /**
     * Default method to find the most likely value. By default the accurate
     * method is used which refines the solution found by using a histogram and
     * a maximum detector.
     */
    public static final MaximumLikelihoodEstimatorMethod DEFAULT_METHOD = 
            MaximumLikelihoodEstimatorMethod.ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR;
    
    /**
     * Array containing input data to be used to find the most likely value.
     */
    protected double[] inputData;
    
    /**
     * Minimum value found on provided input data array.
     */
    protected double minValue;
    
    /**
     * Maximum value found on provided input data array.
     */
    protected double maxValue;
    
    /**
     * Boolean indicating whether minimum and maximum values in array are 
     * already available.
     */
    protected boolean areMinMaxAvailable;
    
    /**
     * Boolean indicating whether this instance is locked because some 
     * computations are being done. While this instance is locked, attempting to
     * change its status or parameters will raise an exception.
     */
    protected boolean locked;
    
    /**
     * Actual Gaussian sigma to be used on each sample when aggregating Gaussian
     * functions centered at each input data sample value.
     */
    protected double gaussianSigma;
    
    /**
     * Constructor.
     * @param gaussianSigma Gaussian sigma to be used on each sample
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */
    public MaximumLikelihoodEstimator(double gaussianSigma) {
        inputData = null;
        minValue = maxValue = 0.0;
        areMinMaxAvailable = false;
        locked = false;
        internalSetGaussianSigma(gaussianSigma);
    }
    
    /**
     * Empty constructor.
     */
    public MaximumLikelihoodEstimator() {
        inputData = null;
        minValue = maxValue = 0.0;
        areMinMaxAvailable = false;
        locked = false;
        gaussianSigma = DEFAULT_GAUSSIAN_SIGMA;
    }
    
    /**
     * Constructor.
     * @param inputData Array containing input data where most likely value must
     * be estimated from.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */
    public MaximumLikelihoodEstimator(double[] inputData, double gaussianSigma) {
        this.inputData = inputData;
        minValue = maxValue = 0.0;
        areMinMaxAvailable = false;
        locked = false;
        internalSetGaussianSigma(gaussianSigma);
    }
    
    /**
     * Constructor.
     * @param minValue Minimum value assumed to be contained within input data
     * array.
     * @param maxValue Maximum value assumed to be contained within input data
     * array.
     * @param inputData Array containing input data where most likely value must
     * be estimated from.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is 
     * negative or zero, or if minValue &lt; maxValue.
     */
    public MaximumLikelihoodEstimator(double minValue, double maxValue,
            double[] inputData, double gaussianSigma) {
        this.inputData = inputData;
        locked = false;
        internalSetMinMaxValues(minValue, maxValue);
        internalSetGaussianSigma(gaussianSigma);
    }
    
    /**
     * Returns method to be used for maximum likelihood estimation on subclasses
     * of this class.
     * @return Method for maximum likelihood estimation.
     */
    public abstract MaximumLikelihoodEstimatorMethod getMethod();
    
    /**
     * Returns minimum value found on provided input data array.
     * @return Minimum value found on provided input data array.
     * @throws NotAvailableException Raised if this value has not yet been set
     * or computed.
     */
    public double getMinValue() throws NotAvailableException {
        if (!areMinMaxValuesAvailable()) {
            throw new NotAvailableException();
        }
        return minValue;
    }
    
    /**
     * Returns maximum value found on provided input data array.
     * @return Maximum value found on provided input data array.
     * @throws NotAvailableException Raised if this value has not yet been set
     * or computed.
     */
    public double getMaxValue() throws NotAvailableException {
        if (!areMinMaxValuesAvailable()) {
            throw new NotAvailableException();
        }
        return maxValue;
    }
    
    /**
     * Sets minimum and maximum value assumed to be found in input data array.
     * @param minValue Minimum value in input data array.
     * @param maxValue Maximum value in input data array.
     * @throws LockedException Exception raised if this instance is locked.
     * This method can only be executed when computations finish and this 
     * instance becomes unlocked.
     * @throws IllegalArgumentException Exception raised if minValue &lt; maxValue.
     */
    public void setMinMaxValues(double minValue, double maxValue)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetMinMaxValues(minValue, maxValue);
    }
    
    /**
     * Returns boolean indicating whether minimum and maximum values in array 
     * are already available.
     * @return Boolean indicating whether minimum and maximum values in array
     * are already available.
     */
    public boolean areMinMaxValuesAvailable() {
        return areMinMaxAvailable;
    }
    
    /**
     * Returns boolean indicating whether this instance is locked because some 
     * computations are being done. While this instance is locked, attempting to
     * change its status or parameters will raise an exception.
     * @return Returns boolean indicating whether this instance is locked.
     */
    public boolean isLocked() {
        return locked;
    }
    
    /**
     * Returns array containing input data to be used to find the most likely 
     * value.
     * @return Returns array containing input data.
     * @throws NotAvailableException Exception raised if input data has not yet
     * been provided.
     */
    public double[] getInputData() throws NotAvailableException {
        if (!isInputDataAvailable()) {
            throw new NotAvailableException();
        }
        return inputData;
    }
    
    /**
     * Sets array containing input data to be used to find the most likely 
     * value.
     * @param inputData Array containing input data.
     * @throws LockedException Exception raised if this instance is locked.
     * This method can only be executed when computations finish and this 
     * instance becomes unlocked.
     */
    public void setInputData(double[] inputData) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.inputData = inputData;
    }
    
    /**
     * Sets array containing input data to be used to find the most likely 
     * value along with the minimum and maximum values assumed to be contained
     * in it.
     * @param inputData Array containing input data.
     * @param minValue Minimum value assumed to be contained in provided input
     * data array.
     * @param maxValue Maximum value assumed to be contained in provided input
     * data array.
     * @throws LockedException Exception raised if this instance is locked.
     * This method can only be executed when computations finish and this
     * instance becomes unlocked.
     * @throws IllegalArgumentException Exception raised if minValue &lt; maxValue.
     */
    public void setInputData(double[] inputData, double minValue, 
            double maxValue) throws LockedException {
        setMinMaxValues(minValue, maxValue);
        this.inputData = inputData;
    }
    
    /**
     * Returns boolean indicating whether input data has already been provided
     * or not.
     * @return True if input data is available and can be retrieved. 
     */
    public boolean isInputDataAvailable() {
        return inputData != null;
    }
    
    /**
     * Returns boolean indicating if enough parameters have been provided in
     * order to start the computation of the maximum likelihood value.
     * Usually providing input data is enough to make this instance ready, but
     * this is dependent of specific implementations of subclasses.
     * @return True if this instance is ready to start the computation of the
     * maximum likelihood value.
     */
    public boolean isReady() {
        return isInputDataAvailable();
    }
    
    /**
     * Returns Gaussian sigma to be used on each sample when aggregating 
     * Gaussian functions centered at each input data sample value.
     * @return Gaussian sigma to be used on each sample.
     */
    public double getGaussianSigma() {
        return gaussianSigma;
    }

    /**
     * Sets Gaussian sigma to be used on each sample when aggregating Gaussian
     * functions centered at each input data sample value.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @throws LockedException Exception raised if this instance is locked.
     * This method can only be executed when computations finish and this
     * instance becomes unlocked
     * @throws IllegalArgumentException Exception raised if provided Gaussian
     * sigma is negative or zero.
     */
    public void setGaussianSigma(double gaussianSigma) 
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetGaussianSigma(gaussianSigma);
    }
    
    /**
     * Starts the estimation of the most likely value contained within provided
     * input data array.
     * @return The most likely value.
     * @throws LockedException Exception raised if this instance is locked.
     * This method can only be executed when computations finish and this
     * instance becomes unlocked.
     * @throws NotReadyException Exception raised if this instance is not yet
     * ready
     * @see #isReady()
     */
    public abstract double estimate() throws LockedException, NotReadyException;
    
    /**
     * Creates an instance of a subclass of this class based on provided method
     * and using provided Gaussian sigma.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @param method Method to estimate maximum likelihood value.
     * @return A maximum likelihood estimator.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */
    public static MaximumLikelihoodEstimator create(double gaussianSigma,
            MaximumLikelihoodEstimatorMethod method) {
        switch (method) {
            case HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR:
                return new HistogramMaximumLikelihoodEstimator(gaussianSigma,
                        HistogramMaximumLikelihoodEstimator.
                        DEFAULT_NUMBER_OF_BINS);
            case ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR:
            default:
                return new AccurateMaximumLikelihoodEstimator(gaussianSigma,
                        AccurateMaximumLikelihoodEstimator.
                        DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION);
                
        }
    }
    
    /**
     * Creates an instance of a subclass of this class using default maximum
     * likelihood estimation method and provided Gaussian sigma.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @return A maximum likelihood estimator.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */
    public static MaximumLikelihoodEstimator create(double gaussianSigma) {
        return create(gaussianSigma, DEFAULT_METHOD);
    }
    
    /**
     * Creates an instance of a subclass of this class using default maximum
     * likelihood estimation method and default Gaussian sigma.
     * @return A maximum likelihood estimator.
     */
    public static MaximumLikelihoodEstimator create() {
        return create(DEFAULT_GAUSSIAN_SIGMA);
    }
    
    /**
     * Creates an instance of a subclass of this class based on provided method
     * and using provided Gaussian sigma and input data array.
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @param method Method to estimate maximum likelihood value
     * @return A maximum likelihood estimator.
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */    
    public static MaximumLikelihoodEstimator create(double[] inputData,
            double gaussianSigma, MaximumLikelihoodEstimatorMethod method) {
        switch (method) {
            case HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR:
                return new HistogramMaximumLikelihoodEstimator(inputData,
                        gaussianSigma, HistogramMaximumLikelihoodEstimator.
                        DEFAULT_NUMBER_OF_BINS);
            case ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR:
            default:
                return new AccurateMaximumLikelihoodEstimator(inputData,
                        gaussianSigma, AccurateMaximumLikelihoodEstimator.
                        DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION);
        }
    }
    
    /**
     * Creates an instance of a subclass of this class using default maximum 
     * likelihood method, provided Gaussian sigma and input data array
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @return A maximum likelihood estimator
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */        
    public static MaximumLikelihoodEstimator create(double[] inputData,
            double gaussianSigma) {
        return create(inputData, gaussianSigma, DEFAULT_METHOD);
    }
    
    /**
     * Creates an instance of a subclass of this class using default maximum 
     * likelihood method and Gaussian sigma, and provided input data array
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @return A maximum likelihood estimator
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero.
     */            
    public static MaximumLikelihoodEstimator create(double[] inputData) {
        return create(inputData, DEFAULT_GAUSSIAN_SIGMA);
    }            
    
    /**
     * Creates an instance of a subclass of this class based on provided method
     * and using provided Gaussian sigma, input data array and minimum/maximum
     * values assumed to be contained in provided array.
     * @param minValue Minimum value assumed to be contained in input data array
     * @param maxValue Maximum value assumed to be contained in input data array
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @param method Method to estimate maximum likelihood value
     * @return A maximum likelihood estimator
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero, or if minValue &lt; maxValue.
     */        
    public static MaximumLikelihoodEstimator create(double minValue,
            double maxValue, double[] inputData, double gaussianSigma,
            MaximumLikelihoodEstimatorMethod method) {
        switch (method) {
            case HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR:
                return new HistogramMaximumLikelihoodEstimator(minValue, 
                        maxValue, inputData, gaussianSigma,
                        HistogramMaximumLikelihoodEstimator.
                        DEFAULT_NUMBER_OF_BINS);
            case ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR:
            default:
                return new AccurateMaximumLikelihoodEstimator(minValue,
                        maxValue, inputData, gaussianSigma, 
                        AccurateMaximumLikelihoodEstimator.
                        DEFAULT_USE_HISTOGRAM_INITIAL_SOLUTION);
        }
    }    
    
    /**
     * Creates an instance of a subclass of this class using default maximum
     * likelihood method and using provided Gaussian sigma, input data array and
     * minimum/maximum values assumed to be contained in provided array.
     * @param minValue Minimum value assumed to be contained in input data array
     * @param maxValue Maximum value assumed to be contained in input data array
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @param gaussianSigma Gaussian sigma to be set for each sample.
     * @return A maximum likelihood estimator
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero, or if minValue &lt; maxValue.
     */            
    public static MaximumLikelihoodEstimator create(double minValue,
            double maxValue, double[] inputData, double gaussianSigma) {
        return create(minValue, maxValue, inputData, gaussianSigma, 
                DEFAULT_METHOD);
    }
    
    /**
     * Creates an instance of a subclass of this class using default maximum
     * likelihood method and default Gaussian sigma, and using provided input 
     * data array and minimum/maximum values assumed to be contained in provided 
     * array.
     * @param minValue Minimum value assumed to be contained in input data array
     * @param maxValue Maximum value assumed to be contained in input data array
     * @param inputData Array containing input data to be used for the 
     * estimation of the maximum likelihood value.
     * @return A maximum likelihood estimator
     * @throws IllegalArgumentException Raised if provided Gaussian sigma is
     * negative or zero, or if minValue &lt; maxValue.
     */                
    public static MaximumLikelihoodEstimator create(double minValue,
            double maxValue, double[] inputData) {
        return create(minValue, maxValue, inputData, DEFAULT_GAUSSIAN_SIGMA);
    }
    
    /**
     * Internal method to compute minimum and maximum values of provided input
     * data array.
     */
    protected void computeMinMaxValues() {
        if (!isInputDataAvailable()) {
            return;
        }

        minValue = Integer.MAX_VALUE;
        maxValue = -Integer.MAX_VALUE;
                
        double value;
        for (double data : inputData) {
            value = data;
            if (value < minValue) minValue = value;
            if (value > maxValue) maxValue = value;
        }
        
        areMinMaxAvailable = true;
    }

    /**
     * Method to set internally minimum and maximum value found in input data
     * array.
     * @param minValue Minimum value in input data array.
     * @param maxValue Maximum value in input data array.
     * @throws IllegalArgumentException Exception raised if minValue &lt; maxValue.
     */
    private void internalSetMinMaxValues(double minValue, double maxValue) {
        if (minValue > maxValue) {
            throw new IllegalArgumentException();
        }

        this.minValue = minValue;
        this.maxValue = maxValue;
        areMinMaxAvailable = true;
    }

    /**
     * Internal method to set Gaussian sigma to be used on each sample when
     * aggregating Gaussian functions centered at each input data sample value.
     * @param gaussianSigma Gaussian sigma to be used on each sample.
     * @throws IllegalArgumentException Exception raised if provided Gaussian
     * sigma is negative or zero.
     */
    private void internalSetGaussianSigma(double gaussianSigma) {
        if (gaussianSigma <= MIN_GAUSSIAN_SIGMA) {
            throw new IllegalArgumentException();
        }
        this.gaussianSigma = gaussianSigma;
    }
}
