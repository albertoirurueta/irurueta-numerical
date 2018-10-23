/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;

/**
 * Base class to fit provided multidimensional data (x1, x2, ..., y1, y2, ...) 
 * to a function made of a linear combination of functions used as a basis
 * (i.e. f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...).
 * Where f0, f1, ... is the function basis which ideally should be formed by
 * orthogonal functions.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public abstract class MultiDimensionLinearFitter extends MultiDimensionFitter {
    
    /**
     * Evaluator of functions
     */
    protected LinearFitterMultiDimensionFunctionEvaluator evaluator;
        
    /**
     * Array where results of function evaluations are stored
     */
    protected double[] afunc;
    
    /**
     * Number of function basis used as a linear combination of functions being
     * fitted
     */
    protected int ma;
    
    /**
     * Constructor
     */
    MultiDimensionLinearFitter() {
        super();
    }
    
    /**
     * Constructor
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public MultiDimensionLinearFitter(Matrix x, double[] y, double[] sig)
            throws IllegalArgumentException {
        super(x, y, sig);
    }

    /**
     * Constructor
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public MultiDimensionLinearFitter(Matrix x, double[] y, double sig)
            throws IllegalArgumentException {
        super(x, y, sig);
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @throws FittingException if evaluation fails
     */
    public MultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator)
            throws FittingException {
        super();
        internalSetFunctionEvaluator(evaluator);
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws FittingException if evaluation fails
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public MultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator,
            Matrix x, double[] y, double[] sig)
            throws FittingException, IllegalArgumentException {
        super(x, y, sig);
        internalSetFunctionEvaluator(evaluator);
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws FittingException if evaluation fails
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */    
    public MultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator,
            Matrix x, double[] y, double sig)
            throws FittingException, IllegalArgumentException {
        super(x, y, sig);
        internalSetFunctionEvaluator(evaluator);
    }     
    
    /**
     * Returns function evaluator to evaluate function at a given point and
     * obtain the evaluation of function basis at such point
     * @return function evaluator
     */    
    public LinearFitterMultiDimensionFunctionEvaluator getFunctionEvaluator() {
        return evaluator;
    }
    
    /**
     * Sets function evaluator to evaluate function at a given point and obtain
     * the evaluation of function basis at such point
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */
    public void setFunctionEvaluator(LinearFitterMultiDimensionFunctionEvaluator evaluator)
            throws FittingException {
        internalSetFunctionEvaluator(evaluator);
    }
    
    /**
     * Internal method to set function evaluator to evaluate function at a given 
     * point and obtain the evaluation of function basis at such point
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */    
    private void internalSetFunctionEvaluator(LinearFitterMultiDimensionFunctionEvaluator evaluator) 
            throws FittingException {
        
        try {
            this.evaluator = evaluator;    
        
            if (evaluator != null) {
                afunc = evaluator.createResultArray();
                ma = afunc.length;
                a = new double[ma];
                covar = new Matrix(ma, ma);
            }
        } catch (AlgebraException e) {
            throw new FittingException(e);
        }
    }
            
    /**
     * Indicates whether provided instance has enough data to start the function
     * fitting.
     * @return true if this instance is ready to start the function fitting, 
     * false otherwise
     */
    @Override
    public boolean isReady() {
        return evaluator != null && x != null && y != null && 
                x.getRows() == y.length && 
                x.getColumns() == evaluator.getNumberOfDimensions();
    }    
}
