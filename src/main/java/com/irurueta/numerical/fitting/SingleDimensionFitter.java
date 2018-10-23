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

import com.irurueta.algebra.Matrix;

import java.util.Arrays;

/**
 * Base class to fit a single dimension function y = f(x) by using provided
 * data (x, y)
 */
@SuppressWarnings("WeakerAccess")
public abstract class SingleDimensionFitter extends Fitter {
    /**
     * Input points x where function f(x) is evaluated.
     */
    protected double[] x;    
    
    /**
     * Result of evaluation of linear single dimensional function f(x) at 
     * provided x points. This is provided as input data along with x array.
     */
    protected double[] y;
    
    /**
     * Standard deviations of each pair of points (x,y). 
     */
    protected double[] sig;    
    
    /**
     * Number of samples (x, y) in provided input data.
     */
    protected int ndat;
    
    /**
     * Estimated parameters of single dimensional function.
     */
    protected double[] a;
    
    /**
     * Covariance of estimated parameters of single dimensional function.
     */
    protected Matrix covar;
    
    /**
     * Estimated chi square value of input data.
     */
    protected double chisq;

    /**
     * Constructor.
     */
    public SingleDimensionFitter() { }
    
    /**
     * Constructor.
     * @param x input points x where function f(x) is evaluated.
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     * length.
     */
    public SingleDimensionFitter(double[] x, double[] y, double[] sig) 
            throws IllegalArgumentException {
        setInputData(x, y, sig);
    }
    
    /**
     * Constructor.
     * @param x input points x where function f(x) is evaluated.
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length.
     */
    public SingleDimensionFitter(double[] x, double[] y, double sig)
            throws IllegalArgumentException {
        setInputData(x, y, sig);
    }
    
    /**
     * Returns input points x where function f(x) is evaluated.
     * @return input points x.
     */
    public double[] getX() {
        return x;
    }
    
    /**
     * Returns result of evaluation of linear single dimensional function f(x) 
     * at provided x points. This is provided as input data along with x array.
     * @return sampled functoin evaluations.
     */
    public double[] getY() {
        return y;
    }
    
    /**
     * Returns standard deviations of each pair of points (x,y).
     * @return standard deviations of each pair of points (x,y).
     */
    public double[] getSig() {
        return sig;
    }    
    
    /**
     * Sets required input data to start function fitting.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points. This is provided as input data along with x array.
     * @param sig standard deviations of each pair of points (x,y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     * size.
     */
    public final void setInputData(double[] x, double [] y, double[] sig) 
            throws IllegalArgumentException {
        if (x.length != y.length || sig.length != x.length) {
            throw new IllegalArgumentException();
        }
        
        this.x = x; 
        this.y = y;
        this.sig = sig;
        
        ndat = x.length;        
    }
    
    /**
     * Sets required input data to start function fitting and assuming constant
     * standard deviation errors in input data.
     * @param x input points x where a linear single dimensional function f(x) =
     * a * f0(x) + b * f1(x) + ...
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points.
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same
     * size.
     */
    public final void setInputData(double[] x, double[] y, double sig)
            throws IllegalArgumentException {
        if (x.length != y.length) {
            throw new IllegalArgumentException();
        }
        
        this.x = x;
        this.y = y;
        
        this.sig = new double[x.length];
        Arrays.fill(this.sig, sig);
        
        ndat = x.length;
    }    
    
    /**
     * Returns estimated parameters of linear single dimensional function.
     * @return estimated parameters.
     */
    public double[] getA() {
        return a;
    }
    
    /**
     * Returns covariance of estimated parameters of linear single dimensional 
     * function.
     * @return covariance of estimated parameters.
     */
    public Matrix getCovar() {
        return covar;
    }
    
    /**
     * Returns estimated chi square value of input data.
     * @return estimated chi square value of input data.
     */
    public double getChisq() {
        return chisq;
    }        
}
