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
 * Base class to fit a multi variate function [y1, y2, ...] = f([x1, x2, ...])
 * by using provided data (x, y).
 */
public abstract class MultiVariateFitter extends Fitter {
    /**
     * Input points x where a multi dimensional function f(x1, x2, ...) is
     * evaluated where each column of the matrix represents each dimension of
     * the point and each row is related to each sample corresponding to
     * provided y pairs of values.
     */
    protected Matrix x;

    /**
     * Result of evaluation of multi variate function f(x1, x2, ...) at
     * provided x points.
     * Each row contains the function evaluation for a given point x, and
     * each column contains values for each output variable f1, f2, ...
     */
    protected Matrix y;

    /**
     * Standard deviations of each pair of points (x, y).
     */
    protected double[] sig;

    /**
     * Number of samples (x, y) in provided input data.
     */
    protected int ndat;

    /**
     * Estimated parameters of linear single dimensional function.
     */
    protected double[] a;

    /**
     * Covariance of estimated parameters of linear single dimensional function.
     */
    protected Matrix covar;

    /**
     * Estimated chi square value of input data.
     */
    protected double chisq;

    /**
     * Constructor.
     */
    protected MultiVariateFitter() {
    }

    /**
     * Constructor.
     *
     * @param x   input points x where a multi variate function f(x1, x2, ...) is
     *            evaluated.
     * @param y   result of evaluation of multi variate function f(x1, x2, ...) at
     *            provided x points.
     * @param sig standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided matrix rows and arrays don't
     *                                  have the same length.
     */
    protected MultiVariateFitter(final Matrix x, final Matrix y, final double[] sig) {
        setInputData(x, y, sig);
    }

    /**
     * Constructor.
     *
     * @param x   input points x where a multi variate function f(x1, x2, ...) is
     *            evaluated.
     * @param y   result of evaluation of multi variate function f(x1, x2, ...) at
     *            provided x points.
     * @param sig standard deviation of all pair of points assuming that
     *            standard deviations are constant.
     * @throws IllegalArgumentException if provided matrix rows and arrays don't
     *                                  have the same length.
     */
    protected MultiVariateFitter(final Matrix x, final Matrix y, final double sig) {
        setInputData(x, y, sig);
    }

    /**
     * Returns input points x where a multi variate function f(x1, x2, ...)
     * is evaluated and where each column of the matrix represents
     * each dimension of the point and each row is related to each sample
     * corresponding to provided y pairs of values.
     *
     * @return input point x.
     */
    public Matrix getX() {
        return x;
    }

    /**
     * Returns result of evaluation of multi variate function f(x) at provided
     * x points. This is provided as input data along with x array.
     *
     * @return result of evaluation.
     */
    public Matrix getY() {
        return y;
    }

    /**
     * Returns standard deviations of each pair of points (x,y).
     *
     * @return standard deviations of each pair of points (x,y).
     */
    public double[] getSig() {
        return sig;
    }

    /**
     * Sets required input data to start function fitting.
     *
     * @param x   input points x where a multi variate function f(x1, x2, ...)
     *            is evaluated and where each column of the matrix represents each
     *            dimension of the point and each row is related to each sample
     *            corresponding to provided y pairs of values.
     * @param y   result of evaluation of multi variate function
     *            f(x1, x2, ...) at provided x points. This is provided as input data along
     *            with x array.
     * @param sig standard deviations of each pair of points (x,y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  size.
     */
    public final void setInputData(final Matrix x, final Matrix y, final double[] sig) {
        if (x.getRows() != y.getRows() || sig.length != y.getRows()) {
            throw new IllegalArgumentException();
        }

        this.x = x;
        this.y = y;
        this.sig = sig;

        ndat = sig.length;
    }

    /**
     * Sets required input data to start function fitting and assuming constant
     * standard deviation errors in input data.
     *
     * @param x   input points x where a multi variate function f(x1, x2, ...)
     *            is evaluated and where each column of the matrix represents each
     *            dimension of the point and each row is related to each sample
     *            corresponding to provided y pairs of values.
     * @param y   result of evaluation of multi variate function
     *            f(x1, x2, ...) at provided x points. This is provided as input data along
     *            with x array.
     * @param sig standard deviations of each pair of points (x,y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  size.
     */
    public final void setInputData(final Matrix x, final Matrix y, final double sig) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException();
        }

        this.x = x;
        this.y = y;

        ndat = y.getRows();
        this.sig = new double[ndat];
        Arrays.fill(this.sig, sig);
    }

    /**
     * Returns estimated parameters of linear single dimensional function.
     *
     * @return estimated parameters.
     */
    public double[] getA() {
        return a;
    }

    /**
     * Returns covariance of estimated parameters of linear single dimensional
     * function.
     *
     * @return covariance of estimated parameters.
     */
    public Matrix getCovar() {
        return covar;
    }

    /**
     * Returns estimated chi square value of input data.
     *
     * @return estimated chi square value of input data.
     */
    public double getChisq() {
        return chisq;
    }
}
