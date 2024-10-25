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
import com.irurueta.algebra.GaussJordanElimination;
import com.irurueta.algebra.Matrix;
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.NotReadyException;

import java.util.Arrays;

/**
 * Fits provided data (x,y) to a function made of a linear combination of
 * functions used as a basis (i.e. f(x) = a * f0(x) + b * f1(x) + ...).
 * Where f0, f1, ... is the function basis which ideally should be formed by
 * orthogonal function.
 * This class is based on the implementation available at Numerical Recipes
 * 3rd Ed, page 791
 */
public class SimpleSingleDimensionLinearFitter extends SingleDimensionLinearFitter {

    /**
     * Determines which parameters can be modified during estimation (if true)
     * and which ones are locked (if false)
     */
    private boolean[] ia;


    /**
     * Constructor.
     */
    public SimpleSingleDimensionLinearFitter() {
        super();
    }

    /**
     * Constructor.
     *
     * @param x   input points x where a linear single dimensional function f(x) =
     *            a * f0(x) + b * f1(x) + ...
     * @param y   result of evaluation of linear single dimensional function f(x)
     *            at provided x points.
     * @param sig standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public SimpleSingleDimensionLinearFitter(final double[] x, final double[] y, final double[] sig) {
        super(x, y, sig);
    }

    /**
     * Constructor.
     *
     * @param x   input points x where a linear single dimensional function f(x) =
     *            a * f0(x) + b * f1(x) + ...
     * @param y   result of evaluation of linear single dimensional function f(x)
     *            at provided x points.
     * @param sig standard deviation of all pair of points assuming that
     *            standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public SimpleSingleDimensionLinearFitter(final double[] x, final double[] y, final double sig) {
        super(x, y, sig);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @throws FittingException if evaluation fails.
     */
    public SimpleSingleDimensionLinearFitter(
            final LinearFitterSingleDimensionFunctionEvaluator evaluator) throws FittingException {
        super();
        setFunctionEvaluator(evaluator);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @param x         input points x where a linear single dimensional function f(x) =
     *                  a * f0(x) + b * f1(x) + ...
     * @param y         result of evaluation of linear single dimensional function f(x)
     *                  at provided x points.
     * @param sig       standard deviation of all pair of points assuming that
     *                  standard deviations are constant.
     * @throws FittingException         if evaluation fails.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length .
     */
    public SimpleSingleDimensionLinearFitter(
            final LinearFitterSingleDimensionFunctionEvaluator evaluator, final double[] x, final double[] y,
            final double[] sig) throws FittingException {
        super(x, y, sig);
        setFunctionEvaluator(evaluator);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @param x         input points x where a linear single dimensional function f(x) =
     *                  a * f0(x) + b * f1(x) + ...
     * @param y         result of evaluation of linear single dimensional function f(x)
     *                  at provided x points.
     * @param sig       standard deviation of all pair of points assuming that
     *                  standard deviations are constant.
     * @throws FittingException         if evaluation fails.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public SimpleSingleDimensionLinearFitter(
            final LinearFitterSingleDimensionFunctionEvaluator evaluator, final double[] x, final double[] y,
            final double sig) throws FittingException {
        super(x, y, sig);
        setFunctionEvaluator(evaluator);
    }

    /**
     * Sets function evaluator to evaluate function at a given point and obtain
     * the evaluation of function basis at such point.
     *
     * @param evaluator function evaluator.
     * @throws FittingException if evaluation fails.
     */
    @Override
    public final void setFunctionEvaluator(
            final LinearFitterSingleDimensionFunctionEvaluator evaluator) throws FittingException {
        super.setFunctionEvaluator(evaluator);
        if (ma > 0) {
            ia = new boolean[ma];
            Arrays.fill(ia, true);
        }
    }

    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value.
     *
     * @throws FittingException  if fitting fails.
     * @throws NotReadyException if enough input data has not yet been provided.
     */
    @Override
    @SuppressWarnings("Duplicates")
    public void fit() throws FittingException, NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }

        try {
            resultAvailable = false;

            int i;
            int j;
            int k;
            int l;
            int m;
            var mfit = 0;
            double ym;
            double wt;
            double sum;
            double sig2i;
            for (j = 0; j < ma; j++) {
                if (ia[j]) {
                    mfit++;
                }
            }
            if (mfit == 0) {
                throw new FittingException("lfit: no parameters to be fitted");
            }
            final var temp = new Matrix(mfit, mfit);
            final var beta = new Matrix(mfit, 1);
            for (i = 0; i < ndat; i++) {
                evaluator.evaluate(x[i], afunc);
                ym = y[i];
                if (mfit < ma) {
                    for (j = 0; j < ma; j++) {
                        if (!ia[j]) ym -= a[j] * afunc[j];
                    }
                }
                sig2i = 1.0 / Math.pow(sig[i], 2.0);
                for (j = 0, l = 0; l < ma; l++) {
                    if (ia[l]) {
                        wt = afunc[l] * sig2i;
                        int index;
                        for (k = 0, m = 0; m <= l; m++) {
                            if (ia[m]) {
                                index = temp.getIndex(j, k++);
                                temp.getBuffer()[index] += wt * afunc[m];
                            }
                        }
                        index = beta.getIndex(j++, 0);
                        beta.getBuffer()[index] += ym * wt;
                    }
                }
            }
            for (j = 1; j < mfit; j++) {
                for (k = 0; k < j; k++) {
                    temp.setElementAt(k, j, temp.getElementAt(j, k));
                }
            }
            GaussJordanElimination.process(temp, beta);
            for (j = 0, l = 0; l < ma; l++) {
                if (ia[l]) {
                    a[l] = beta.getElementAt(j++, 0);
                }
            }
            chisq = 0.0;
            for (i = 0; i < ndat; i++) {
                evaluator.evaluate(x[i], afunc);
                sum = 0.0;
                for (j = 0; j < ma; j++) {
                    sum += a[j] * afunc[j];
                }
                chisq += Math.pow((y[i] - sum) / sig[i], 2.0);
            }
            for (j = 0; j < mfit; j++) {
                for (k = 0; k < mfit; k++) {
                    covar.setElementAt(j, k, temp.getElementAt(j, k));
                }
            }
            for (i = mfit; i < ma; i++) {
                for (j = 0; j < i + 1; j++) {
                    covar.setElementAt(i, j, 0.0);
                    covar.setElementAt(j, i, 0.0);
                }
            }
            k = mfit - 1;
            for (j = ma - 1; j >= 0; j--) {
                if (ia[j]) {
                    for (i = 0; i < ma; i++) {
                        swap(covar.getBuffer(), covar.getBuffer(),
                                covar.getIndex(i, k), covar.getIndex(i, j));
                    }
                    for (i = 0; i < ma; i++) {
                        swap(covar.getBuffer(), covar.getBuffer(),
                                covar.getIndex(k, i), covar.getIndex(j, i));
                    }
                    k--;
                }
            }

            resultAvailable = true;

        } catch (final AlgebraException | EvaluationException e) {
            throw new FittingException(e);
        }
    }

    /**
     * Prevents parameter at position i of linear combination of basis functions
     * to be modified during function fitting.
     *
     * @param i   position of parameter to be retained.
     * @param val value to be set for parameter at position i.
     */
    public void hold(final int i, final double val) {
        ia[i] = false;
        a[i] = val;
    }

    /**
     * Releases parameter at position i of linear combination of basis functions,
     * so it can be modified again if needed.
     *
     * @param i position of parameter to be released.
     */
    public void free(final int i) {
        ia[i] = true;
    }

    /**
     * Swaps values of arrays at provided positions.
     *
     * @param array1 1st array.
     * @param array2 2nd array.
     * @param pos1   1st position.
     * @param pos2   2nd position.
     */
    private static void swap(final double[] array1, final double[] array2, final int pos1, final int pos2) {
        final var value1 = array1[pos1];
        final var value2 = array2[pos2];
        array1[pos1] = value2;
        array2[pos2] = value1;
    }
}
