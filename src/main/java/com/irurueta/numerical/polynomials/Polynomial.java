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
package com.irurueta.numerical.polynomials;

import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Complex;
import com.irurueta.numerical.NumericalException;
import com.irurueta.numerical.roots.FirstDegreePolynomialRootsEstimator;
import com.irurueta.numerical.roots.LaguerrePolynomialRootsEstimator;
import com.irurueta.numerical.roots.PolynomialRootsEstimator;
import com.irurueta.numerical.roots.SecondDegreePolynomialRootsEstimator;
import com.irurueta.numerical.roots.ThirdDegreePolynomialRootsEstimator;
import com.irurueta.numerical.signal.processing.Convolver1D;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Contains a polynomial and common operations done with polynomials.
 * This implementation only supports polynomials with real parameters.
 */
public class Polynomial implements Serializable {

    /**
     * Minimum derivative / integration order.
     */
    private static final int MIN_ORDER = 1;

    /**
     * Minimum allowed length in polynomial parameters.
     */
    public static final int MIN_VALID_POLY_PARAMS_LENGTH = 1;

    /**
     * Constant defining machine precision
     */
    public static final double EPS = 1e-10;

    /**
     * Array containing parameters defining a polynomial.
     * For a polynomial having the expression p(x) = a + b*x + c*x^2 + ...
     * provided array must be [a, b, c, ...]
     */
    private double[] mPolyParams;

    /**
     * Constructor.
     * Creates a polynomial initialized to zero.
     */
    public Polynomial() {
        mPolyParams = new double[MIN_VALID_POLY_PARAMS_LENGTH];
    }

    /**
     * Constructor.
     *
     * @param numberOfParameters number of parameters of polynomial to create.
     * @throws IllegalArgumentException if number of parameters is less than 1.
     */
    public Polynomial(final int numberOfParameters) {
        if (numberOfParameters < MIN_VALID_POLY_PARAMS_LENGTH) {
            throw new IllegalArgumentException(
                    "at least 1 parameter is required");
        }
        mPolyParams = new double[numberOfParameters];
    }

    /**
     * Constructor.
     * For a polynomial having the expression p(x) = a + b*x + c*x^2 + ...
     * provided array must be [a, b, c, ...]
     *
     * @param polyParams parameters defining a polynomial.
     * @throws IllegalArgumentException if provided array does not have at least
     *                                  length 2.
     */
    public Polynomial(final double... polyParams) {
        setPolyParams(polyParams);
    }

    /**
     * Gets array defining parameters of polynomial.
     * A polynomial having the expression p(x) = a + b*x + c*x^2 + ...
     * has an array of the form [a, b, c, ...].
     *
     * @return parameters defining a polynomial.
     */
    public double[] getPolyParams() {
        return mPolyParams;
    }

    /**
     * Sets array defining parameters of polynomial.
     * A polynomial having the expression p(x) = a + b*x + c*x^2 + ...
     * has an array of the form [a, b, c, ...].
     *
     * @param polyParams array defining parameters of polynomial. Must have at
     *                   least length 2.
     * @throws IllegalArgumentException if provided array does not have at least
     *                                  length 2.
     */
    public final void setPolyParams(final double... polyParams) {
        if (polyParams.length < MIN_VALID_POLY_PARAMS_LENGTH) {
            throw new IllegalArgumentException("must have at least length 1");
        }

        mPolyParams = polyParams;
    }

    /**
     * Gets degree of polynomial.
     *
     * @return degree of polynomial.
     */
    public int getDegree() {
        for (int i = mPolyParams.length - 1; i >= 1; i--) {
            if (Math.abs(mPolyParams[i]) > EPS) {
                return i;
            }
        }

        return 0;
    }

    /**
     * Adds this polynomial to another one and stores the result into provided
     * instance.
     *
     * @param other  other polynomial to be added.
     * @param result instance where result will be stored.
     */
    @SuppressWarnings("Duplicates")
    public void add(final Polynomial other, final Polynomial result) {
        final int maxLength = Math.max(mPolyParams.length, other.mPolyParams.length);
        final int minLength = Math.min(mPolyParams.length, other.mPolyParams.length);

        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != maxLength) {
            resultPolyParams = new double[maxLength];
        }

        for (int i = 0; i < minLength; i++) {
            resultPolyParams[i] = mPolyParams[i] + other.mPolyParams[i];
        }

        if (mPolyParams.length > other.mPolyParams.length) {
            // this is longer than other
            System.arraycopy(mPolyParams, minLength, resultPolyParams,
                    minLength, maxLength - minLength);
        } else {
            // other is longer than this
            System.arraycopy(other.mPolyParams, minLength, resultPolyParams,
                    minLength, maxLength - minLength);
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Adds another polynomial to this polynomial.
     *
     * @param other other polynomial to be added.
     */
    public void add(final Polynomial other) {
        add(other, this);
    }

    /**
     * Adds this polynomial to another one and returns a new polynomial as a
     * result.
     *
     * @param other other polynomial to be added.
     * @return a new polynomial containing the sum.
     */
    public Polynomial addAndReturnNew(final Polynomial other) {
        final int length = Math.max(mPolyParams.length, other.mPolyParams.length);
        final Polynomial result = new Polynomial(length);
        add(other, result);

        return result;
    }

    /**
     * Subtract other polynomial from this one and stores the result into
     * provided instance.
     *
     * @param other  other polynomial to be subtracted from this one.
     * @param result instance where result will be stored.
     */
    @SuppressWarnings("Duplicates")
    public void subtract(final Polynomial other, final Polynomial result) {
        final int maxLength = Math.max(mPolyParams.length, other.mPolyParams.length);
        final int minLength = Math.min(mPolyParams.length, other.mPolyParams.length);

        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != maxLength) {
            resultPolyParams = new double[maxLength];
        }

        for (int i = 0; i < minLength; i++) {
            resultPolyParams[i] = mPolyParams[i] - other.mPolyParams[i];
        }

        if (mPolyParams.length > other.mPolyParams.length) {
            // this is longer than other
            System.arraycopy(mPolyParams, minLength, resultPolyParams,
                    minLength, maxLength - minLength);
        } else {
            // other is longer than this
            for (int i = minLength; i < maxLength; i++) {
                resultPolyParams[i] = -other.mPolyParams[i];
            }
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Subtracts another polynomial form this one.
     *
     * @param other other polynomial to be subtracted from this one.
     */
    public void subtract(final Polynomial other) {
        subtract(other, this);
    }

    /**
     * Subtract other polynomial from this one and returns a new polynomial as a
     * result.
     *
     * @param other other polynomial to be subtracted from this one.
     * @return a new polynomial containing result of subtraction.
     */
    public Polynomial subtractAndReturnNew(final Polynomial other) {
        final int length = Math.max(mPolyParams.length, other.mPolyParams.length);
        final Polynomial result = new Polynomial(length);
        subtract(other, result);

        return result;
    }

    /**
     * Multiplies two polynomials.
     *
     * @param other  other polynomial to multiply with.
     * @param result instance where resulting polynomial will be stored.
     */
    public void multiply(final Polynomial other, final Polynomial result) {
        final int thisLength = mPolyParams.length;
        final int otherLength = other.mPolyParams.length;
        final int resultLength = thisLength + otherLength - 1;
        if (result.mPolyParams.length != resultLength || result == this) {
            // if length does not match or result is stored in this polynomial,
            // create new polynomial array of parameters
            result.mPolyParams = Convolver1D.convolve(mPolyParams,
                    other.mPolyParams);
        } else {
            // if length is the same, overwrite values
            Convolver1D.convolve(mPolyParams, other.mPolyParams,
                    result.mPolyParams);
        }
    }

    /**
     * Multiplies this polynomial with another one.
     *
     * @param other other polynomial to multiply with.
     */
    public void multiply(final Polynomial other) {
        multiply(other, this);
    }

    /**
     * Multiplies two polynomials and returns a new instance containing result.
     *
     * @param other other polynomial to multiply with.
     * @return a new polynomial containing result of multiplication.
     */
    public Polynomial multiplyAndReturnNew(final Polynomial other) {
        final int thisLength = mPolyParams.length;
        final int otherLength = other.mPolyParams.length;
        final int resultLength = thisLength + otherLength - 1;
        final Polynomial result = new Polynomial(resultLength);
        Convolver1D.convolve(mPolyParams, other.mPolyParams,
                result.mPolyParams);

        return result;
    }

    /**
     * Multiplies all parameters of this polynomial by an scalar and stores the
     * result into provided polynomial instance.
     *
     * @param scalar scalar to multiply parameters with.
     * @param result instance where result will be stored.
     */
    public void multiplyByScalar(final double scalar, final Polynomial result) {
        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != mPolyParams.length || result == this) {
            resultPolyParams = new double[mPolyParams.length];
        }
        ArrayUtils.multiplyByScalar(mPolyParams, scalar, resultPolyParams);
        result.mPolyParams = resultPolyParams;
    }

    /**
     * Multiplies all parameters of this polynomial by provided scalar.
     *
     * @param scalar scalar to multiply parameters with.
     */
    public void multiplyByScalar(final double scalar) {
        multiplyByScalar(scalar, this);
    }

    /**
     * Multiplies all parameters of this polynomial by an scalar and returns a
     * new polynomial containing the result.
     *
     * @param scalar scalar to multiply parameters with.
     * @return a new polynomial containing the result of the operation.
     */
    public Polynomial multiplyByScalarAndReturnNew(final double scalar) {
        final Polynomial result = new Polynomial(mPolyParams.length);
        multiplyByScalar(scalar, result);
        return result;
    }

    /**
     * Gets roots of polynomial.
     *
     * @return estimated roots of this polynomial
     * @throws NumericalException if roots estimation fails.
     */
    public Complex[] getRoots() throws NumericalException {
        final int degree = getDegree();

        final PolynomialRootsEstimator estimator;
        switch (degree) {
            case 0:
                // no roots
                return null;
            case 1:
                // first degree
                estimator = new FirstDegreePolynomialRootsEstimator(
                        mPolyParams);
                break;
            case 2:
                // second degree
                estimator = new SecondDegreePolynomialRootsEstimator(
                        mPolyParams);
                break;
            case 3:
                // third degree
                estimator = new ThirdDegreePolynomialRootsEstimator(
                        mPolyParams);
                break;
            default:
                // greater degree

                // copy real parameters into complex values
                final Complex[] polyParams = new Complex[mPolyParams.length];
                for (int i = 0; i < mPolyParams.length; i++) {
                    polyParams[i] = new Complex(mPolyParams[i]);
                }
                estimator = new LaguerrePolynomialRootsEstimator(polyParams);
                break;
        }

        estimator.estimate();
        return estimator.getRoots();
    }

    /**
     * Evaluates polynomial at provided value.
     *
     * @param x value to evaluate polynomial at.
     * @return result of polynomial evaluation.
     */
    public double evaluate(final double x) {
        double result = 0.0;
        double powX = 1.0;
        for (double polyParam : mPolyParams) {
            result += polyParam * powX;
            powX *= x;
        }

        return result;
    }

    /**
     * Computes derivative of polynomial.
     *
     * @param result instance where derivative will be stored.
     */
    @SuppressWarnings("Duplicates")
    public void derivative(final Polynomial result) {
        final int resultLength = mPolyParams.length - 1;
        final int resultLength2 = Math.max(resultLength, 1);

        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != resultLength2 || result == this) {
            resultPolyParams = new double[resultLength2];
        }
        if (resultLength == 0) {
            resultPolyParams[0] = 0.0;
        }

        for (int i = 0, j = 1; i < resultLength; i++, j++) {
            resultPolyParams[i] = j * mPolyParams[j];
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Replaces this instance by its derivative.
     */
    public void derivative() {
        derivative(this);
    }

    /**
     * Computes derivative of polynomial.
     *
     * @return a new instance containing derivative.
     */
    public Polynomial derivativeAndReturnNew() {
        final int resultLength = Math.max(mPolyParams.length - 1, 1);
        final Polynomial result = new Polynomial(resultLength);
        derivative(result);
        return result;
    }

    /**
     * Evaluates derivative of polynomial at provided value.
     *
     * @param x value to evaluate derivative of polynomial at.
     * @return result of evaluation of derivative.
     */
    public double evaluateDerivative(final double x) {
        double result = 0.0;
        double powX = 1.0;
        for (int j = 1; j < mPolyParams.length; j++) {
            result += j * mPolyParams[j] * powX;
            powX *= x;
        }

        return result;
    }

    /**
     * Computes second derivative of polynomial.
     *
     * @param result instance where second derivative will be stored.
     */
    @SuppressWarnings("Duplicates")
    public void secondDerivative(final Polynomial result) {
        final int resultLength = mPolyParams.length - 2;
        final int resultLength2 = Math.max(resultLength, 1);

        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != resultLength2 || result == this) {
            resultPolyParams = new double[resultLength2];
        }
        if (resultLength == 0) {
            resultPolyParams[0] = 0.0;
        }

        for (int i = 0, j = 2, k = 1; i < resultLength; i++, j++, k++) {
            resultPolyParams[i] = j * k * mPolyParams[j];
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Replaces this instance by its second derivative.
     */
    public void secondDerivative() {
        secondDerivative(this);
    }

    /**
     * Computes second derivative of polynomial.
     *
     * @return a new instance containing second derivative.
     */
    public Polynomial secondDerivativeAndReturnNew() {
        final int resultLength = Math.max(mPolyParams.length - 2, 1);
        final Polynomial result = new Polynomial(resultLength);
        secondDerivative(result);
        return result;
    }

    /**
     * Evaluates second derivative of polynomial at provided value.
     *
     * @param x value to evaluate second derivative of polynomial at.
     * @return result of evaluation of second derivative.
     */
    public double evaluateSecondDerivative(final double x) {
        double result = 0.0;
        double powX = 1.0;
        for (int j = 2, k = 1; j < mPolyParams.length; j++, k++) {
            result += j * k * mPolyParams[j] * powX;
            powX *= x;
        }

        return result;
    }

    /**
     * Computes nth-order derivative of polynomial.
     *
     * @param order  order of derivative to compute. Must be at least 1.
     * @param result instance where nth-order derivative will be stored.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    @SuppressWarnings("Duplicates")
    public void nthDerivative(final int order, final Polynomial result) {
        if (order < MIN_ORDER) {
            throw new IllegalArgumentException();
        }

        final int resultLength = mPolyParams.length - order;
        final int resultLength2 = Math.max(resultLength, 1);

        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != resultLength2 || result == this) {
            resultPolyParams = new double[resultLength2];
        }
        if (resultLength == 0) {
            resultPolyParams[0] = 0.0;
        }

        for (int i = 0, j = order; i < resultLength; i++, j++) {
            int param = j;
            for (int k = 1; k < order; k++) {
                param *= j - k;
            }
            resultPolyParams[i] = param * mPolyParams[j];
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Replaces this instance by its nth-order derivative.
     *
     * @param order order of derivative to compute. Must be at least 1.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    public void nthDerivative(final int order) {
        nthDerivative(order, this);
    }

    /**
     * Computes nth-order derivative of polynomial.
     *
     * @param order order of derivative to compute. Must be at least 1.
     * @return a new instance containing nth-order derivative.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    public Polynomial nthDerivativeAndReturnNew(final int order) {
        final int resultLength = Math.max(mPolyParams.length - order, 1);
        final Polynomial result = new Polynomial(resultLength);
        nthDerivative(order, result);
        return result;
    }

    /**
     * Evaluates nth-derivative of polynomial at provided value.
     *
     * @param x     value to evaluate nth-derivative of polynomial at.
     * @param order order of derivative to evaluate. Must be at least 1.
     * @return result of evaluation of nth-derivative.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    public double evaluateNthDerivative(final double x, final int order) {
        if (order < MIN_ORDER) {
            throw new IllegalArgumentException("order must be at least 1");
        }

        double result = 0.0;
        double powX = 1.0;
        for (int i = order; i < mPolyParams.length; i++) {
            int param = i;
            for (int j = 1; j < order; j++) {
                param *= i - j;
            }
            result += param * mPolyParams[i] * powX;
            powX *= x;
        }

        return result;
    }

    /**
     * Computes polynomial containing the integration of current one.
     * Because infinite polynomials exist with different constant values,
     * constant term can be provided as well.
     *
     * @param result   instance where resulting polynomial will be stored.
     * @param constant constant term.
     */
    public void integration(final Polynomial result, final double constant) {
        final int resultLength = mPolyParams.length + 1;
        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != resultLength || result == this) {
            resultPolyParams = new double[resultLength];
        }

        resultPolyParams[0] = constant;
        for (int i = 0, j = 1; i < mPolyParams.length; i++, j++) {
            resultPolyParams[j] = mPolyParams[i] / j;
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Computes polynomial containing the integration of current one and
     * assuming a zero constant term.
     *
     * @param result instance where resulting polynomial will be stored.
     */
    public void integration(final Polynomial result) {
        integration(result, 0.0);
    }

    /**
     * Updates this instance to contain its integration.
     *
     * @param constant constant term.
     */
    public void integration(final double constant) {
        integration(this, constant);
    }

    /**
     * Updates this instance to contain its integration using a zero constant
     * term.
     */
    public void integration() {
        integration(this);
    }

    /**
     * Computes polynomial containing the integration of current one.
     * Because infinite polynomials exist with different constant values,
     * constant term can be provided as well.
     *
     * @param constant constant term.
     * @return a new instance containing integration polynomial.
     */
    public Polynomial integrationAndReturnNew(final double constant) {
        final Polynomial result = new Polynomial(mPolyParams.length + 1);
        integration(result, constant);
        return result;
    }

    /**
     * Computes polynomial containing the integration of current one and
     * assuming a zero constant term.
     *
     * @return a new instance containing integration polynomial.
     */
    public Polynomial integrationAndReturnNew() {
        return integrationAndReturnNew(0.0);
    }

    /**
     * Integrate polynomial within provided interval.
     *
     * @param startX start of integration interval.
     * @param endX   end of integration interval.
     * @return result of integration.
     */
    public double integrateInterval(final double startX, final double endX) {

        double resultStart = 0.0;
        double resultEnd = 0.0;
        double powStartX = startX;
        double powEndX = endX;
        double polyParam;
        for (int i = 0, j = 1; i < mPolyParams.length; i++, j++) {
            polyParam = mPolyParams[i] / j;
            resultStart += polyParam * powStartX;
            powStartX *= startX;

            resultEnd += polyParam * powEndX;
            powEndX *= endX;
        }

        return resultEnd - resultStart;
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     * Because infinite polynomials exist with different constant values,
     * constant terms for each integration order can be provided as well.
     *
     * @param order     order of integration to compute. Must be at least 1.
     * @param result    instance where resulting polynomial will be stored.
     * @param constants constant terms for each integration order. Must have a
     *                  length equal to order if provided.
     * @throws IllegalArgumentException if provided order is less than 1 or if
     *                                  constants does not have length equal to order.
     */
    public void nthIntegration(final int order, final Polynomial result, final double[] constants) {
        if (order < MIN_ORDER) {
            throw new IllegalArgumentException("order must be at least 1");
        }
        if (constants != null && constants.length != order) {
            throw new IllegalArgumentException(
                    "length of constants must be order");
        }
        final int resultLength = mPolyParams.length + order;
        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != resultLength || result == this) {
            resultPolyParams = new double[resultLength];
        }

        for (int i = 0; i < order; i++) {
            if (constants != null) {
                int param = 1;
                for (int k = 1; k <= i; k++) {
                    param *= k;
                }
                resultPolyParams[i] = constants[i] / param;
            } else {
                resultPolyParams[i] = 0.0;
            }
        }
        for (int i = 0, j = order; i < mPolyParams.length; i++, j++) {
            int param = j;
            for (int k = 1; k < order; k++) {
                param *= j - k;
            }
            resultPolyParams[j] = mPolyParams[i] / param;
        }

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     *
     * @param order  order of integration to compute. Must be at least 1.
     * @param result instance where resulting polynomial will be stored.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    public void nthIntegration(final int order, final Polynomial result) {
        nthIntegration(order, result, null);
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     * Because infinite polynomials exist with different constant values,
     * constant terms for each integration order can be provided as well.
     *
     * @param order     order of integration to compute. Must be at least 1.
     * @param constants constant terms for each integration order. Must have a
     *                  length equal to order if provided.
     * @throws IllegalArgumentException if provided order is less than 1 or if
     *                                  constants does not have length equal to order.
     */
    public void nthIntegration(final int order, final double[] constants) {
        nthIntegration(order, this, constants);
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     *
     * @param order order of integration to compute. Must be at least 1.
     */
    public void nthIntegration(final int order) {
        nthIntegration(order, (double[]) null);
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     * Because infinite polynomials exist with different constant values,
     * constant terms for each integration order can be provided as well.
     *
     * @param order     order of integration to compute. Must be at least 1.
     * @param constants constant terms for each integration order. Must have a
     *                  length equal to order if provided.
     * @return a new polynomial containing the nth-order integration.
     * @throws IllegalArgumentException if provided order is less than 1 or if
     *                                  constants does not have length equal to order.
     */
    public Polynomial nthIntegrationAndReturnNew(final int order, final double[] constants) {
        final Polynomial result = new Polynomial();
        nthIntegration(order, result, constants);
        return result;
    }

    /**
     * Computes polynomial containing the nth-order integration of current one.
     *
     * @param order order of integration to compute. Must be at least 1.
     * @return a new polynomial containing the nth-order integration.
     * @throws IllegalArgumentException if provided order is less than 1 or if
     *                                  constants does not have length equal to order.
     */
    public Polynomial nthIntegrationAndReturnNew(final int order) {
        return nthIntegrationAndReturnNew(order, null);
    }

    /**
     * Computes nth-integration over provided interval.
     *
     * @param startX    start of integration interval.
     * @param endX      end of integration interval.
     * @param order     order of integration. Must be at least 1.
     * @param constants constant terms for each integration order. Must have a
     *                  length equal to order if provided.
     * @return result of integration.
     * @throws IllegalArgumentException if provided order is less than 1 or if
     *                                  constants does not have length equal to order.
     */
    public double nthOrderIntegrateInterval(
            final double startX, final double endX,
            final int order, final double[] constants) {
        if (order < MIN_ORDER) {
            throw new IllegalArgumentException();
        }
        if (constants != null && constants.length != order) {
            throw new IllegalArgumentException();
        }

        double resultStart = 0.0;
        double resultEnd = 0.0;
        double powStartX = 1.0;
        double powEndX = 1.0;
        double polyParam;
        for (int i = 0; i < order; i++) {
            if (constants != null) {
                int param = 1;
                for (int k = 1; k <= i; k++) {
                    param *= k;
                }
                polyParam = constants[i] / param;
                resultStart += polyParam * powStartX;
                resultEnd += polyParam * powEndX;
            }
            powStartX *= startX;
            powEndX *= endX;
        }

        for (int i = 0, j = order; i < mPolyParams.length; i++, j++) {
            int param = j;
            for (int k = 1; k < order; k++) {
                param *= j - k;
            }
            polyParam = mPolyParams[i] / param;
            resultStart += polyParam * powStartX;
            powStartX *= startX;

            resultEnd += polyParam * powEndX;
            powEndX *= endX;
        }

        return resultEnd - resultStart;
    }

    /**
     * Computes nth-integration over provided interval.
     *
     * @param startX start of integration interval.
     * @param endX   end of integration interval.
     * @param order  order of integration. Must be at least 1.
     * @return result of integration.
     * @throws IllegalArgumentException if provided order is less than 1.
     */
    public double nthOrderIntegrateInterval(
            final double startX, final double endX, final int order) {
        return nthOrderIntegrateInterval(startX, endX, order, null);
    }

    /**
     * Trims polynomial to remove all terms above degree that can be neglected.
     *
     * @param result instance where result will be stored.
     */
    public void trim(final Polynomial result) {
        final int degree = getDegree();
        final int resultLength = degree + 1;

        final double[] resultPolyParams;
        if (result.mPolyParams.length != resultLength) {
            resultPolyParams = new double[resultLength];
        } else {
            resultPolyParams = result.mPolyParams;
        }
        System.arraycopy(mPolyParams, 0, resultPolyParams, 0, resultLength);

        result.mPolyParams = resultPolyParams;
    }

    /**
     * Trims this polynomial to remove all terms above degree that can be
     * neglected.
     */
    public void trim() {
        trim(this);
    }

    /**
     * Trims this polynomial to remove all terms above degree that can be
     * neglected and returns the result as a new polynomial.
     *
     * @return a new trimmed polynomial.
     */
    public Polynomial trimAndReturnNew() {
        final Polynomial result = new Polynomial();
        trim(result);
        return result;
    }

    /**
     * Normalizes parameters of this polynomial so that the array of parameters
     * has unitary norm and stores result into provided instance.
     * Normalization keeps location of real roots, but other roots or
     * properties of polynomials might change.
     *
     * @param result instance where normalized polynomial will be stored.
     */
    public void normalize(final Polynomial result) {
        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != mPolyParams.length) {
            resultPolyParams = new double[mPolyParams.length];
        }
        ArrayUtils.normalize(mPolyParams, resultPolyParams);
        result.mPolyParams = resultPolyParams;
    }

    /**
     * Normalizes this polynomial so that the array of parameters has unitary
     * norm.
     * Normalization keeps location of real roots, but other roots or
     * properties of polynomials might change.
     */
    public void normalize() {
        normalize(this);
    }

    /**
     * Normalizes parameters of this polynomial so that the array of parameters
     * has unitary norm and returns result as a new polynomial instance.
     * Normalization keeps location of real roots, but other roots or
     * properties of polynomials might change.
     *
     * @return a new normalized polynomial instance.
     */
    public Polynomial normalizeAndReturnNew() {
        final Polynomial result = new Polynomial(mPolyParams.length);
        normalize(result);
        return result;
    }

    /**
     * Normalizes parameters of this polynomial so that the highest degree term
     * becomes 1.0 and stores result into provided instance.
     *
     * @param result instance where result of normalization will be stored.
     */
    public void normalizeHighestDegreeTerm(final Polynomial result) {
        final int degree = getDegree();
        final double term = mPolyParams[degree];
        double[] resultPolyParams = result.mPolyParams;
        if (resultPolyParams.length != mPolyParams.length) {
            resultPolyParams = new double[mPolyParams.length];
        }
        ArrayUtils.multiplyByScalar(mPolyParams, 1.0 / term, resultPolyParams);
        result.mPolyParams = resultPolyParams;
    }

    /**
     * Normalizes parameters of this polynomial so that the highest degree term
     * becomes 1.0.
     */
    public void normalizeHighestDegreeTerm() {
        normalizeHighestDegreeTerm(this);
    }

    /**
     * Normalizes parameters of this polynomial so that the highest degree term
     * becomes 1.0 and returns the result as a new instance.
     *
     * @return a new normalized polynomial.
     */
    public Polynomial normalizeHighestDegreeTermAndReturnNew() {
        final Polynomial result = new Polynomial(mPolyParams.length);
        normalizeHighestDegreeTerm(result);
        return result;
    }

    /**
     * Gets location of maxima in this polynomial.
     *
     * @return location of maxima or null if polynomial has no maxima.
     * @throws NumericalException if maxima cannot be determined due to
     *                            numerical instabilities.
     */
    public double[] getMaxima() throws NumericalException {
        return getMaxima(EPS);
    }

    /**
     * Gets location of maxima in this polynomial.
     *
     * @param threshold threshold to allow possible small deviations in first
     *                  derivative respect to pure real roots. This should be a very small
     *                  positive value.
     * @return location of maxima or null if polynomial has no maxima.
     * @throws NumericalException       if maxima cannot be determined due to
     *                                  numerical instabilities.
     * @throws IllegalArgumentException if provided threshold is negative.
     */
    @SuppressWarnings("Duplicates")
    public double[] getMaxima(final double threshold) throws NumericalException {
        if (threshold < 0.0) {
            throw new IllegalArgumentException();
        }

        final Polynomial derivative = derivativeAndReturnNew();

        // roots of derivative contains either minima or maxima.
        final Complex[] derivativeRoots = derivative.getRoots();
        final List<Complex> maxima = new ArrayList<>();
        if (derivativeRoots != null) {
            for (Complex derivativeRoot : derivativeRoots) {
                if (Math.abs(derivativeRoot.getImaginary()) > threshold) {
                    // root is imaginary (not allowed)
                    continue;
                }

                final double x = derivativeRoot.getReal();
                final double secondDerivativeEval = evaluateSecondDerivative(x);
                if (secondDerivativeEval < 0.0) {
                    // is maxima
                    maxima.add(derivativeRoot);
                }
            }
        }

        // return real parts of maxima, since we only allow real roots of first
        // derivative
        if (maxima.isEmpty()) {
            return null;
        }

        final double[] result = new double[maxima.size()];
        int i = 0;
        for (final Complex m : maxima) {
            result[i] = m.getReal();
            i++;
        }

        return result;
    }

    /**
     * Gets location of minima in this polynomial.
     *
     * @return location of minima or null if polynomial has no minima.
     * @throws NumericalException if minima cannot be determined due to
     *                            numerical instabilities.
     */
    public double[] getMinima() throws NumericalException {
        return getMinima(EPS);
    }

    /**
     * Gets location of minima in this polynomial.
     *
     * @param threshold threshold to allow possible small deviations in first
     *                  derivative respect to pure real roots. This should be a very small
     *                  positive value.
     * @return location of minima or null if polynomial has no minima.
     * @throws NumericalException       if minima cannot be determined due to
     *                                  numerical instabilities.
     * @throws IllegalArgumentException if provided threshold is negative.
     */
    @SuppressWarnings("Duplicates")
    public double[] getMinima(double threshold) throws NumericalException {
        if (threshold < 0.0) {
            throw new IllegalArgumentException();
        }

        final Polynomial derivative = derivativeAndReturnNew();

        // roots of derivative contains either minima or maxima.
        final Complex[] derivativeRoots = derivative.getRoots();
        final List<Complex> minima = new ArrayList<>();
        if (derivativeRoots != null) {
            for (final Complex derivativeRoot : derivativeRoots) {
                if (Math.abs(derivativeRoot.getImaginary()) > threshold) {
                    //root is imaginary (not allowed)
                    continue;
                }

                final double x = derivativeRoot.getReal();
                final double secondDerivativeEval = evaluateSecondDerivative(x);
                if (secondDerivativeEval >= 0.0) {
                    // is minima
                    minima.add(derivativeRoot);
                }
            }
        }

        // return real parts of minima, since we only allow real roots of first
        // derivative
        if (minima.isEmpty()) {
            return null;
        }

        final double[] result = new double[minima.size()];
        int i = 0;
        for (final Complex m : minima) {
            result[i] = m.getReal();
            i++;
        }

        return result;
    }

    /**
     * Gets location of minima or maxima (i.e. extrema) in this polynomial.
     *
     * @return location of minima or maxima, or null if polynomial has no
     * minima or maxima.
     * @throws NumericalException if minima or maxima cannot be determined due
     *                            to numerical instabilities.
     */
    public double[] getExtrema() throws NumericalException {
        return getExtrema(EPS);
    }

    /**
     * Gets location of minima or maxima (i.e. extrema) in this polynomial.
     *
     * @param threshold threshold to allow possible small deviations in first
     *                  derivative respect to pure real roots. This should be a very small
     *                  positive value.
     * @return location of minima or maxima, or null if polynomial has no minima
     * or maxima.
     * @throws NumericalException       if minima or maxima cannot be determined due
     *                                  to numerical instabilities.
     * @throws IllegalArgumentException if provided threshold is negative.
     */
    @SuppressWarnings("DuplicatedCode")
    public double[] getExtrema(final double threshold)
            throws NumericalException {
        if (threshold < 0.0) {
            throw new IllegalArgumentException("threshold must be positive");
        }

        final Polynomial derivative = derivativeAndReturnNew();

        // roots of derivative contains either minima or maxima.
        final Complex[] derivativeRoots = derivative.getRoots();
        final List<Complex> minimaOrMaxima = new ArrayList<>();
        if (derivativeRoots != null) {
            for (final Complex derivativeRoot : derivativeRoots) {
                if (Math.abs(derivativeRoot.getImaginary()) > threshold) {
                    // root is imaginary (not allowed)
                    continue;
                }

                minimaOrMaxima.add(derivativeRoot);
            }
        }

        // return real parts of roots, since we only allow real roots of first
        // derivative
        if (minimaOrMaxima.isEmpty()) {
            return null;
        }

        final double[] result = new double[minimaOrMaxima.size()];
        int i = 0;
        for (final Complex m : minimaOrMaxima) {
            result[i] = m.getReal();
            i++;
        }

        return result;
    }
}
