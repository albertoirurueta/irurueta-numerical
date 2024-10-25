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
package com.irurueta.numerical.roots;

import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;

import java.util.Arrays;

/**
 * This class estimates the roots of a polynomial of degree n.
 * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
 * then the array of parameters is [an, a(n-1), ... a1, a0]
 * This class supports polynomials having either real or complex parameters.
 */
public class LaguerrePolynomialRootsEstimator extends PolynomialRootsEstimator {

    // In this implementation we have increased MR and MT to increase accuracy
    // by iterating a larger but finite number of times

    /**
     * Constant that affects the number of iterations.
     */
    public static final int MR = 80;

    /**
     * Constant that affects the number of iterations.
     */
    public static final int MT = 100;

    /**
     * Maximum number of iterations.
     */
    public static final int MAXIT = MT * MR;

    /**
     * Constant considered as machine precision for Laguerre method.
     */
    public static final double LAGUER_EPS = 1e-10;

    /**
     * Constant considered as machine precision.
     */
    public static final double EPS = 1e-14;

    /**
     * Constant indicating whether roots will be refined.
     */
    public static final boolean DEFAULT_POLISH_ROOTS = true;

    /**
     * Minimum allowed length in polynomial parameters.
     */
    public static final int MIN_VALID_POLY_PARAMS_LENGTH = 2;

    /**
     * Array containing values for Laguerre method.
     */
    private static final double[] frac =
            {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};


    /**
     * Indicates if roots should be refined.
     */
    private boolean polishRoots;

    /**
     * Constructor.
     *
     * @param polishRoots Boolean to determine whether roots should be refined.
     */
    public LaguerrePolynomialRootsEstimator(final boolean polishRoots) {
        super();
        this.polishRoots = polishRoots;
    }

    /**
     * Empty constructor.
     */
    public LaguerrePolynomialRootsEstimator() {
        super();
        this.polishRoots = DEFAULT_POLISH_ROOTS;
    }

    /**
     * Constructor.
     *
     * @param polyParams  Array containing polynomial parameters.
     * @param polishRoots Boolean indicating whether roots will be refined.
     * @throws IllegalArgumentException Raised if length of provided parameters
     *                                  is not valid. It has to be greater or equal than 2.
     */
    public LaguerrePolynomialRootsEstimator(final Complex[] polyParams, final boolean polishRoots) {
        super();
        this.polishRoots = polishRoots;
        internalSetPolynomialParameters(polyParams);
    }

    /**
     * Constructor.
     *
     * @param polyParams Array containing polynomial parameters.
     * @throws IllegalArgumentException Raised if length of provided parameters
     *                                  is not valid. It has to be greater or equal than 2.
     */
    public LaguerrePolynomialRootsEstimator(final Complex[] polyParams) {
        super();
        this.polishRoots = DEFAULT_POLISH_ROOTS;
        internalSetPolynomialParameters(polyParams);
    }

    /**
     * Estimates the roots of provided polynomial.
     *
     * @throws LockedException         Raised if this instance is locked estimating a
     *                                 root.
     * @throws NotReadyException       Raised if this instance is not ready because
     *                                 polynomial parameters have not been provided.
     * @throws RootEstimationException Raised if roots cannot be estimated for
     *                                 some reason (lack of convergence, etc.).
     */
    @Override
    public void estimate() throws LockedException, NotReadyException, RootEstimationException {

        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        // polynomial must be at least degree 1
        if (polyParams.length < MIN_VALID_POLY_PARAMS_LENGTH) {
            throw new RootEstimationException();
        }

        locked = true;

        final var a = polyParams;
        roots = new Complex[a.length - 1];

        int i;
        final var its = new int[1];
        var x = new Complex();
        Complex b;
        Complex c;
        final var m = a.length - 1;

        final var ad = Arrays.copyOf(a, a.length);
        for (var j = m - 1; j >= 0; j--) {
            x.setRealAndImaginary(0.0, 0.0);
            final var adV = Arrays.copyOf(ad, j + 2);
            internalLaguer(adV, x, its);
            if (Math.abs(x.getImaginary()) <= 2.0 * EPS * Math.abs(x.getReal())) {
                x = new Complex(x.getReal(), 0.0);
            }
            roots[j] = new Complex(x);
            b = new Complex(ad[j + 1]);
            for (var jj = j; jj >= 0; jj--) {
                c = new Complex(ad[jj]);
                ad[jj] = new Complex(b);
                b.multiply(x);
                b.add(c);
            }
        }
        if (polishRoots) {
            for (var j = 0; j < m; j++) {
                internalLaguer(a, roots[j], its);
            }
        }
        for (var j = 1; j < m; j++) {
            x = new Complex(roots[j]);
            for (i = j - 1; i >= 0; i--) {
                if (roots[i].getReal() <= x.getReal()) {
                    break;
                }
                roots[i + 1] = new Complex(roots[i]);
            }
            roots[i + 1] = new Complex(x);
        }

        locked = false;
    }

    /**
     * Returns boolean indicating whether roots are refined after an initial
     * estimation.
     *
     * @return True if roots are refined, false otherwise.
     */
    public boolean areRootsPolished() {
        return polishRoots;
    }

    /**
     * Sets boolean indicating whether roots will be refined after an initial
     * estimation.
     *
     * @param enable True if roots will be refined, false otherwise.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setPolishRootsEnabled(final boolean enable) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        polishRoots = enable;
    }

    /**
     * Internal method to set parameters of a polynomial, taking into account
     * that a polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * then the array of parameters is [an, a(n - 1), ... a1, a0]
     * Polynomial parameters can be either real or complex values
     * This method does not check if this class is locked.
     *
     * @param polyParams Polynomial parameters.
     * @throws IllegalArgumentException Raised if the length of the array is not
     *                                  valid.
     */
    @Override
    protected final void internalSetPolynomialParameters(final Complex[] polyParams) {
        if (polyParams.length < MIN_VALID_POLY_PARAMS_LENGTH) {
            throw new IllegalArgumentException();
        }
        this.polyParams = polyParams;
    }

    /**
     * Internal method to compute a root after decomposing and decreasing the
     * degree of the polynomial.
     *
     * @param a   Remaining polynomial parameters (on 1st iteration, the whole
     *            polynomial is provided, on subsequent iterations, the polynomial is
     *            deflated and the degree is reduced).
     * @param x   Estimated root.
     * @param its number of iterations needed to achieve the estimation.
     * @throws RootEstimationException Raised if root couldn't be estimated
     *                                 because of lack of convergence.
     */
    private void internalLaguer(final Complex[] a, final Complex x, final int[] its) throws RootEstimationException {

        Complex x1;
        Complex b;
        Complex g;
        Complex g2;
        final var dx = new Complex();
        final var d = new Complex();
        final var f = new Complex();
        final var h = new Complex();
        final var sq = new Complex();
        var gp = new Complex();
        final var gm = new Complex();
        final var m = a.length - 1;
        for (var iter = 1; iter <= MAXIT; iter++) {
            its[0] = iter;
            b = new Complex(a[m]);
            var err = b.getModulus();
            d.setRealAndImaginary(0.0, 0.0);
            f.setRealAndImaginary(0.0, 0.0);
            final var abx = x.getModulus();
            for (var j = m - 1; j >= 0; j--) {
                f.multiply(x);
                f.add(d);

                d.multiply(x);
                d.add(b);

                b.multiply(x);
                b.add(a[j]);

                err = b.getModulus() + abx * err;
            }
            err *= LAGUER_EPS;
            if (b.getModulus() <= err) {
                return;
            }
            g = d.divideAndReturnNew(b);
            g2 = g.powAndReturnNew(2.0);
            f.divide(b, h);
            h.multiplyByScalar(-2.0);
            h.add(g2);

            h.multiplyByScalar(m, sq);
            sq.subtract(g2);
            sq.multiplyByScalar(m - 1.0);
            sq.sqrt();

            g.add(sq, gp);
            g.subtract(sq, gm);

            final var abp = gp.getModulus();
            final var abm = gm.getModulus();
            if (abp < abm) {
                gp = gm;
            }
            if (Math.max(abp, abm) > 0.0) {
                dx.setRealAndImaginary(m, 0.0);
                dx.divide(gp);
            } else {
                dx.setModulusAndPhase(1.0 + abx, iter);
            }
            x1 = x.subtractAndReturnNew(dx);
            if (x.equals(x1)) {
                return;
            }
            if (iter % MT != 0) {
                x.copyFrom(x1);
            } else {
                int pos = Math.min(iter / MT, frac.length - 1);
                x.subtract(dx.multiplyByScalarAndReturnNew(frac[pos]));
            }
        }
        // too many iterations in Laguerre
        locked = false;
        throw new RootEstimationException();
    }
}
