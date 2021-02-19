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
import com.irurueta.algebra.Utils;
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.ChiSqDist;
import com.irurueta.statistics.MaxIterationsExceededException;

import java.util.Arrays;

/**
 * Fits provided data (x,y) to a generic non-linear function using
 * Levenberg-Marquardt iterative algorithm.
 * This class is based on the implementation available at Numerical Recipes 3rd
 * Ed, page 801.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class LevenbergMarquardtSingleDimensionFitter
        extends SingleDimensionFitter {

    /**
     * Default convergence parameter. Number of times that tolerance is assumed
     * to be reached to consider that algorithm has finished iterating.
     */
    public static final int DEFAULT_NDONE = 4;

    /**
     * Default maximum number of iterations.
     */
    public static final int DEFAULT_ITMAX = 5000;

    /**
     * Default tolerance to reach convergence.
     */
    public static final double DEFAULT_TOL = 1e-3;

    /**
     * Indicates whether covariance must be adjusted or not after fitting is finished.
     */
    public static final boolean DEFAULT_ADJUST_COVARIANCE = true;

    /**
     * Convergence parameter.
     */
    private int ndone = DEFAULT_NDONE;

    /**
     * Maximum number of iterations.
     */
    private int itmax = DEFAULT_ITMAX;

    /**
     * Tolerance to reach convergence.
     */
    private double tol = DEFAULT_TOL;

    /**
     * Evaluator of functions.
     */
    private LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator;

    /**
     * Number of function parameters to be estimated.
     */
    private int ma;

    /**
     * Determines which parameters can be modified during estimation (if true)
     * and which ones are locked (if false).
     */
    private boolean[] ia;

    /**
     * Curvature matrix.
     */
    private Matrix alpha;

    /**
     * Number of parameters to be fitted.
     */
    private int mfit = 0;

    /**
     * Mean square error.
     */
    private double mse = 0.0;

    /**
     * Indicates whether covariance must be adjusted or not.
     * When covariance adjustment is enabled, then covariance is recomputed taking
     * into account input samples, input standard deviations of the samples and
     * jacobians of the model function over estimated parameters using the following
     * expression: Cov = (J'*W*J)^-1 where:
     * Cov is the covariance of estimated parameters
     * J is a matrix containing the Jacobians of the function over estimated parameters
     * for each input parameter x. Each row of J matrix contains an evaluation of
     * the model function Jacobian for i-th input parameter x. Each column of J matrix
     * contains the partial derivative of model function over j-th estimated parameter.
     * W is the inverse of input variances. It's a diagonal matrix containing the
     * reciprocal of the input variances (squared input standard deviations). That is:
     * W = diag(w) where k element of w is wk = 1 / sigmak^2, which corresponds to
     * the k-th standard deviation of input sample k.
     * By default covariance is adjusted after fitting finishes.
     */
    private boolean adjustCovariance = DEFAULT_ADJUST_COVARIANCE;

    /**
     * Constructor.
     */
    public LevenbergMarquardtSingleDimensionFitter() {
        super();
    }

    /**
     * Constructor.
     *
     * @param x   input points x where function f(x) is evaluated.
     * @param y   result of evaluation of linear single dimensional function f(x)
     *            at provided x points.
     * @param sig standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public LevenbergMarquardtSingleDimensionFitter(
            final double[] x, final double[] y, final double[] sig) {
        super(x, y, sig);
    }

    /**
     * Constructor.
     *
     * @param x   input points x where function f(x) is evaluated.
     * @param y   result of evaluation of linear single dimensional function f(x)
     *            at provided x points.
     * @param sig standard deviation of all pair of points assuming that
     *            standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public LevenbergMarquardtSingleDimensionFitter(
            final double[] x, final double[] y, final double sig) {
        super(x, y, sig);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @throws FittingException if evaluation fails.
     */
    public LevenbergMarquardtSingleDimensionFitter(
            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator)
            throws FittingException {
        this();
        setFunctionEvaluator(evaluator);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @param x         input points x where function f(x) is evaluated.
     * @param y         result of evaluation of linear single dimensional function f(x)
     *                  at provided x points.
     * @param sig       standard deviations of each pair of points (x, y).
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     * @throws FittingException         if evaluation fails.
     */
    public LevenbergMarquardtSingleDimensionFitter(
            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator,
            final double[] x, final double[] y, final double[] sig)
            throws FittingException {
        this(x, y, sig);
        setFunctionEvaluator(evaluator);
    }

    /**
     * Constructor.
     *
     * @param evaluator evaluator to evaluate function at provided point and
     *                  obtain the evaluation of function basis at such point.
     * @param x         input points x where function f(x) is evaluated.
     * @param y         result of evaluation of linear single dimensional function f(x)
     *                  at provided x points.
     * @param sig       standard deviation of all pair of points assuming that
     *                  standard deviations are constant.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     * @throws FittingException         if evaluation fails.
     */
    public LevenbergMarquardtSingleDimensionFitter(
            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator,
            final double[] x, final double[] y, final double sig) throws FittingException {
        this(x, y, sig);
        setFunctionEvaluator(evaluator);
    }

    /**
     * Returns convergence parameter.
     *
     * @return convergence parameter.
     */
    public int getNdone() {
        return ndone;
    }

    /**
     * Sets convergence parameter.
     *
     * @param ndone convergence parameter.
     * @throws IllegalArgumentException if provided value is less than 1.
     */
    public void setNdone(final int ndone) {
        if (ndone < 1) {
            throw new IllegalArgumentException();
        }
        this.ndone = ndone;
    }

    /**
     * Returns maximum number of iterations.
     *
     * @return maximum number of iterations.
     */
    public int getItmax() {
        return itmax;
    }

    /**
     * Sets maximum number of iterations.
     *
     * @param itmax maximum number of iterations.
     * @throws IllegalArgumentException if provided value is zero or negative.
     */
    public void setItmax(final int itmax) {
        if (itmax <= 0) {
            throw new IllegalArgumentException();
        }
        this.itmax = itmax;
    }

    /**
     * Returns tolerance to reach convergence.
     *
     * @return tolerance to reach convergence.
     */
    public double getTol() {
        return tol;
    }

    /**
     * Sets tolerance to reach convergence.
     *
     * @param tol tolerance to reach convergence.
     * @throws IllegalArgumentException if provided value is zero or negative.
     */
    public void setTol(final double tol) {
        if (tol <= 0.0) {
            throw new IllegalArgumentException();
        }
        this.tol = tol;
    }

    /**
     * Returns function evaluator to evaluate function at a given point and
     * obtain function derivatives respect to each provided parameter
     *
     * @return function evaluator
     */
    public LevenbergMarquardtSingleDimensionFunctionEvaluator getFunctionEvaluator() {
        return evaluator;
    }

    /**
     * Sets function evaluator to evaluate function at a given point and obtain
     * function derivatives respect to each provided parameter
     *
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */
    public final void setFunctionEvaluator(
            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator)
            throws FittingException {
        internalSetFunctionEvaluator(evaluator);
    }

    /**
     * Indicates whether provided instance has enough data to start the function
     * fitting.
     *
     * @return true if this instance is ready to start the function fitting,
     * false otherwise
     */
    @Override
    public boolean isReady() {
        return evaluator != null && x != null && y != null &&
                x.length == y.length;
    }

    /**
     * Returns curvature matrix.
     * Curvature matrix is an indicatiion o fht curvature of the error of the
     * function being fitted on parameters dimensions.
     * The larger the curvatures are, the more likely that parameters can be correctly
     * fitted because a deep enough valley has been found to converge to an optimal
     * solution.
     * Typically curvature is proportional to the inverse of the covariance matrix.
     *
     * @return curvature matrix.
     */
    public Matrix getAlpha() {
        return alpha;
    }

    /**
     * Returns degrees of freedom of computed chi square value.
     * Degrees of fredom is equal to the number of sampled data minus the
     * number of estimated parameters.
     *
     * @return degrees of freedom of computed chi square value.
     */
    public int getChisqDegreesOfFreedom() {
        return ndat - ma;
    }

    /**
     * Gets mean square error produced by estimated parameters respect to
     * provided sample data.
     *
     * @return mean square error.
     */
    public double getMse() {
        return mse;
    }

    /**
     * Gets the probability of finding a smaller chi square value.
     * The smaller the found chi square value is, the better the fit of the estimated
     * parameters to the actual parameter.
     * Thus, the smaller the chance of finding a smaller chi square value, then the
     * better the estimated fit is.
     *
     * @return probability of finding a smaller chi square value (better fit), expressed
     * as a value between 0.0 and 1.0.
     * @throws MaxIterationsExceededException if convergence of incomplete
     *                                        gamma function cannot be reached. This is rarely thrown and happens
     *                                        usually for numerically unstable input values.
     */
    public double getP() throws MaxIterationsExceededException {
        return ChiSqDist.cdf(getChisq(), getChisqDegreesOfFreedom());
    }

    /**
     * Gets a measure of quality of estimated fit as a value between 0.0 and 1.0.
     * The larger the quality value is, the better the fit that has been estimated.
     *
     * @return measure of quality of estimated fit.
     * @throws MaxIterationsExceededException if convergence of incomplete
     *                                        gamma function cannot be reached. This is rarely thrown and happens
     *                                        usually for numerically unstable input values.
     */
    public double getQ() throws MaxIterationsExceededException {
        return 1.0 - getP();
    }

    /**
     * Indicates whether covariance must be adjusted or not.
     * When covariance adjustment is enabled, then covariance is recomputed taking
     * into account input samples, input standard deviations of the samples and
     * jacobians of the model function over estimated parameters using the following
     * expression: Cov = (J'*W*J)^-1 where:
     * Cov is the covariance of estimated parameters
     * J is a matrix containing the Jacobians of the function over estimated parameters
     * for each input parameter x. Each row of J matrix contains an evaluation of
     * the model function Jacobian for i-th input parameter x. Each column of J matrix
     * contains the partial derivative of model function over j-th estimated parameter.
     * W is the inverse of input variances. It's a diagonal matrix containing the
     * reciprocal of the input variances (squared input standard deviations). That is:
     * W = diag(w) where k element of w is wk = 1 / sigmak^2, which corresponds to
     * the k-th standard deviation of input sample k.
     * By default covariance is adjusted after fitting finishes.
     * More info about confidence os estimated parameters can be found here:
     * http://people.duke.edu/~hpgavin/ce281/lm.pdf
     * https://www8.cs.umu.se/kurser/5DA001/HT07/lectures/lsq-handouts.pdf
     * Numerical Recipes 3rd Ed, page 812
     *
     * @return true if covariance must be adjusted, false otherwise.
     */
    public boolean isCovarianceAdjusted() {
        return adjustCovariance;
    }

    /**
     * Specifies whether covariance must be adjusted or not.
     * When covariance adjustment is enabled, then covariance is recomputed taking
     * into account input samples, input standard deviations of the samples and
     * jacobians of the model function over estimated parameters using the following
     * expression: Cov = (J'*W*J)^-1 where:
     * Cov is the covariance of estimated parameters
     * J is a matrix containing the Jacobians of the function over estimated parameters
     * for each input parameter x. Each row of J matrix contains an evaluation of
     * the model function Jacobian for i-th input parameter x. Each column of J matrix
     * contains the partial derivative of model function over j-th estimated parameter.
     * W is the inverse of input variances. It's a diagonal matrix containing the
     * reciprocal of the input variances (squared input standard deviations). That is:
     * W = diag(w) where k element of w is wk = 1 / sigmak^2, which corresponds to
     * the k-th standard deviation of input sample k.
     * By default covariance is adjusted after fitting finishes.
     *
     * @param adjustCovariance true if covariance must be adjusted, false otherwise.
     */
    public void setCovarianceAdjusted(final boolean adjustCovariance) {
        this.adjustCovariance = adjustCovariance;
    }

    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value.
     * If chi square value is close to 1, the fit is usually good.
     * If it is much larger, then error cannot be properly fitted.
     * If it is close to zero, then the model overfits the error.
     * Methods {@link #getP()} and {@link #getQ()} can also be used to determine
     * the quality of the fit.
     *
     * @throws FittingException  if fitting fails.
     * @throws NotReadyException if enough input data has not yet been provided.
     */
    @Override
    public void fit() throws FittingException, NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }

        try {
            resultAvailable = false;

            int j;
            int k;
            int l;
            int iter;
            int done = 0;
            double alamda = 0.001;
            double ochisq;
            final double[] atry = new double[ma];
            final double[] beta = new double[ma];
            final double[] da = new double[ma];

            // number of parameters to be fitted
            mfit = 0;
            for (j = 0; j < ma; j++) {
                if (ia[j]) {
                    mfit++;
                }
            }

            final Matrix oneda = new Matrix(mfit, 1);
            final Matrix temp = new Matrix(mfit, mfit);

            // initialization
            mrqcof(a, alpha, beta);
            for (j = 0; j < ma; j++) {
                atry[j] = a[j];
            }

            ochisq = chisq;
            for (iter = 0; iter < itmax; iter++) {

                if (done == ndone) {
                    // last pass. Use zero alamda
                    alamda = 0.0;
                }

                for (j = 0; j < mfit; j++) {
                    // alter linearized fitting matrix, by augmenting diagonal
                    // elements
                    for (k = 0; k < mfit; k++) {
                        covar.setElementAt(j, k, alpha.getElementAt(j, k));
                    }
                    covar.setElementAt(j, j, alpha.getElementAt(j, j) * (1.0 + alamda));
                    for (k = 0; k < mfit; k++) {
                        temp.setElementAt(j, k, covar.getElementAt(j, k));
                    }
                    oneda.setElementAt(j, 0, beta[j]);
                }

                // matrix solution
                GaussJordanElimination.process(temp, oneda);

                for (j = 0; j < mfit; j++) {
                    for (k = 0; k < mfit; k++) {
                        covar.setElementAt(j, k, temp.getElementAt(j, k));
                    }
                    da[j] = oneda.getElementAt(j, 0);
                }

                if (done == ndone) {
                    // Converged. Clean up and return
                    covsrt(covar);
                    covsrt(alpha);

                    if (adjustCovariance) {
                        adjustCovariance();
                    }

                    resultAvailable = true;

                    return;
                }

                // did the trial succeed?
                for (j = 0, l = 0; l < ma; l++) {
                    if (ia[l]) {
                        atry[l] = a[l] + da[j++];
                    }
                }

                mrqcof(atry, covar, da);
                if (Math.abs(chisq - ochisq) < Math.max(tol, tol * chisq)) {
                    done++;
                }

                if (chisq < ochisq) {
                    // success, accept the new solution
                    alamda *= 0.1;
                    ochisq = chisq;
                    for (j = 0; j < mfit; j++) {
                        for (k = 0; k < mfit; k++) {
                            alpha.setElementAt(j, k, covar.getElementAt(j, k));
                        }
                        beta[j] = da[j];
                    }
                    for (l = 0; l < ma; l++) {
                        a[l] = atry[l];
                    }
                } else {
                    // failure, increase alamda
                    alamda *= 10.0;
                    chisq = ochisq;
                }
            }

            // too many iterations
            throw new FittingException("too many iterations");

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
     * Releases parameter at position i of linear combination of basis functions
     * so it can be modified again if needed.
     *
     * @param i position of parameter to be released.
     */
    public void free(final int i) {
        ia[i] = true;
    }

    /**
     * Adjusts covariance.
     * Covariance must be adjusted to produce more real results close to the scale
     * of problem, otherwise estimated covariance will just be a measure of
     * goodness similar to chi square value because it will be the inverse of
     * the curvature matrix, which is just a solution of the covariance up to scale.
     * <p>
     * Covariance is adjusted taking into account input samples, input standard
     * deviations of the samples and jacobians of the model function over estimated
     * parameters using the following expression: Cov = (J'*W*J)^-1 where:
     * Cov is the covariance of estimated parameters
     * J is a matrix containing the Jacobians of the function over estimated parameters
     * for each input parameter x. Each row of J matrix contains an evaluation of
     * the model function Jacobian for i-th input parameter x. Each column of J matrix
     * contains the partial derivative of model function over j-th estimated parameter.
     * W is the inverse of input variances. It's a diagonal matrix containing the
     * reciprocal of the input variances (squared input standard deviations). That is:
     * W = diag(w) where k element of w is wk = 1 / sigmak^2, which corresponds to
     * the k-th standard deviation of input sample k.
     *
     * @throws AlgebraException    if there are numerical instabilities.
     * @throws EvaluationException if function evaluation fails.
     */
    private void adjustCovariance() throws AlgebraException, EvaluationException {

        final Matrix invCov = new Matrix(a.length, a.length);
        final Matrix tmp1 = new Matrix(a.length, 1);
        final Matrix tmp2 = new Matrix(1, a.length);
        final Matrix tmpInvCov = new Matrix(a.length, a.length);
        final double[] derivatives = new double[a.length];
        final int chiSqrDegreesOfFreedom = getChisqDegreesOfFreedom();
        for (int i = 0; i < ndat; i++) {
            evaluator.evaluate(i, x[i], a, derivatives);

            tmp1.fromArray(derivatives);
            tmp2.fromArray(derivatives);

            tmp1.multiply(tmp2, tmpInvCov);

            final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sig[i] * sig[i]);
            tmpInvCov.multiplyByScalar(w);
            invCov.add(tmpInvCov);
        }

        covar = Utils.inverse(invCov);
    }

    /**
     * Internal method to set function evaluator to evaluate function at a given
     * point and obtain function derivatives respect to each provided parameter
     *
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */
    private void internalSetFunctionEvaluator(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator)
            throws FittingException {

        try {
            this.evaluator = evaluator;

            if (evaluator != null) {
                a = evaluator.createInitialParametersArray();
                ma = a.length;
                covar = new Matrix(ma, ma);
                alpha = new Matrix(ma, ma);
                ia = new boolean[ma];
                Arrays.fill(ia, true);
            }
        } catch (final AlgebraException e) {
            throw new FittingException(e);
        }
    }

    /**
     * Used by {@link #fit()} to evaluate the linearized fitting matrix alpha, and
     * vector beta to calculate chi square.
     *
     * @param a     estimated parameters so far.
     * @param alpha curvature (i.e. fitting) matrix.
     * @param beta  array where derivative increments for each parameter are
     *              stored.
     * @throws EvaluationException if function evaluation fails.
     */
    private void mrqcof(final double[] a, final Matrix alpha, final double[] beta)
            throws EvaluationException {

        int i;
        int j;
        int k;
        int l;
        int m;
        double ymod;
        double wt;
        double sig2i;
        double dy;
        final double[] dyda = new double[ma];

        // initialize (symmetric) alpha, beta
        for (j = 0; j < mfit; j++) {
            for (k = 0; k <= j; k++) {
                alpha.setElementAt(j, k, 0.0);
            }
            beta[j] = 0.0;
        }

        chisq = 0.0;
        mse = 0.0;
        final int degreesOfFreedom = getChisqDegreesOfFreedom();
        for (i = 0; i < ndat; i++) {
            // summation loop over all data
            ymod = evaluator.evaluate(i, x[i], a, dyda);

            sig2i = 1.0 / (sig[i] * sig[i]);
            dy = y[i] - ymod;
            for (j = 0, l = 0; l < ma; l++) {
                if (ia[l]) {
                    wt = dyda[l] * sig2i;
                    final double[] alphaBuffer = alpha.getBuffer();
                    for (k = 0, m = 0; m < l + 1; m++) {
                        if (ia[m]) {
                            final int index = alpha.getIndex(j, k++);
                            alphaBuffer[index] += wt * dyda[m];
                        }
                    }
                    beta[j++] += dy * wt;
                }
            }

            // add to mse
            mse += dy * dy / (double) Math.abs(degreesOfFreedom);

            // and find chi square
            chisq += dy * dy * sig2i / (double) degreesOfFreedom;
        }

        // fill in the symmetric side
        for (j = 1; j < mfit; j++) {
            for (k = 0; k < j; k++) {
                alpha.setElementAt(k, j, alpha.getElementAt(j, k));
            }
        }
    }

    /**
     * Expand in storage the covariance matrix covar, so as to take into account
     * parameters that are being held fixed. (For the latter, return zero
     * covariances).
     *
     * @param covar covariance matrix.
     */
    private void covsrt(Matrix covar) {
        int i;
        int j;
        int k;
        for (i = mfit; i < ma; i++) {
            for (j = 0; j < i + 1; j++) {
                covar.setElementAt(i, j, 0.0);
                covar.setElementAt(j, i, 0.0);
            }
        }

        k = mfit - 1;
        for (j = ma - 1; j >= 0; j--) {
            if (ia[j]) {
                final double[] buffer = covar.getBuffer();
                for (i = 0; i < ma; i++) {
                    final int pos1 = covar.getIndex(i, k);
                    final int pos2 = covar.getIndex(i, j);
                    swap(buffer, buffer, pos1, pos2);
                }

                for (i = 0; i < ma; i++) {
                    final int pos1 = covar.getIndex(k, i);
                    final int pos2 = covar.getIndex(j, i);
                    swap(buffer, buffer, pos1, pos2);
                }

                k--;
            }
        }
    }

    /**
     * Swaps values of arrays at provided positions.
     *
     * @param array1 1st array.
     * @param array2 2nd array.
     * @param pos1   1st position.
     * @param pos2   2nd position.
     */
    private void swap(final double[] array1, final double[] array2, final int pos1, final int pos2) {
        final double value1 = array1[pos1];
        final double value2 = array2[pos2];
        array1[pos1] = value2;
        array2[pos2] = value1;
    }
}
