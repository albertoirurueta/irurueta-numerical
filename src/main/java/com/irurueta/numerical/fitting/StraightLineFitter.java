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

import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.Gamma;
import com.irurueta.statistics.MaxIterationsExceededException;

/**
 * Fits provided data (x,y) to a straight line following equation y = a + b*x,
 * estimates parameters a and b their variances, covariance and their chi square
 * value.
 * This class is based on the implementation available at Numerical Recipes
 * 3rd Ed, page 784.
 */
public class StraightLineFitter extends Fitter {

    /**
     * Array containing x coordinates of input data to be fitted to a straight
     * line.
     */
    private double[] x;

    /**
     * Array containing y coordinates of input data to be fitted to a straight
     * line.
     */
    private double[] y;

    /**
     * Standard deviations of each pair of points (x,y). This is optional, if
     * not provided, variances of a and b will be estimated assuming equal
     * error for all input points.
     */
    private double[] sig;

    /**
     * Estimated "a" parameter of line following equation y = a + b*x
     */
    private double a;

    /**
     * Estimated "b" parameter of line following equation y = a + b*X
     */
    private double b;

    /**
     * Estimated standard deviation of parameter "a".
     */
    private double siga;

    /**
     * Estimated standard deviation of parameter "b".
     */
    private double sigb;

    /**
     * Estimated chi square value.
     */
    private double chi2;

    /**
     * Estimated goodness-of-fit probability (i.e. that the fit would have a
     * chi square value equal or larger than the estimated one).
     */
    private double q;

    /**
     * Estimated standard deviation of provided input data. This is only
     * estimated if array of standard deviations of input points is not provided.
     */
    private double sigdat;

    /**
     * Constructor.
     */
    public StraightLineFitter() {
        q = 1.0;
        chi2 = sigdat = 0.0;
    }

    /**
     * Constructor.
     *
     * @param x x coordinates of input data to be fitted to a straight line.
     * @param y y coordinates of input data to be fitted to a straight line.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public StraightLineFitter(final double[] x, final double[] y) {
        this();
        setInputData(x, y);
    }

    /**
     * Constructor.
     *
     * @param x   x coordinates of input data to be fitted to a straight line.
     * @param y   y coordinates of input data to be fitted to a straight line.
     * @param sig standard deviation (i.e. errors) of provided data. This is
     *            optional, if not provided, variances of a and b will be estimated
     *            assuming equal error for all input points.
     * @throws IllegalArgumentException if provided arrays don't have the same
     *                                  length.
     */
    public StraightLineFitter(final double[] x, final double[] y, final double[] sig) {
        this();
        setInputDataAndStandardDeviations(x, y, sig);
    }

    /**
     * Returns array containing x coordinates of input data to be fitted to a
     * straight line.
     *
     * @return array containing x coordinates of input data to be fitted to a
     * straight line.
     */
    public double[] getX() {
        return x;
    }

    /**
     * Returns array containing y coordinates of input data to be fitted to a
     * straight line.
     *
     * @return array containing y coordinates of input data to be fitted to a
     * straight line.
     */
    public double[] getY() {
        return y;
    }

    /**
     * Returns standard deviations of each pair of points (x,y). This is
     * optional, if not provided, variances of a and b will be estimated
     * assuming equal error for all input points.
     *
     * @return standard deviations of each pair of points (x,y).
     */
    public double[] getSig() {
        return sig;
    }

    /**
     * Sets input data to fit a straight line to.
     *
     * @param x x coordinates.
     * @param y y coordinates.
     * @throws IllegalArgumentException if arrays don't have the same length.
     */
    public final void setInputData(final double[] x, final double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException();
        }

        this.x = x;
        this.y = y;
        this.sig = null;
    }

    /**
     * Sets input data and standard deviations of input data to fit a straight
     * line to.
     *
     * @param x   x coordinates.
     * @param y   y coordinates.
     * @param sig standard deviations of each pair of points (x,y). This is
     *            optional, if not provided, variances of a and b will be estimated
     *            assuming equal error for all input points.
     * @throws IllegalArgumentException if arrays don't have the same length.
     */
    public final void setInputDataAndStandardDeviations(
            final double[] x, final double[] y, final double[] sig) {
        if (sig != null) {
            if (x.length != y.length || y.length != sig.length) {
                throw new IllegalArgumentException();
            }

            this.x = x;
            this.y = y;
            this.sig = sig;
        } else {
            setInputData(x, y);
        }
    }


    /**
     * Indicates whether this instance is ready because enough input data has
     * been provided to start the fitting process.
     *
     * @return true if this fitter is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return x != null && y != null && x.length == y.length && (sig == null || sig.length == y.length);
    }

    /**
     * Returns estimated "a" parameter of line following equation y = a + b*x
     *
     * @return estimated "a" parameter.
     */
    public double getA() {
        return a;
    }

    /**
     * Returns estimated "b" parameter of line following equation y = a + b*x
     *
     * @return estimated "b" parameter
     */
    public double getB() {
        return b;
    }

    /**
     * Returns estimated standard deviation of parameter "a".
     *
     * @return estimated standard deviation of parameter "a".
     */
    public double getSigA() {
        return siga;
    }

    /**
     * Returns estimated standard deviation of parameter "b".
     *
     * @return estimated standard deviation of parameter "b".
     */
    public double getSigB() {
        return sigb;
    }

    /**
     * Returns estimated chi square value.
     *
     * @return estimated chi square value.
     */
    public double getChi2() {
        return chi2;
    }

    /**
     * Returns estimated goodness-of-fit probability (i.e. that the fit would
     * have a chi square value equal or larger than the estimated one).
     *
     * @return estimated goodness-of-fit probability.
     */
    public double getQ() {
        return q;
    }

    /**
     * Returns estimated standard deviation of provided input data. This is only
     * estimated if array of standard deviations of input points is not provided.
     *
     * @return estimated standard deviation of provided input data.
     */
    public double getSigdat() {
        return sigdat;
    }

    /**
     * Fits a straight line following equation y = a + b*x to provided data
     * (x, y) so that parameters associated a, b can be estimated along with
     * their variances, covariance and chi square value.
     *
     * @throws FittingException  if fitting fails.
     * @throws NotReadyException if enough input data has not yet been provided.
     */
    @Override
    public void fit() throws FittingException, NotReadyException {
        if (!isReady()) {
            throw new NotReadyException();
        }

        resultAvailable = false;

        if (sig != null) {
            fitWithSig();
        } else {
            fitWithoutSig();
        }

        resultAvailable = true;
    }

    /**
     * Fits data when standard deviations of input data is provided.
     *
     * @throws FittingException if fitting fails.
     */
    private void fitWithSig() throws FittingException {
        final var gam = new Gamma();
        int i;
        double ss = 0.0;
        double sx = 0.0;
        double sy = 0.0;
        double st2 = 0.0;
        double t;
        double wt;
        final double sxoss;
        final var ndata = x.length;
        b = 0.0;
        for (i = 0; i < ndata; i++) {
            wt = 1.0 / Math.pow(sig[i], 2.0);
            ss += wt;
            sx += x[i] * wt;
            sy += y[i] * wt;
        }
        sxoss = sx / ss;
        for (i = 0; i < ndata; i++) {
            t = (x[i] - sxoss) / sig[i];
            st2 += t * t;
            b += t * y[i] / sig[i];
        }
        b /= st2;
        a = (sy - sx * b) / ss;
        siga = Math.sqrt((1.0 + sx * sx / (ss * st2)) / ss);
        sigb = Math.sqrt(1.0 / st2);
        for (i = 0; i < ndata; i++) {
            chi2 += Math.pow((y[i] - a - b * x[i]) / sig[i], 2.0);
        }
        try {
            if (ndata > 2) {
                q = gam.gammq(0.5 * (ndata - 2), 0.5 * chi2);
            }
        } catch (final MaxIterationsExceededException e) {
            throw new FittingException(e);
        }
    }

    /**
     * Fits data when standard deviations of input data is not provided.
     */
    private void fitWithoutSig() {
        int i;
        final double ss;
        var sx = 0.0;
        var sy = 0.0;
        var st2 = 0.0;
        double t;
        final double sxoss;
        final var ndata = x.length;
        b = 0.0;
        for (i = 0; i < ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata;
        sxoss = sx / ss;
        for (i = 0; i < ndata; i++) {
            t = x[i] - sxoss;
            st2 += t * t;
            b += t * y[i];
        }
        b /= st2;
        a = (sy - sx * b) / ss;
        siga = Math.sqrt((1.0 + sx * sx / (ss * st2)) / ss);
        sigb = Math.sqrt(1.0 / st2);
        for (i = 0; i < ndata; i++) {
            chi2 += Math.pow(y[i] - a - b * x[i], 2.0);
        }
        if (ndata > 2) {
            sigdat = Math.sqrt(chi2 / (ndata - 2));
        }
        siga *= sigdat;
        sigb *= sigdat;
    }
}
