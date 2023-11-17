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
 * This class evaluates a multidimensional function along a line, such line is
 * defined by an input point and a given direction.
 * Using provided input point and direction, the multidimensional function's
 * input parameters are determined, so they all lay on a line.
 */
public class DirectionalEvaluator {

    /**
     * Listener to evaluate a multidimensional function.
     */
    protected MultiDimensionFunctionEvaluatorListener listener;

    /**
     * Point used as a reference to determine the function's input parameters
     * along a line.
     */
    protected double[] point;

    /**
     * Vector indicating the direction of the line where the function is
     * evaluated.
     */
    protected double[] direction;

    /**
     * Point currently being evaluated in the multidimensional function.
     * This is used internally.
     */
    protected final double[] p;

    /**
     * Constructor.
     *
     * @param listener  Listener to evaluate a multidimensional function.
     * @param point     Point used as a reference to determine the function's input
     *                  parameters along a line.
     * @param direction Vector indicating the direction of the line where the
     *                  function is evaluated.
     * @throws IllegalArgumentException Raised if point and direction don't have
     *                                  the same length.
     */
    public DirectionalEvaluator(
            final MultiDimensionFunctionEvaluatorListener listener, final double[] point,
            final double[] direction) {

        setPointAndDirection(point, direction);
        this.listener = listener;
        p = new double[point.length];
    }

    /**
     * Returns listener to evaluate a multidimensional function.
     *
     * @return Listener to evaluate a multidimensional function.
     */
    public MultiDimensionFunctionEvaluatorListener getListener() {
        return listener;
    }

    /**
     * Sets listener to evaluate a multidimensional function.
     *
     * @param listener Listener to evaluate a multidimensional function.
     */
    public void setListener(final MultiDimensionFunctionEvaluatorListener listener) {
        this.listener = listener;
    }

    /**
     * Returns point used as a reference to determine the function's input
     * parameters along a line.
     *
     * @return Point used as a reference to determine the function's input
     * parameters along a line
     */
    public double[] getPoint() {
        return point;
    }

    /**
     * Returns array indicating the direction of the line where the function is
     * evaluated.
     *
     * @return Array indicating the direction of the line where the function is
     * evaluated.
     */
    public double[] getDirection() {
        return direction;
    }

    /**
     * Sets point used as a reference to determine the function's input
     * parameters along a line and the direction of the line where the
     * function is evaluated.
     *
     * @param point     Point used as a reference to determine the function's input
     *                  parameters along a line.
     * @param direction Array indicating the direction of the line where the
     *                  function is evaluated.
     * @throws IllegalArgumentException Raised if point and direction don't have
     *                                  the same length.
     */
    public final void setPointAndDirection(final double[] point, final double[] direction) {
        if (point.length != direction.length) {
            throw new IllegalArgumentException();
        }

        this.point = point;
        this.direction = direction;
    }

    /**
     * Evaluates a function using current listener at a distance x from current
     * point using current direction
     *
     * @param x Distance from provided point using provided direction where
     *          function is being evaluated
     * @return Result of evaluating the function
     * @throws EvaluationException Thrown if function evaluation fails
     */
    public double evaluateAt(final double x) throws EvaluationException {
        for (int i = 0; i < point.length; i++) {
            p[i] = point[i] + x * direction[i];
        }
        return listener.evaluate(p);
    }
}
