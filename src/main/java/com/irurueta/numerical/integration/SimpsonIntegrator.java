/*
 * Copyright (C) 2023 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.integration;

import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Base integrator for implementations based on Simpson's method.
 * Simpson's method is an optimization of Trapezoidal quadrature integrator.
 * Implementations of this class will in general be more efficient than
 * Trapezoidal quadrature integrators (i.e., require fewer function evaluations) when the function
 * to be integrated has a finite fourth derivative (i.e., a continuous third derivative).
 *
 * <p>
 * The extrapolation step combines the current quadrature estimate {@code st} and the previous one {@code ost} as
 * {@code s = (ratio * st - ost) / (ratio - 1)}. The ratio is the factor by which the quadrature error is reduced
 * on each refinement: {@code 2^2 = 4} for quadratures halving their step (e.g., trapezoidal), which yields
 * {@code (4 * st - ost) / 3}, and {@code 3^2 = 9} for quadratures tripling their step (e.g., mid-point), which
 * yields {@code (9 * st - ost) / 8}.
 * See Numerical Recipes, 3rd ed., section 4.2 (p. 165) and section 4.4 (p. 169).
 *
 * @param <T> a quadrature.
 */
public abstract class SimpsonIntegrator<T extends Quadrature> extends Integrator {

    /**
     * Default accuracy.
     */
    public static final double EPS = 1e-10;

    /**
     * Extrapolation ratio for quadratures that halve their step on each refinement (2^2).
     * Leads to extrapolation {@code (4 * st - ost) / 3}.
     */
    protected static final double TRAPEZOIDAL_RATIO = 4.0;

    /**
     * Extrapolation ratio for quadratures that triple their step on each refinement (3^2).
     * Leads to extrapolation {@code (9 * st - ost) / 8}.
     */
    protected static final double TRIPLING_RATIO = 9.0;

    /**
     * Default maximum number of allowed steps.
     */
    protected static final int DEFAULT_JMAX = 20;

    /**
     * Maximum number of allowed steps for quadratures tripling their step on each refinement.
     */
    protected static final int TRIPLING_JMAX = 14;

    /**
     * Minimum required number of steps.
     */
    private static final int JMIN = 5;

    /**
     * Quadrature used for integration.
     */
    private final T q;

    /**
     * Required accuracy.
     */
    private final double eps;

    /**
     * Extrapolation ratio.
     */
    private final double ratio;

    /**
     * Maximum number of allowed steps.
     */
    private final int jmax;

    /**
     * Constructor using halving step extrapolation ratio and default maximum number of steps.
     *
     * @param q   Quadrature used for integration.
     * @param eps Required accuracy.
     */
    protected SimpsonIntegrator(final T q, final double eps) {
        this(q, eps, TRAPEZOIDAL_RATIO, DEFAULT_JMAX);
    }

    /**
     * Constructor.
     *
     * @param q     Quadrature used for integration.
     * @param eps   Required accuracy.
     * @param ratio Extrapolation ratio (4.0 for quadratures halving their step, 9.0 for tripling their step).
     *              Must be greater than 1; otherwise the extrapolation divides by zero (or flips sign) and
     *              integration never converges.
     * @param jmax  Maximum number of allowed steps. Must be greater than 6, since convergence is only checked
     *              after the minimum number of steps; otherwise integration always fails.
     */
    protected SimpsonIntegrator(final T q, final double eps, final double ratio, final int jmax) {
        this.q = q;
        this.eps = eps;
        this.ratio = ratio;
        this.jmax = jmax;
    }

    /**
     * Integrates function between provided lower and upper limits.
     *
     * @return result of integration.
     * @throws IntegrationException if integration fails for numerical reasons.
     */
    @Override
    public double integrate() throws IntegrationException {
        try {
            double s;
            double st;
            var ost = 0.0;
            var os = 0.0;
            for (var j = 0; j < jmax; j++) {
                st = q.next();
                s = (ratio * st - ost) / (ratio - 1.0);
                if (j > JMIN && (Math.abs(s - os) < eps * Math.abs(os) || (s == 0.0 && os == 0.0))) {
                    return s;
                }
                os = s;
                ost = st;
            }
        } catch (final EvaluationException e) {
            throw new IntegrationException(e);
        }

        // Too many steps
        throw new IntegrationException();
    }

    /**
     * Gets type of integrator.
     *
     * @return type of integrator.
     */
    @Override
    public IntegratorType getIntegratorType() {
        return IntegratorType.SIMPSON;
    }

    /**
     * Creates an integrator using Simpson's method.
     *
     * @param a              Lower limit of integration.
     * @param b              Upper limit of integration.
     * @param listener       listener to evaluate a single dimension function at required points.
     * @param eps            required accuracy.
     * @param quadratureType quadrature type.
     * @return created integrator.
     * @throws IllegalArgumentException if provided quadrature type is not supported.
     */
    public static SimpsonIntegrator<Quadrature> create(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener, final double eps,
            final QuadratureType quadratureType) {
        return switch (quadratureType) {
            case TRAPEZOIDAL -> cast(new SimpsonTrapezoidalQuadratureIntegrator(a, b, listener, eps));
            case MID_POINT -> cast(new SimpsonMidPointQuadratureIntegrator(a, b, listener, eps));
            case INFINITY_MID_POINT -> cast(new SimpsonInfinityMidPointQuadratureIntegrator(a, b, listener, eps));
            case LOWER_SQUARE_ROOT_MID_POINT ->
                    cast(new SimpsonLowerSquareRootMidPointQuadratureIntegrator(a, b, listener, eps));
            case UPPER_SQUARE_ROOT_MID_POINT ->
                    cast(new SimpsonUpperSquareRootMidPointQuadratureIntegrator(a, b, listener, eps));
            case DOUBLE_EXPONENTIAL_RULE ->
                    cast(new SimpsonDoubleExponentialRuleQuadratureIntegrator(a, b, listener, eps));
            default -> throw new IllegalArgumentException();
        };
    }

    /**
     * Creates an integrator using Simpson's method and having default accuracy.
     *
     * @param a              Lower limit of integration.
     * @param b              Upper limit of integration.
     * @param listener       listener to evaluate a single dimension function at required points.
     * @param quadratureType quadrature type.
     * @return created integrator.
     * @throws IllegalArgumentException if provided quadrature type is not supported.
     */
    public static SimpsonIntegrator<Quadrature> create(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener,
            final QuadratureType quadratureType) {
        return switch (quadratureType) {
            case TRAPEZOIDAL -> cast(new SimpsonTrapezoidalQuadratureIntegrator(a, b, listener));
            case MID_POINT -> cast(new SimpsonMidPointQuadratureIntegrator(a, b, listener));
            case INFINITY_MID_POINT -> cast(new SimpsonInfinityMidPointQuadratureIntegrator(a, b, listener));
            case LOWER_SQUARE_ROOT_MID_POINT ->
                    cast(new SimpsonLowerSquareRootMidPointQuadratureIntegrator(a, b, listener));
            case UPPER_SQUARE_ROOT_MID_POINT ->
                    cast(new SimpsonUpperSquareRootMidPointQuadratureIntegrator(a, b, listener));
            case DOUBLE_EXPONENTIAL_RULE -> cast(new SimpsonDoubleExponentialRuleQuadratureIntegrator(a, b, listener));
            default -> throw new IllegalArgumentException();
        };
    }

    /**
     * Creates an integrator using Simpson's method and default quadrature type.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     * @param eps      required accuracy.
     * @return created integrator.
     */
    public static SimpsonIntegrator<Quadrature> create(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener, final double eps) {
        return create(a, b, listener, eps, DEFAULT_QUADRATURE_TYPE);
    }

    /**
     * Creates an integrator using Simpson's method and having default accuracy and quadrature type.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     * @return created integrator.
     */
    public static SimpsonIntegrator<Quadrature> create(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener) {
        return create(a, b, listener, DEFAULT_QUADRATURE_TYPE);
    }

    /**
     * Casts integrator to a quadrature integrator without wildcard parameter.
     *
     * @param integrator integrator to be cast.
     * @return cast integrator.
     */
    private static SimpsonIntegrator<Quadrature> cast(final SimpsonIntegrator<?> integrator) {
        //noinspection unchecked
        return (SimpsonIntegrator<Quadrature>) integrator;
    }
}
