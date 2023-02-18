/**************************************************************************/
/*  minimization.cpp                                                      */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#include "minimization.h"

#include "core/error/error_macros.h"
#include "core/typedefs.h"
#include "math_funcs.h"

// Define the golden ratio
const real_t GOLDEN_RATIO = (1.0 + Math::sqrt(5.0)) / 2.0;
const real_t GOLDEN_SECTION = 1.0 / GOLDEN_RATIO;
const real_t TINY = 1.0e-20;
const real_t TOL = 1e-6;
const int MAX_ITERATIONS = 100;

void Minimization::bracketing_triplet_from_interval(void *data, real_function *f, real_t *ax, real_t *bx, real_t *cx, real_t *fa, real_t *fb, real_t *fc) {
	real_t a = *ax, b = *bx, c = *cx;
	real_t x, w, v, fx, fw, fv, u, fu;
	real_t d, e;
	real_t p, q, r, tol1, tol2;

	fx = f(data, b);
	fv = fw = fx;
	v = w = b - a;
	x = c;
	u = d = e = 0.0;
	fu = *fb;

	for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
		real_t xm = 0.5 * (a + c);
		tol1 = TOL * Math::abs(x) + TINY;
		tol2 = tol1 * 2.0;

		// If the difference between x and xm is within tolerance, we're done
		if (Math::abs(x - xm) <= (tol2 - 0.5 * (c - a))) {
			*ax = a;
			*bx = x;
			*cx = c;
			*fa = f(data, a);
			*fb = fx;
			*fc = f(data, c);
			return;
		}

		if (Math::abs(e) > tol1) {
			// Use Brent's method to find a new point
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);

			if (q > 0.0) {
				p = -p;
			} else {
				q = -q;
			}

			r = e;
			e = d;
			// Check if the new point is within the bracketing interval and satisfies the conditions for success
			if (Math::abs(p) < Math::abs(0.5 * q * r) && p > q * (a - x) && p < q * (c - x)) {
				d = p / q;
				u = x + d;

				// If the new point is too close to a or c, adjust it slightly
				if (u - a < tol2 || c - u < tol2) {
					d = tol1 * ((x < xm) ? 1.0 : -1.0);
				}
			} else {
				// If Brent's method fails, use golden section search
				if (x < xm) {
					e = c - x;
				} else {
					e = a - x;
				}

				d = GOLDEN_SECTION * e;
			}
		} else {
			// If the previous step was not successful, use golden section search
			if (x < xm) {
				e = c - x;
			} else {
				e = a - x;
			}

			d = GOLDEN_SECTION * e;
		}
		// Compute the new point and its function value
		if (Math::abs(d) >= tol1) {
			u = x + d;
		} else {
			u = x + ((d > 0.0) ? tol1 : -tol1);
		}

		fu = f(data, u);

		// Update the bracketing triplet and their function values
		if (fu <= fx) {
			if (u >= x) {
				a = x;
			} else {
				c = x;
			}

			v = w;
			fv = fw;
			w = x;
			fw = fx;
			x = u;
			fx = fu;
		} else {
			if (u >= x) {
				c = u;
			} else {
				a = u;
			}

			if (fu <= fw || w == x) {
				v = w;
				fv = fw;
				w = u;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

	// Maximum iterations reached without finding a bracket
	*ax = a;
	*bx = x;
	*cx = c;
	*fa = f(data, a);
	*fb = fx;
	*fc = f(data, c);
	ERR_FAIL_MSG("bracketing_triplet_from_interval failed to find a bracket");
}

real_t Minimization::get_local_minimum(void *data, real_function *f, real_function *dw, real_t ax, real_t bx, real_t cx, real_t tol, real_t *xmin) {
	real_t a = (ax < cx) ? ax : cx;
	real_t b = (ax > cx) ? ax : cx;
	real_t x = bx;
	real_t w = bx;
	real_t v = bx;
	real_t fx = f(data, x);
	real_t fx_prev = fx;
	real_t fw = fx;
	real_t fv = fx;
	real_t e = 0.0;
	real_t d = 0.0;

	// Pre-compute constants
	const real_t tol1 = tol * Math::abs(x) + TOL;
	const real_t tol2 = 2.0 * tol1;
	const real_t tol3 = tol1 / 10.0;

	for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
		const real_t xm = 0.5 * (a + b);

		// Check if we have converged
		if (Math::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
			*xmin = x;
			return fx;
		}

		if (b - a < tol3) {
			// Bracketing interval is small enough
			*xmin = x;
			return fx;
		}

		if (Math::abs(fx - fx_prev) < tol) {
			// Change in function value is small enough
			*xmin = x;
			return fx;
		}

		fx_prev = fx;
		real_t u = x;

		// Use the secant method to calculate the next step size
		if (Math::abs(e) > tol1) {
			const real_t r = (x - w) * (fx - fv);
			real_t q = (x - v) * (fx - fw);
			real_t p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) {
				p = -p;
			}
			q = Math::abs(q);
			const real_t etemp = e;
			e = d;
			if (Math::abs(p) >= Math::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
				// Perform golden section step with adaptive step size
				const real_t e_sign = signbit(x >= xm ? a - x : b - x);
				e = a - x + b - x;
				const real_t step_size = e * e_sign / (q + TINY);
				d = step_size;
			} else {
				// Use parabolic interpolation with adaptive step size
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) {
					const real_t e_sign = signbit(x >= xm ? a - x : b - x);
					e = a - x + b - x;
					const real_t step_size = e * e_sign / (q + TINY);
					d = step_size;
				}
			}
		} else {
			// Perform golden section step with adaptive step size
			const real_t e_sign = signbit(x >= xm ? a - x : b - x);
			e = a - x + b - x;
			const real_t step_size = e * e_sign / (GOLDEN_RATIO + TINY);
			d = step_size;
		}

		// Use the Newton-Raphson method to refine the estimate of the minimum
		const real_t dw_u = dw ? dw(data, u) : 0.0;
		real_t f_u = f(data, u);
		if (dw_u != 0.0) {
			const real_t u_newton = u - f_u / dw_u;
			if (u_newton > a && u_newton < b && Math::abs(u_newton - u) < tol1) {
				const real_t f_u_newton = f(data, u_newton);
				if (f_u_newton < f_u) {
					u = u_newton;
					f_u = f_u_newton;
				}
			}
		}

		// Check if we have converged
		if (Math::abs(x - u) < tol) {
			*xmin = x;
			return fx;
		}

		// Update bracketing interval
		if (f_u <= fx) {
			if (u >= x) {
				a = x;
			} else {
				b = x;
			}
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = f_u;
		} else {
			if (u < x) {
				a = u;
			} else {
				b = u;
			}
			if (f_u <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = f_u;
			} else if (f_u <= fv || v == x || v == w) {
				v = u;
				fv = f_u;
			}
		}
	}

	// Maximum number of iterations exceeded
	*xmin = x;
	ERR_FAIL_V_MSG(fx, "get_local_minimum failed to converge.");
}
