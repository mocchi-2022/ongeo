#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <cstdio>

#include "ONGEO_NumericalIntegration.hpp"

// Reference : B. Guenter, R. Parent,
//   "Computing the arc length of parametric curves".
//    IEEE Computer Graphics and Applications, Vol 10, No.3, May. 1990, pp.72-78.

namespace{
	struct IntegralFunc{
		ON_BezierCurve bc;
		double operator()(double t) const{
			ON_3dPoint pt;
			ON_3dVector dt;
			bc.Ev1Der(t, pt, dt);
			return dt.Length();
		}
	};
}

double ONGEO_LengthBezierCurve_NumericalIntegration(const ON_BezierCurve &bc, const double *subdivided_prms, int cnt, double *estimated_deviation){
	ONGEO_NumericalIntegration_GaussKronrod<5> integ(0, 0);

	IntegralFunc func;
	func.bc = bc;

	double t_prv = *subdivided_prms;
	double sum_l = 0, sum_err = 0;
	for (int i = 1; i < cnt; ++i){
		double t_cur = subdivided_prms[i];
		integ.SetRange(t_prv, t_cur);
		double err, l = integ.Calc(func, err);
		sum_l += l;
		sum_err += err;
		t_prv = t_cur;
	}
	if (estimated_deviation) *estimated_deviation = sum_err;
	return sum_l;
}

double ONGEO_LengthBezierCurve_AdaptiveQuadrate(const ON_BezierCurve &bc, const double *subdivided_prms, int cnt, double tolerance, double *estimated_deviation){
	IntegralFunc func;
	func.bc = bc;

	double t_prv = *subdivided_prms;
	double sum_l = 0, sum_err = 0;
	for (int i = 1; i < cnt; ++i){
		double t_cur = subdivided_prms[i];
		double err, l = ONGEO_AdaptiveQuadrature<ONGEO_NumericalIntegration_GaussKronrod<5>, IntegralFunc>(func, t_prv, t_cur, tolerance * (t_cur - t_prv), err);
		sum_l += l;
		sum_err += err;
		t_prv = t_cur;
	}
	if (estimated_deviation) *estimated_deviation = sum_err;
	return sum_l;
}
