/*
 * FitNurbsToPoints
 * Copylight (C) 2020 mocchi
 * mocchi_2003@yahoo.co.jp
 * License: Boost ver.1
 */

#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>

// Reference : Piegl, L. and Tiller, W.
//   "The NURBS Book", Second edition.
//    Springer-Verlag Berlin, Heidelberg, 1997.

// pp.229 - 234 : Chapter 6.1, "Point Inversion and Projection for Curves and Surfaces"

bool ONGEO_NearestPointNurbsCurve_Newton(const ON_NurbsCurve &nc, const ON_3dPoint &pt_query, double &t, ON_Interval *range_, int iteration_max){
	ON_3dPoint pt;
	ON_3dVector der1, der2;
	ON_Interval range;
	if (range_) range = *range_;
	else range = nc.Domain();
	for(int i = 0; i < iteration_max; ++i){
		nc.Ev2Der(t, pt, der1, der2);
		double der1_len = der1.Length();
		ON_3dVector diff = pt - pt_query;
		double diff_len = diff.Length();
		double numerator = ON_DotProduct(der1, diff);

		bool cond_1 = diff_len <= ON_ZERO_TOLERANCE;
		bool cond_2 = std::abs(numerator / (der1_len * diff_len)) <= ON_ZERO_TOLERANCE;
		bool cond_4 = false;

		if (!cond_1 || !cond_2){
			double tn = t - numerator / (ON_DotProduct(der2, diff) + der1_len * der1_len);
			if (!nc.IsClosed()){
				if (tn < range.Min()) tn = range.Min();
				else if (tn > range.Max()) tn = range.Max();
			}else{
				if (tn < range.Min()) tn = range.Max() - (range.Min() - tn);
				else if (tn > range.Max()) tn = range.Min() - (tn - range.Max());
			}
			cond_4 = std::abs((tn - t) * der1_len) <= ON_ZERO_TOLERANCE;
			t = tn;
		}
		if (cond_1 || cond_2 || cond_4){
			return true;
		}
	}
	return false;
}

