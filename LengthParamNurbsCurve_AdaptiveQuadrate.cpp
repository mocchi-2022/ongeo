#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "ONGEO_NumericalIntegration.hpp"
#include "ONGEO_RootFinding.hpp"

// Reference : B. Guenter, R. Parent,
//   "Computing the arc length of parametric curves".
//    IEEE Computer Graphics and Applications, Vol 10, No.3, May. 1990, pp.72-78.

namespace{
	struct IntegralFunc{
		ON_NurbsCurve nc;
		double operator()(double t) const{
			ON_3dPoint pt;
			ON_3dVector dt;
			nc.Ev1Der(t, pt, dt);
			return dt.Length();
		}
	};
}

ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate::ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate(const ON_NurbsCurve &nc, double tolerance){
	IntegralFunc func;
	func.nc = nc;
	this->nc = nc;
	this->tolerance = tolerance;

	if (!nc.IsValid()) return;

	// Nurbs 曲線を区間分割
	{
		ON_SimpleArray<double> span_vect(nc.SpanCount()+1);
		if (span_vect.Capacity() == 1) return;

		span_vect.SetCount(span_vect.Capacity());
		nc.GetSpanVector(span_vect.First());

		double cpl = nc.ControlPolygonLength();

		nc_prms.AppendNew();
		nc.GetDomain(&nc_prms[0], 0);
		
		for (int j = 0; j <= nc.CVCount() - nc.Order(); ++j){
			ON_BezierCurve bc;
			if (!nc.ConvertSpanToBezier(j, bc)) continue;
			ON_Interval p_int(*nc_prms.Last(), nc.Knot(j+nc.Order()-1));

			ON_SimpleArray<double> tess_prm;
			ONGEO_TessellateBezierCurve_QuasiInterpolating(bc, cpl / 1000.0, tess_prm);

			for (int i = 1; i < tess_prm.Count(); ++i){
				nc_prms.Append(p_int.ParameterAt(tess_prm[i]));
			}
		}
	}

	double prm_range = *nc_prms.Last() - *nc_prms.First();
	
	// 適応求積法で区間ごとの曲線長を算出し、区間ごとの距離の積算値を nc_lengths に保持
	double t_prv = nc_prms[0];
	double sum_l = 0, sum_err = 0;
	nc_lengths.Append(0);
	for (int i = 1; i < nc_prms.Count(); ++i){
		double t_cur = nc_prms[i];
		double err, l = ONGEO_AdaptiveQuadrature<ONGEO_NumericalIntegration_GaussKronrod<5>, IntegralFunc>(func, t_prv, t_cur, tolerance * (t_cur - t_prv) / prm_range, err);
		sum_l += l;
		sum_err += err;
		t_prv = t_cur;

		nc_lengths.Append(sum_l);
	}
	estimated_deviation = sum_err;
}

double ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate::Length() const{
	if (nc_lengths.Count() == 0) return 0;
	return *nc_lengths.Last();
}

double ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate::ParamToLength(double nc_prm) const{
	IntegralFunc func;
	func.nc = nc;
	double prm_range = *nc_prms.Last() - *nc_prms.First();

	// 指定されたパラメータ値を超えない最も大きなパラメータ値を保持した配列から探し、
	// そのインデックスを取得する。
	const double *p = std::upper_bound(nc_prms.First(), nc_prms.Last()+1, nc_prm);
	int p_idx = p - nc_prms.First();
	if (p_idx == 0) return -1;

	// 指定されたパラメータ値と同じ値が配列に入っている場合は、対応する距離そのままを返す。
	if (*(p-1) == nc_prm) return nc_lengths[p_idx-1];

	// 配列から探したパラメータ値と指定されたパラメータ値で決まる区間の曲線長を算出し、足して返す。
	double err;
	double l = ONGEO_AdaptiveQuadrature<ONGEO_NumericalIntegration_GaussKronrod<5>, IntegralFunc>(func, *(p-1), nc_prm, tolerance * (nc_prm - *(p-1)) / prm_range, err);

	return nc_lengths[p_idx-1] + l;
}
double ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate::LengthToParam(double nc_length) const{
	IntegralFunc func;
	func.nc = nc;
	double prm_range = *nc_prms.Last() - *nc_prms.First();

	const double *l = std::upper_bound(nc_lengths.First(), nc_lengths.Last()+1, nc_length);
	int l_idx = l - nc_lengths.First();
	if (l_idx == 0) return -1;

	// 指定された距離と同じ値が配列に入っている場合は、対応するパラメータそのままを返す。
	if (*(l-1) == nc_length) return nc_prms[l_idx-1];


	ON_Interval p_range(nc_prms[l_idx-1], nc_prms[l_idx]);
	ON_Interval l_range(nc_lengths[l_idx-1], nc_lengths[l_idx]);
	double diff_l = nc_length - *(l-1);
#if 0
	// 配列から見つけた距離から区間を特定し、区間内の距離からパラメータをはさみうち法で算出し、返す。
	for(;;){
		double prm_estimate = p_range.ParameterAt(l_range.NormalizedParameterAt(nc_length));
		double err;
		double itg_l = ONGEO_AdaptiveQuadrature<ONGEO_NumericalIntegration_GaussKronrod<5>, IntegralFunc>(func, p_range.Min(), prm_estimate, tolerance * (prm_estimate - p_range.Min()) / prm_range, err);
		if (diff_l >= itg_l){
			if (diff_l - itg_l <= tolerance) return prm_estimate;
			l_range.m_t[0] += itg_l;
			p_range.m_t[0] = prm_estimate;
			diff_l = nc_length - l_range.m_t[0];
		}else{
			if (itg_l - diff_l < tolerance) return prm_estimate;
			l_range.m_t[1] = l_range.m_t[0] + itg_l;
			p_range.m_t[1] = prm_estimate;
		}
	}
#else
	return ONGEO_FindRootByBrentMethod(p_range.m_t[0], p_range.m_t[1], 
		[&] (double t) {
			double err, itg_l = 0;
			if (std::abs(t - p_range.m_t[0]) > tolerance){
				itg_l = ONGEO_AdaptiveQuadrature<ONGEO_NumericalIntegration_GaussKronrod<5>, IntegralFunc>(func, p_range.m_t[0], t, tolerance * (t - p_range.m_t[0]) / prm_range, err);
			}
			return l_range.m_t[0] +itg_l - nc_length;
		}
		, tolerance);
#endif
}

