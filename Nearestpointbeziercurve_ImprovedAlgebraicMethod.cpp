// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <stack>
#include <map>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <cstdio>

// Reference : Xiao-Diao Chen, Yin Zhou, Zhenyu Shu, Hua Su, Jean-Claude Paul,
//   "Improved Algebraic Algorithm On Point Projection For Bezier Curves".
//    IMSCCS '07 Proceedings of the Second International Multi-Symposiums on Computer and Computational Sciences
//    Pages 158-163

namespace {
inline void assign_coef3(std::vector<double> &coefs, int order, double *coef[3]){
	coefs.resize(order * 3);
	coef[0] = &coefs[0], coef[1] = coef[0] + order, coef[2] = coef[1] + order;
}
// Bni(t) = nCi * t^i * (1-t)^(n-i)
void ExpandBernsteinCoef(int n, int i, double *coefs){
	int *te = ONGEO_Polynomial_CalcPascalTriangle(n-i-1);
	double ni = static_cast<double>(ONGEO_Polynomial_CalcBinomialCoef(n-1, i));
	for (int j = 0; j < n-i; ++j){
		coefs[j] = ni * te[j];
		if (((n-i-j) & 1) == 0) coefs[j] *= -1;
	}
	for (int j = n-i; j < n; ++j) coefs[j] = 0;
}
}

double ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod(const ON_BezierCurve *bcs, int num_bcs, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest){
	double dist2_min = std::numeric_limits<double>::max();
	double zero2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;

//	std::vector<double> q_coefs;
	std::vector<double> v_coefs, dv_coefs, w_coef, dw_coef, 
		wdv_coefs, dwv_coefs, wp_coefs, wdv_dwv_coefs, wp_v_coefs, wdv_dwv_dot_wp_v_coefs, g_coef, sturm;
	std::vector<double> dv_dot_p_v_coefs;
	std::vector<double> dg_coefs;

	for (int k = 0; k < num_bcs; ++k){
		int dim = bcs[k].Dimension();
		int order = bcs[k].m_order;
#if 0
		ON_PolynomialCurve pc(bcs[k]);
		double p0[4], p1[4];
		pc.Evaluate(0, 0, 4, p0);
		pc.Evaluate(1, 0, 4, p1);
#endif
//		double *q_coef[3];
		double *v_coef[3], *dv_coef[3];
//		assign_coef3(q_coefs, order, q_coef);
		assign_coef3(v_coefs, order, v_coef);
		assign_coef3(dv_coefs, order-1, dv_coef);

		int order_g;
		if (bcs[k].m_is_rat){
			w_coef.resize(order);
			w_coef.assign(w_coef.size(), 0);
#if 0
			// ON_PolynomialCurveから、 q(u)、V(u)それぞれの係数列を読み込み
			for (int j = 0, jinv = order-1; j < order; ++j, --jinv){
				double w = w_coef[j] = pc.m_cv[jinv][3];
				for (int i = 0; i < dim; ++i){
					v_coef[i][j] = pc.m_cv[jinv][i];
//					v_coef[i][j] = q_coef[i][j] / w;
				}
			}
#else
			std::vector<double> coefs(order);
			v_coefs.assign(v_coefs.size(), 0);
			for (int j = 0; j < order; ++j){
				ExpandBernsteinCoef(order, j, &coefs[0]);
				double cv[4];
				bcs[k].GetCV(j, ON::homogeneous_rational, cv);
				for (int i = 0; i < order; ++i){
					for (int h = 0; h < dim; ++h){
						v_coef[h][i] += coefs[i] * cv[h];
					}
					w_coef[i] += coefs[i] * cv[dim];
				}
			}
#endif
			// V(u) から V'(u)を計算
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Differential(v_coef[i], order, dv_coef[i]);
			}
			dw_coef.resize(order-1);
			ONGEO_Polynomial_Differential(&w_coef[0], order, &dw_coef[0]);

			// w(u)V'(u) と w'(u)V'(u)を計算
			double *wdv_coef[3], *dwv_coef[3];
			int order2_2 = order*2-2;
			assign_coef3(wdv_coefs, order2_2, wdv_coef);
			assign_coef3(dwv_coefs, order2_2, dwv_coef);
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Multiply(&w_coef[0], order, dv_coef[i], order-1, wdv_coef[i]);
				ONGEO_Polynomial_Multiply(&dw_coef[0], order-1, v_coef[i], order, dwv_coef[i]);
			}

			// w(u)V'(u) - w'(u)V(u)を計算
			double *wdv_dwv_coef[3];
			assign_coef3(wdv_dwv_coefs, order2_2, wdv_dwv_coef);
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Subtract(wdv_coef[i], order2_2, dwv_coef[i], order2_2, wdv_dwv_coef[i], order2_2);
			}

			// w(u)p - V(u)を計算
			double *wp_coef[3], *wp_v_coef[3];
			assign_coef3(wp_coefs, order, wp_coef);
			assign_coef3(wp_v_coefs, order, wp_v_coef);
			for (int j = 0; j < dim; ++j){
				for (int i = 0; i < order; ++i) wp_coef[j][i] = w_coef[i] * pt_query[j];
				ONGEO_Polynomial_Subtract(wp_coef[j], order, v_coef[j], order, wp_v_coef[j], order);
			}

			// x, y, z各要素について、 (w(u)V'(u) - w'(u)V(u)) x (w(u)p - V(u))を計算
			double *wdv_dwv_dot_wp_v_coef[3];
			int order3_3 = order*3-3;
			assign_coef3(wdv_dwv_dot_wp_v_coefs, order3_3, wdv_dwv_dot_wp_v_coef);
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Multiply(wdv_dwv_coef[i], order2_2, wp_v_coef[i], order, wdv_dwv_dot_wp_v_coef[i]);
			}

			// 上記の係数列のx, y, z各要素について総和を取り(内積)、 g(u)とする。 
			g_coef.resize(order3_3);
			if (dim == 1){
				for (int h = 0; h < order3_3; ++h){
					g_coef[h] = wdv_dwv_dot_wp_v_coef[0][h];
				}
			}
			if (dim >= 2){
				ONGEO_Polynomial_Add(wdv_dwv_dot_wp_v_coef[0], order3_3, wdv_dwv_dot_wp_v_coef[1], order3_3, &g_coef[0], order3_3);
			}
			if (dim == 3){
				ONGEO_Polynomial_Add(wdv_dwv_dot_wp_v_coef[2], order3_3, &g_coef[0], order3_3, &g_coef[0], order3_3);
			}
			order_g = order3_3;
		}else{
			std::vector<double> coefs(order);
			v_coefs.assign(v_coefs.size(), 0);
			for (int j = 0; j < order; ++j){
				ExpandBernsteinCoef(order, j, &coefs[0]);
				double cv[3];
				bcs[k].GetCV(j, ON::not_rational, cv);
				for (int i = 0; i < order; ++i){
					for (int h = 0; h < dim; ++h){
						v_coef[h][i] += coefs[i] * cv[h];
					}
				}
			}

			// V(u) から V'(u)を計算
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Differential(v_coef[i], order, dv_coef[i]);
			}

			// p - V(u)を計算
			double *p_v_coef[3];
			assign_coef3(wp_v_coefs, order, p_v_coef);
			for (int j = 0; j < dim; ++j){
				ONGEO_Polynomial_Subtract(&static_cast<const double *>(pt_query)[j], 1, v_coef[j], order, p_v_coef[j], order);
			}

			// x, y, z各要素について、 V'(u) x (p - V(u))を計算
			double *dv_dot_p_v_coef[3];
			int order2_2 = order*2-2;
			assign_coef3(wdv_dwv_dot_wp_v_coefs, order2_2, dv_dot_p_v_coef);
			for (int i = 0; i < dim; ++i){
				ONGEO_Polynomial_Multiply(dv_coef[i], order-1, p_v_coef[i], order, dv_dot_p_v_coef[i]);
			}

			// 上記の係数列のx, y, z各要素について総和を取り(内積)、 g(u)とする。 
			g_coef.resize(order2_2);
			if (dim == 1){
				for (int h = 0; h < order2_2; ++h){
					g_coef[h] = dv_dot_p_v_coef[0][h];
				}
			}
			if (dim >= 2){
				ONGEO_Polynomial_Add(dv_dot_p_v_coef[0], order2_2, dv_dot_p_v_coef[1], order2_2, &g_coef[0], order2_2);
			}
			if (dim == 3){
				ONGEO_Polynomial_Add(dv_dot_p_v_coef[2], order2_2, &g_coef[0], order2_2, &g_coef[0], order2_2);
			}
			order_g = order2_2;
		}

		dg_coefs.resize(order_g-1);
		ONGEO_Polynomial_Differential(&g_coef[0], order_g, &dg_coefs[0]);

		// g(u)の Sturm Sequence を求める。 
		int num_strum = order_g * order_g;
		sturm.resize(num_strum);
		ONGEO_Polynomial_CreateSturmSequence(&g_coef[0], order_g, &sturm[0]);

		struct Interval{
			double a, b;
			int cnta, cntb;
			Interval(double a_, int cnta_, double b_, int cntb_) : a(a_), cnta(cnta_), b(b_), cntb(cntb_){
			}
		};
		std::stack<Interval> ints;
		ints.push(Interval(
			0.0, ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(0.0, &sturm[0], order_g),
			1.0, ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(1.0, &sturm[0], order_g)
		));
		std::vector<std::pair<double, double> > candidate;
		std::vector<int> candidate_cnt_first;
		while(ints.size()){
			Interval interval = ints.top();
			ints.pop();
			if (interval.cnta == interval.cntb) continue;
			if (interval.cnta - interval.cntb >= 2){
				double m = (interval.a + interval.b) * 0.5;
				double dd = interval.b - interval.a;
				int cntm = ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(m, &sturm[0], order_g);
				if (dd < DBL_EPSILON || interval.a == m || interval.b == m){
					candidate.push_back(std::make_pair(m, m));
					candidate_cnt_first.push_back(cntm);
				}else{
					ints.push(Interval(interval.a, interval.cnta, m, cntm));
					ints.push(Interval(m, cntm, interval.b, interval.cntb));
				}
			}else if (interval.cnta - interval.cntb == 1){
				double ag = ONGEO_Polynomial_Evaluate(interval.a, &g_coef[0], order_g);
				double bg = ONGEO_Polynomial_Evaluate(interval.b, &g_coef[0], order_g);
				if (ag == 0){
					candidate.push_back(std::make_pair(interval.a, interval.a));
					candidate_cnt_first.push_back(interval.cnta);
				}
				if (bg == 0){
					candidate.push_back(std::make_pair(interval.b, interval.b));
					candidate_cnt_first.push_back(interval.cntb);
				}
				if (ag * bg < 0){
					candidate.push_back(std::make_pair(interval.a, interval.b));
					candidate_cnt_first.push_back(interval.cnta);
				}
			}
		}
		for (size_t i = 0; i < candidate.size(); ++i){
			std::pair<double, double> crange = candidate[i];
			double dist2_prev = std::numeric_limits<double>::max();
			double tc, gcv;
			double gv[2] = {
				ONGEO_Polynomial_Evaluate(crange.first, &g_coef[0], order_g),
				ONGEO_Polynomial_Evaluate(crange.second, &g_coef[0], order_g)
			};

			ON_3dPoint ptc;
			if (crange.first < crange.second){
				// 2分法
#if 0
				int cnta = candidate_cnt_first[i];
				for(;;){
					tc = (crange.first + crange.second) * 0.5;
					gcv = ONGEO_Polynomial_Evaluate(tc, &g_coef[0], order_g);
					ptc = bcs[k].PointAt(tc);
					double dist2 = (ptc-pt_query).LengthSquared();
					double subd2 = std::abs(dist2_prev - dist2);
					dist2_prev = dist2;
					if (subd2 < ON_ZERO_TOLERANCE) break;
					if (gv[0] * gcv < 0){
						gv[1] = gcv;
						crange.second = tc;
					}else{
						gv[0] = gcv;
						crange.first = tc;
					}
//					int cntm = ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(tc, &sturm[0], order3_3);
//					if (cnta == cntm) crange.first = tc;
//					else crange.second = tc;
				}
#elif 1
				// Brent法
				tc = ONGEO_Polynomial_FindRootByBrentMethod(crange.first, crange.second, &g_coef[0], order_g);
				gcv = ONGEO_Polynomial_Evaluate(tc, &g_coef[0], order_g);
				ptc = bcs[k].PointAt(tc);
				dist2_prev = (ptc-pt_query).LengthSquared();
#elif 0
				// ニュートン法
				for(;;){
					double dt = ONGEO_Polynomial_Evaluate(tc, &g_coef[0], order_g) / ONGEO_Polynomial_Evaluate(tc, &dg_coefs[0], order_g-1);
					if (std::abs(dt) < ON_ZERO_TOLERANCE) break;
					tc -= dt;
				}
				ptc = bcs[k].PointAt(tc);
				dist2_prev = (ptc-pt_query).LengthSquared();
#endif
			}else{
				tc = crange.first;
				ptc = bcs[k].PointAt(tc);
				dist2_prev = (ptc-pt_query).LengthSquared();
			}
			// 今までの最短距離よりも短い場合は更新
			if (dist2_min > dist2_prev){
				dist2_min = dist2_prev;
				t = tc;
				pt_nearest = ptc;
				bc_nearest = &bcs[k];
			}
		}
		for (double tt = 0; tt <= 1.0; tt += 1.0){
			double dist2;
			ON_3dPoint pt = bcs[k].PointAt(tt);
			dist2 = (pt-pt_query).LengthSquared();
			if (dist2_min > dist2){
				dist2_min = dist2;
				t = tt;
				pt_nearest = pt;
				bc_nearest = &bcs[k];
			}
		}
	}
	return std::sqrt(dist2_min);
}
