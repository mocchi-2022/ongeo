// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <stack>
#include <map>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>

//#define DEBUGOUT

#include "Profile.h"

#ifdef DEBUGOUT
#define DPRINTF(X, ...) std::fprintf(X, __VA_ARGS__), std::fflush(fp);
#else
#define DPRINTF(X, ...) 
#endif

namespace {
void ProjectToRayPlane(const ON_Plane &rnPlane, const ON_BezierSurface &bez, ON_SimpleArray<ON_2dPoint> &cp2d){
	PROF("ProjectToRayPlane(polynomial)");
	cp2d.SetCapacity(bez.Order(0)*bez.Order(1));
	cp2d.SetCount(cp2d.Capacity());
	for (int j = 0, jj = 0; j < bez.Order(1); ++j, jj += bez.Order(0)){
		for (int i = 0; i < bez.Order(0); ++i){
			ON_3dPoint pt;
			bez.GetCV(i, j, pt);
			rnPlane.ClosestPointTo(pt, &cp2d[jj+i][0], &cp2d[jj+i][1]);
		}
	}
}

void ProjectToRayPlane(const ON_Plane &rnPlane, const ON_BezierSurface &bez, ON_SimpleArray<ON_3dPoint> &cp3d){
	PROF("ProjectToRayPlane(rational)");
	cp3d.SetCapacity(bez.Order(0)*bez.Order(1));
	cp3d.SetCount(cp3d.Capacity());
	for (int j = 0, jj = 0; j < bez.Order(1); ++j, jj += bez.Order(0)){
		for (int i = 0; i < bez.Order(0); ++i){
			ON_3dPoint pt;
			bez.GetCV(i, j, pt);
			rnPlane.ClosestPointTo(pt, &cp3d[jj+i][0], &cp3d[jj+i][1]);
			cp3d[jj+i][2] = bez.Weight(i, j);
		}
	}
}

template <typename T> bool IsInside(ON_SimpleArray<T> &cppln, int order_u, int order_v){
	PROF("IsInside");
	double minv[2], maxv[2];
	minv[0] = maxv[0] = cppln[0].x;
	minv[1] = maxv[1] = cppln[0].y;
	for (int i = 1; i < cppln.Count(); ++i){
		maxv[0] = std::max(cppln[i].x, maxv[0]);
		minv[0] = std::min(cppln[i].x, minv[0]);
		maxv[1] = std::max(cppln[i].y, maxv[0]);
		minv[1] = std::min(cppln[i].y, minv[0]);
		if ((minv[0] * maxv[0]) < 0 && (minv[1] * maxv[1]) < 0) return true;
	}
	return false;
}

template <typename T> bool IsMonotonic(ON_SimpleArray<T> &cppln, int order_u, int order_v){
	PROF("IsMonotonic");
	bool sign[2][2];
	sign[0][0] = ((cppln[1].x - cppln[0].x) > 0) ? true : false;
	sign[0][1] = ((cppln[1].y - cppln[0].y) > 0) ? true : false;
	sign[1][0] = ((cppln[order_u].x - cppln[0].x) > 0) ? true : false;
	sign[1][1] = ((cppln[order_u].y - cppln[0].y) > 0) ? true : false;
	for (int j = 1, jj = order_u; j < order_v - 1; ++j, jj += order_u){
		for (int i = 1; i < order_u - 1; ++i){
			if (sign[0][0] != (((cppln[jj+i+1].x - cppln[jj+i].x) > 0) ? true : false)) return false;
			if (sign[0][1] != (((cppln[jj+i+1].y - cppln[jj+i].y) > 0) ? true : false)) return false;
			if (sign[1][0] != (((cppln[jj+i+order_u].x - cppln[jj+i].x) > 0) ? true : false)) return false;
			if (sign[1][1] != (((cppln[jj+i+order_u].y - cppln[jj+i].y) > 0) ? true : false)) return false;
		}
	}
	return true;
}

bool GuessNewton(ON_Plane &pln, ON_BezierSurface &bez, ON_SimpleArray<ON_2dPoint> &cp, double uv[2]){
	PROF("GuessNewton(polynomial)");
	ON_BezierSurface bez_t(3, false, bez.Order(0), bez.Order(1));
	for (int j = 0, jj = 0; j < bez.Order(1); ++j, jj += bez.Order(0)){
		for (int i = 0; i < bez.Order(0); ++i){
			ON_3dPoint pt;
			bez.GetCV(i, j, pt);
			ON_3dPoint ptt(cp[i+jj].x, cp[i+jj].y, pln.DistanceTo(pt));
			bez_t.SetCV(i, j, ptt);
		}
	}
	ON_3dPoint p00, p0n, pn0;
	bez_t.GetCV(0, 0, p00);
	bez_t.GetCV(bez.Order(0), 0, p0n);
	bez_t.GetCV(0, bez.Order(1), pn0);
	ON_3dVector vec[2] = {p0n-p00, pn0-p00};
	double length[2] = {vec[0].Length(), vec[1].Length()};
	vec[0] /= length[0], vec[1] /= length[1];

	double ztole2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;

	uv[0] = -ON_DotProduct(ON_3dVector(p00), vec[0]);
	uv[1] = -ON_DotProduct(ON_3dVector(p00), vec[1]);
	for(;;){
		if (uv[0] < 0 || uv[1] < 0 || uv[0] > 1 || uv[1] > 1) return false;
		ON_3dVector pt = bez_t.PointAt(uv[0], uv[1]);
		if (pt.LengthSquared() < ztole2) return true;
		uv[0] -= ON_DotProduct(pt, vec[0]) / length[0];
		uv[1] -= ON_DotProduct(pt, vec[1]) / length[1];
	}
	return false;
}

bool GuessNewton(ON_Plane &pln, ON_BezierSurface &bez, ON_SimpleArray<ON_3dPoint> &cp, double uv[2]){
	PROF("GuessNewton(rational)");
	ON_BezierSurface bez_t(3, true, bez.Order(0), bez.Order(1));
	for (int j = 0, jj = 0; j < bez.Order(1); ++j, jj += bez.Order(0)){
		for (int i = 0; i < bez.Order(0); ++i){
			ON_4dPoint pt;
			bez.GetCV(i, j, pt);
			ON_4dPoint ptt(cp[i+jj].x, cp[i+jj].y, pln.DistanceTo(pt), pt.w);
			bez_t.SetCV(i, j, ptt);
		}
	}
	ON_3dPoint p00, p0n, pn0;
	bez_t.GetCV(0, 0, p00);
	bez_t.GetCV(bez.Order(0)-1, 0, p0n);
	bez_t.GetCV(0, bez.Order(1)-1, pn0);
	p00.z = p0n.z = pn0.z = 0;
	ON_3dVector vec[2] = {p0n-p00, pn0-p00};
	double length[2] = {vec[0].Length(), vec[1].Length()};
	vec[0] /= length[0], vec[1] /= length[1];

	double ztole2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;

	uv[0] = -ON_DotProduct(ON_3dVector(p00), vec[0]) / length[0];
	uv[1] = -ON_DotProduct(ON_3dVector(p00), vec[1]) / length[1];
	for(;;){
		if (uv[0] < 0 || uv[1] < 0 || uv[0] > 1 || uv[1] > 1) return false;
		ON_3dVector pt = bez_t.PointAt(uv[0], uv[1]);
		pt.z = 0;
		if (pt.LengthSquared() < ztole2) return true;
		uv[0] -= ON_DotProduct(pt, vec[0]) / length[0];
		uv[1] -= ON_DotProduct(pt, vec[1]) / length[1];
	}
	return false;
}

template <typename T> struct type_traits{
};
template <> struct type_traits<ON_3dPoint>{
	typedef ON_2dPoint BelowPoint_T;
};
template <> struct type_traits<ON_4dPoint>{
	typedef ON_3dPoint BelowPoint_T;
};
}

template <typename T> int IntersectRayBezier_Subdiv_Newton_(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance){
	PROF("IntersectRayBezier_Subdiv_Newton_");

	std::stack<std::pair<ON_Interval, ON_Interval> > stack_ints;
	stack_ints.push(std::make_pair(ON_Interval(0,1), ON_Interval(0,1)));
	ON_Plane rayplane(ray.m_P, ray.m_V);
	std::vector<double> du(bez.m_order[0]), dv(bez.m_order[1]);
	std::vector<int> jis(bez.m_order[1]);
	std::vector<double> err2_vertices(bez.m_order[0]*bez.m_order[1]);
	ON_SimpleArray<type_traits<T>::BelowPoint_T> cppln;
	struct Result{
		double t, u, v;
		ON_3dPoint ptsrf;
		ON_3dPoint ptlin;
	};
	std::vector<Result> results;
	for (int j = 0, ji = 0; j < bez.m_order[1]; ++j, ji += bez.m_order[0]) jis[j] = ji;
	DPRINTF(fp, "\n\n\n'==== start ==== %d %d %s\n", bez.Order(0), bez.Order(1), bez.IsRational() ? "Rational" : "Polynomial");
	while(stack_ints.size()){
		ON_BezierSurface bez_t(bez);
		ON_Interval intu = stack_ints.top().first;
		ON_Interval intv = stack_ints.top().second;
		stack_ints.pop();
		bez_t.Trim(0, intu);
		bez_t.Trim(1, intv);

		{
			ON_3dPoint p0, pn;
			bez_t.GetCV(0, 0, p0);
			bez_t.GetCV(bez_t.Order(0)-1, 0, pn);
			ON_3dVector yaxis1 = ON_CrossProduct(pn - p0, ray.m_V);
			ON_3dVector xaxis1 = ON_CrossProduct(ray.m_V, yaxis1);
			xaxis1.Unitize(), yaxis1.Unitize();
			ON_Plane rayplane1(ray.m_P, xaxis1, yaxis1);
			ProjectToRayPlane(rayplane1, bez_t, cppln);
			if (!IsInside(cppln, bez_t.Order(0), bez_t.Order(1))) continue;
		}

		ON_3dPoint p0, pn;
		bez_t.GetCV(0, 0, p0);
		bez_t.GetCV(bez_t.Order(0)-1, bez_t.Order(1)-1, pn);
		ON_3dVector yaxis2 = ON_CrossProduct(pn - p0, ray.m_V);
		ON_3dVector xaxis2 = ON_CrossProduct(ray.m_V, yaxis2);
		xaxis2.Unitize(), yaxis2.Unitize();
		ON_Plane rayplane2(ray.m_P, xaxis2, yaxis2);
		ProjectToRayPlane(rayplane2, bez_t, cppln);

		if (IsMonotonic(cppln, bez_t.Order(0), bez_t.Order(1))){
			double uv[2];
			if (GuessNewton(rayplane2, bez_t, cppln, uv)){
				Result r;
				r.u = intu.m_t[0]+(intu.m_t[1]-intu.m_t[0])*uv[0];
				r.v = intv.m_t[0]+(intv.m_t[1]-intv.m_t[0])*uv[1];
				results.push_back(r);
			}

		}else{
			for (int i = 0; i < bez.m_order[0]; ++i){
				du[i] = (intu.m_t[1] * static_cast<double>(i) +
						 intu.m_t[0] * static_cast<double>(bez.m_order[0] - 1 - i)) /
						 static_cast<double>(bez.m_order[0]-1);
			}
			for (int j = 0; j < bez.m_order[1]; ++j){
				dv[j] = (intv.m_t[1] * static_cast<double>(j) +
						 intv.m_t[0] * static_cast<double>(bez.m_order[1] - 1 - j)) /
						 static_cast<double>(bez.m_order[1]-1);
			}
			for (int j = 0, ji = 0; j < bez.m_order[1] - 1; ++j, ji += bez.m_order[0]){
				for (int i = 0; i < bez.m_order[0] - 1; ++i){
					stack_ints.push(std::make_pair(ON_Interval(du[i],du[i+1]), ON_Interval(dv[j],dv[j+1])));
				}
			}
		}
	}
#ifdef DEBUGOUT
	std::fclose(fp);
#endif
	if (results.size() == 0) return 0;
	// tを昇順ソートして、隣りとのtの差がトレランス以上かどうかを判断する。
	// トレランス以下の場合は外す。
	std::vector<int> indices(results.size());
	for (int i = 0; i < static_cast<int>(indices.size()); ++i) indices[i] = i;
	struct Sort{
		const std::vector<Result> &v;
		Sort(const std::vector<Result> &v_) : v(v_){}
		bool operator ()(const int &lhs, const int &rhs) const{
			return v[lhs].t < v[rhs].t;
		}
	};
	std::sort(indices.begin(), indices.end(), Sort(results));
	for (size_t i = indices.size()-1; i > 0; --i){
		if (results[indices[i]].t - results[indices[i-1]].t < tolerance) indices.erase(indices.begin()+i-1);
	}
	for (size_t i = 0; i < indices.size(); ++i){
		Result &r = results[indices[i]];
		tuvints.Append(ON_3dPoint(r.t, r.u, r.v));
		ptsrfs.Append(r.ptsrf);
		ptlins.Append(r.ptlin);
	}
	return 0;
}

int ONGEO_IntersectRayBezier_Subdiv_Newton(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance){
	if (bez.IsRational()){
		return IntersectRayBezier_Subdiv_Newton_<ON_4dPoint>(ray, bez, tuvints, ptsrfs, ptlins, tolerance);
	}else{
		return IntersectRayBezier_Subdiv_Newton_<ON_3dPoint>(ray, bez, tuvints, ptsrfs, ptlins, tolerance);
	}
}
