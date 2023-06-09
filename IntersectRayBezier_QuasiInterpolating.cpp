﻿// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
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

#include "ONGEO_QuasiInterpolating.hpp"

//#define DEBUGOUT

#include "Profile.h"

#ifdef DEBUGOUT
#define DPRINTF(X, ...) std::fprintf(X, __VA_ARGS__), std::fflush(fp);
#else
#define DPRINTF(X, ...) 
#endif

#define ZT2 (ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE)

namespace {
template <typename T> void EvaluateBiLinearSurf(const ON_3dRay &ray, T cods[4], double u, double v, double &t, ON_3dPoint &ptsrf, ON_3dPoint &ptlin){
	PROF("EvaluateBiLinearSurf");
	T pt;
	double uv = u * v;
	for (int i = 0; i < ONGEO_QI_PointDimTraits<T, ON_BezierSurface>::DIM; ++i){
		pt[i] = (cods[3][i] - cods[2][i] - cods[1][i] + cods[0][i]) * uv + (cods[1][i] - cods[0][i]) * u + (cods[2][i] - cods[0][i]) * v + cods[0][i];
	}
	ptsrf = pt;
	ON_Line lin(ray.m_P, ray.m_P+ray.m_V);
	lin.ClosestPointTo(ptsrf, &t);
	ptlin = lin.PointAt(t);
}

// 0:解が得られた。1:解が得られない。-1:入力点数が4個ではない。
// inrange 0:両方範囲外、1:[0]の方が範囲内、2:[1]の方が範囲内、3:両方範囲内
int IntersectZeroPointBiLinear(const ON_2dPoint qicp2d[4], double u[2], double v[2], int &inrange, double err2, double tolerance){
	PROF("IntersectZeroPointBiLinear");
	const double tole_z2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;
	inrange = 0;
	bool inv_u = false ,inv_v = false;
	ON_2dVector ds, du, dv;
	ON_2dPoint p11 = qicp2d[0], p12 = qicp2d[1], p21 = qicp2d[2], p22 = qicp2d[3];
	du = p12 - p11, dv = p21 - p11;
	double dul2 = du.LengthSquared(), dvl2 = dv.LengthSquared();

	// 縮退しているときは反対側を使う
	if (dul2 < tole_z2){
		inv_v = true;
		ON_2dPoint pw;
		pw = p11, p11 = p21, p21 = pw;
		pw = p12, p12 = p22, p22 = pw;
		du = p12 - p11, dv *= -1;
		dul2 = du.LengthSquared();
		if (dul2 < tole_z2){
			// 面とZeroPointの線が平行
			return 1;
		}
	}else if (dvl2 < tole_z2){
		inv_u = true;
		ON_2dPoint pw;
		pw = p11, p11 = p12, p12 = pw;
		pw = p21, p21 = p22, p22 = pw;
		du *= -1, dv = p21 - p11;
		dvl2 = du.LengthSquared();
		if (dvl2 < tole_z2){
			// 面とZeroPointの線が平行
			return 1;
		}
	}
	ON_2dPoint p = p11;
	ds = p22 - p21 - p12 + p11;

	double denom_u = 2*(ds[0]*du[1]-ds[1]*du[0]);
	double denom_v = 2*(ds[0]*dv[1]-ds[1]*dv[0]);

	// 平行な辺が存在するとき
	bool zero_u = std::abs(denom_u) < tolerance;
	bool zero_v = std::abs(denom_v) < tolerance;
	if (zero_u && zero_v) { // 平行四辺形、菱形、長方形、正方形
		u[0] = (dv[1]*p[0]-dv[0]*p[1])/(du[1]*dv[0]-du[0]*dv[1]);
		v[0] = (du[0]*p[1]-du[1]*p[0])/(du[1]*dv[0]-du[0]*dv[1]);
		u[1] = -1, v[1] = -1;
	}else{
		// (0,0) = u*v*ds+u*du+v*dv+p を u,vについて解くと下記のようになる。
		double du1dv0 = du[1]*dv[0], du0dv1 = du[0]*dv[1];
		double du1dv0_du0dv1 = du1dv0-du0dv1;
		double ds0p1_ds1p0 = ds[0]*p[1]-ds[1]*p[0];
		double R =
			du1dv0_du0dv1*du1dv0_du0dv1 + ds0p1_ds1p0*ds0p1_ds1p0
			+ 2*( 2*(du[1]*dv[1]*ds[0]*p[0] + du[0]*dv[0]*ds[1]*p[1]) - (ds[0]*p[1]+ds[1]*p[0])*(du0dv1+du1dv0) );
		if (R < 0) return 1;
		double R_2 = std::sqrt(R);
		double nume_u = -du1dv0_du0dv1-ds0p1_ds1p0;
		double nume_v = -du1dv0_du0dv1+ds0p1_ds1p0;

		// 解は2つできる。関数を使う側でu,vともに0-1の中に入っている方を採用することで、双一次曲線と直線との交点を
		// 求めることができる。

		if (!zero_u){
			u[0] = (nume_u + R_2) / denom_u;
			u[1] = (nume_u - R_2) / denom_u;
		}

		if (!zero_v){
			v[0] = -(nume_v + R_2) / denom_v;
			v[1] = -(nume_v - R_2) / denom_v;
		}

		// denom_u や denom_v が小さすぎて除数に使えない場合は下記の計算で求める。
		if (zero_v) {
			for (int i = 0; i < 2; ++i){
				ON_2dPoint denom = ds * u[i] + dv, nume = du * u[i] + p11;
				if (std::abs(denom.x) > ON_ZERO_TOLERANCE){
					v[i] = - nume.x / denom.x;
				}else{
					v[i] = - nume.y / denom.y;
				}
			}
		} else if (zero_u) {
			for (int i = 0; i < 2; ++i){
				ON_2dPoint denom = ds * v[i] + du, nume = dv * v[i] + p11;
				if (std::abs(denom.x) > ON_ZERO_TOLERANCE){
					u[i] = - nume.x / denom.x;
				}else{
					u[i] = - nume.y / denom.y;
				}
			}
		}

		if (inv_u) u[0] = 1 - u[0], u[1] = 1 - u[1];
		if (inv_v) v[0] = 1 - v[0], v[1] = 1 - v[1];
	}
//	double margin = err2 ? std::sqrt(err2 / (dul2 + dvl2)) : 0;
	double margin = 0.25;

	if (u[0] >= -margin && u[0] <= 1+margin && v[0] >= -margin && v[0] <= 1+margin) inrange |= 1;
	if (u[1] >= -margin && u[1] <= 1+margin && v[1] >= -margin && v[1] <= 1+margin) inrange |= 2;

	return 0;
}

void ProjectToRayPlane(const ON_Plane &rnPlane, const ON_SimpleArray<ON_3dPoint> &qicp, ON_2dPointArray &qicp2d){
	PROF("ProjectToRayPlane<ON_3dPoint>");
	qicp2d.SetCapacity(qicp.Count());
	qicp2d.SetCount(qicp.Count());
	for (int i = 0; i < qicp.Count(); ++i){
		rnPlane.ClosestPointTo(qicp[i], &qicp2d[i][0], &qicp2d[i][1]);
	}
}

void ProjectToRayPlane(const ON_Plane &rnPlane, const ON_SimpleArray<ON_4dPoint> &qicp, ON_2dPointArray &qicp2d){
	PROF("ProjectToRayPlane<ON_4dPoint>");
	qicp2d.SetCapacity(qicp.Count());
	qicp2d.SetCount(qicp.Count());
	for (int i = 0; i < qicp.Count(); ++i){
		rnPlane.ClosestPointTo(qicp[i], &qicp2d[i][0], &qicp2d[i][1]);
		qicp2d[i] *= qicp[i].w;
	}
}
}
/// 0  : 線分に接していなく、かつ、閉じた領域の外
/// 1  : 線分に接している、または閉じた領域の中
__declspec(dllexport) int IntersectionTest(const ON_2dPointArray &qicp2d, int iu, int iv, int stride, double error2, double *err2_vertices){
	PROF("IntersectionTest");
	const int index = iv * stride + iu;
	const int indices[] = {index + stride, index, index + 1, index + 1 + stride, index + stride, index};

	ON_2dVector pts[6];
	for (int i = 1; i < 5; ++i){
		pts[i] = qicp2d[indices[i]];
	}
	pts[0] = pts[4];
	pts[5] = pts[1];

	// 4点のBoundingBoxを計算し、その距離^2がerror2よりも大きければ閉じた領域の外で線分にも接していない。
	ON_2dPoint bbmin = pts[0], bbmax = pts[0];
	for (int i = 0; i < 4; ++i){
		if (bbmin.x > pts[i].x) bbmin.x = pts[i].x;
		else if (bbmax.x < pts[i].x) bbmax.x = pts[i].x;
		if (bbmin.y > pts[i].y) bbmin.y = pts[i].y;
		else if (bbmax.y < pts[i].y) bbmax.y = pts[i].y;
	}
	bool bx = (bbmin.x * bbmax.x < 0), by = (bbmin.y * bbmax.y < 0);
	if (!bx && !by){ // 4隅が全て同じ象限に固まっている
		ON_2dPoint ptnear = bbmin;
		if (std::abs(ptnear.x) > std::abs(bbmax.x)) ptnear.x = bbmax.x;
		if (std::abs(ptnear.y) > std::abs(bbmax.y)) ptnear.y = bbmax.y;
		if (error2 < ptnear.x * ptnear.x + ptnear.y * ptnear.y) return 0;
	}else if (!bx){  // 1,4象限、または2,3象限にまたがっている
		if (std::abs(bbmin.x) < std::abs(bbmax.x)){
			if (error2 < bbmin.x * bbmin.x) return 0;
		}else{
			if (error2 < bbmax.x * bbmax.x) return 0;
		}
	}else if (!by){  // 1,2象限、または3,4象限にまたがっている
		if (std::abs(bbmin.y) < std::abs(bbmax.y)){
			if (error2 < bbmin.y * bbmin.y) return 0;
		}else{
			if (error2 < bbmax.y * bbmax.y) return 0;
		}
	}

	// 頂点との距離がerror以内だったら交差あり
	for (int i = 1; i < 5; ++i){
		double &err2_vertex = err2_vertices[indices[i]];
		if (err2_vertex < 0) err2_vertex = (pts[i].x * pts[i].x + pts[i].y * pts[i].y);
		if (err2_vertex < error2) return 1;
	}

	// エッジに接していたら交差あり
	double vy[5], cc[5];
	int ints_count = 0;
	for (int i = 1; i < 5; ++i){
		// 陰関数で計算することで高速化
		ON_2dVector p1 = pts[i], p2 = pts[i+1];
		ON_2dVector v = p2 - p1;
		double c = p1.x*p2.y - p1.y*p2.x;
		cc[i] = c;
		// 線分p1-p2と重なる直線をl1とすると、l1 : -v.y * x + v.x * y + c = 0 となる。
		// 原点を通り、l1と垂直な直線l2は、   l2 :  v.x * x + v.y * y = 0 となる。
		double l2inv = 1.0 / (v.x*v.x + v.y*v.y); // v.LengthSquared();
		// l1 と l2 の交点Iは、 (c*v.y*l2inv, -c*v.x*l2inv) となる。
		// また、直線l1と原点との距離は |c*l2inv * (v.y, -v.x)| = c*l2inv*sqrt(v.LengthSqured()) = c/v.Length()
		// よって、距離の2乗は c*c*l2invとなる。
		vy[i] = v.y;
		double cl2inv = c * l2inv;
		// 距離計算
		if (c*cl2inv < error2){
			// 線分と交点との内外判定
			// 交点が原点に来るように線分を平行移動して判定
			if (std::abs(v.x) > std::abs(v.y)){
				double xi = -v.y*cl2inv;
				p1.x += xi, p2.x += xi;
				if (p1.x * p2.x < 0) return 1;
			}else{
				double yi = v.x*cl2inv;
				p1.y += yi, p2.y += yi;
				if (p1.y * p2.y < 0) return 1;
			}
		}
	}

	// 内外判定
	// y=0の水平線を横切り、かつそのときのxが正である回数
	for (int i = 1; i < 5; ++i){
		// 下記の2条件のどちらかを満たせばカウント
		// 条件1 : pts[i].yとpts[i+1].yの符号が逆で、その線分のy=0のときのxの値が正
		// 条件2 : pts[i].yがちょうど0で、かつ、前後のyの符号が逆(=V字型になっていない)
		if ((pts[i].y * pts[i+1].y < 0 && cc[i] * vy[i] >= 0) || 
			(pts[i].x > 0 && pts[i].y == 0 && pts[i-1].y * pts[i+1].y < 0)) ++ints_count;
	}

	return ints_count & 1;
}

template <typename T> int IntersectRayBezier_QuasiInterpolating_(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance){
	PROF("IntersectRayBezier_QuasiInterpolating_");

	std::stack<std::pair<ON_Interval, ON_Interval> > stack_ints;
	stack_ints.push(std::make_pair(ON_Interval(0,1), ON_Interval(0,1)));
	ON_2dPointArray qicp2d;
	ON_SimpleArray<T> qicp;
	ON_Plane rayplane(ray.m_P, ray.m_V);
	std::vector<double> du(bez.m_order[0]), dv(bez.m_order[1]);
	std::vector<int> jis(bez.m_order[1]);
	std::vector<double> err2_vertices(bez.m_order[0]*bez.m_order[1]);
	struct Result{
		double t, u, v;
		ON_3dPoint ptsrf;
		ON_3dPoint ptlin;
	};
	std::vector<Result> results;
	for (int j = 0, ji = 0; j < bez.m_order[1]; ++j, ji += bez.m_order[0]) jis[j] = ji;
#ifdef DEBUGOUT
	FILE *fp = std::fopen("d:/pnttest.txt", "a"); // dbg
#endif
	DPRINTF(fp, "\n\n\n'==== start ==== %d %d %s\n", bez.Order(0), bez.Order(1), bez.IsRational() ? "Rational" : "Polynomial");
	while(stack_ints.size()){
		ON_BezierSurface bez_t(bez);
		ON_Interval intu = stack_ints.top().first;
		ON_Interval intv = stack_ints.top().second;
		stack_ints.pop();
		{
			PROF("Trim2");
			bez_t.Trim(0, intu);
			bez_t.Trim(1, intv);
		}
		double err = ONGEO_QI_GetError<T>(bez_t), err2 = err * err;
		if (err2 < ZT2) err2 = ZT2;
		ONGEO_QI_GetQuasiInterpControlPoints(bez_t, qicp);
#ifdef DEBUGOUT
		double wmax, wmin;
		ONGEO_CalculateMinMaxWeight(bez_t, wmin, wmax);
		DPRINTF(fp, "MaxWeight:%f\tMinWeight:%f\n", wmax, wmin);
		ON_BoundingBox bb = bez_t.BoundingBox();
		DPRINTF(fp, "Error:%f\tBBSize:%f\n", err, bb.m_max.DistanceTo(bb.m_min));
		// dbg
		{
			std::fprintf(fp, "' Quasi-\n");
			{
				std::fprintf(fp, "Rhino.AddNurbsSurface Array(%d,%d), Array( _\n  ", bez.Order(0), bez.Order(1));
				for (int h = 0; h < qicp.Count(); ++h){
					std::fprintf(fp, "Array(%f,%f,%f), ", qicp[h].x, qicp[h].y, qicp[h].z);
				}
				std::string k[2];
				for (int h = 0; h < 2; ++h){
					for (int g = 0; g <= bez.Degree(h); ++g){
						char buf[10]; std::sprintf(buf, "%d%c", g, (g < bez.Degree(h)) ? ',' : ' ');
						k[h] += buf;
					}
				}
				std::fprintf(fp, "Array(0,0,0)), _\n Array(%s), Array(%s), Array(1,1)\n", k[0].c_str(), k[1].c_str());
			}
			if (PointDimTraits<T>::DIM == 4){
				std::fprintf(fp, "Rhino.AddNurbsSurface Array(%d,%d), Array( _\n  ", bez.Order(0), bez.Order(1));
				for (int h = 0; h < qicp.Count(); ++h){
					ON_3dPoint p = qicp[h];
					std::fprintf(fp, "Array(%f,%f,%f), ", p.x, p.y, p.z);
				}
				std::string k[2];
				for (int h = 0; h < 2; ++h){
					for (int g = 0; g <= bez.Degree(h); ++g){
						char buf[10]; std::sprintf(buf, "%d%c", g, (g < bez.Degree(h)) ? ',' : ' ');
						k[h] += buf;
					}
				}
				std::fprintf(fp, "Array(0,0,0)), _\n Array(%s), Array(%s), Array(1,1)\n", k[0].c_str(), k[1].c_str());
			}
		}
		// dbg
#endif
		ProjectToRayPlane(rayplane, qicp, qicp2d);

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

		for (int j = 0, ji = 0; j < bez.m_order[1]; ++j, ji += bez.m_order[0]){
			for (int i = 0; i < bez.m_order[0]; ++i){
				err2_vertices[ji+i] = -1;
			}
		}

		for (int j = 0, ji = 0; j < bez.m_order[1] - 1; ++j, ji += bez.m_order[0]){
			int jim = ji + j;
			for (int i = 0; i < bez.m_order[0] - 1; ++i){
#ifdef DEBUGOUT
				// dbg
				std::fprintf(fp, "' Test u:%f - %f, v:%f - %f\n", du[i], du[i+1], dv[j], dv[j+1]);
				{
					ON_3dPoint p[] = {qicp[ji+i], qicp[ji+i+1], qicp[ji+bez.m_order[0]+i], qicp[ji+bez.m_order[0]+i+1]};
					std::fprintf(fp, "Rhino.AddNurbsSurface Array(2,2), Array( _\n  ");
					for (int h = 0; h < 4; ++h){
						std::fprintf(fp, "Array(%f,%f,%f), ", p[h].x, p[h].y, p[h].z);
					}
					std::fprintf(fp, "Array(0,0,0)) _\n Array(0,1), Array(0,1), Array(1,1)\n");
					ON_2dPoint p2[] = {qicp2d[ji+i], qicp2d[ji+i+1], qicp2d[ji+bez.m_order[0]+i], qicp2d[ji+bez.m_order[0]+i+1]};
					for (int h = 0; h < 4; ++h){
						std::fprintf(fp, "'2D  %f,%f\n", p2[h].x, p2[h].y);
					}
				}
				// dbg
#endif
				if (IntersectionTest(qicp2d, i, j, bez.m_order[0], err2, &err2_vertices[0])){
					// dbg
					DPRINTF(fp, "'... True\n");
					// dbg
					if (err < tolerance){
						ON_2dPoint pts[] = {qicp2d[ji+i], qicp2d[ji+i+1], qicp2d[ji+bez.m_order[0]+i], qicp2d[ji+bez.m_order[0]+i+1]};
						T pts3d[] = {qicp[ji+i], qicp[ji+i+1], qicp[ji+bez.m_order[0]+i], qicp[ji+bez.m_order[0]+i+1]};
						ON_3dPoint ptsrf, ptlin;
						double u[2], v[2];
						int inrange;
						IntersectZeroPointBiLinear(pts, u, v, inrange, err2, tolerance);
						DPRINTF(fp, "'IntersectZeroPointBiLinear => (%f %f), (%f %f),  %d\n", u[0], v[0], u[1], v[1], inrange);
						for (int h = 0; h < 2; ++h){
							if ((inrange & (h+1)) == 0) continue;
							Result r;
							EvaluateBiLinearSurf(ray, pts3d, u[h], v[h], r.t, r.ptsrf, r.ptlin);
							if (r.t < 0) continue;
							r.u = du[i]+(du[i+1]-du[i])*u[h];
							r.v = dv[j]+(dv[j+1]-dv[j])*v[h];
							results.push_back(r);
#ifdef DEBUGOUT
							// dbg
							fprintf(fp, "' ===== Detect%d =====\n", h + 1);
							fprintf(fp, "' t:%f u:%f v:%f ptsrf:(%f %f %f) ptlins:(%f %f %f)\n",
								tuvints.Last()->x, tuvints.Last()->y, tuvints.Last()->z,
								ptsrf.x, ptsrf.y, ptsrf.z,
								ptlin.x, ptlin.y, ptlin.z);
							// dbg
#endif
						}
					}else{
						stack_ints.push(std::make_pair(ON_Interval(du[i],du[i+1]), ON_Interval(dv[j],dv[j+1])));
					}
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

int ONGEO_IntersectRayBezier_QuasiInterpolating(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance){
	if (bez.IsRational()){
		return IntersectRayBezier_QuasiInterpolating_<ON_4dPoint>(ray, bez, tuvints, ptsrfs, ptlins, tolerance);
	}else{
		return IntersectRayBezier_QuasiInterpolating_<ON_3dPoint>(ray, bez, tuvints, ptsrfs, ptlins, tolerance);
	}
}
