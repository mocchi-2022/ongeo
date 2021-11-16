// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#include "opennurbs.h"
#include "ONGEO.h"

#include <cmath>

#include <vector>
#include <limits>
#include <list>
#include <set>
#include <map>
#include <algorithm>

#include <cstdio>

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

namespace{
	template <int N> struct rotator{
	};
	template <> struct rotator<3>{
		double rad_s, rad_c;
		ON_3dVector crvec;
		rotator(const ON_3dPoint *ctrlpts, size_t order){
			ON_3dVector cptaxis = ((ctrlpts[1] - ctrlpts[0]) + (ctrlpts[order-1] - ctrlpts[order-2]));
			double length2 = cptaxis.LengthSquared();
			if (!length2){
				cptaxis = (ctrlpts[order-1] - ctrlpts[0]);
				length2 = cptaxis.LengthSquared();
			}
			if (length2){
				cptaxis /= std::sqrt(length2);
				ON_3dVector xaxis(1,0,0);
				crvec = ON_CrossProduct(cptaxis, xaxis);
				double length = crvec.Length();
				rad_s = length;
				rad_c = ON_DotProduct(cptaxis, xaxis);
				crvec /= length;
			}else{
				rad_s = 0, rad_c = 1;
				crvec.Set(1,0,0);
			}
		}
		void rotate(ON_3dVector &vec){
			vec.Rotate(rad_s, rad_c, crvec);
		}
	};
	template <> struct rotator<2>{
		double rad_s, rad_c;
		rotator(const ON_3dPoint *ctrlpts, size_t order){
			ON_2dVector cptaxis = ((ctrlpts[1] - ctrlpts[0]) + (ctrlpts[order-1] - ctrlpts[order-2]));
			double length2 = cptaxis.LengthSquared();
			if (!length2){
				cptaxis = (ctrlpts[order-1] - ctrlpts[0]);
				length2 = cptaxis.LengthSquared();
			}
			if (length2){
				cptaxis /= std::sqrt(length2);
				ON_2dVector xaxis(1,0);
				rad_s = ON_CrossProduct(cptaxis, xaxis)[2];
				rad_c = ON_DotProduct(cptaxis, xaxis);
			}else{
				rad_s = 0, rad_c = 1;
			}
		}
		void rotate(ON_3dVector &vec){
			ON_2dVector v(vec);
			v.Rotate(rad_s, rad_c);
			vec = v;
		}
	};

	template <int N> struct type_traits{
	};
	template <int N> double distance2(const ON_BoundingBox &bb, const ON_3dPoint &query){
		const ON_3dPoint &min_pt = bb.Min();
		const ON_3dPoint &max_pt = bb.Max();

		double d, sum = 0;
		if (min_pt.x > query.x) sum += ((d = (min_pt.x - query.x)), d * d); 
		else if (max_pt.x < query.x) sum += ((d = (query.x - max_pt.x)), d * d);

		if (min_pt.y > query.y) sum += ((d = (min_pt.y - query.y)), d * d); 
		else if (max_pt.y < query.y) sum += ((d = (query.y - max_pt.y)), d * d);

		if (min_pt.z > query.z) sum += (d = (min_pt.z - query.z)), d * d; 
		else if (max_pt.z < query.z) sum += ((d = (query.z - max_pt.z)), d * d);

		return sum;
	}
}

typedef int RESULT;
const RESULT INPUT_GEOMETRY_INVALID = -1;
const RESULT TOLERANCE_INVALID = -3;
template <int N> double Nearest_BBClipping(const ON_BezierCurve *bc_begin, const ON_BezierCurve *bc_end, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest){
	if (bc_begin > bc_end || !bc_begin) return INPUT_GEOMETRY_INVALID;
	if (tolerance <= 0) return TOLERANCE_INVALID;
	// first : tの範囲, second : 分割後曲線
	typedef std::pair<std::pair<double, double>, std::pair<ON_BezierCurve, const ON_BezierCurve *> > data;
	std::list<data> data_buffer;
	std::vector<data *> bcs1, bcs2, *bcs, *bcs_next, bcs_last;
	bcs = &bcs1; bcs_next = &bcs2;

	ON_BezierCurve bc1, bc2;
	double ts[3];
	double certain_dist2 = std::numeric_limits<double>::max();
	// certain_dist2 : クエリと曲線上の適当に見つけたとある点との距離の二乗の中で最小の値
	double min_dist2 = 0;

	const int num_ptsearch = 1;
	const double ztole2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;
	double tolerance2 = tolerance * tolerance;

	ts[0] = 0, ts[1] = 1.0;
	int order_max = 0;
	for (const ON_BezierCurve *iter = bc_begin; iter != bc_end; ++iter){
		if (iter->Order() < 2) continue;
		data_buffer.push_back(std::make_pair(std::make_pair(ts[0], ts[1]), std::make_pair(*iter,iter)));
		bcs->push_back(&data_buffer.back());
		if (order_max < iter->Order()) order_max = iter->Order();
		int o2 = iter->Order() * iter->Order();
		double td = 1.0 / static_cast<double>(o2), t = 0;
		for (int i = 0; i <= o2; ++i){
			certain_dist2 = std::min(certain_dist2, (iter->PointAt(t) - pt_query).LengthSquared());
			t += td;
		}
	}
	if (bcs->size() == 0 || order_max < 2) return INPUT_GEOMETRY_INVALID;
	ON_SimpleArray<ON_3dPoint> ctrlpts(order_max);
		std::vector<double> mindist(bcs->size());
	std::vector<bool> push_candidate(bcs->size());
	while(bcs->size()){
		min_dist2 = std::numeric_limits<double>::max();
		mindist.resize(bcs->size());
		push_candidate.resize(bcs->size());
		for (size_t j = 0; j < bcs->size(); ++j){
			data &d = *(*bcs)[j];
			int order = d.second.first.Order();
			for (int i = 0; i < order; ++i){
				d.second.first.GetCV(i, ctrlpts[i]);
			}
			certain_dist2 = std::min(certain_dist2, std::min(
				(ctrlpts[0] - pt_query).LengthSquared(),
				(ctrlpts[order - 1] - pt_query).LengthSquared()
			));

			ON_3dVector pt_query_rotate = pt_query;
			{
				rotator<N-1> r(&ctrlpts[0], order);
				r.rotate(pt_query_rotate);
				for (int i = 0; i < order; ++i){
					ON_3dVector v(ctrlpts[i]);
					r.rotate(v);
					ctrlpts[i] = v;
				}
			}

			ON_BoundingBox bb;
			for (int i= 0; i < order; ++i){
				bb.Set(ctrlpts[i], 1);
			}
			double &min_dist2_cur = mindist[j];

			min_dist2_cur = distance2<N>(bb, ON_3dPoint(pt_query_rotate));
			if (min_dist2 > min_dist2_cur) min_dist2 = min_dist2_cur;
			if (min_dist2 > certain_dist2 + ON_ZERO_TOLERANCE) continue;
			// min_dist2 : 曲線セグメント群とクエリとの距離の最小値の二乗よりも必ず小さい値

			// 各制御点間距離がトレランス以下の場合分割打ち切り
			bool candidate = true;
			for (int i = 0; i < order-1; ++i){
				if ((ctrlpts[i]-ctrlpts[i+1]).LengthSquared() > tolerance2){
					candidate = false;
					break;
				}
			}
			push_candidate[j] = candidate;
		}
		bcs_next->clear();
		for (size_t i = 0; i < bcs->size(); ++i){
			if (certain_dist2 + ztole2 < mindist[i]) continue;
			data &d = *(*bcs)[i];
			if (push_candidate[i]){
				bcs_last.push_back(&d);
				continue;
			}
			ts[0] = d.first.first;
			ts[2] = d.first.second;
			ts[1] = (ts[0] + ts[2]) * 0.5;
			d.second.first.Split(0.5, bc1, bc2);
			data_buffer.push_back(std::make_pair(std::make_pair(ts[0], ts[1]), std::make_pair(bc1,d.second.second)));
			bcs_next->push_back(&data_buffer.back());
			data_buffer.push_back(std::make_pair(std::make_pair(ts[1], ts[2]), std::make_pair(bc2,d.second.second)));
			bcs_next->push_back(&data_buffer.back());
		}
		std::swap(bcs_next, bcs);
	}

	pt_nearest.Set(0,0,0);
	bc_nearest = 0;
	t = 0;
	double distance2_nearest = std::numeric_limits<double>::max();
	for (size_t i = 0; i < bcs_last.size(); ++i){
		data &d = *bcs_last[i];
		for (int h = 0; h < d.second.first.Order(); ++h){
			ON_3dPoint pt;
			d.second.first.GetCV(h, pt);
			double dist = (pt-pt_query).LengthSquared();
			if (distance2_nearest >= dist){
				distance2_nearest = dist;
				double hh = (double)h / (double)(d.second.first.Order()-1);
				t = (1-hh)*d.first.first + hh*d.first.second;
				pt_nearest = pt;
				bc_nearest = d.second.second;
			}
		}
	}
	return std::sqrt(distance2_nearest);
}

double ONGEO_NearestPointBezierCurve_BBClipping(const ON_BezierCurve *bc_begin, const ON_BezierCurve *bc_end, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest){
	bool z_is_0 = true;
	double pt[4];
	if (!bc_begin) return INPUT_GEOMETRY_INVALID;
	for (const ON_BezierCurve *iter = bc_begin; iter != bc_end; ++iter){
		const ON_BezierCurve &cur = *iter;
		for (int i = 0; i < cur.CVCount(); ++i){
			cur.GetCV(i, ON::homogeneous_rational, pt);
			if (pt[2]){
				z_is_0 = false;
				break;
			}
		}
		if (z_is_0) break;
	}
	if (z_is_0){
		return Nearest_BBClipping<3>(bc_begin, bc_end, tolerance, pt_query, bc_nearest, t, pt_nearest);
	}else{
		return Nearest_BBClipping<4>(bc_begin, bc_end, tolerance, pt_query, bc_nearest, t, pt_nearest);
	}
}
