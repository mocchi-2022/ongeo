// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <limits>

void ONGEO_GetBezierLoops(const ON_BrepFace &face, ON_SimpleArray<ON_BezierCurve> &loop_crvs, ON_SimpleArray<int> &num_crvs_in_a_loop){
	num_crvs_in_a_loop.Empty();
	loop_crvs.Empty();
	for (int k = 0; k < face.LoopCount(); ++k){
		const ON_BrepLoop *l = face.Loop(k);
		int num_crvs_prev = loop_crvs.Count();
		for (int j = 0; j < l->TrimCount(); ++j){
			const ON_BrepTrim *trim = l->Trim(j);
			const ON_Curve *curv = trim->TrimCurveOf();

			ON_NurbsCurve nbc;
			curv->GetNurbForm(nbc);
			for (int i = 0; i <= nbc.m_cv_count - nbc.m_order; ++i){
				ON_BezierCurve bc;
				if (!nbc.ConvertSpanToBezier(i, bc)) continue;
				if (bc.IsValid()) loop_crvs.Append(bc);
			}
		}
		num_crvs_in_a_loop.Append(loop_crvs.Count() - num_crvs_prev);
	}
}

bool ONGEO_UVPointIsInside(const ON_SimpleArray<ON_BezierCurve> &loop_crvs, const ON_SimpleArray<int> &num_crvs_in_a_loop, const ON_2dPoint &uv, double tolerance){
	double t;
	const ON_BezierCurve *bcn;
	ON_3dPoint pt, ptdmy;
	double dist = ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod(loop_crvs.First(), loop_crvs.Count(), tolerance, ON_3dPoint(uv), bcn, t, pt);
	ON_3dVector dir_t, dir_tf, dir_tp;
	bcn->EvTangent(t, pt, dir_t);
	if (t <= ON_ZERO_TOLERANCE){
		const ON_BezierCurve *bcn_fwd = (bcn != loop_crvs.First()) ? bcn - 1 : bcn + loop_crvs.Count() - 1;
		bcn_fwd->EvTangent(1.0, ptdmy, dir_tf);
		dir_t.Unitize(), dir_tf.Unitize();
		dir_t += dir_tf;
	}else if(t >= 1 - ON_ZERO_TOLERANCE){
		const ON_BezierCurve *bcn_prev = (bcn != loop_crvs.First() + loop_crvs.Count() - 1) ? bcn + 1 : loop_crvs.First();
		bcn_prev->EvTangent(1.0, ptdmy, dir_tp);
		dir_t.Unitize(), dir_tp.Unitize();
		dir_t += dir_tp;
	}
	ON_2dVector dir_ct = ON_2dPoint(pt) - uv;
	return (dir_ct.Length() < tolerance || ON_CrossProduct(dir_ct, ON_2dVector(dir_t))[2] >= 0);
}

bool ONGEO_UVPointIsInside(const ON_Polyline *loop_pols, int num_loop_pols, const ON_2dPoint &uv, double tolerance){
	double t;
	ON_3dPoint pt, ptdmy, uv3d(uv);
	bool rc = false;
	double dist2_min = std::numeric_limits<double>::max();
	int imin = -1;
	for (int i = 0; i < num_loop_pols; ++i){
		ON_BoundingBox bb = loop_pols[i].BoundingBox();
		double dist2 = (bb.ClosestPoint(uv3d)-uv3d).LengthSquared();
		if (dist2_min < dist2) continue;

		bool ret = loop_pols[i].ClosestPointTo(uv3d, &t);
		if (!ret) continue;
		rc |= ret;
		pt = loop_pols[i].PointAt(t);
		dist2 = (uv-ON_2dPoint(pt)).LengthSquared();
		if (dist2_min > dist2){
			dist2_min = dist2;
			imin = i;
		}
	}
	if (!rc) return false;
	const ON_Polyline &loop_pol = loop_pols[imin];

	ON_3dVector dir_t = loop_pol.TangentAt(t), dir_tf, dir_tp;

	double ti;
	double tf = std::modf(t, &ti);
	if (tf <= ON_ZERO_TOLERANCE){
		dir_tf = loop_pol.TangentAt(t-0.5);
		dir_t += dir_tf;
	}else if(tf >= 1 - ON_ZERO_TOLERANCE){
		dir_tp = loop_pol.TangentAt(t+0.5);
		dir_t += dir_tp;
	}
	ON_2dVector dir_ct = ON_2dPoint(loop_pol.PointAt(t)) - uv;
	return (dir_ct.Length() < tolerance || ON_CrossProduct(dir_ct, ON_2dVector(dir_t))[2] >= 0);
}
