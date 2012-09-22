// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"

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
	ON_3dPoint pt;
	double dist = ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod(loop_crvs.First(), loop_crvs.Count(), tolerance, ON_3dPoint(uv), bcn, t, pt);
	ON_3dVector dir_t;
	bcn->EvTangent(t, pt, dir_t);
	ON_2dVector dir_ct = ON_2dPoint(pt) - uv;
	return (dir_ct.Length() < tolerance || ON_CrossProduct(dir_ct, ON_2dVector(dir_t))[2] >= 0);
}
