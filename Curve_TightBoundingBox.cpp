// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "ONGEO.h"
#include <stack>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "Profile.h"

int ONGEO_CalculateTightBoundingBox(const ON_BezierCurve *bcs, int num_bcs, double tolerance, ON_BoundingBox &bb){
	ON_BoundingBox rbb, tbb;
	for (int i = 0; i < num_bcs; ++i){
		rbb.Union(bcs[i].BoundingBox());
	}

	ON_SimpleArray<ON_BezierCurve> crvs;
	crvs.SetCapacity(num_bcs);
	crvs.SetCount(num_bcs);
	for (int j = 0; j < 3; ++j){
		if (rbb.m_min[j] == rbb.m_max[j]){
			tbb.m_min[j] = tbb.m_max[j] = rbb.m_min[j];
		}else{
			for (int i = 0; i < num_bcs; ++i){
				crvs[i].Create(1, bcs[i].m_is_rat, bcs[i].m_order);
				for (int h = 0; h < bcs[i].CVCount(); ++h){
					double *dcv = crvs[i].CV(h), *scv = bcs[i].CV(h);
					dcv[0] = scv[j];
					if (bcs[i].m_is_rat) dcv[1] = scv[bcs[i].m_dim];
				}
			}
			ON_3dPoint ptqmin(0,0,0), ptqmax(0,0,0), ptn;
			ptqmin[0] = rbb.m_min[j];
			ptqmax[0] = rbb.m_max[j];
			const ON_BezierCurve *bcn;
			double t;
			ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod
				(crvs.First(), num_bcs, tolerance, ptqmin, bcn, t, ptn);
			tbb.m_min[j] = ptn[0] - tolerance;
			ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod
				(crvs.First(), num_bcs, tolerance, ptqmax, bcn, t, ptn);
			tbb.m_max[j] = ptn[0] + tolerance;
		}
	}
	bb = tbb;

	return 0;
}

int ONGEO_CalculateTightBoundingBox(const ONGEO_NurbsCurveBezierCache &nbc_c, double tolerance, ON_BoundingBox &bb){
	return ONGEO_CalculateTightBoundingBox(nbc_c.bcs, nbc_c.bcs.Count(), tolerance, bb);
}
