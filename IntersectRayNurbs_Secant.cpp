// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <algorithm>

bool ONGEO_IntersectRayNurbs_Secant(const ON_3dRay &ray, const ON_NurbsSurface &srf, const ON_2dPoint &uv0, const ON_2dPoint &uv1, ON_3dPoint &tuv, double tolerance, int max_iter){
	double tol2 = tolerance * tolerance;
	ON_2dPoint uv[3] = {uv1, uv0, ON_2dPoint()};
	for(int i = 0; i < max_iter; ++i){
		ON_3dPoint pt = srf.PointAt(uv[1].x, uv[1].y);
		ON_3dVector du = (pt - srf.PointAt(uv[0].x, uv[1].y));
		ON_3dVector dv = (pt - srf.PointAt(uv[1].x, uv[0].y));

		ON_Plane pln;
		pln.origin = pt, pln.xaxis = du, pln.yaxis = dv;
		double dulen = pln.xaxis.LengthAndUnitize();
		double dvlen = pln.yaxis.Length();

		pln.yaxis -= ON_DotProduct(pln.yaxis, pln.xaxis)*pln.xaxis;
		pln.yaxis.Unitize();
		pln.zaxis = ON_CrossProduct( pln.xaxis, pln.yaxis ), pln.zaxis.Unitize();
		pln.UpdateEquation();

		ON_Intersect(ON_Line(ray.m_P, ray.m_P+ray.m_V), pln, &tuv.x);
		ON_3dPoint ptint = ray.m_P+ray.m_V * tuv.x;
		ON_3dVector gap = ptint - pt;
		double dist2 = gap.LengthSquared();

		if (dist2 < tol2){
			if (!srf.Domain(0).Includes(uv[1].x) || !srf.Domain(1).Includes(uv[1].y)){
				return false;
			}
			ON_3dVector dt = ptint - ray.m_P;
			tuv.y = uv[1].x, tuv.z = uv[1].y;
			return true;
		}

		uv[2].x = ON_DotProduct(gap, pln.xaxis) * (uv[1].x - uv[0].x) / dulen + uv[1].x;
		uv[2].y = ON_DotProduct(gap, pln.yaxis) * (uv[1].y - uv[0].y) / dvlen + uv[1].y;
		uv[0] = uv[1], uv[1] = uv[2];
	}
	return false;
}
