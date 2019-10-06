// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.

#include "ONGEO.h"

ONGEO_NurbsCurveBezierCache::ONGEO_NurbsCurveBezierCache(const ON_NurbsCurve &nc){
	if (nc.SpanCount() == 0) return;

	int order = nc.Order();
	prms.Append(nc.Knot(order-2));

	for (int j = 0; j <= nc.CVCount() - order; ++j){
		double prm_end = nc.Knot(j+order-1);
		if (*prms.Last() == prm_end) continue;
		prms.Append(prm_end);

		ON_BezierCurve &bc = bcs.AppendNew();
		nc.ConvertSpanToBezier(j, bc);
	}
}

double ONGEO_NurbsCurveBezierCache::GetNurbsParameterFromBezierParameter(int bezier_idx, double bezier_t) const{
	if (bezier_idx < 0 || bezier_idx >= prms.Count() - 1) return ON_UNSET_VALUE;
	return ON_Interval(prms[bezier_idx], prms[bezier_idx+1]).ParameterAt(bezier_t);
}
