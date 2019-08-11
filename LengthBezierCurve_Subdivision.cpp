#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <cstdio>

#include "ONGEO_CurveSubdivision.hpp"

// Reference : Maharavo Randrianarivony,
//   "Length Estimation of Rational Bezier Curves and Application to CAD Parametrization".

inline void GetLengthBounds(const ON_BezierCurve &bc_cur, double &lb, double &ub){
	ON_3dPoint pt_start, pt_end;
	bc_cur.GetCV(0, pt_start);
	bc_cur.GetCV(bc_cur.CVCount()-1, pt_end);
	lb = pt_end.DistanceTo(pt_start);
	ub = bc_cur.ControlPolygonLength();
}

struct LengthOperator{
	// input
	double tolerance;

	// output
	double sum_length;
	LengthOperator(double tolerance_) : tolerance(tolerance_), sum_length(0){
	}
	bool TestConverge(const ON_BezierCurve &bc_cur, const ON_Interval &t_int){
		double lb, ub;
		GetLengthBounds(bc_cur, lb, ub);
		bool converged = ((ub - lb) < tolerance * t_int.Length());
		if (converged){
			sum_length += 0.5 * (ub + lb);
		}
		return converged;
	}
};

double ONGEO_LengthBezierCurve_SimpleSubdivision(const ON_BezierCurve &bc, double tolerance){
	double lb_0, ub_0;
	GetLengthBounds(bc, lb_0, ub_0);

	// ‚Ù‚Ú’¼ü
	if ((ub_0 - lb_0) < tolerance){
		return 0.5 * (ub_0 + lb_0);
	}

	LengthOperator lo(tolerance);
	ONGEO_CurveSubdivision(lo, bc);
	return lo.sum_length;
}

