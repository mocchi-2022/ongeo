#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <cstdio>

#include "ONGEO_CurveSubdivision.hpp"
#include "ONGEO_QuasiInterpolating.hpp"

inline double GetRoughSag(const ON_BezierCurve &bc_cur){
	ON_Line l;
	bc_cur.GetCV(0, l.from);
	bc_cur.GetCV(bc_cur.CVCount()-1, l.to);
	double sag_max = 0;
	for (int i = 1; i < bc_cur.CVCount() - 1; ++i){
		ON_3dPoint pt;
		bc_cur.GetCV(i, pt);
		double dist = l.MinimumDistanceTo(pt);
		if (sag_max < dist) sag_max = dist;
	}
	return sag_max;
}

struct SimpleTessellateOperator{
	// input
	double tolerance;

	// output
	ON_SimpleArray<double> prms;
	SimpleTessellateOperator(double tolerance_) : tolerance(tolerance_){}
	void Initialize(double tolerance_){
		tolerance = tolerance_;
		prms.Empty();
	}
	bool TestConverge(const ON_BezierCurve &bc_cur, const ON_Interval &t_int){
		double rsag = GetRoughSag(bc_cur);
		bool converged = rsag < tolerance;
		if (converged){
			prms.Append(t_int.m_t[0]);
		}
		return converged;
	}
	void EndConverge(){
		prms.Append(1.0);
	}
};

void ONGEO_TessellateBezierCurve_Simple(const ON_BezierCurve &bc, double tolerance, ON_SimpleArray<double> &prms){
	SimpleTessellateOperator to(tolerance);
	ONGEO_CurveSubdivision(to, bc);
	to.EndConverge();

	ONGEO_Swap(prms, to.prms);
}

template<typename T> struct QuasiInterpolatingTessellateOperator{
	// input
	double tolerance;

	// output
	ON_SimpleArray<double> prms;
	QuasiInterpolatingTessellateOperator(double tolerance_) : tolerance(tolerance_){}
	void Initialize(double tolerance_){
		tolerance = tolerance_;
		prms.Empty();
	}
	bool TestConverge(const ON_BezierCurve &bc_cur, const ON_Interval &t_int){
		double qi_sag = ONGEO_QI_GetError<T>(bc_cur);
		bool qi_converged = qi_sag < tolerance;
		if (qi_converged){
#if 1
			int idx_c = bc_cur.CVCount() / 2;
			double t0 = static_cast<double>(idx_c) / (bc_cur.CVCount()-1);
			double t1 = static_cast<double>(idx_c+1) / (bc_cur.CVCount()-1);
			ON_3dPoint ptcrv = bc_cur.PointAt((t0+t1)*0.5);
			ON_3dPoint ptlin = (bc_cur.PointAt(t0) + bc_cur.PointAt(t1)) * 0.5;
			if (ptcrv.DistanceTo(ptlin) > tolerance) return false;
#endif

			for (int i = 0; i < bc_cur.CVCount() - 1; ++i){
				double t = static_cast<double>(i) / (bc_cur.CVCount() - 1);
				prms.Append(t_int.ParameterAt(t));
			}
			return true;
		}
		return false;
	}
	void EndConverge(){
		prms.Append(1.0);
	}
};

void ONGEO_TessellateBezierCurve_QuasiInterpolating(const ON_BezierCurve &bc, double tolerance, ON_SimpleArray<double> &prms){
	if (bc.IsRational()){
		QuasiInterpolatingTessellateOperator<ON_4dPoint> qo(tolerance);
		ONGEO_CurveSubdivision(qo, bc);
		qo.EndConverge();

		ONGEO_Swap(prms, qo.prms);
	}else{
		QuasiInterpolatingTessellateOperator<ON_3dPoint> qo(tolerance);
		ONGEO_CurveSubdivision(qo, bc);
		qo.EndConverge();

		ONGEO_Swap(prms, qo.prms);
	}
}

