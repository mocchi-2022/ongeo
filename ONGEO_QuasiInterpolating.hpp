#define NOMINMAX
#include "ONGEO.h"
#include <cmath>
#include <stack>
#include <map>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>


// Reference : Yohan D.Fougerolle, Sandrine Lanquetin, Marc Neveu, and Thierry Lauthelier,
//   "A Geometric Algorithm for Ray/Bezier Surfaces Intersection using Quasi-interpolating Control Net".
//    Signal Image Technology and Internet Based Systems, 2008.


// ======
// === 曲線、曲面に共通な処理 ===
// ======
inline double ONGEO_QI_GetZInf(int dim){
	if (dim == 1) return 0;
	if (dim == 2) return 0.0625; // 1/16
	return static_cast<double>(((dim + 1) / 2) * (dim / 2))/(2.0*static_cast<double>(dim)) - 0.25;
}

inline void ONGEO_QI_GetQuasiInterpCP(const ON_3dPoint &pp, const ON_3dPoint &pc, const ON_3dPoint &pn, ON_3dPoint &po){
	po.x = (pp.x + pc.x * 2.0 + pn.x) * 0.25;
	po.y = (pp.y + pc.y * 2.0 + pn.y) * 0.25;
	po.z = (pp.z + pc.z * 2.0 + pn.z) * 0.25;
}

inline void ONGEO_QI_GetQuasiInterpCP(const ON_4dPoint &pp, const ON_4dPoint &pc, const ON_4dPoint &pn, ON_4dPoint &po){
	po.x = (pp.x + pc.x * 2.0 + pn.x) * 0.25;
	po.y = (pp.y + pc.y * 2.0 + pn.y) * 0.25;
	po.z = (pp.z + pc.z * 2.0 + pn.z) * 0.25;
	po.w = (pp.w + pc.w * 2.0 + pn.w) * 0.25;
}

inline double ONGEO_QI_GetDelta2LengthSquared(const ON_3dPoint &pp, const ON_3dPoint &pc, const ON_3dPoint &pn){
	double d, dd = 0;
	d = pp.x - pc.x * 2.0 + pn.x, d *= d, dd += d;
	d = pp.y - pc.y * 2.0 + pn.y, d *= d, dd += d;
	d = pp.z - pc.z * 2.0 + pn.z, d *= d, dd += d;
	return dd;
//	return ON_3dVector(pp - pc * 2.0 + pn).LengthSquared();
}

inline double ONGEO_QI_GetDelta2LengthSquared(const ON_4dPoint &pp, const ON_4dPoint &pc, const ON_4dPoint &pn){
	double d, dd = 0;
	d = pp.x - pc.x * 2.0 + pn.x, d *= d, dd += d;
	d = pp.y - pc.y * 2.0 + pn.y, d *= d, dd += d;
	d = pp.z - pc.z * 2.0 + pn.z, d *= d, dd += d;
	d = pp.w - pc.w * 2.0 + pn.w, d *= d, dd += d;
	return dd;
}

template <typename T, typename S> struct ONGEO_QI_PointDimTraits{
};
// ======
// === 曲線用の処理 ===
// ======
template <> struct ONGEO_QI_PointDimTraits<ON_3dPoint, ON_BezierCurve>{
	static const int DIM = 3;
	static double GetErrFromErr2(const ON_BezierCurve &bez, double errmax2){
		return ONGEO_QI_GetZInf(bez.m_order-1) * std::sqrt(errmax2);
	}
};

template <> struct ONGEO_QI_PointDimTraits<ON_4dPoint, ON_BezierCurve>{
	static const int DIM = 4;
	static double GetErrFromErr2(const ON_BezierCurve &bez, double errmax2){
		double wmin = bez.Weight(0);
		for (int i = 0; i < bez.m_order; ++i){
			if (wmin > bez.Weight(i)) wmin = bez.Weight(i);
		}
		return ONGEO_QI_PointDimTraits<ON_3dPoint, ON_BezierCurve>::GetErrFromErr2(bez, errmax2) / wmin;
	}
};

template <typename T> void ONGEO_QI_GetQuasiInterpControlPoints(const ON_BezierCurve &bez, ON_SimpleArray<T> &qicp){
	ON_SimpleArray<T> work_cp;
    work_cp.SetCapacity(bez.m_order);
	work_cp.SetCount(work_cp.Capacity());

	qicp.SetCapacity(work_cp.Count());
	qicp.SetCount(qicp.Capacity());
	for (int i = 0; i < bez.m_order; ++i){
		bez.GetCV(i, j, work_cp[i]);
	}
	for (int i = 1; i < bez.m_order-1; ++i){
		ONGEO_QI_GetQuasiInterpCP(work_cp[i-1], work_cp[i], work_cp[i+1], qicp[i]);
	}
	qicp[0] = work_cp[0];
	qicp[bez.m_order-1] = work_cp[bez.m_order-1];
}

template <typename T> double ONGEO_QI_GetError(const ON_BezierCurve &bez){
	double errmax2 = 0;
	const int ord = bez.m_order;
	ON_SimpleArray<T> cp(ord);
	for (int i = 0; i < ord; ++i){
		ONGEO_GetCVHomogeneous(bez, i, cp[i]);
	}
	for (int i = 1; i < ord-1; ++i){
		double d = ONGEO_QI_GetDelta2LengthSquared(cp[i-1], cp[i], cp[i+1]);
		if (errmax2 < d) errmax2 = d;
	}

	return ONGEO_QI_PointDimTraits<T, ON_BezierCurve>::GetErrFromErr2(bez, errmax2);
}

// ======
// === 曲面用の処理 ===
// ======
template <> struct ONGEO_QI_PointDimTraits<ON_3dPoint, ON_BezierSurface>{
	static const int DIM = 3;
	static double GetErrFromErr2(const ON_BezierSurface &bez, double errmax2u, double errmax2v){
		return ONGEO_QI_GetZInf(bez.m_order[0]-1) * std::sqrt(errmax2u) + ONGEO_QI_GetZInf(bez.m_order[1]-1) * std::sqrt(errmax2v);
	}
};

template <> struct ONGEO_QI_PointDimTraits<ON_4dPoint, ON_BezierSurface>{
	static const int DIM = 4;
	static double GetErrFromErr2(const ON_BezierSurface &bez, double errmax2u, double errmax2v){
		double wmin = bez.Weight(0, 0);
		for (int j = 0; j < bez.m_order[1]; ++j){
			for (int i = 0; i < bez.m_order[0]; ++i){
				if (wmin > bez.Weight(i, j)) wmin = bez.Weight(i, j);
			}
		}
		return ONGEO_QI_PointDimTraits<ON_3dPoint, ON_BezierSurface>::GetErrFromErr2(bez, errmax2u, errmax2v) / wmin;
	}
};

template <typename T> void ONGEO_QI_GetQuasiInterpControlPoints(const ON_BezierSurface &bez, ON_SimpleArray<T> &qicp){
	ON_SimpleArray<T> work_cp;
    work_cp.SetCapacity(bez.m_order[0] * bez.m_order[1]);
	work_cp.SetCount(work_cp.Capacity());
	// U方向
	ON_SimpleArray<T> work_u(bez.m_order[0]);
	for (int j = 0, ji = 0; j < bez.m_order[1]; ++j, ji += bez.m_order[0]){
		for (int i = 0; i < bez.m_order[0]; ++i){
			bez.GetCV(i, j, work_u[i]);
		}
		for (int i = 1; i < bez.m_order[0]-1; ++i){
			ONGEO_QI_GetQuasiInterpCP(work_u[i-1], work_u[i], work_u[i+1], work_cp[ji+i]);
		}
		work_cp[ji] = work_u[0], work_cp[ji+bez.m_order[0]-1] = work_u[bez.m_order[0]-1];
	}
//	ON_SimpleArray<T> work_v(bez.m_order[1]);

    qicp.SetCapacity(work_cp.Count());
	qicp.SetCount(qicp.Capacity());
	// V方向
	for (int i = 0; i < bez.m_order[0]; ++i){
		int j, ji;
		for (j = 1, ji = bez.m_order[0]; j < bez.m_order[1]-1; ++j, ji += bez.m_order[0]){
			ONGEO_QI_GetQuasiInterpCP(work_cp[ji+i-bez.m_order[0]], work_cp[ji+i], work_cp[ji+i+bez.m_order[0]], qicp[ji+i]);
		}
		qicp[i] = work_cp[i], qicp[ji+i] = work_cp[ji+i];
	}
}

template <typename T> double ONGEO_QI_GetError(const ON_BezierSurface &bez){
	// U方向
	double errmax2u = 0;
	const int ord_u = bez.m_order[0], ord_v = bez.m_order[1];
	ON_SimpleArray<T> cp(ord_u * ord_v);
	for (int j = 0, ji = 0; j < ord_v; ++j, ji += ord_u){
		for (int i = 0; i < ord_u; ++i){
			ONGEO_GetCVHomogeneous(bez, i, j, cp[ji+i]);
		}
		for (int i = 1; i < ord_u-1; ++i){
			double d = ONGEO_QI_GetDelta2LengthSquared(cp[ji+i-1], cp[ji+i], cp[ji+i+1]);
			if (errmax2u < d) errmax2u = d;
		}
	}
	// V方向
	double errmax2v = 0;
	for (int i = 0; i < ord_u; ++i){
		for (int j = 1, ji = ord_u; j < ord_v-1; ++j, ji += ord_u){
			double d = ONGEO_QI_GetDelta2LengthSquared(cp[ji-ord_u], cp[ji], cp[ji+ord_u]);
			if (errmax2v < d) errmax2v = d;
		}
	}

	return ONGEO_QI_PointDimTraits<T, ON_BezierSurface>::GetErrFromErr2(bez, errmax2u, errmax2v);
}
