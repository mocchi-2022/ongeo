/*
 * FitNurbsToPoints
 * Copylight (C) 2019 mocchi
 * mocchi_2003@yahoo.co.jp
 * License: Boost ver.1
 */

#include "ONGEO.h"
#include <cmath>
#include <numeric>
#ifdef max
#undef max
#endif

// Reference : Piegl, L. and Tiller, W.
//   "The NURBS Book", Second edition.
//    Springer-Verlag Berlin, Heidelberg, 1997.

// pp.361 - 453 : Chapter 9, "Curve and Surface Fitting"

namespace{
bool CreateParameter(const ON_3dPoint *pts, int pt_cnt, int pt_stride, ONGEO_NI_ParameterMethod pmethod, double *prm){
	int n = pt_cnt - 1;
	switch(pmethod){
		case ONGEO_NI_EquallySpaced: {
			for (int i = 0; i <= n; ++i){
				prm[i] = static_cast<double>(i) / static_cast<double>(n);
			}
			break;
		}
		case ONGEO_NI_ChordLength: {
			double chord_length = 0;
			ON_3dPoint pt_prv = pts[0];
			prm[0] = 0;
			for (int i = 1; i <= n; ++i){
				ON_3dPoint pt_cur = pts[i*pt_stride];
				chord_length += pt_cur.DistanceTo(pt_prv);
				prm[i] = chord_length;
				pt_prv = pt_cur;
			}
			if (chord_length > 0){
				for (int i = 1; i <= n; ++i){
					prm[i] /= chord_length;
				}
			}
			break;
		}
		case ONGEO_NI_Centripetal: {
			double sq_chord_length = 0;
			ON_3dPoint pt_prv = pts[0];
			prm[0] = 0;
			for (int i = 1; i <= n; ++i){
				ON_3dPoint pt_cur = pts[i*pt_stride];
				sq_chord_length += std::sqrt(pt_cur.DistanceTo(pt_prv));
				prm[i] = sq_chord_length;
				pt_prv = pt_cur;
			}
			if (sq_chord_length > 0){
				for (int i = 1; i <= n; ++i){
					prm[i] /= sq_chord_length;
				}
			}
		}
	}
	return true;
}

struct IEquationSolver{
	virtual void Initialize(const ON_Matrix &N) = 0;
	virtual void Solve(int pt_cnt, const ON_3dPoint *pts, ON_3dPoint *Q, int stride) const = 0;
};

// 通常は LU 分解を使うところであるが、 まずは有り物(OpenNurbs の逆行列を求める関数)で実現。
// 性能的に不足であるようなら LU 分解の追加も検討する。
struct Solver_ON_Matrix_Inverse : public IEquationSolver{
	ON_Matrix Ninv;
	virtual void Initialize(const ON_Matrix &N){
		Ninv = N;
		Ninv.Invert(ON_ZERO_TOLERANCE);
	}
	virtual void Solve(int pt_cnt, const ON_3dPoint *pts, ON_3dPoint *Q, int stride) const{
		int n = pt_cnt - 1;
		ON_Matrix P(n+1, 1), Qm(n+1, 1);
		for (int dim = 0; dim < 3; ++dim){
			for (int i = 0; i <= n; ++i){
				P[i][0] = pts[i][dim];
			}
			Qm.Multiply(Ninv, P);
			for (int i = 0; i <= n; ++i){
				Q[i*stride][dim] = Qm[i][0];
			}
		}
	}
};

void InitializeSolver(int pt_cnt, int order, double *prm, double *knot, int knot_count, IEquationSolver *solver){
	int n = pt_cnt - 1;
	int m = knot_count;
	int p = order - 1;
	for (int i = 0; i < p; ++i){
		knot[i] = 0;
		knot[m-i-1] = 1;
	}
	for (int j = 1; j <= n - p; ++j){
		double s = 0;
		for (int i = j; i < j + p; ++i){
			s += prm[i];
		}
		knot[j+p-1] = s / static_cast<double>(p);
	}
	ON_Matrix N(n+1, n+1);
	for (int j = 0; j < n + 1; ++j){
		for (int i = 0; i < n + 1; ++i){
			N[j][i] = 0;
		}
	}
	N[0][0] = 1; N[n][n] = 1;

	ON_SimpleArray<double> basis(order * order); basis.SetCount(basis.Capacity());
	for (int j = 1; j < n; ++j){
		double t = prm[j];
		int span_index = ON_NurbsSpanIndex(order, n+1, knot, t, 0, 0);
		for (int i = 0; i < span_index; ++i){
			N[j][i] = 0;
		}
		ON_EvaluateNurbsBasis(order, knot + span_index, t, basis.First());
		for (int i = span_index; i < span_index + order; ++i){
			N[j][i] = basis[i-span_index];
		}
		for (int i = span_index + order; i <= n; ++i){
			N[j][i] = 0;
		}
	}

	solver->Initialize(N);
}
}

bool ONGEO_FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, ONGEO_NI_Solver solver_id){
	if (pt_cnt < order || !pts) return false;
	int n = pt_cnt - 1;

	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	nc.Create(3, false, order, n+1);

	ON_SimpleArray<double> prm(n+1); prm.SetCount(n+1);
	CreateParameter(pts, pt_cnt, 1, pmethod, prm);


	InitializeSolver(pt_cnt, order, prm.First(), nc.m_knot, nc.KnotCount(), solver);

	ON_ClassArray<ON_3dPoint> Q(n+1); Q.SetCount(n+1);
	solver->Solve(n+1, pts, Q.First(), 1);
	for (int i = 0; i <= n; ++i){
		nc.SetCV(i, Q[i]);
	}

	return true;
}

bool ONGEO_FitNurbsSurfaceToPointGrid(const ON_3dPoint *pts, int u_count, int v_count, int u_order, int v_order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsSurface &nf, ONGEO_NI_Solver solver_id){
	if (u_count < u_order || v_count < v_order || !pts) return false;
	int n = u_count - 1, m = v_count - 1;

	ON_NurbsCurve nc;
	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	nf.Create(3, false, u_order, v_order, n+1, m+1);

	ON_SimpleArray<double> prm_buf(std::max(n,m)+1); prm_buf.SetCount(prm_buf.Capacity());

	// U方向
	ON_SimpleArray<double> prm_u(n+1); prm_u.SetCount(n+1);
	CreateParameter(pts, n+1, 1, pmethod, prm_u.First());
	for (int j = 1; j <= m; ++j){
		CreateParameter(pts + j * (n+1), n+1, 1, pmethod, prm_buf.First());
		for (int i = 0; i <= n; ++i) prm_u[i] += prm_buf[i];
	}
	for (int i = 0; i <= n; ++i) prm_u[i] /= static_cast<double>(m+1);

	InitializeSolver(n + 1, u_order, prm_u.First(), nf.m_knot[0], nf.KnotCount(0), solver);

	ON_ClassArray<ON_3dPoint> R((n+1) * (m+1)); R.SetCount(R.Capacity());
	for (int j = 0; j <= m; ++j){
		solver->Solve(n+1, pts + j * (n+1), R.First() + j, m+1);
	}

	// V方向
	ON_SimpleArray<double> prm_v(m+1); prm_v.SetCount(m+1);
	CreateParameter(pts, m+1, n+1, pmethod, prm_v.First());
	for (int i = 1; i <= n; ++i){
		CreateParameter(pts + i, m+1, n+1, pmethod, prm_buf.First());
		for (int j = 0; j <= m; ++j) prm_v[j] += prm_buf[j];
	}
	for (int j = 0; j <= m; ++j) prm_v[j] /= static_cast<double>(n+1);

	InitializeSolver(m + 1, v_order, prm_v.First(), nf.m_knot[1], nf.KnotCount(1), solver);

	ON_ClassArray<ON_3dPoint> Q(m+1); Q.SetCount(m+1);
	for (int i = 0; i <= n; ++i){
		solver->Solve(m+1, R + i * (m+1), Q.First(), 1);
		for (int j = 0; j <= m; ++j){
			nf.SetCV(i, j, Q[j]);
		}
	}
	return true;
}

