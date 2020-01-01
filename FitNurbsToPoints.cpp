/*
 * FitNurbsToPoints
 * Copylight (C) 2019,2020 mocchi
 * mocchi_2003@yahoo.co.jp
 * License: Boost ver.1
 */

#include "ONGEO.h"
#include <cmath>
#include <numeric>
#include <algorithm>
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

void InitializeSolver_ExactlyPassThroughPoints(int pt_cnt, int order, double *prm, double *knot, int knot_count, IEquationSolver *solver){
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
	if (order > 2){
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

// R のサイズは cpt_cnt - 2
void InitializeSolver_PassThroughPoints_LSQ_EndPointsConstraint(int pt_cnt, int order, int cpt_cnt, double *prm, double *knot, int knot_count, const ON_3dPoint *pts, ON_3dPoint *R, IEquationSolver *solver){
	int n = cpt_cnt - 1;
	int m = pt_cnt - 1;
	int p = order - 1;
	for (int i = 0; i < p; ++i){
		knot[i] = 0;
		knot[knot_count-i-1] = 1;
	}

	double d = static_cast<double>(m + 1) / static_cast<double>(n - p + 1);
	for (double j = 1.0; j <= n - p; j += 1.0){
		double i = std::floor(j * d);
		int ii = static_cast<int>(i);
		double a = j * d - i;
		knot[static_cast<int>(j)+p-1] = (1.0 - a) * prm[ii-1] + a * prm[ii];
	}

	ON_Matrix N(m-1, n-1);
	ON_SimpleArray<double> basis(order * order); basis.SetCount(basis.Capacity());
	ON_ClassArray<ON_3dPoint> rpt(m-1);

	for (int k = 1; k <= m-1; ++k){
		int r = k;
		double t = prm[k];

		int span_index = ON_NurbsSpanIndex(order, n+1, knot, t, 0, 0);
		ON_EvaluateNurbsBasis(order, knot + span_index, t, basis.First());
		double basis_0 = 0, basis_n = 0;
		for (int c = 0; c <= n; ++c){
			double &dest = (c == 0) ? basis_0 : ((c == n) ? basis_n : N[r-1][c-1]);
			int basis_idx = c - span_index;
			if (basis_idx < 0 || basis_idx > order){
				dest = 0;
			}else{
				dest = basis[basis_idx];
			}
		}
		rpt.Append(pts[k] - pts[0] * basis_0 - pts[m] * basis_n);
	}

	for (int c = 1; c <= n-1; ++c){
		ON_3dPoint Rc(0,0,0);
		for (int r = 1; r <= m-1; ++r){
			Rc += rpt[r-1] * N[r-1][c-1];
		}
		R[c-1] = Rc;
	}

	ON_Matrix NT = N;
	NT.Transpose();
	ON_Matrix NN;
	NN.Multiply(NT, N);
	solver->Initialize(NN);
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

	InitializeSolver_ExactlyPassThroughPoints(pt_cnt, order, prm.First(), nc.m_knot, nc.KnotCount(), solver);

	ON_ClassArray<ON_3dPoint> Q(n+1); Q.SetCount(n+1);
	if (order == 2){
		for (int i = 0; i < pt_cnt; ++i) Q[i] = pts[i];
	}else{
		solver->Solve(n+1, pts, Q.First(), 1);
	}
	for (int i = 0; i <= n; ++i){
		nc.SetCV(i, Q[i]);
	}

	return true;
}

double ONGEO_FitNurbsCurveToPointArray_LeastSquare(const ON_3dPoint *pts, int pt_cnt, int order, int cpt_cnt, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, double *prm_, double *err_, ONGEO_NI_Solver solver_id){
	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	ON_NurbsCurve pol;
	ONGEO_FitNurbsCurveToPointArray(pts, pt_cnt, 2, pmethod, pol, solver_id);

	int n = pt_cnt - 1;
	ON_ClassArray<double> prm_buf;
	double *prm = prm_;
	if (!prm){
		prm_buf.SetCapacity(pt_cnt);
		prm_buf.SetCount(pt_cnt);
		prm = prm_buf.First();
	}
	for (int i = 0; i < n+1; ++i){
		prm[i] = pol.Knot(i);
	}

	nc.Create(3, false, order, cpt_cnt);
	ON_ClassArray<ON_3dPoint> R(nc.CVCount()-2); R.SetCount(R.Capacity());
	InitializeSolver_PassThroughPoints_LSQ_EndPointsConstraint(pt_cnt, order, nc.CVCount(), prm, nc.m_knot, nc.KnotCount(), pts, R, solver);
	ON_ClassArray<ON_3dPoint> P(nc.CVCount()); P.SetCount(P.Capacity());
	P[0] = pts[0];
	P[P.Count()-1] = pts[pt_cnt-1];
	solver->Solve(R.Count(), R.First(), P.First()+1, 1);
	for (int i = 0; i < nc.CVCount(); ++i) nc.SetCV(i, P[i]);

	double err_max = 0;
	for (int i = 1; i < pt_cnt - 1; ++i){
		double err = -1;
		double t = prm[i];
		if (ONGEO_NearestPointNurbsCurve_Newton(nc, pts[i], t)){
			if (prm_) prm_[i] = t;
			err = pts[i].DistanceTo(nc.PointAt(t));
			if (err_) err_[i] = err;
		}
		if (err_max < err) err_max = err;
	}
	return err_max;
}

#if 0
double CalculateBr(const ON_NurbsCurve &nc, int r, int s){
	int ord = nc.Order();
	int p = ord - 1;
	int last = r - s;
	int first = r - p;
	int off = first - 1;
	const double *knot = nc.Knot();

	ON_ClassArray<ON_3dPoint> temp(last + 2 - off);
	temp.SetCount(temp.Capacity());
	nc.GetCV(off, temp[0]);
	nc.GetCV(last+1, temp[last + 1 - off]);
	int i = first, j = last;
	int ii = 1, jj = last - off;
	while(j-i > 0){
		double alfi = (knot[r] - knot[i]) / (knot[i+ord] - knot[i]);
		double alfj = (knot[r] - knot[j]) / (knot[j+ord] - knot[j]);
		ON_3dPoint Pi, Pj;
		nc.GetCV(i, Pi); nc.GetCV(j, Pj);
		temp[ii] = (Pi - temp[ii-1] * (1.0 - alfi)) / alfi;
		temp[jj] = (Pj - temp[jj+1] * alfj) / (1.0 - alfj);
		i = i + 1, ii = ii + 1;
		j = j - 1, jj = jj - 1;
	}
	if (j - 1 < 0){
		return temp[ii-1].DistanceTo(temp[jj+1]);
	}else{
		double alfi = (knot[r] - knot[i]) / (knot[i+ord] - knot[i]);
		ON_3dPoint Pi;
		nc.GetCV(i, Pi);
		return Pi.DistanceTo(temp[ii+1] * alfi + (1.0 - alfi) * temp[ii-1]);
	}
}

bool ONGEO_FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, double tolerance, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, ONGEO_NI_Solver solver_id){
	if (tolerance <= 0) return ONGEO_FitNurbsCurveToPointArray(pts, pt_cnt, order, pmethod, nc, solver_id);

	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	ON_SimpleArray<double> basis(order * order); basis.SetCount(basis.Capacity());

	ON_NurbsCurve nc_tmp;
	ONGEO_FitNurbsCurveToPointArray(pts, pt_cnt, 2, pmethod, nc_tmp, solver_id);

	int n = pt_cnt - 1;
	ON_SimpleArray<double> prm(n+1); prm.SetCount(n+1);
	ON_SimpleArray<double> err(n+1); err.SetCount(n+1);
	for (int i = 0; i < n+1; ++i){
		prm[i] = nc_tmp.Knot(i);
		err[i] = 0;
	}

	int p = order - 1;
	// ノット除去
	for (int deg = 1; deg <= p; ++deg){
		ON_NurbsCurve nc_tmp2;
		nc_tmp2.Create(3, false, deg+1, nc_tmp.CVCount());
		// ノット除去
		while(true){
			// 最小の Br 値と、そのBrに対応するノット番号を取得
			double Br_min = -1, err_min = -1;
			const double *knot = nc_tmp.Knot();
			int s = nc_tmp.KnotMultiplicity(0);
			int r_min = -1, r = s-1;
			s = nc_tmp.KnotMultiplicity(r+1);
			r += s;
			while(s > 0){
				double Br = CalculateBr(nc_tmp, r, s);

				int p_s = deg + s;
				s = nc_tmp.KnotMultiplicity(r+1);
				if (r + s >= nc_tmp.KnotCount()-1) break;
#if 0
				nc_tmp2 = nc_tmp;
				nc_tmp2.RemoveSpan(r-(p-1));
				double Br;
				if ((p_s & 1) == 0){
					int k = p_s / 2;
					double a = (knot[r] - knot[r-k]) / (knot[r-k+p+1] - knot[r-k]);
					ON_3dPoint P, P1r, P1l;
					nc_tmp.GetCV(r-k, P);
					nc_tmp2.GetCV(r-k+1, P1r);
					int rk_1 = r-k-1;
					nc_tmp2.GetCV(rk_1 >= 0 ? rk_1 : 0, P1l);
					Br = P.DistanceTo(P1r * a + P1l * (1.0-a));
				}else{
					int k = (p_s + 1) / 2;
					ON_3dPoint P1r, P1l;
					nc_tmp2.GetCV(r-k+1, P1r);
					nc_tmp2.GetCV(r-k, P1l);
					Br = P1r.DistanceTo(P1l);
				}
#endif
				if (Br_min < 0 || Br_min > Br) Br_min = Br, r_min = r;
				r += s;
			}
			if (Br_min == -1) break;
			{
				int span_index = ON_NurbsSpanIndex(order, nc_tmp.CVCount(), knot, knot[r_min], 0, 0);
				ON_EvaluateNurbsBasis(order, knot + span_index, knot[r_min], basis.First());
				int s_min = nc_tmp.KnotMultiplicity(r_min);
				int p_s = deg + s_min;
				if ((p_s & 1) == 0){
					int k = p_s / 2;
					err_min = basis[r_min - k - span_index+1] * Br_min;
				}else{
					int k = (p_s + 1) / 2;
					double a = (knot[r_min] - knot[r_min-k+1])/(knot[r_min-k+deg+2] - knot[r_min-k+1]);
					err_min = (1.0 - a) * basis[r_min - k + 1 - span_index+1] * Br_min;
				}
			}
			if (err_min > tolerance) break;
			int idx_s = ON_NurbsSpanIndex(order, nc_tmp.CVCount(), knot, knot[r_min], 0, 0);
			int idx_n = idx_s + 2 * deg - 1;

			ptrdiff_t diff_s = std::lower_bound(prm.First(), prm.Last()+1, knot[idx_s]) - prm.First();
			ptrdiff_t diff_e = std::upper_bound(prm.First(), prm.Last()+1, knot[idx_n]) - prm.First();
			double *es = err.First() + diff_s, *ee = err.First() + diff_e;
	
			bool knot_removable = true;
			for (double *e = es; e != ee; ++e){
				if (*e + err_min < tolerance) continue;
				knot_removable = false;
				break;
			}
			if (knot_removable){
				for (double *e = es; e != ee; ++e) *e += err_min;
				nc_tmp.RemoveSpan(r_min-(order-2));
			}else break;
		}
		if (deg == p){
			nc = nc_tmp;
			break;
		}

		// 次数を1上げて最小二乗法で曲線を再生成
		nc_tmp.Create(3, false, deg+2, nc_tmp.CVCount());
		ON_ClassArray<ON_3dPoint> R(nc_tmp.CVCount()-2); R.SetCount(R.Capacity());
		InitializeSolver_PassThroughPoints_LSQ_EndPointsConstraint(pt_cnt, deg+2, nc_tmp.CVCount(), prm, nc_tmp.m_knot, nc_tmp.KnotCount(), pts, R, solver);
		ON_ClassArray<ON_3dPoint> P(nc_tmp.CVCount()); P.SetCount(P.Capacity());
		P[0] = pts[0];
		P[P.Count()-1] = pts[pt_cnt-1];
		solver->Solve(R.Count(), R.First(), P.First()+1, 1);
		for (int i = 0; i < nc_tmp.CVCount(); ++i) nc_tmp.SetCV(i, P[i]);

		double err_max = 0;
		for (int i = 1; i < pt_cnt - 1; ++i){
			double t = prm[i];
			if (ONGEO_NearestPointNurbsCurve_Newton(nc_tmp, pts[i], t)){
				prm[i] = t;
				err[i] = pts[i].DistanceTo(nc_tmp.PointAt(t));
			}
			if (err_max < err[i]) err_max = err[i];
		}
	}

	return true;
}
#endif


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

	InitializeSolver_ExactlyPassThroughPoints(n + 1, u_order, prm_u.First(), nf.m_knot[0], nf.KnotCount(0), solver);

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

	InitializeSolver_ExactlyPassThroughPoints(m + 1, v_order, prm_v.First(), nf.m_knot[1], nf.KnotCount(1), solver);

	ON_ClassArray<ON_3dPoint> Q(m+1); Q.SetCount(m+1);
	for (int i = 0; i <= n; ++i){
		solver->Solve(m+1, R + i * (m+1), Q.First(), 1);
		for (int j = 0; j <= m; ++j){
			nf.SetCV(i, j, Q[j]);
		}
	}
	return true;
}

