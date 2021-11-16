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
			break;
		}
		case ONGEO_NI_Manual: {
			break;
		}
	}
	return true;
}

struct IEquationSolver{
	virtual void Initialize(const ON_Matrix &N) = 0;
	void Solve(int pt_cnt, const ON_3dPoint *pts, ON_3dPoint *Q, int stride) const {
		Solve(pt_cnt, pts, 0, 0, Q, stride);
	}
	virtual void Solve(int pt_cnt, const ON_3dPoint *pts, int *constraint_indices, int constraint_cnt, ON_3dPoint *Q, int stride) const = 0;
};

// 通常は LU 分解を使うところであるが、 まずは有り物(OpenNurbs の逆行列を求める関数)で実現。
// 性能的に不足であるようなら LU 分解の追加も検討する。
struct Solver_ON_Matrix_Inverse : public IEquationSolver{
	ON_Matrix Ninv;
	virtual void Initialize(const ON_Matrix &N){
		Ninv = N;
		Ninv.Invert(ON_ZERO_TOLERANCE);
	}
	virtual void Solve(int pt_cnt, const ON_3dPoint *pts, int *constraint_indices, int constraint_cnt, ON_3dPoint *Q, int stride) const{
		int n = pt_cnt - 1;
		ON_Matrix P(n+1, 1), Qm(n+1, 1);
		for (int dim = 0; dim < 3; ++dim){
			for (int i = 0; i <= n; ++i){
				int ii = constraint_indices ? constraint_indices[i] : i;
				P[i][0] = pts[ii][dim];
			}
			Qm.Multiply(Ninv, P);
			for (int i = 0; i <= n; ++i){
				Q[i*stride][dim] = Qm[i][0];
			}
		}
	}
};

void InitializeSolver_ExactlyPassThroughPoints(int pt_cnt, int order, double *prm, int *constraint_indices, int constraint_cnt, double *knot, int knot_count, IEquationSolver *solver){
	int n = (constraint_indices) ? constraint_cnt - 1 : pt_cnt - 1;
	int m = knot_count;
	int p = order - 1;
	int prm_first_idx = 0, prm_last_idx = pt_cnt - 1;
	if (constraint_indices) prm_first_idx = constraint_indices[0], prm_last_idx = constraint_indices[n];
	for (int i = 0; i < p; ++i){
		knot[i] = prm[prm_first_idx];
		knot[m-i-1] = prm[prm_last_idx];
	}
	for (int j = 1; j <= n - p; ++j){
		double s = 0;
		for (int i = j; i < j + p; ++i){
			s += prm[constraint_indices ? constraint_indices[i] : i];
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
			double t = prm[constraint_indices ? constraint_indices[j] : j];
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
		double jd = j * d;
		double i = std::floor(jd);
		int ii = static_cast<int>(i);
		double a = jd - i;
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

bool FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, int *constraint_indices, int constraint_cnt, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, double *prm_, ONGEO_NI_Solver solver_id){
	if (pt_cnt < order || !pts) return false;
	int n = (constraint_indices) ? constraint_cnt - 1 : pt_cnt - 1;

	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	nc.Create(3, false, order, n+1);

	ON_SimpleArray<double> prm;
	if (prm_ == 0){
		if (pmethod == ONGEO_NI_Manual) return false;
		prm.SetCapacity(pt_cnt); prm.SetCount(pt_cnt);
		prm_ = prm;
	}
	CreateParameter(pts, pt_cnt, 1, pmethod, prm_);

	InitializeSolver_ExactlyPassThroughPoints(pt_cnt, order, prm_, constraint_indices, constraint_cnt, nc.m_knot, nc.KnotCount(), solver);

	ON_ClassArray<ON_3dPoint> Q(n+1); Q.SetCount(n+1);
	if (order == 2){
		if (constraint_indices){
			for (int i = 0; i < constraint_cnt; ++i) Q[i] = pts[constraint_indices[i]];
		}else{
			for (int i = 0; i < pt_cnt; ++i) Q[i] = pts[i];
		}
	}else{
		solver->Solve(n+1, pts, constraint_indices, constraint_cnt, Q.First(), 1);
	}
	for (int i = 0; i <= n; ++i){
		nc.SetCV(i, Q[i]);
	}

	return true;
}

}

bool ONGEO_FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, double *prm_, ONGEO_NI_Solver solver_id){
	return FitNurbsCurveToPointArray(pts, pt_cnt, order, 0, 0, pmethod, nc, prm_, solver_id);
}

#if 0
bool ONGEO_FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, ONGEO_NI_Solver solver_id){
	if (pt_cnt < order || !pts) return false;
	int n = pt_cnt - 1;

	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	nc.Create(3, false, order, n+1);

	ON_SimpleArray<double> prm;
	if (prm_ == 0){
		prm.SetCapacity(n+1); prm.SetCount(n+1);
		prm_ = prm;
	}
	CreateParameter(pts, pt_cnt, 1, pmethod, prm_);

	InitializeSolver_ExactlyPassThroughPoints(pt_cnt, order, prm_, nc.m_knot, nc.KnotCount(), solver);

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
#endif

double ONGEO_FitNurbsCurveToPointArray_LeastSquare(const ON_3dPoint *pts, int pt_cnt, int order, int cpt_cnt, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, double *prm_, double *err_, ONGEO_NI_Solver solver_id){
	Solver_ON_Matrix_Inverse solver_omi;
	IEquationSolver *solver = 0;
	if (solver_id == ONGEO_NI_Solver_ON_Matrix_Inverse){
		solver = &solver_omi;
	}

	ON_ClassArray<double> prm_buf;
	double *prm = prm_;
	if (!prm){
		prm_buf.SetCapacity(pt_cnt);
		prm_buf.SetCount(pt_cnt);
		prm = prm_buf.First();
	}
	int n = pt_cnt - 1;
	ON_NurbsCurve pol;
	ONGEO_FitNurbsCurveToPointArray(pts, pt_cnt, 2, pmethod, pol, 0, solver_id);

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

// Reference : Razdan, A.
//    Knot Placement for B-Spline Curve Approximation
//    Arizona State University, 1999.
namespace{
// pt_indices_marked は呼び出し元で pt_cnt 分の領域を確保すること。戻り値の数だけ書き込まれる。
int EstimateNumberOfPoints(const ON_3dPoint *pts, int pt_cnt, double sigma, int *pt_indices_marked){
	if (pt_cnt < 2) return 0;
	double arc_length = 0;
	int n = 0;
	if (pt_indices_marked) pt_indices_marked[n] = 0;
	n++;
	ON_3dPoint pt_marked = pts[0], pt_prev = pts[0];
	for (int i = 1; i < pt_cnt - 1; ++i){
		ON_3dPoint pt_cur = pts[i];
		double chord_length = pt_cur.DistanceTo(pt_marked);
		arc_length += pt_cur.DistanceTo(pt_prev);
		pt_prev = pt_cur;
		double ratio = arc_length / chord_length;
		if (sigma < ratio){
			if (pt_indices_marked) pt_indices_marked[n] =i;
			arc_length = 0;
			pt_marked = pt_cur;
			n++;
		}
	}
	if (pt_indices_marked) pt_indices_marked[n] = pt_cnt - 1;
	n++;
	return n;
}

// arc_integral 及び kappa_integral は呼び出し元で pt_cnt - 1 の領域を確保すること。
void ArcKappaIntegral(const ON_3dPoint *pts, int pt_cnt, double *arc_integral, double *kappa_integral){
	if (pt_cnt <= 2) return;
	ON_3dPoint pt1 = pts[0], pt2 = pts[1];
	kappa_integral[0] = 0;
	arc_integral[0] = pt2.DistanceTo(pt1);
	for (int i = 1; i < pt_cnt - 1; ++i){
		ON_3dPoint pt3 = pts[i+1];
		double r = ON_Circle(pt1, pt2, pt3).Radius();
		double arc = pt3.DistanceTo(pt2);
		arc_integral[i] = arc + arc_integral[i-1];

		// pt1、pt2、pt3 が同一線上にある場合は曲率 0 だが、逆関数を求めづらいため、トレランスを足す。
		double kappa = (r == 0) ? ON_ZERO_CURVATURE_TOLERANCE : 1.0 / r;
		pt1 = pt2, pt2 = pt3;
		kappa_integral[i] = kappa + kappa_integral[i-1];
	}
}

// kappa_integral は KappaIntegral の計算結果を渡す。 pt_cnt - 1 個の値が入っているはず。
// pt_indices_selected は呼び出し元で pt_cnt 分の領域を確保すること。戻り値の数だけ書き込まれる。
int KappaParametrization(const ON_3dPoint *pts, int pt_cnt, int enp, const double *kappa_integral, int *pt_indices_selected){
	if (pt_cnt < 2) return 0;
	if (enp == 2){
		pt_indices_selected[0] = 0;
		pt_indices_selected[1] = pt_cnt - 1;
		return 2;
	}
	int ki_count = pt_cnt - 1;
	int n = 0;

	ON_Interval ki_range(kappa_integral[0], kappa_integral[pt_cnt - 2]);

	int ki_cursor = 0;
	pt_indices_selected[n++] = ki_cursor;
	for (int j = 1; j < enp; ++j){
		double normalized_interval = static_cast<double>(j) / static_cast<double>(enp - 1);
		double ki = ki_range.ParameterAt(normalized_interval);
		int ki_cnt_search = ki_count - ki_cursor;
		int ki_inv = ON_SearchMonotoneArray(kappa_integral + ki_cursor, ki_cnt_search, ki);
		if (ki_inv < 0 || ki_inv == ki_cnt_search) return -1;
		if (ki_inv > 0) ki_cursor += ki_inv, pt_indices_selected[n++] = ki_cursor;
	}
	if (pt_indices_selected[n-1] != pt_cnt - 1){
		pt_indices_selected[n-1] = pt_cnt - 1;
	}
	return n;
}

// pt_indices_added は呼び出し元で pt_cnt - pt_indices_cnt_already_selected だけ確保しておくこと。戻り値の数だけ書き込まれる。
int AdaptiveKnotSequenceGeneration(const ON_3dPoint *pts, int pt_cnt, const double *prm, const int *pt_indices_already_selected, int pt_indices_cnt_already_selected, double rho, int *pt_indices_added){
	if (pt_indices_cnt_already_selected == 0) return 0;
	double rho_inv = 1.0 / rho;
	double m = pt_indices_cnt_already_selected - 1;
	int n = 0;

	ON_3dPoint pt1 = pts[pt_indices_already_selected[0]];
	ON_3dPoint pt2 = pts[pt_indices_already_selected[1]];
	double delta1 = pt2.DistanceTo(pt1);
	int ridx_prev = -1;
	for (int i = 2; i <= m; ++i){
		ON_3dPoint pt3 = pts[pt_indices_already_selected[i]];
		double delta2 = pt3.DistanceTo(pt2);
		double ratio_delta = delta1 < ON_ZERO_TOLERANCE ? 0 : delta2 / delta1;
		int ridx = delta2 > delta1 ? i - 1 : i - 2;
		delta1 = delta2;
		pt2 = pt3;

		if (ratio_delta >= rho_inv && ratio_delta <= rho) continue;

		// 既に分割済のセグメントを選択した場合は飛ばす。
		// 短い、長い、短い と並んだセグメント列に到達した場合に、同一セグメントを2回分割してしまうのを防止するための措置。
		if (ridx == ridx_prev) continue;
		ridx_prev = ridx;

		// 2つの長さのバランスが悪いところを見つけたら、
		// 長い方のセグメントの両端のパラメータの中央値に最も近いインデックス点を探し、
		// 未選択の点であれば追加。選択済の点の場合は追加しない。
		double prm1 = prm[pt_indices_already_selected[ridx]];
		double prm2 = prm[pt_indices_already_selected[ridx+1]];
		double prm_mid = (prm1 + prm2) * 0.5;
		int prm_mid_idx[2];
		prm_mid_idx[0] = ON_SearchMonotoneArray(prm, pt_cnt, prm_mid);
		bool last_idx = prm_mid_idx[0] == pt_cnt - 1;

		if (!last_idx){
			prm_mid_idx[1] = (prm_mid_idx[0] < pt_cnt - 1) ? prm_mid_idx[0] + 1 : prm_mid_idx[0];
			// prm_mid に近い順に並べ替え
			double d1 = prm[prm_mid_idx[0]] - prm_mid;
			double d2 = prm_mid - prm[prm_mid_idx[1]];
			if (d1 > d2){
				std::swap(prm_mid_idx[0], prm_mid_idx[1]);
			}
		}

		const int *pidx[2];
		pidx[0] = ON_BinarySearchIntArray(prm_mid_idx[0], pt_indices_already_selected, pt_indices_cnt_already_selected);
		if (pidx[0]){
			pidx[1] = last_idx ? pidx[0] : ON_BinarySearchIntArray(prm_mid_idx[1], pt_indices_already_selected, pt_indices_cnt_already_selected);
		}

		int prm_idx_toadd;
		if (!pidx[0]) prm_idx_toadd = prm_mid_idx[0];
		else if (!pidx[1]) prm_idx_toadd = prm_mid_idx[1];
		else continue;

		if (n > 0 && ON_BinarySearchIntArray(prm_idx_toadd, pt_indices_added, n)) continue;

		pt_indices_added[n++] = prm_idx_toadd;
	}

	return n;
}
}

double ONGEO_FitNurbsCurveToPointArray_AdaptiveKnotSelection(const ON_3dPoint *pts, int pt_cnt, int order_want, double tolerance, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, double *prm_, double *err_, ON_SimpleArray<int> *selected_indices_, ONGEO_NI_Solver solver_id){
	ON_SimpleArray<double> prm;
	if (prm_ == 0){
		if (pmethod == ONGEO_NI_Manual) return -1;
		prm.SetCapacity(pt_cnt); prm.SetCount(pt_cnt);
		prm_ = prm;
	}
	CreateParameter(pts, pt_cnt, 1, pmethod, prm_);

	ON_SimpleArray<double> arc_integral(pt_cnt - 1), kappa_integral(pt_cnt - 1);
	arc_integral.SetCount(pt_cnt - 1);
	kappa_integral.SetCount(pt_cnt - 1);
	ArcKappaIntegral(pts, pt_cnt, arc_integral.First(), kappa_integral.First());

	double ai_last = (*arc_integral.Last());

	ON_SimpleArray<int> selected_indices;
	if (selected_indices_ == 0){
		selected_indices_ = &selected_indices;
	}
	selected_indices_->SetCapacity(pt_cnt); selected_indices_->SetCount(pt_cnt);

	// トレランスから大まかな sigma を求める。
	double sigma;
	{
		double seg_len = ai_last / static_cast<double>(pt_cnt - 1);
		ON_3dPoint pta(0, 0, 0), ptb(seg_len * 0.5, tolerance, 0), ptc(seg_len, 0, 0);
		sigma = ON_Arc(pta, ptb, ptc).Length() / seg_len;
	}
	// トレランスで制御したいため、最初は少なめに見積もる。
	int enp = EstimateNumberOfPoints(pts, pt_cnt, sigma, 0);
	if (enp > 50) enp /= order_want * 2;
	else if (enp > 10) enp /= order_want;

	if (enp < 2) enp = 2;
	
	int kp_cnt = KappaParametrization(pts, pt_cnt, enp, kappa_integral.First(), selected_indices_->First());
	int kp_aksg_cnt = kp_cnt;
	for(;;){
		int aksg_cnt = AdaptiveKnotSequenceGeneration(pts, pt_cnt, prm_, selected_indices_->First(), kp_aksg_cnt, 3.0, selected_indices_->First() + kp_aksg_cnt);
		if (aksg_cnt == 0) break;
		kp_aksg_cnt += aksg_cnt;
		ON_SortIntArray(ON::heap_sort, selected_indices_->First(), kp_aksg_cnt);
	}

	// 収束計算
	double err_max;
	for (;;){
		int order = (order_want > kp_aksg_cnt) ? kp_aksg_cnt : order_want;
		FitNurbsCurveToPointArray(pts, pt_cnt, order, selected_indices_->First(), kp_aksg_cnt, ONGEO_NI_Manual, nc, prm_, solver_id);

		ON_SimpleArray<double> segment_over_tol(kp_aksg_cnt-1); segment_over_tol.SetCount(kp_aksg_cnt-1);
		for (int i = 0; i < segment_over_tol.Count(); ++i) segment_over_tol[i] = 0;

		int segment_idx = 0;
		err_max = 0;
		bool need_converge = false;
		for (int i = 1; i < pt_cnt - 1; ++i){
			while((*selected_indices_)[segment_idx+1] <= i && segment_idx < kp_aksg_cnt - 1){
				++segment_idx;
			}
			double err = -1;
			double t = prm_[i];
			if (ONGEO_NearestPointNurbsCurve_Newton(nc, pts[i], t)){
				prm_[i] = t;
				err = pts[i].DistanceTo(nc.PointAt(t));
				if (err > tolerance){
					need_converge = true;
					if (segment_over_tol[segment_idx] < err) segment_over_tol[segment_idx] = err;
				}
				if (err_) err_[i] = err;
			}
			if (err_max < err) err_max = err;
		}

		if (!need_converge) break;

		int *selected_indices_to_add = selected_indices_->First() + kp_aksg_cnt;
		int add_cnt = 0;
		for (int i = 0; i < segment_over_tol.Count(); ++i){
			if (segment_over_tol[i] == 0) continue;
			ON_Interval idx_range((*selected_indices_)[i], (*selected_indices_)[i+1]);
			double divide = 0.5;
			if (segment_over_tol[i] / tolerance > 100.0){
				divide = 0.25;
			}
			int idx_prev = static_cast<int>(idx_range.m_t[0]);
			for (double nrm_p = divide; nrm_p < 1.0; nrm_p += divide){
				int midx = idx_range.ParameterAt(nrm_p);
				if (idx_prev == midx) continue;
				selected_indices_to_add[add_cnt++] = midx;
				idx_prev = midx;
			}
		}
		if (add_cnt == 0) break;
		kp_aksg_cnt += add_cnt;
		ON_SortIntArray(ON::heap_sort, selected_indices_->First(), kp_aksg_cnt);
	}
	selected_indices_->SetCount(kp_aksg_cnt);

	return err_max;
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

	InitializeSolver_ExactlyPassThroughPoints(n + 1, u_order, prm_u.First(), 0, 0, nf.m_knot[0], nf.KnotCount(0), solver);

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

	InitializeSolver_ExactlyPassThroughPoints(m + 1, v_order, prm_v.First(), 0, 0, nf.m_knot[1], nf.KnotCount(1), solver);

	ON_ClassArray<ON_3dPoint> Q(m+1); Q.SetCount(m+1);
	for (int i = 0; i <= n; ++i){
		solver->Solve(m+1, R + i * (m+1), Q.First(), 1);
		for (int j = 0; j <= m; ++j){
			nf.SetCV(i, j, Q[j]);
		}
	}
	return true;
}

