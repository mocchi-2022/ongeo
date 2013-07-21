// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "opennurbs.h"
#include "ONGEO.h"
#include <cmath>
#include <algorithm>

// Reference : H.Akima,
//   "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures".
//    Journal of the Association for Computing Machinary, Vol 17, No.4, Oct. 1970, pp.589-602.

struct ONGEO_Interpolation_Akima_Univariate::Impl{
	ON_SimpleArray<double> xa, ya, ta; // ta.Count() = ya.Count() - 4 * ydim  (外挿で求めた傾きを除外)
	ON_Interval intx;
	int ydim;
	int xiprev;
	void CalcSlopeT(){
		intx.m_t[0] = xa[2], intx.m_t[1] = *(xa.Last()-2);

		xiprev = -1;
		ON_SimpleArray<double> slopes;

		// 点間の傾きを求める。
		slopes.SetCapacity((xa.Count()-1)*ydim);
		slopes.SetCount(slopes.Capacity());
		for (int j = 2, ji = 2 * ydim; j < xa.Count() - 3; ++j, ji += ydim){
			double dx = xa[j+1] - xa[j];
			for (int i = 0; i < ydim; ++i){
				slopes[ji+i] = (ya[ydim+ji+i] - ya[ji+i]) / dx;
			}
		}

		// 外挿点、及び外挿点に対応するslopeを求める。
		int ydim2 = ydim + ydim, ydim3 = ydim * 3;

		if (xa.Count() < 7){
			ta.SetCapacity(ydim);
			ta.SetCount(ta.Capacity());
			double dx = xa[3]-xa[2];
			for (int i = 0; i < ydim; ++i){
				ta[i] = (ya[ydim3+i] - ya[ydim2+i]) / dx;
			}
			return; // 6個のときは直線近似 taに傾きを入れる。
		}

		double dx42, dx21, dx10;
		double *xal = xa.Last(), *yal = ya.Last(), *sll = slopes.Last();

		dx42 = (xa[4] - xa[2]), xa[1] = xa[3] - dx42, xa[0] = xa[2] - dx42;
		dx21 = xa[2] - xa[1], dx10 = xa[1] - xa[0];

		for (int i = 0; i < ydim; ++i){
			ya[ydim+i] = ya[ydim2+i] - dx21 * (slopes[ydim3+i] - 2.0 * slopes[ydim2+i]);
			slopes[ydim+i] = (ya[ydim2+i] - ya[ydim+i]) / dx21;

			ya[i] = ya[ydim+i] - dx10 * (slopes[ydim2+i] - 2.0 * slopes[ydim+i]);
			slopes[i] = (ya[ydim+i] - ya[i]) / dx10;
		}

		dx42 = (xal[-2] - xal[-4]), xal[-1] = xal[-3] + dx42, xal[0] = xal[-2] + dx42;
		dx21 = xal[-1] - xal[-2], dx10 = xal[0] - xal[-1];

		for (int i = 0; i < ydim; ++i){
			yal[-ydim-i] = yal[-ydim2-i] - dx21 * (sll[-ydim3-i] - 2.0 * sll[-ydim2-i]);
			sll[-ydim-i] = (yal[-ydim-i] - yal[-ydim2-i]) / dx21;

			yal[-i] = yal[-ydim-i] - dx10 * (sll[-ydim2-i] - 2.0 * sll[-ydim-i]);
			sll[-i] = (yal[-i] - yal[-ydim-i]) / dx10;
		}

		// slope tを求める。
		ta.SetCapacity(ya.Count() - 4 * ydim);
		ta.SetCount(ta.Capacity());
		for (int j = 0, ji = 0; j < xa.Count() - 4; ++j, ji += ydim){
			double *m1 = slopes.First() + ji;
			double *m2 = m1 + ydim, *m3 = m2 + ydim, *m4 = m3 + ydim;
			for (int i = 0; i < ydim; ++i){
				if (m1[i] == m2[i] && m3[i] != m4[i]) ta[ji+i] = m2[i];
				else if (m3[i] == m4[i] && m1[i] != m2[i]) ta[ji+i] = m3[i];
				else if (m2[i] == m3[i]) ta[ji+i] = m2[i];
				else{
					double m4_3 = std::abs(m4[i]-m3[i]), m2_1 = std::abs(m2[i]-m1[i]);
					ta[ji+i] = (m4_3 * m2[i] + m2_1 * m3[i]) / (m4_3 + m2_1);
				}
			}
		}
	}
};

ONGEO_Interpolation_Akima_Univariate::ONGEO_Interpolation_Akima_Univariate(){
	pimpl = new Impl();
}
ONGEO_Interpolation_Akima_Univariate::~ONGEO_Interpolation_Akima_Univariate(){
	delete pimpl;
}

bool ONGEO_Interpolation_Akima_Univariate::SetPoints(double *xa, double *ya, int num, int ydim){
	// 昇順チェック
	for (int i = 0; i < num - 1; ++i) if (xa[i+1]<=xa[i]) return false;
	if (num < 2) return false;

	pimpl->xa.SetCapacity(num+4);
	pimpl->xa.SetCount(2);
	pimpl->xa.Append(num, xa);
	pimpl->xa.SetCount(pimpl->xa.Capacity());

	pimpl->ya.SetCapacity(pimpl->xa.Capacity()*ydim);
	pimpl->ya.SetCount(2*ydim);
	pimpl->ya.Append(num*ydim, ya);
	pimpl->ya.SetCount(pimpl->ya.Capacity());

	pimpl->ydim = ydim;

	// taを計算
	pimpl->CalcSlopeT();
	return true;
}

double *ONGEO_Interpolation_Akima_Univariate::Evaluate(double x, double y[]){
	if (!y) return 0;
	int ydim = pimpl->ydim, ydim2 = pimpl->ydim * 2;
	if (pimpl->xa.Count() == 6){
		double dx = (x - pimpl->xa[2]);
		for (int i = 0; i < ydim; ++i){
			y[i] = pimpl->ta[i] * dx + pimpl->ya[ydim2+i];
		}
		return y;
	}else{
		double *x1 = 0, *y1 = 0, *t1 = 0;
		if (x <= pimpl->intx.m_t[0]){
			for (int i = 0; i < ydim; ++i){
				y[i] = pimpl->ya[ydim2 + i];
			}
			return y;
		}
		if (pimpl->intx.m_t[1] <= x){
			for (int i = 0; i < ydim; ++i){
				y[i] = pimpl->ya[pimpl->ya.Count() - 3 * ydim + i];
			}
			return y;
		}
		if (pimpl->xiprev >= 2){
			double *xs = pimpl->xa.First() + pimpl->xiprev;
			if (*xs < x && x < *(xs+1))
				x1 = xs,
				y1 = pimpl->ya.First() + pimpl->xiprev * ydim,
				t1 = pimpl->ta.First() + (pimpl->xiprev - 2) * ydim;
		}
		if (!x1){
			x1 = std::upper_bound(pimpl->xa.First(), pimpl->xa.Last() + 1, x) - 1;
			pimpl->xiprev = static_cast<int>(x1 - pimpl->xa.First());
			y1 = pimpl->ya.First() + pimpl->xiprev * ydim;
			t1 = pimpl->ta.First() + (pimpl->xiprev - 2) * ydim;
		}
		double *y2 = y1 + ydim, *t2 = t1 + ydim;
		double dx21 = *(x1+1) - *x1;
		double dx = x - *x1;
		double dx2 = dx * dx;
		double dx3 = dx2 * dx;
		for (int i = 0; i < ydim; ++i){
			double p0 = y1[i], p1 = t1[i];
			double ds = (y2[i] - y1[i]) / dx21;
			double p2 = (3.0 * ds - 2 * t1[i] - t2[i]) / dx21;
			double p3 = (t1[i] + t2[i] - 2.0 * ds) / (dx21 * dx21);
			y[i] = p0 + p1 * dx + p2 * dx2 + p3 * dx3;
		}
		return y;
	}
}

ONGEO_DECL ONGEO_Interpolation_Akima_Univariate *ONGEO_Interpolation_Akima_Univariate_New(){
	return new ONGEO_Interpolation_Akima_Univariate();
}

ONGEO_DECL void ONGEO_Interpolation_Akima_Univariate_Delete(ONGEO_Interpolation_Akima_Univariate *ths){
	delete ths;
}

ONGEO_DECL bool ONGEO_Interpolation_Akima_Univariate_SetPoints(ONGEO_Interpolation_Akima_Univariate *ths, double *xa, double *ya, int num, int ydim){
	if (!ths) return false;
	return ths->SetPoints(xa, ya, num, ydim);
}
ONGEO_DECL double *ONGEO_Interpolation_Akima_Univariate_Evaluate(ONGEO_Interpolation_Akima_Univariate *ths, double x, double y[]){
	if (!ths) return 0;
	return ths->Evaluate(x, y);
}
