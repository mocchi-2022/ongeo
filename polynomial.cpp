// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.

#define NOMINMAX
#include <cmath>
#include <xutility>
#include <limits>
#include <vector>
#include <windows.h>
#include "ONGEO.h"

// 多項式関数に値を代入して計算する。
// @param [in] t 代入する値
// @param [in] coef 多項式の係数列
// @param [in] num 係数の個数(=次数+1)
// @return 計算結果
double ONGEO_Polynomial_Evaluate(double t, const double *coef, int num){
	double s = coef[num-1];
	double tt = t;
	for (int i = num - 2; i >= 0; --i){
		s += coef[i] * tt;
		tt *= t;
	}
	return s;
}

// 2つの多項式の積を多項式として展開する。
// @param [in] coef1 1つめの多項式の係数列
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_mul 2つの多項式の積を展開した多項式の係数列
//	([num1+num2-2]次の多項式となるため、入力側でnum1+num2-1個分の係数の領域を確保する必要がある。)
void ONGEO_Polynomial_Multiply(const double *coef1, int num1,const  double *coef2, int num2, double *coef_mul){
	int num_mul = num1 + num2 - 1;
	for (int j = 0; j < num_mul; ++j){
		coef_mul[j] = 0;
		int is = std::max(0, j - num2 + 1);
		int ie = std::min(j, num1 - 1);
		for (int i = is; i <= ie; ++i){
			coef_mul[j] += coef1[i] * coef2[j-i];
		}
	}
}

// 2つの多項式の和を求める。
// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_add 多項式1 + 多項式2を示す多項式の係数列
// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
// @return true:成功、 false:失敗
bool ONGEO_Polynomial_Add(const double *coef1, int num1, const double *coef2, int num2, double *coef_add, int num_add){
	if (num1 > num_add || num2 > num_add) return false;
	if ((!coef1 && num1) || (!coef2 && num2)) return false;
	int num_min = std::min(num1, num2);
	for (int i = 0; i < num_min; ++i){
		coef_add[num_add - i - 1] = coef1[num1 - i - 1] + coef2[num2 - i - 1];
	}
	int num_ic;
	if (num1 < num2){
		for (int i = num1; i < num2; ++i){
			coef_add[num_add - i - 1] = coef2[num2 - i - 1];
		}
		num_ic = num2;
	}else{
		for (int i = num2; i < num1; ++i){
			coef_add[num_add - i - 1] = coef1[num1 - i - 1];
		}
		num_ic = num1;
	}
	for (int i = num_ic; i < num_add; ++i) coef_add[num_add - i - 1] = 0;
	return true;
}

// 2つの多項式の差を求める。
// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_sub 多項式1 - 多項式2を示す多項式の係数列
// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
// @return true:成功、 false:失敗
bool ONGEO_Polynomial_Subtract(const double *coef1, int num1, const double *coef2, int num2, double *coef_sub, int num_sub){
	if (num1 > num_sub || num2 > num_sub) return false;
	if ((!coef1 && num1) || (!coef2 && num2)) return false;
	int num_min = std::min(num1, num2);
	for (int i = 0; i < num_min; ++i){
		coef_sub[num_sub - i - 1] = coef1[num1 - i - 1] - coef2[num2 - i - 1];
	}
	int num_ic;
	if (num1 < num2){
		for (int i = num1; i < num2; ++i){
			coef_sub[num_sub - i - 1] = -coef2[num2 - i - 1];
		}
		num_ic = num2;
	}else{
		for (int i = num2; i < num1; ++i){
			coef_sub[num_sub - i - 1] = coef1[num1 - i - 1];
		}
		num_ic = num1;
	}
	for (int i = num_ic; i < num_sub; ++i) coef_sub[num_sub - i - 1] = 0;
	return true;
}

// 多項式を一階微分する。
// @param [in] coef 多項式の係数列
// @param [in] num 多項式の係数の個数(=次数+1)
// @param [out] coef_dif 微分した多項式の係数列
//    (呼び出し元でnum-1個分の領域を確保すること
//     coefと同じ場所に結果を上書きしたい場合は、この引数として&coef[1]を渡し、本関数呼出し後にcoef[0] = 0.0とすること)
void ONGEO_Polynomial_Differential(const double *coef, int num, double *coef_dif){
	for (int i = num - 2; i >= 0; --i){
		coef_dif[i] = (num - 1 - i) * coef[i];
	}
}

namespace {
	// 係数列の絶対値の最大値で係数列を正規化したとき、高次側から順に、
	// 絶対値がON_ZERO_TOLERANCE以下となる係数を探し続け、見つかる限り0を入れる。
	int ZeroRoundSturmSequence(double *strum_line, int num){
		// 絶対値の最大値
		double absmax = 0;
		for (int i = 0; i < num; ++i){
			double a = std::abs(strum_line[i]);
			absmax = std::max(absmax, a);
		}
		int num_zero = 0;
		for (int i = 0; i < num; ++i){
			if (std::abs(strum_line[i]) / absmax < ON_ZERO_TOLERANCE){
				++num_zero;
				strum_line[i] = 0;
			}else break;
		}
		return num_zero;
	}
}

// 多項式の係数列からSturm列を生成する。
// @param coef [in] 多項式の係数列 (次数の大きい順)
// @param num [in] 多項式の係数の数
// @param s [out] Sturm列 (呼び出し元でnf*nf個分の領域を確保すること)
void ONGEO_Polynomial_CreateSturmSequence(const double *coef, int num, double *strum){
	int max_nz = 0;
//	for (max_nz = 0; coef[max_nz] == 0 && max_nz < num; ++max_nz);
//	for (int i = 0; i < num; ++i) strum[max_nz * num + i] = coef[i];
	for (int i = num; i < num * num; ++i) strum[i] = 0;
	for (int i = 0; i < num; ++i) strum[i] = coef[i];
	int num_zero = ZeroRoundSturmSequence(strum, num);
	if (num_zero){
		if (num_zero == num) return;
		for (int i = 0; i < num; ++i) std::swap(strum[num_zero * num + i], strum[i]);
		max_nz = num_zero;
	}
	if (max_nz == num - 1) return;

	ONGEO_Polynomial_Differential(strum + max_nz * num, num, strum + (num * (1 + max_nz)) + 1);

	for (int j = 2+max_nz, jj = (2+max_nz) * num; j < num; ++j, jj += num){
		double qx = strum[j-2+jj-2*num] / strum[j-1+jj-num];
		double qc = (strum[j-1+jj-2*num] - strum[j  +jj-num]*qx) / strum[j-1+jj-num];
		strum[jj+j-2] = -(strum[j-2+jj-2*num] - qx*strum[j-1+jj-num]);
		for (int i = j-1; i < num-1; ++i){
			strum[jj+i] = -(strum[i+jj-2*num]-(qx*strum[i+1+jj-num]+qc*strum[i+jj-num]));
		}
		strum[jj+num-1] = -(strum[jj-num-1]-qc*strum[jj-1]);
		num_zero = ZeroRoundSturmSequence(&strum[jj+j-2], num - (j-2));
		if (num_zero > 2){
			if (num_zero + j-2 >= num) return;
			for (int i = j-2; i < num; ++i) std::swap(strum[jj+i], strum[jj+(num_zero-2)*num+i]);
			j += num_zero - 2;
		}
	}
}

// Sturm列を用いて、指定した値での符号反転回数を計算する。
// @param t [in] 値
// @param sturm [in] Sturm列
// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
// @return 符号反転回数
int ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(double t, const double *strum, int num){
	int count = 0;
	double fp;
	int is = 0;
	do{
		fp = ONGEO_Polynomial_Evaluate(t, strum + is * num, num);
		++is;
	}while(fp == 0 && is < num);
	for (int i = is; i < num; ++i){
		double fi = ONGEO_Polynomial_Evaluate(t, strum + i * num, num);
		if (fi == 0) continue;
		if (fi * fp < 0) ++count;
		fp = fi;
	}
	return count;
}

// Sturm列を用いて、指定した区間内の多項式の根の数を計算する。
// @param t1 [in] 区間の下限
// @param t2 [in] 区間の上限(t1 < t2 であること)
// @param sturm [in] Sturm列
// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
// @return 解の個数(エラーの場合、負数)
int ONGEO_Polynomial_CalculateNumRoot(double t1, double t2, const double *strum, int num){
	if (t1 > t2) return -1;
	return
		ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(t1, strum, num) - 
		ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(t2, strum, num);
}

struct PascalTriangle{
	int maxdim;
	std::vector<int *> triangles;
	CRITICAL_SECTION cs;
	PascalTriangle(){
		::InitializeCriticalSection(&cs);
		Get(5);
	}
	int *Get(int dim){
		if (dim < 0) return 0;
		if (triangles.size() > static_cast<size_t>(dim)){
			return triangles[dim];
		}

		::EnterCriticalSection(&cs);
		for (size_t i = triangles.size(); i <= static_cast<size_t>(dim); ++i){
			triangles.push_back(new int[i+1]);
			int *cur = triangles.back();
			if (i == 0){
				cur[0] = 1;
				continue;
			}
			int *prev = triangles[triangles.size()-2];
			cur[0] = prev[0];
			cur[i] = prev[i-1];
			for (size_t h = 1; h < i; ++h){
				cur[h] = prev[h-1]+prev[h];
			}
		}
		::LeaveCriticalSection(&cs);
		return triangles[dim];
	}
}pt;
int *ONGEO_Polynomial_CalcPascalTriangle(int dim){
	return pt.Get(dim);
}
int ONGEO_Polynomial_CalcBinomialCoef(int m, int n){
	if (m < 0 || n < 0) return 0;
	if (m < n) return 0;
	return pt.Get(m)[n];
}


// Brent法を用いて根を求める。区間 ti1 - ti2 の中で単調増加、または単調減少であること。
// @param ti1 [in] 狭めたい区間の下限
// @param ti2 [in] 狭めたい区間の上限 ti1 < ti2であること。
// @param coef [in] 多項式の係数
// @param num [in] 多項式の係数の数
// @return 根となる値
// Reference : R.P.Brent,
//   "An algorithm with guaranteed convergence for finding a zero of a function".
//    The Computer Journal, 14(1971), pp.422-425.
double ONGEO_Polynomial_FindRootByBrentMethod(double ti1, double ti2, const double *coef, int num) {
	double a = ti1, b = ti2;
	double fa = ONGEO_Polynomial_Evaluate(a, coef, num);
	double fb = ONGEO_Polynomial_Evaluate(b, coef, num);
	double c = a, fc = fa, d, dprev;
	d = dprev = b - a;
	// d : 収束処理でbを動かす量
	for(;;){
		// 根に近い方をb,fb、遠い方をa=c,fa=fcとする。
		if (std::abs(fc) < std::abs(fb)){
			a = b, b = c, c = a;
			fa = fb, fb = fc, fc = fa;
		}
		double tol = 2.0 * DBL_EPSILON * std::abs(b) + ON_ZERO_TOLERANCE;
		double c_b = c - b;
		double m = 0.5 * c_b;
		if (std::abs(m) <= tol || fb == 0) break;

		// bisectionの場合のd、dprev
		d = dprev = m;
		if (std::abs(dprev) >= tol && std::abs(fa) > std::abs(fb)){
			double r3 = fb / fa;
			double p, q;
			if (a == c){
				// linear interpolation
				p = c_b * r3;
				q = 1.0 - r3;
			}else{
				// inverse quadratic interpolation
				double ifc = 1.0 / fc;
				double r1 = fa * ifc, r2 = fb * ifc;
				p = r3 * (c_b * r1 * (r1 - r2) - (b - a) * (r2 - 1));
				q = (r1 - 1) * (r2 - 1) * (r3 - 1);
			}
			if (p > 0) q = -q;
			else p = -p;
			if (2.0 * p < 1.5 * c_b * q - std::abs(tol * q) && p < std::abs(0.5 * dprev * q)){
				// 上記条件に合う場合は、dとdprevをlinear interpolation、
				// またはinverse quadratic interpolation由来のものに書き換える。
				dprev = d;
				d = p / q;
			}
		}
		// a, fa, b, fbを更新
		a = b, fa = fb;
		b += (std::abs(d) > tol) ? d : ((m > 0) ? tol : -tol);
		fb = ONGEO_Polynomial_Evaluate(b, coef, num);
		if (fb * fc > 0){
			c = a, fc = fa;
			dprev = b - a;
		}
	}
	return b;
}
#if 0
// Smaleの基準を用いて根の範囲を狭める。区間 ti1 - ti2 の中で単調増加、または単調減少であること。
// @param ti1 [in] 狭めたい区間の下限
// @param ti2 [in] 狭めたい区間の上限 ti1 < ti2であること。
// @param coef [in] 多項式の係数
// @param num [in] 多項式の係数の数
// @param to1 [out] 狭めた区間の下限
// @param to2 [out] 狭めた区間の上限
// @return true:成功、 false:失敗
bool ONGEO_Polynomial_TightenRootRangeBySmale(double ti1, double ti2, double *coef, int num, double &to1, double &to2) {
	/// 参考文献
	/// Yinyu Ye,
	/// Combining Binary Search and Newton's Method to Compute Real Roots for a Class of Real Functions.
	/// Journal Of Complexity 10, pp 271-280 (1994)

	// 最大次数を求める。
	int d = num;
	for (int i = 0; i < num; ++i){
		if (std::abs(coef[i]) < ON_ZERO_TOLERANCE) --d;
		else break;
	}
	double alpha = static_cast<double>(d - 1) * 0.5;
	double beta = 
		(ONGEO_Polynomial_Evaluate(ti1, coef, num) < ONGEO_Polynomial_Evaluate(ti2, coef, num)) ?
		1.0 / (1.0 - 1.0 / (8.0 * alpha)) : (1.0 + 1.0 / (8.0 * alpha);
	std::vector<double> b;
	double beta_pw = 1, th = ti2 / ti1;
	for(;;){
		if (beta_pw
		b.push_back(beta_pw);
	}
	


}
#endif