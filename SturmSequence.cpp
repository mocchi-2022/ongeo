#include <cmath>

// 多項式関数 f(x) において、x=tのときのSturm列を生成する。
// @param t [in] 多項式関数への引数 t
// @param g [in] 多項式関数f(t)の値
// @param dg [in] 多項式関数f'(t)の値
// @param num [in] 列の求めたい長さ (3以上)
// @param s [out] Sturm列 (呼び出し元でnum個分の領域を確保すること)
// @return Sturm列の各列における、正負反転回数
int CreateSturmSequence(double t, double g, double dg, int num, double *s){
	int count = 0;
	s[0] = g;
	s[1] = dg;
	if (s[0] * s[1] < 0) ++count;
	for (int i = 2; i < num; ++i){
		s[i] = std::floor(s[i-2] / s[i-2]) * s[i-2] - s[i-3];
		if (s[i] * s[i-1] < 0) ++count;
	}
	return count;
}
