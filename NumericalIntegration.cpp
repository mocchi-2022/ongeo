#include "ONGEO.h"
#include "ONGEO_NumericalIntegration.hpp"

// 数値積分用の係数列

// Gauss-Kronrod 公式用の係数
// 参考
// http://www.ktech.biz/jp/archives/1009
// http://www.ktech.biz/jp/archives/1045

template<>
ONGEO_DECL void ONGEO_GaussKronrod_Coefs<3>(const double *&t_i_, const double *&wg_i_, const double *&wk_i_){
	static const double t_i[]  = { 0, 0.434243749346802, 0.774596669241483, 0.96049126870802 };
	static const double wg_i[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };
	static const double wk_i[] = { 0.450916538658474, 0.401397414775962, 0.268488089868333, 0.104656226026467 };
	t_i_ = t_i, wg_i_ = wg_i, wk_i_ = wk_i;
}

template<>
ONGEO_DECL void ONGEO_GaussKronrod_Coefs<5>(const double *&t_i_, const double *&wg_i_, const double *&wk_i_){
	static const double t_i[]  = { 0, 0.279630413161783, 0.538469310105683, 0.754166726570849, 0.906179845938664, 0.984085360094842 };
	static const double wg_i[] = { 0.236926885056189, 0.478628670499367, 0.568888888888889, 0.478628670499367, 0.236926885056189 };
	static const double wk_i[] = { 0.282987417857491, 0.272849801912559, 0.241040339228648, 0.186800796556493, 0.115233316622474, 0.0425820367510818 };
	t_i_ = t_i, wg_i_ = wg_i, wk_i_ = wk_i;
}

template<>
ONGEO_DECL void ONGEO_GaussKronrod_Coefs<7>(const double *&t_i_, const double *&wg_i_, const double *&wk_i_){
	static const double t_i[]  = { 0, 0.207784955007898, 0.405845151377397, 0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342758, 0.991455371120813 };
	static const double wg_i[] = { 0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870 };
	static const double wk_i[] = { 0.209482141084728, 0.204432940075299, 0.190350578064785, 0.169004726639268, 0.140653259715526, 0.104790010322250, 0.0630920926299790, 0.0229353220105292 };
	t_i_ = t_i, wg_i_ = wg_i, wk_i_ = wk_i;
}
