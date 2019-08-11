#define NOMINMAX
#include "opennurbs.h"

template <int N> void ONGEO_GaussKronrod_Coefs(const double *&t_i_, const double *&wg_i_, const double *&wk_i_);

template <int N> struct ONGEO_NumericalIntegration_GaussKronrod{
	double t[N*2+1], wg[N], wk[N*2+1];
	ONGEO_NumericalIntegration_GaussKronrod(double a, double b){
		SetRange(a, b);
	}
	void SetRange(double a, double b){
		if (a == b) return;
		const double *t_i, *wg_i, *wk_i;
		ONGEO_GaussKronrod_Coefs<N>(t_i, wg_i, wk_i);

		// ãÊä‘ïœçX
		double ba_diff = (b - a) * 0.5, ba_mean = (a + b) * 0.5;
		for (int i = 1; i <= N; ++i) t[N+i] = ba_diff * t_i[i] + ba_mean, t[N-i] = - ba_diff * t_i[i] + ba_mean;
		t[N] = ba_diff * t_i[0] + ba_mean;
		for (int i = 0; i < N; ++i) wg[i] = ba_diff * wg_i[i];
		for (int i = 1; i <= N; ++i) wk[N+i] = wk[N-i] = ba_diff * wk_i[i];
		wk[N] = ba_diff * wk_i[0];
	}

	template <typename Func>
	double Calc(const Func &f, double &err){
		double fv[N*2+1];
		for (int i = 0; i < N*2+1; ++i) fv[i] = f(t[i]);
		double g = 0;
		for (int i = 0; i < N; ++i){
			g += fv[i*2+1] * wg[i];
		}
		double k = 0;
		for (int i = 0; i < N*2+1; ++i) k += fv[i] * wk[i];
		err = std::pow(200.0 * std::abs(k - g), 1.5);
		return k;
	}
};

template <typename Integrator, typename Func> double ONGEO_AdaptiveQuadrature(const Func &f, double a, double b, double tolerance, double &err_result){
	double range = b - a;
	Integrator integ(0, 0);
	ON_SimpleArray<double> stack;

	double sum = 0;
	err_result = 0;

	stack.Append(b);
	stack.Append(a);

	while(stack.Count() > 1){
		double t_s = *stack.Last();
		stack.SetCount(stack.Count()-1);

		double s_a = t_s, s_b = *stack.Last();

		integ.SetRange(s_a, s_b);

		double sect_tolerance = tolerance * ((s_b - s_a) / range);

		double err;
		double r = integ.Calc(f, err);
		if (err < sect_tolerance){
			sum += r;
			err_result += err;
		}else{
			double t_m = (s_a + s_b) * 0.5;
			stack.Append(t_m);
			stack.Append(t_s);
		}
	}
	return sum;
};
