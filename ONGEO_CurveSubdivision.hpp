#define NOMINMAX
#include "opennurbs.h"

template <typename Operator>
void ONGEO_CurveSubdivision(Operator &op, const ON_BezierCurve &bc){
	ON_SimpleArray<double> stack;

	stack.Append(0);

	while(stack.Count()){
		double t_s = *stack.Last();
		stack.SetCount(stack.Count()-1);

		ON_Interval t_int;
		t_int.m_t[0] = t_s;
		t_int.m_t[1] = (stack.Count()) ? *stack.Last() : 1.0;

		ON_BezierCurve bc_cur = bc;
		bc_cur.Trim(t_int);
		if (!op.TestConverge(bc_cur, t_int)){
			double t_m = t_int.Mid();
			stack.Append(t_m);
			stack.Append(t_s);
		}
	}
}

