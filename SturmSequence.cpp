#include <cmath>

// �������֐� f(x) �ɂ����āAx=t�̂Ƃ���Sturm��𐶐�����B
// @param t [in] �������֐��ւ̈��� t
// @param g [in] �������֐�f(t)�̒l
// @param dg [in] �������֐�f'(t)�̒l
// @param num [in] ��̋��߂������� (3�ȏ�)
// @param s [out] Sturm�� (�Ăяo������num���̗̈���m�ۂ��邱��)
// @return Sturm��̊e��ɂ�����A�������]��
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
