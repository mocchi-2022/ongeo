/*
 * eigen3_cardano
 * Copylight (C) 2012 mocchi
 * mocchi_2003@yahoo.co.jp
 * License: Boost ver.1
 */
// �J���_�m�̌����ŌŗL�l�����߁A�A��1���������𒼐ډ����ČŗL�x�N�g�������߂�v���O����

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <complex>

#include "ONGEO.h"
#include "Profile.h"

namespace {
typedef std::complex<double> dcomp_t;

struct polar_rep{
	double abs;
	double arg;
	polar_rep(double abs_, double arg_) : abs(abs_), arg(arg_){
	}
	polar_rep conj(){
		return polar_rep(abs, -arg);
	}
};

// http://people.freebsd.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf
double cubert_real(double a){
	// a ���A b * 2^n �ɕ��� (0.125 <= b < 1.0�A n��3�̔{��)
	int n;
	double b = std::frexp(a, &n);
	int na = (n < 0) ? ((-n) % 3) : ((3 - (n % 3)) % 3);
	n += na;
	b /= static_cast<double>(1 << na);

	// ���̎��_��n�͕K��3�̔{��
	// �w������n/3�ŗ���������������ƂɂȂ�B
	n /= 3;

	// �������̗�����
	b = -0.46946116*b*b + 1.072302*b + 0.3812513; // ���x 6bit

	double r = std::ldexp(b, n);

	r = 2.0/3.0 * r + 1.0/3.0 * a / (r * r); // 12bit
	r = 2.0/3.0 * r + 1.0/3.0 * a / (r * r); // 24bit
	r = 2.0/3.0 * r + 1.0/3.0 * a / (r * r); // 48bit
	r = 2.0/3.0 * r + 1.0/3.0 * a / (r * r); // 96bit
	return r;
}

polar_rep cubert_polar_rep(dcomp_t &x){
	double xabs = std::sqrt(x.real()*x.real()+x.imag()*x.imag());
	return polar_rep(cubert_real(xabs), std::atan2(x.imag(), x.real()) / 3.0);
}

const static dcomp_t w1(-0.5,  0.86602540378443864676372317075294);
const static dcomp_t w2(-0.5, -0.86602540378443864676372317075294);
}

bool ONGEO_EigenValue3_Cardano(double *A[], double l[6]){
	PROF("EigenValue3_Cardano");
	const double a11 = A[0][0], a12 = A[0][1], a13 = A[0][2];
	const double a21 = A[1][0], a22 = A[1][1], a23 = A[1][2];
	const double a31 = A[2][0], a32 = A[2][1], a33 = A[2][2];
	// �ŗL������ |A - lI| = 0����l�̕����������B
	const double a = -1.0;
	const double b = a11 + a22 + a33;
	const double c = -a11*a22 - a11*a33 + a12*a21 + a13*a31 - a22*a33 + a23*a32;
	const double d =
		  a11*a22*a33 - a11*a23*a32 - a12*a21*a33
		+ a12*a23*a31 + a13*a21*a32 - a13*a22*a31;

	// a*l^3 + b*l^2 + c*l + d = 0 ���O���������̉��̌����ŉ���
	//  R1 : �������̒��̕������̒��g
	dcomp_t R1sq;
	double R1 = 3*(27*a*a*d*d - 18*a*b*c*d + 4*a*c*c*c + 4*b*b*b*d - b*b*c*c);
	if (R1 >= 0) R1sq = std::sqrt(R1)*12.0*a;
	else R1sq = dcomp_t(0, std::sqrt(-R1)*12.0*a);
	//  R2 : �������̒��g���畽�����̍�������������
	dcomp_t R2(4*(-27*a*a*d + 9*a*b*c - 2*b*b*b), 0);

	dcomp_t R2_plus_R1sq = R2 + R1sq;
	dcomp_t R2_minus_R1sq = R2 - R1sq;
	polar_rep R2_1cb = cubert_polar_rep(R2_plus_R1sq);
	// R2_plus_R1sq��R2_minus_R1sq�́A�Ƃ��Ɏ���(R1�����̂Ƃ�)�A�܂��͕��f�����̊֌W(R1�����̂Ƃ�)
	// ���f����(R2_minus_R1sq = std::conj(R2_plus_R1sq))�̂Ƃ��A R2_2cb = std::conj(R2_1cb)
	polar_rep R2_2cb = (R2_1cb.arg != 0) ? R2_1cb.conj() : polar_rep(cubert_real(R2_minus_R1sq.real()), 0);

	double a3i = 1.0 / (3.0*a);
	if (R1 < 0){
		double cos_R2_1cb = std::cos(R2_1cb.arg);
		double sin_R2_1cb = std::sqrt(1 - cos_R2_1cb * cos_R2_1cb);
		l[0] = (R2_1cb.abs*cos_R2_1cb - b)*a3i;
		l[2] = (R2_1cb.abs*(cos_R2_1cb*w1.real() - sin_R2_1cb*w1.imag()) - b)*a3i;
		l[4] = (R2_1cb.abs*(cos_R2_1cb*w1.real() + sin_R2_1cb*w1.imag()) - b)*a3i;
		l[1] = l[3] = l[5] = 0;
	}else{
		double a6i = a3i * 0.5;
		l[0] = (-2*b+R2_1cb.abs+R2_2cb.abs)*a6i, l[1] = 0;
		dcomp_t s2 = (-2*b+w1*R2_1cb.abs+w2*R2_2cb.abs)*a6i;
		l[2] = s2.real(), l[3] = s2.imag();
		dcomp_t s3 = (-2*b+w2*R2_1cb.abs+w1*R2_2cb.abs)*a6i;
		l[4] = s3.real(), l[5] = s3.imag();
	}
	// �e���̐��x���j���[�g���@�ŏグ��B
	for (int i = 0; i < 6; ++i){
		if (l[i] == 0) continue;
		double ll = l[i];
		double l2 = ll * ll, l3 = ll * l2;
		l[i] = ll - (a*l3 + b*l2 + c*ll + d) / (3*a*l2 + 2*b*ll + c);
	}
	return true;
}

bool ONGEO_EigenVector3(double *A[], double l, double v[3]){
	PROF("EigenVector3");
	const double a11 = A[0][0]-l, a12 = A[0][1], a13 = A[0][2];
	const double a21 = A[1][0], a22 = A[1][1]-l, a23 = A[1][2];
	const double a31 = A[2][0], a32 = A[2][1], a33 = A[2][2]-l;
	double denom_i[] = {
		a11*a22-a12*a21, a11*a23-a13*a21, a12*a23-a13*a22,
		a11*a32-a12*a31, a11*a33-a13*a31, a12*a33-a13*a32,
		a21*a32-a22*a31, a21*a33-a23*a31, a22*a33-a23*a32
	};
	int index_max = 0;
	double amax = std::abs(denom_i[0]);
	for (int i = 1; i < 9; ++i){
		double a = std::abs(denom_i[i]);
		if (amax < a) amax = a, index_max = i;
	}
	if (amax < DBL_EPSILON * 7) return false;
	double denom = denom_i[index_max];
	switch(index_max){
		case 0:
			v[0] = (a12*a23-a13*a22) / denom, v[1] = (a13*a21-a11*a23) / denom, v[2] = 1.0;
			return true;
		case 1:
			v[0] = (a13*a22-a12*a23) / denom, v[1] = 1.0, v[2] = (a12*a21-a11*a22) / denom;
			return true;
		case 2:
			v[0] = 1.0, v[1] = (a13*a21-a11*a23) / denom, v[2] = (a11*a22-a12*a21) / denom;
			return true;
		case 3:
			v[0] = (a12*a33-a13*a32) / denom, v[1] = (a13*a31-a11*a33) / denom, v[2] = 1.0;
			return true;
		case 4:
			v[0] = (a13*a32-a12*a33) / denom, v[1] = 1.0, v[2] = (a12*a31-a11*a32) / denom;
			return true;
		case 5:
			v[0] = 1.0, v[1] = (a13*a31-a11*a33) / denom, v[2] = (a11*a32-a12*a31) / denom;
			return true;
		case 6:
			v[0] = (a22*a33-a23*a32) / denom, v[1] = (a23*a31-a21*a33) / denom, v[2] = 1.0;
			return true;
		case 7:
			v[0] = (a23*a32-a22*a33) / denom, v[1] = 1.0, v[2] = (a22*a31-a21*a32) / denom;
			return true;
		case 8:
			v[0] = 1.0, v[1] = (a23*a31-a21*a33) / denom, v[2] = (a21*a32-a22*a31) / denom;
			return true;
	}
	return false;
}

void ONGEO_Create_Covariance_Matrix3x3(double *vecs, int num_vecs, bool calc_corr_coef_matrix, double *A[3], double *meanvec){
	// ���U�����U�s������߂�B
	ON_3dVector mean = ON_3dVector(0,0,0);
	// �܂��͊e�����̕���
	for (int i = 0, i3 = 0; i < num_vecs; ++i, i3 += 3){
		mean += ON_3dVector(vecs + i3);
	}
	double denom = 1.0 / static_cast<double>(num_vecs);
	mean *= denom;
	A[0][0] = A[0][1] = A[0][2] = 0;
	A[1][0] = A[1][1] = A[1][2] = 0;
	A[2][0] = A[2][1] = A[2][2] = 0;

	// ���U�E�����U�s��쐬
	for (int i = 0, i3 = 0; i < num_vecs; ++i, i3 += 3){
		ON_3dVector dif = ON_3dVector(vecs + i3) - mean;
		A[0][0] += dif.x * dif.x;
		A[0][1] += dif.x * dif.y;
		A[0][2] += dif.x * dif.z;
		A[1][1] += dif.y * dif.y;
		A[1][2] += dif.y * dif.z;
		A[2][2] += dif.z * dif.z;
	}
	A[0][0] *= denom, A[0][1] *= denom, A[0][2] *= denom;
	A[1][1] *= denom, A[1][2] *= denom, A[2][2] *= denom;

	if (calc_corr_coef_matrix){
		// ���֍s��ɕϊ�
		double s[3] = {std::sqrt(A[0][0]), std::sqrt(A[1][1]), std::sqrt(A[2][2])};
		A[0][0] = (s[0] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
		A[1][1] = (s[1] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
		A[2][2] = (s[2] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
		double ss;
		A[0][1] = ((ss = s[0]*s[1]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][1] /= ss;
		A[0][2] = ((ss = s[0]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][2] /= ss;
		A[1][2] = ((ss = s[1]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[1][2] /= ss;
	}

	A[1][0] = A[0][1], A[2][0] = A[0][2], A[2][1] = A[1][2];
	if (meanvec){
		for (int h = 0; h < 3; ++h) meanvec[h] = mean[h];
	}
}

// ���U�E�����U�s�񂩂瑊�֌W���s��𐶐�
void ONGEO_Covariance2CorrCoef_Matrix3x3(double *A[3], double *B[3]){
	double s[3] = {std::sqrt(A[0][0]), std::sqrt(A[1][1]), std::sqrt(A[2][2])};
	B[0][0] = (s[0] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
	B[1][1] = (s[1] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
	B[2][2] = (s[2] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
	double ss;
	B[0][1] = ((ss = s[0]*s[1]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][1] /= ss;
	B[0][2] = ((ss = s[0]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][2] /= ss;
	B[1][2] = ((ss = s[1]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[1][2] /= ss;
	B[1][0] = B[0][1], B[2][0] = B[0][2], B[2][1] = B[1][2];
}
