// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#include "opennurbs.h"

#ifdef _MSC_VER
#ifdef ONGEO_EXPORTS
#define ONGEO_CLASS __declspec(dllexport)
#define ONGEO_DECL __declspec(dllexport)
#else
#define ONGEO_CLASS
#define ONGEO_DECL
#endif
#else
#define ONGEO_CLASS
#define ONGEO_DECL
#endif
/// �N�G���_�ɍł��߂��L��Bezier�Ȑ���̓_�����߂�(BBClipping��������)�B
/// @param [in] bcs �x�W�G�Ȑ��̔z��̐擪�v�f�̃|�C���^
/// @param [in] num_bcs �x�W�G�Ȑ��̐�
/// @param [in] tolerance �����g�������X
/// @param [in] pt_query �N�G���_
/// @param [out] bc_nearest �ł��߂��_�������Ă���x�W�G�Ȑ�
/// @param [out] t �ł��߂��_�������Ă���x�W�G�Ȑ���ł̃p�����[�^�l(�e�x�W�G�Ȑ��̃p�����[�^�͈̔͂�0 - 1)
/// @param [out] pt_nearest �ł��߂��_�̎O�������W
/// @return �N�G���_�ƍŋߓ_�Ƃ̋���
ONGEO_DECL double ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod(const ON_BezierCurve *bcs, int num_bcs, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest);

/// �w�肵���_�ɍł��߂��L��Bezier�Ȑ���̓_�����߂�(���x���x�����ߔ񐄏�)�B
/// @param [in] bcs �x�W�G�Ȑ��̔z��̐擪�v�f�̃|�C���^
/// @param [in] num_bcs �x�W�G�Ȑ��̐�
/// @param [in] tolerance �����g�������X
/// @param [in] pt_query �N�G���_
/// @param [out] bc_nearest �ł��߂��_�������Ă���x�W�G�Ȑ�
/// @param [out] t �ł��߂��_�������Ă���x�W�G�Ȑ���ł̃p�����[�^�l(�e�x�W�G�Ȑ��̃p�����[�^�͈̔͂�0 - 1)
/// @param [out] pt_nearest �ł��߂��_�̎O�������W
/// @return �N�G���_�ƍŋߓ_�Ƃ̋���
ONGEO_DECL double ONGEO_NearestPointBezierCurve_BBClipping(const ON_BezierCurve *bc_begin, const ON_BezierCurve *bc_end, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest);

/// ������(ON_3dRay)�ƗL��Bezier�ȖʂƂ̌�_�����߂�B
/// @param [in] ray ������
/// @param [in] bez �x�W�G�Ȗ�
/// @param [out] tuvints ��_�̔��������̃p�����[�^�A�Ȗʑ���u,v�p�����[�^��3�̒l�ō\�������x�N�g���̔z��
/// @param [out] ptsrfs ��_�̋Ȗʑ��̎O�������W�̔z��
/// @param [out] ptlins ��_�̔��������̎O�������W�̔z��
/// @param [in] tolerance �����_����g�������X
/// @return 0:�����A����ȊO:���s
ONGEO_DECL int ONGEO_IntersectRayBezier_QuasiInterpolating(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance);

ONGEO_DECL void ONGEO_CalculateMinMaxWeight(const ON_BezierSurface &src, double &wmin, double &wmax);

/// �x�W�G�Ȗʂ��ޑ�܂���BoundingSphere�𐶐�����B
ONGEO_DECL int ONGEO_CalculateRoughBoundingSphere(const ON_BezierSurface &src, ON_3dPoint &center, double &radius);

/// �J���_�m�̌�����3x3�s��̌ŗL�l�����߂�B
/// @param [in] A 3x3�̍s��
/// @param [out] l 3�̌ŗL�l (l[0] : 1�ڂ̉��̎����� l[1] : 1�ڂ̉��̋������A l[2] : 2�ڂ̉��̎������A ...)
/// @return true : �v�Z����, false : ���s
ONGEO_DECL bool ONGEO_EigenValue3_Cardano(double *A[], double l[6]);

/// 3x3�s��ƌŗL�l����A�ŗL�l�ɑΉ������ŗL�x�N�g�������߂�B�ŗL�l�����f���̏ꍇ�͔�Ή�
/// @param [in] A 3x3�̍s��
/// @param [in] l �ŗL�l (����)
/// @param [out] �ŗL�x�N�g�� (x, y, z�����ꂩ�̐�����1�ƂȂ��Ă���B�K�v�ɉ����ĒP�ʉ����邱��)
/// @return true : �v�Z����, false : ���s
ONGEO_DECL bool ONGEO_EigenVector3(double *A[], double l, double v[3]);

// 3�����x�N�g���z�񂩂�3x3�̕��U�����U�s��A�܂��͑��֌W���s������߂�B
// @param vecs [in] �x�N�g���z��
// @param num_vecs [in] �x�N�g���̐�
// @param calc_corr_coef_matrix [in] false:���U�����U�����߂� true:���֌W���s������߂�B
// @param A [out] 3x3�̍s�� �̈�͌Ăяo�����Ŋm�ۂ��邱�ƁB
// @param meanvec [out] ��NULL�̂Ƃ��A���σx�N�g����Ԃ��B
ONGEO_DECL void ONGEO_Create_Covariance_Matrix3x3(double *vecs, int num_vecs, bool calc_corr_coef_matrix, double *A[3], double *meanvec = 0);

// 3x3�̕��U�����U�s�񂩂瑊�֌W���s������߂�B
// @param A [in] 3x3�̕��U�����U�s��
// @param B [out] 3x3�̑��֌W���s��
ONGEO_DECL void ONGEO_Covariance2CorrCoef_Matrix3x3(double *A[3], double *B[3]);

/// u�A�܂���v�ŕ������J��Ԃ��A���ꂼ��̕������ꂽ�Ȗʃp�b�`���m�[�h�ɂ���Sphere2����
/// �e�m�[�h�ɂ͒��S�_�Ɣ��a�A�����������i�[�����B
struct ONGEO_CLASS ONGEO_SphereTree{
	struct Node{
		ON_3dPoint center; /// ���̒��_
		double radius2;    /// ���a��2��
		int bez_index;     /// ���Ԗڂ̃x�W�G�Ȗʂ� �����Ώۂ�Nurbs�ȖʌQ�̏ꍇ��-2
		int direction;     /// -3 : �G���[�m�[�h�A-2 : Nurbs�ȖʌQ�̕���, -1 : �����Ȃ�, 0 : split at u, 1 : split at v
		int childLeft, childRight; // nodes�̃C���f�b�N�X�B-1�̏ꍇ�͎q����
	};
	ON_SimpleArray<int> nurbs_index_to_bez_first_index;
	int num_root_bezsurfs;
	ON_SimpleArray<ON_BezierSurface> bezs;
	ON_SimpleArray<Node> nodes;
	ONGEO_SphereTree(int num_surfs, const ON_NurbsSurface *nbsurf);
	~ONGEO_SphereTree();

	int GetRootNodeIndex() const;

	/// �؂𐶐�����B �����I���������ݒ肳��Ă��Ȃ��ꍇ�͖؂𐶐����Ȃ��B
	/// @param radius [in] ���a��radius�ȉ��ɂȂ�܂ŕ������J��Ԃ��B0�̏ꍇ�͔��a�𕪊��I�������Ƃ��Ȃ��B
	/// @param weight_rate [in]  �d�݂̍ő�l�ƍŏ��l�̔�(�ő�l/�ŏ��l)�����̒l�ȉ��ɂȂ�܂ŕ������J��Ԃ��B0�̏ꍇ�͕����I�������Ƃ��Ȃ��B
	/// @param level [in] �ő啪�����x�� -1�̏ꍇ�͕����I�������Ƃ��Ȃ��B
	/// @return 0:�����A����ȊO:���s
	int CreateTree(double radius2, double weight_rate, int level);

	/// �x�W�G�Ȗʂ̃C���f�b�N�X����A(�R���X�g���N�^�ŗ^����)���ƂȂ�Nurbs�Ȗʂ̃C���f�b�N�X���擾����B 
	/// @param bez_index [in] �x�W�G�Ȗʂ̃C���f�b�N�X
	/// @return Nurbs�Ȗʂ̃C���f�b�N�X
	int GetNurbsIndexFromBezIndex(int bez_index) const;

	/// �x�W�G�Ȗʂ̃C���f�b�N�X����A����Nurbs�Ȗʂ��番�����ꂽ�ʂ̂����A�擪�̃x�W�G�Ȗʂ̃C���f�b�N�X���擾����B 
	/// @param bez_index [in] �x�W�G�Ȗʂ̃C���f�b�N�X
	/// @return ����Nurbs�Ȗʂ��番�����ꂽ�ʂ̒��Ő擪�̃x�W�G�Ȗʂ̃C���f�b�N�X
	int GetFirstBezIndexFromBezIndex(int bez_index) const;

	/// Nurbs�ȖʌQ�ƃx�W�G�Ȗʂ̃C���f�b�N�X����A���̃x�W�G�ȖʂɑΉ�����Nurbs�Ȗʏ�̃p�����[�^�͈͂����߂�B 
	/// @param nbsurf [in] Nurbs�Ȗʂ̔z��̃|�C���^�B�R���X�g���N�^�ŗ^�������̂Ɠ���ł��邱��
	/// @param bez_index [in] �x�W�G�Ȗʂ̃C���f�b�N�X
	/// @param range [out] �p�����[�^�͈� [0] : u�����A [1] : v����
	/// @return bez_index��������Nurbs�Ȗʂ̃C���f�b�N�X
	int GetNurbsIntervalFromBezIndex(const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]) const;

	struct Result{
		int bez_index;
		int node_index;
		ON_Interval uint, vint;
	};
	/// ray�Ɩ؂̌����e�X�g�����{����B�t�m�[�h��H�������ʂ�results�ɏ����o���B
	/// ray��m_V�͒P�ʉ�����Ă��邱��
	/// �Ώۂ�bezs���o�ꎟ�Ȗʂ̂Ƃ��͖�������results�ɏ����o���B
	void RayIntersectTest(const ON_3dRay &ray, ON_SimpleArray<Result> &results) const;
};

// Face������Loop�Q��UV�J�[�u���x�W�G�Ȑ��Q�ɕ�������B
// @param [in] face �Ώۂ�Face�v�f
// @param [out] loop_crvs Loop��UV�J�[�u�Q
// @param [out] num_crvs_in_a_loop �eLoop������BezierUV�J�[�u�������Ă��邩�������B num_crvs_in_a_loop[0]���O�����[�v���\������Ȑ��̐��A[1]�ȍ~���������[�v
//   �Ⴆ�΁Aloop_crvs�� num_crvs_in_a_loop[0]�Ԗځ`num_crvs_in_a_loop[0]+num_crvs_in_a_loop[1]-1�Ԗڂ܂ł���ڂ̓����̃��[�v���\������Bezier�Ȑ��������B
ONGEO_DECL void ONGEO_GetBezierLoops(const ON_BrepFace &face, ON_SimpleArray<ON_BezierCurve> &loop_crvs, ON_SimpleArray<int> &num_crvs_in_a_loop);

// UV�J�[�u�������x�W�G�Ȑ��Q����Auv�_�̓��O��������{����B
// @param [in] loop_crvs Loop��UV�J�[�u�Q
// @param [in] num_crvs_in_a_loop �eLoop������BezierUV�J�[�u�������Ă��邩
// @param [in] uv uv�_
// @param [in] tolerance uv��ԓ��ł̈ʒu�g�������X
// @return true:uv�_�̓��[�v�̒��Afalse:uv�_�̓��[�v�̊O
ONGEO_DECL bool ONGEO_UVPointIsInside(const ON_SimpleArray<ON_BezierCurve> &loop_crvs, const ON_SimpleArray<int> &num_crvs_in_a_loop, const ON_2dPoint &uv, double tolerance);

ONGEO_DECL ONGEO_SphereTree *ONGEO_NewSphereTree(int num, const ON_NurbsSurface *nbsurfs);
ONGEO_DECL void ONGEO_DeleteSphereTree(ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_CreateTree(ONGEO_SphereTree *, double radius2, double weight_rate, int level);
ONGEO_DECL void ONGEO_SphereTree_RayIntersectTest(const ONGEO_SphereTree *, const ON_3dRay &ray, ON_SimpleArray<ONGEO_SphereTree::Result> &results);
ONGEO_DECL int ONGEO_SphereTree_GetRootNodeIndex(const ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetFirstBezIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIntervalFromBezIndex(const ONGEO_SphereTree *, const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]);

// �������֐��ɒl�������Čv�Z����B
// @param [in] t �������l
// @param [in] coef �������̌W����
// @param [in] num �W���̌�(=����+1)
// @return �v�Z����
ONGEO_DECL double ONGEO_Polynomial_Evaluate(double t, const double *coef, int num);

// 2�̑������̘a�����߂�B
// @param [in] coef1 1�߂̑������̌W���� num1 == 0 �̂Ƃ��Ɍ���Acoef1�̓k���|�C���^�ł���
// @param [in] num1 1�߂̑������̌W���̌�(=����+1)
// @param [in] coef2 2�߂̑������̌W���� num2 == 0 �̂Ƃ��Ɍ���Acoef2�̓k���|�C���^�ł���
// @param [in] num2 2�߂̑������̌W���̌�(=����+1)
// @param [out] coef_add ������1 + ������2�������������̌W����
// @param [in] num �o�͑������̌W���̌� (num1 <= num ���� num2 <= num�ł��邱��  coef1 �܂��� coef2�Ɠ����ꏊ���w���Ă��Ă��K�؂Ɍv�Z�\)
// @return true:�����A false:���s
ONGEO_DECL bool ONGEO_Polynomial_Add(const double *coef1, int num1, const double *coef2, int num2, double *coef_add, int num_add);

// 2�̑������̍������߂�B
// @param [in] coef1 1�߂̑������̌W���� num1 == 0 �̂Ƃ��Ɍ���Acoef1�̓k���|�C���^�ł���
// @param [in] num1 1�߂̑������̌W���̌�(=����+1)
// @param [in] coef2 2�߂̑������̌W���� num2 == 0 �̂Ƃ��Ɍ���Acoef2�̓k���|�C���^�ł���
// @param [in] num2 2�߂̑������̌W���̌�(=����+1)
// @param [out] coef_sub ������1 - ������2�������������̌W����
// @param [in] num �o�͑������̌W���̌� (num1 <= num ���� num2 <= num�ł��邱��  coef1 �܂��� coef2�Ɠ����ꏊ���w���Ă��Ă��K�؂Ɍv�Z�\)
// @return true:�����A false:���s
ONGEO_DECL bool ONGEO_Polynomial_Subtract(const double *coef1, int num1, const double *coef2, int num2, double *coef_sub, int num_sub);

	// 2�̑������̐ς𑽍����Ƃ��ēW�J����B
// @param [in] coef1 1�߂̑������̌W����
// @param [in] num1 1�߂̑������̌W���̌�(=����+1)
// @param [in] coef2 2�߂̑������̌W����
// @param [in] num2 2�߂̑������̌W���̌�(=����+1)
// @param [out] coef_mul 2�̑������̐ς�W�J�����������̌W����
//	([num1+num2-2]���̑������ƂȂ邽�߁A���͑���num1+num2-1���̌W���̗̈���m�ۂ���K�v������B)
ONGEO_DECL void ONGEO_Polynomial_Multiply(const double *coef1, int num1, const double *coef2, int num2, double *coef_mul);

// ����������K��������B
// @param [in] coef �������̌W����
// @param [in] num �������̌W���̌�(=����+1)
// @param [out] coef_dif ���������������̌W����
//    (�Ăяo������num-1���̗̈���m�ۂ��邱��
//     coef�Ɠ����ꏊ�Ɍ��ʂ��㏑���������ꍇ�́A���̈����Ƃ���&coef[1]��n���A�{�֐��ďo�����coef[0] = 0.0�Ƃ��邱��)
ONGEO_DECL void ONGEO_Polynomial_Differential(const double *coef, int num, double *coef_dif);

// �������̌W���񂩂�Sturm��𐶐�����B
// @param coef [in] �������̌W���� (�����̑傫����)
// @param num [in] �������̌W���̐�
// @param s [out] Sturm�� (�Ăяo������nf*nf���̗̈���m�ۂ��邱��)
ONGEO_DECL void ONGEO_Polynomial_CreateSturmSequence(const double *coef, int num, double *strum);

// Sturm���p���āA�w�肵���l�ł̕������]�񐔂��v�Z����B
// @param t [in] �l
// @param sturm [in] Sturm��
// @param num [in] �������̌W���̐� (Sturm��̔z�񒷂�num*num)
// @return �������]��
ONGEO_DECL int ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(double t, const double *strum, int num);

// Sturm���p���āA�w�肵����ԓ��̑������̍��̐����v�Z����B
// @param t1 [in] ��Ԃ̉���
// @param t2 [in] ��Ԃ̏��(t1 < t2 �ł��邱��)
// @param sturm [in] Sturm��
// @param num [in] �������̌W���̐� (Sturm��̔z�񒷂�num*num)
// @return ���̌�(�G���[�̏ꍇ�A����)
ONGEO_DECL int ONGEO_Polynomial_CalculateNumRoot(double t1, double t2, const double *strum, int num);

// �p�X�J���̎O�p�`��dim+1�s�ڂ��v�Z���A�֐����ŗp�ӂ����z���Ԃ��B
// @param dim [in] �v�Z�Ώۂ̎���
// @return �p�X�J���̎O�p�`��dim+1�s�� ���� dim+1�̔z���Ԃ��B�z��͊֐����Ǘ����邽�߁Adelete[]���ŉ�����Ă͂Ȃ�Ȃ��B
ONGEO_DECL int *ONGEO_Polynomial_CalcPascalTriangle(int dim);

// �p�X�J���̎O�p�`���g�p���ē񍀌W�� mCn ���v�Z����B
// @param m, n [in] �W��
// @return �񍀌W��
ONGEO_DECL int ONGEO_Polynomial_CalcBinomialCoef(int m, int n);

// Brent�@��p���č������߂�B��� ti1 - ti2 �̒��ŒP�������A�܂��͒P�������ł��邱�ƁB
// @param ti1 [in] ���߂�����Ԃ̉���
// @param ti2 [in] ���߂�����Ԃ̏�� ti1 < ti2�ł��邱�ƁB
// @param coef [in] �������̌W��
// @param num [in] �������̌W���̐�
// @return ���ƂȂ�l
ONGEO_DECL double ONGEO_Polynomial_FindRootByBrentMethod(double ti1, double ti2, const double *coef, int num);
