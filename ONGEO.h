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
/// クエリ点に最も近い有理Bezier曲線上の点を求める(BBClippingよりも高速)。
/// @param [in] bcs ベジエ曲線の配列の先頭要素のポインタ
/// @param [in] num_bcs ベジエ曲線の数
/// @param [in] tolerance 距離トレランス
/// @param [in] pt_query クエリ点
/// @param [out] bc_nearest 最も近い点が属しているベジエ曲線
/// @param [out] t 最も近い点が属しているベジエ曲線上でのパラメータ値(各ベジエ曲線のパラメータの範囲は0 - 1)
/// @param [out] pt_nearest 最も近い点の三次元座標
/// @return クエリ点と最近点との距離
ONGEO_DECL double ONGEO_NearestPointBezierCurve_ImprovedAlgebraicMethod(const ON_BezierCurve *bcs, int num_bcs, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest);

/// 指定した点に最も近い有理Bezier曲線上の点を求める(速度が遅いため非推奨)。
/// @param [in] bcs ベジエ曲線の配列の先頭要素のポインタ
/// @param [in] num_bcs ベジエ曲線の数
/// @param [in] tolerance 距離トレランス
/// @param [in] pt_query クエリ点
/// @param [out] bc_nearest 最も近い点が属しているベジエ曲線
/// @param [out] t 最も近い点が属しているベジエ曲線上でのパラメータ値(各ベジエ曲線のパラメータの範囲は0 - 1)
/// @param [out] pt_nearest 最も近い点の三次元座標
/// @return クエリ点と最近点との距離
ONGEO_DECL double ONGEO_NearestPointBezierCurve_BBClipping(const ON_BezierCurve *bc_begin, const ON_BezierCurve *bc_end, double tolerance, const ON_3dPoint &pt_query, const ON_BezierCurve *&bc_nearest, double &t, ON_3dPoint &pt_nearest);

/// 半直線(ON_3dRay)と有理Bezier曲面との交点を求める。
/// @param [in] ray 半直線
/// @param [in] bez ベジエ曲面
/// @param [out] tuvints 交点の半直線側のパラメータ、曲面側のu,vパラメータの3つの値で構成されるベクトルの配列
/// @param [out] ptsrfs 交点の曲面側の三次元座標の配列
/// @param [out] ptlins 交点の半直線側の三次元座標の配列
/// @param [in] tolerance 同一交点判定トレランス
/// @return 0:成功、それ以外:失敗
ONGEO_DECL int ONGEO_IntersectRayBezier_QuasiInterpolating(const ON_3dRay &ray, const ON_BezierSurface &bez, ON_3dPointArray &tuvints, ON_3dPointArray &ptsrfs, ON_3dPointArray &ptlins, double tolerance);

ONGEO_DECL void ONGEO_CalculateMinMaxWeight(const ON_BezierSurface &src, double &wmin, double &wmax);

/// ベジエ曲面を包む大まかなBoundingSphereを生成する。
ONGEO_DECL int ONGEO_CalculateRoughBoundingSphere(const ON_BezierSurface &src, ON_3dPoint &center, double &radius);

/// カルダノの公式で3x3行列の固有値を求める。
/// @param [in] A 3x3の行列
/// @param [out] l 3つの固有値 (l[0] : 1つ目の解の実数部 l[1] : 1つ目の解の虚数部、 l[2] : 2つ目の解の実数部、 ...)
/// @return true : 計算成功, false : 失敗
ONGEO_DECL bool ONGEO_EigenValue3_Cardano(double *A[], double l[6]);

/// 3x3行列と固有値から、固有値に対応した固有ベクトルを求める。固有値が複素数の場合は非対応
/// @param [in] A 3x3の行列
/// @param [in] l 固有値 (実数)
/// @param [out] 固有ベクトル (x, y, zいずれかの成分が1となっている。必要に応じて単位化すること)
/// @return true : 計算成功, false : 失敗
ONGEO_DECL bool ONGEO_EigenVector3(double *A[], double l, double v[3]);

// 3次元ベクトル配列から3x3の分散共分散行列、または相関係数行列を求める。
// @param vecs [in] ベクトル配列
// @param num_vecs [in] ベクトルの数
// @param calc_corr_coef_matrix [in] false:分散共分散を求める true:相関係数行列を求める。
// @param A [out] 3x3の行列 領域は呼び出し元で確保すること。
// @param meanvec [out] 非NULLのとき、平均ベクトルを返す。
ONGEO_DECL void ONGEO_Create_Covariance_Matrix3x3(double *vecs, int num_vecs, bool calc_corr_coef_matrix, double *A[3], double *meanvec = 0);

// 3x3の分散共分散行列から相関係数行列を求める。
// @param A [in] 3x3の分散共分散行列
// @param B [out] 3x3の相関係数行列
ONGEO_DECL void ONGEO_Covariance2CorrCoef_Matrix3x3(double *A[3], double *B[3]);

/// u、またはvで分割を繰り返し、それぞれの分割された曲面パッチをノードにしたSphere2分木
/// 各ノードには中心点と半径、分割方向が格納される。
struct ONGEO_CLASS ONGEO_SphereTree{
	struct Node{
		ON_3dPoint center; /// 球の中点
		double radius2;    /// 半径の2乗
		int bez_index;     /// 何番目のベジエ曲面か 分割対象がNurbs曲面群の場合は-2
		int direction;     /// -3 : エラーノード、-2 : Nurbs曲面群の分割, -1 : 分割なし, 0 : split at u, 1 : split at v
		int childLeft, childRight; // nodesのインデックス。-1の場合は子無し
	};
	ON_SimpleArray<int> nurbs_index_to_bez_first_index;
	int num_root_bezsurfs;
	ON_SimpleArray<ON_BezierSurface> bezs;
	ON_SimpleArray<Node> nodes;
	ONGEO_SphereTree(int num_surfs, const ON_NurbsSurface *nbsurf);
	~ONGEO_SphereTree();

	int GetRootNodeIndex() const;

	/// 木を生成する。 分割終了条件が設定されていない場合は木を生成しない。
	/// @param radius [in] 半径がradius以下になるまで分割を繰り返す。0の場合は半径を分割終了条件としない。
	/// @param weight_rate [in]  重みの最大値と最小値の比(最大値/最小値)がこの値以下になるまで分割を繰り返す。0の場合は分割終了条件としない。
	/// @param level [in] 最大分割レベル -1の場合は分割終了条件としない。
	/// @return 0:成功、それ以外:失敗
	int CreateTree(double radius2, double weight_rate, int level);

	/// ベジエ曲面のインデックスから、(コンストラクタで与えた)元となるNurbs曲面のインデックスを取得する。 
	/// @param bez_index [in] ベジエ曲面のインデックス
	/// @return Nurbs曲面のインデックス
	int GetNurbsIndexFromBezIndex(int bez_index) const;

	/// ベジエ曲面のインデックスから、同一Nurbs曲面から分割された面のうち、先頭のベジエ曲面のインデックスを取得する。 
	/// @param bez_index [in] ベジエ曲面のインデックス
	/// @return 同一Nurbs曲面から分割された面の中で先頭のベジエ曲面のインデックス
	int GetFirstBezIndexFromBezIndex(int bez_index) const;

	/// Nurbs曲面群とベジエ曲面のインデックスから、そのベジエ曲面に対応するNurbs曲面上のパラメータ範囲を求める。 
	/// @param nbsurf [in] Nurbs曲面の配列のポインタ。コンストラクタで与えたものと同一であること
	/// @param bez_index [in] ベジエ曲面のインデックス
	/// @param range [out] パラメータ範囲 [0] : u方向、 [1] : v方向
	/// @return bez_indexが属するNurbs曲面のインデックス
	int GetNurbsIntervalFromBezIndex(const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]) const;

	struct Result{
		int bez_index;
		int node_index;
		ON_Interval uint, vint;
	};
	/// rayと木の交差テストを実施する。葉ノードを辿った結果をresultsに書き出す。
	/// rayのm_Vは単位化されていること
	/// 対象のbezsが双一次曲面のときは無条件にresultsに書き出す。
	void RayIntersectTest(const ON_3dRay &ray, ON_SimpleArray<Result> &results) const;
};

// Faceが持つLoop群のUVカーブをベジエ曲線群に分解する。
// @param [in] face 対象のFace要素
// @param [out] loop_crvs LoopのUVカーブ群
// @param [out] num_crvs_in_a_loop 各Loopが何個のBezierUVカーブを持っているかを示す。 num_crvs_in_a_loop[0]が外側ループを構成する曲線の数、[1]以降が内側ループ
//   例えば、loop_crvsの num_crvs_in_a_loop[0]番目～num_crvs_in_a_loop[0]+num_crvs_in_a_loop[1]-1番目までが一つ目の内側のループを構成するBezier曲線を示す。
ONGEO_DECL void ONGEO_GetBezierLoops(const ON_BrepFace &face, ON_SimpleArray<ON_BezierCurve> &loop_crvs, ON_SimpleArray<int> &num_crvs_in_a_loop);

// UVカーブを示すベジエ曲線群から、uv点の内外判定を実施する。
// @param [in] loop_crvs LoopのUVカーブ群
// @param [in] num_crvs_in_a_loop 各Loopが何個のBezierUVカーブを持っているか
// @param [in] uv uv点
// @param [in] tolerance uv空間内での位置トレランス
// @return true:uv点はループの中、false:uv点はループの外
ONGEO_DECL bool ONGEO_UVPointIsInside(const ON_SimpleArray<ON_BezierCurve> &loop_crvs, const ON_SimpleArray<int> &num_crvs_in_a_loop, const ON_2dPoint &uv, double tolerance);

ONGEO_DECL ONGEO_SphereTree *ONGEO_NewSphereTree(int num, const ON_NurbsSurface *nbsurfs);
ONGEO_DECL void ONGEO_DeleteSphereTree(ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_CreateTree(ONGEO_SphereTree *, double radius2, double weight_rate, int level);
ONGEO_DECL void ONGEO_SphereTree_RayIntersectTest(const ONGEO_SphereTree *, const ON_3dRay &ray, ON_SimpleArray<ONGEO_SphereTree::Result> &results);
ONGEO_DECL int ONGEO_SphereTree_GetRootNodeIndex(const ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetFirstBezIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIntervalFromBezIndex(const ONGEO_SphereTree *, const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]);

// 多項式関数に値を代入して計算する。
// @param [in] t 代入する値
// @param [in] coef 多項式の係数列
// @param [in] num 係数の個数(=次数+1)
// @return 計算結果
ONGEO_DECL double ONGEO_Polynomial_Evaluate(double t, const double *coef, int num);

// 2つの多項式の和を求める。
// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_add 多項式1 + 多項式2を示す多項式の係数列
// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
// @return true:成功、 false:失敗
ONGEO_DECL bool ONGEO_Polynomial_Add(const double *coef1, int num1, const double *coef2, int num2, double *coef_add, int num_add);

// 2つの多項式の差を求める。
// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_sub 多項式1 - 多項式2を示す多項式の係数列
// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
// @return true:成功、 false:失敗
ONGEO_DECL bool ONGEO_Polynomial_Subtract(const double *coef1, int num1, const double *coef2, int num2, double *coef_sub, int num_sub);

	// 2つの多項式の積を多項式として展開する。
// @param [in] coef1 1つめの多項式の係数列
// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
// @param [in] coef2 2つめの多項式の係数列
// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
// @param [out] coef_mul 2つの多項式の積を展開した多項式の係数列
//	([num1+num2-2]次の多項式となるため、入力側でnum1+num2-1個分の係数の領域を確保する必要がある。)
ONGEO_DECL void ONGEO_Polynomial_Multiply(const double *coef1, int num1, const double *coef2, int num2, double *coef_mul);

// 多項式を一階微分する。
// @param [in] coef 多項式の係数列
// @param [in] num 多項式の係数の個数(=次数+1)
// @param [out] coef_dif 微分した多項式の係数列
//    (呼び出し元でnum-1個分の領域を確保すること
//     coefと同じ場所に結果を上書きしたい場合は、この引数として&coef[1]を渡し、本関数呼出し後にcoef[0] = 0.0とすること)
ONGEO_DECL void ONGEO_Polynomial_Differential(const double *coef, int num, double *coef_dif);

// 多項式の係数列からSturm列を生成する。
// @param coef [in] 多項式の係数列 (次数の大きい順)
// @param num [in] 多項式の係数の数
// @param s [out] Sturm列 (呼び出し元でnf*nf個分の領域を確保すること)
ONGEO_DECL void ONGEO_Polynomial_CreateSturmSequence(const double *coef, int num, double *strum);

// Sturm列を用いて、指定した値での符号反転回数を計算する。
// @param t [in] 値
// @param sturm [in] Sturm列
// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
// @return 符号反転回数
ONGEO_DECL int ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(double t, const double *strum, int num);

// Sturm列を用いて、指定した区間内の多項式の根の数を計算する。
// @param t1 [in] 区間の下限
// @param t2 [in] 区間の上限(t1 < t2 であること)
// @param sturm [in] Sturm列
// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
// @return 解の個数(エラーの場合、負数)
ONGEO_DECL int ONGEO_Polynomial_CalculateNumRoot(double t1, double t2, const double *strum, int num);

// パスカルの三角形のdim+1行目を計算し、関数側で用意した配列を返す。
// @param dim [in] 計算対象の次数
// @return パスカルの三角形のdim+1行目 長さ dim+1の配列を返す。配列は関数が管理するため、delete[]等で解放してはならない。
ONGEO_DECL int *ONGEO_Polynomial_CalcPascalTriangle(int dim);

// パスカルの三角形を使用して二項係数 mCn を計算する。
// @param m, n [in] 係数
// @return 二項係数
ONGEO_DECL int ONGEO_Polynomial_CalcBinomialCoef(int m, int n);

// Brent法を用いて根を求める。区間 ti1 - ti2 の中で単調増加、または単調減少であること。
// @param ti1 [in] 狭めたい区間の下限
// @param ti2 [in] 狭めたい区間の上限 ti1 < ti2であること。
// @param coef [in] 多項式の係数
// @param num [in] 多項式の係数の数
// @return 根となる値
ONGEO_DECL double ONGEO_Polynomial_FindRootByBrentMethod(double ti1, double ti2, const double *coef, int num);
