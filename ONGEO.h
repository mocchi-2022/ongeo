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

template <typename T> void ONGEO_Swap(ON_SimpleArray<T> &lhs, ON_SimpleArray<T> &rhs){
	int sz_l[] = {lhs.Count(), lhs.Capacity()};
	int sz_r[] = {rhs.Count(), rhs.Capacity()};
	T *ka_l = lhs.KeepArray();
	T *ka_r = rhs.KeepArray();
	lhs.SetArray(ka_r, sz_r[0], sz_r[1]);
	rhs.SetArray(ka_l, sz_l[0], sz_l[1]);
}

inline void ONGEO_GetCVHomogeneous(const ON_BezierSurface &bez, int i, int j, ON_3dPoint &pt){
	bez.GetCV(i, j, pt);
}

inline void ONGEO_GetCVHomogeneous(const ON_BezierSurface &bez, int i, int j, ON_4dPoint &pt){
	bez.GetCV(i, j, ON::homogeneous_rational, pt);
}

inline void ONGEO_GetCVHomogeneous(const ON_BezierCurve &bez, int i, ON_3dPoint &pt){
	bez.GetCV(i, pt);
}

inline void ONGEO_GetCVHomogeneous(const ON_BezierCurve &bez, int i, ON_4dPoint &pt){
	bez.GetCV(i, ON::homogeneous_rational, pt);
}


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

/// 有理ベジエ曲線の曲線長を再分割による収束計算で求める(制御点列長と始点-終点距離の差がトレランス以下になるまで)。
/// @param [in] bc ベジエ曲線
/// @param [in] tolerance トレランス
/// @return 曲線長
ONGEO_DECL double ONGEO_LengthBezierCurve_SimpleSubdivision(const ON_BezierCurve &bc, double tolerance);

/// 有理ベジエ曲線の曲線長を、引数で指定したパラメータで区間分割→数値積分で求める。
/// @param [in] bc ベジエ曲線
/// @param [in] 再分割パラメータ
/// @param [out] estimate_deviation 推定誤差
/// @return 曲線長
ONGEO_DECL double ONGEO_LengthBezierCurve_NumericalIntegration(const ON_BezierCurve &bc, const double *subdivided_prms, int cnt, double *estimated_deviation = 0);

/// 有理ベジエ曲線の曲線長を、引数で指定したパラメータで区間分割→適応求積法で求める。
/// @param [in] bc ベジエ曲線
/// @param [in] 再分割パラメータ
/// @param [in] トレランス
/// @param [out] estimate_deviation 推定誤差
/// @return 曲線長
ONGEO_DECL double ONGEO_LengthBezierCurve_AdaptiveQuadrate(const ON_BezierCurve &bc, const double *subdivided_prms, int cnt, double tolerance, double *estimated_deviation = 0);

/// Nurbs曲線の曲線長と曲線パラメータを相互変換する。
struct ONGEO_CLASS ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate{
	ON_NurbsCurve nc;

	ON_SimpleArray<double> nc_prms;
	ON_SimpleArray<double> nc_lengths;

	double tolerance;
	double estimated_deviation;

	ONGEO_LengthParamNurbsCurve_AdaptiveQuadrate(const ON_NurbsCurve &nc, double tolerance);
	double Length() const;
	double ParamToLength(double nc_prm) const;
	double LengthToParam(double nc_length) const;
};

/// 有理ベジエ曲線をテセレーションする(始点-終点を結ぶ線分と、各制御点との距離の最大値がトレランス以下になるまで)。
/// @param [in] bc ベジエ曲線
/// @param [in] tolerance トレランス
/// @param [out] prms テセレーション結果を示すパラメータ配列
ONGEO_DECL void ONGEO_TessellateBezierCurve_Simple(const ON_BezierCurve &bc, double tolerance, ON_SimpleArray<double> &prms);

/// 有理ベジエ曲線をテセレーションする(Quasi-Interpolation Control Net を生成したときの、曲線との距離の最大値がトレランス以下になるまで)。
/// @param [in] bc ベジエ曲線
/// @param [in] tolerance トレランス
/// @param [out] prms テセレーション結果を示すパラメータ配列
ONGEO_DECL void ONGEO_TessellateBezierCurve_QuasiInterpolating(const ON_BezierCurve &bc, double tolerance, ON_SimpleArray<double> &prms);

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

/// ベジエ曲線群を包む正確なBoundingBoxを生成する。
ONGEO_DECL int ONGEO_CalculateTightBoundingBox(const ON_BezierCurve *bcs, int num_bcs, double tolerance, ON_BoundingBox &bb);

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

/// 3次元ベクトル配列から3x3の分散共分散行列、または相関係数行列を求める。
/// @param vecs [in] ベクトル配列
/// @param num_vecs [in] ベクトルの数
/// @param calc_corr_coef_matrix [in] false:分散共分散を求める true:相関係数行列を求める。
/// @param A [out] 3x3の行列 領域は呼び出し元で確保すること。
/// @param meanvec [out] 非NULLのとき、平均ベクトルを返す。
ONGEO_DECL void ONGEO_Create_Covariance_Matrix3x3(double *vecs, int num_vecs, bool calc_corr_coef_matrix, double *A[3], double *meanvec = 0);

/// 3x3の分散共分散行列から相関係数行列を求める。
/// @param A [in] 3x3の分散共分散行列
/// @param B [out] 3x3の相関係数行列
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
	ON_ClassArray<ON_BezierSurface> bezs;
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

/// 複数のBrepに対して交点計算するためのオペレータ。
struct ONGEO_CLASS ONGEO_BrepsRayIntersect{
	struct LoopsInAFace{
		ON_ClassArray<ON_BezierCurve> loop_crvs;
		ON_SimpleArray<int> num_crvs_in_a_loop;
	};

	ON_ClassArray<ON_NurbsSurface> nbsurfs;
	ON_ClassArray<LoopsInAFace> loops_faces;
	ON_SimpleArray<int> nbs2brep; // nbsurfのインデックスからbrepのインデックスを求める配列
	ONGEO_SphereTree *st;

	ONGEO_BrepsRayIntersect(const ON_Brep **breps, int num_breps);
	ONGEO_BrepsRayIntersect(const ONGEO_BrepsRayIntersect &rhs);
	ONGEO_BrepsRayIntersect &operator =(const ONGEO_BrepsRayIntersect &rhs);
	~ONGEO_BrepsRayIntersect();

	struct Result{
		int nbs_index;
		ON_2dPoint uv;
		double t;
		ON_3dPoint codsrf, codlin;
	};

	enum TestStage{
		BeforeIntersectRayBezier, ///< ベジエ曲面との交点計算直前
		BeforeUVPointIsInside,    ///< UV空間における内外判定の直前。この直前のタイミングでresultsに交点情報が追記される。falseを返すとresultsに今回追記された交点情報が削除される。
		AfterUVPointIsInside      ///< UV空間における内外判定後。falseを返すとresultsに今回追記された交点情報が削除される。
	};

	/// 登録してあるBrepに対して交点計算を実施する。
	/// Test関数を引数として与えると、ベジエ曲面との交点計算、UV内外判定前、及び後の3ステージで関数が呼ばれる。
	/// Test関数がfalseを返した場合、交点計算対象から外れる。全交点のうち、特定のもののみが必要となるケース
	/// (例えば、原点に最も近い/遠い交点のみが必要といった場合)で、計算効率を向上できる。
	/// Test関数には、Run呼び出し時に引数として与えられたresult配列の参照が渡される。
	/// Test関数内でresultの内容を変えることができる。
	void Run(const ON_3dRay &ray, ON_SimpleArray<Result> &results, double tolerance = 1e-5, bool (*Test)(void *data, TestStage stage, int nbsurf_index, ON_SimpleArray<Result> &results) = 0, void *data = 0) const;
};

ONGEO_DECL ONGEO_BrepsRayIntersect *ONGEO_New_BrepsRayIntersect(const ON_Brep **breps, int num_breps);

/// 登録してあるBrepに対して交点計算を実施する。
/// Test関数を引数として与えると、ベジエ曲面との交点計算、UV内外判定前、及び後の3ステージで関数が呼ばれる。
/// Test関数がfalseを返した場合、交点計算対象から外れる。全交点のうち、特定のもののみが必要となるケース
/// (例えば、原点に最も近い/遠い交点のみが必要といった場合)で、計算効率を向上できる。
ONGEO_DECL void ONGEO_BrepsRayIntersect_Run(const ONGEO_BrepsRayIntersect *bri, const ON_3dRay &ray, double tolerance, ON_SimpleArray<ONGEO_BrepsRayIntersect::Result> &result, bool (*Test)(void *data, ONGEO_BrepsRayIntersect::TestStage stage, int nbsurf_index, ON_SimpleArray<ONGEO_BrepsRayIntersect::Result> &results) = 0, void *data = 0);

/// Faceが持つLoop群のUVカーブをベジエ曲線群に分解する。
/// @param [in] face 対象のFace要素
/// @param [out] loop_crvs LoopのUVカーブ群
/// @param [out] num_crvs_in_a_loop 各Loopが何個のBezierUVカーブを持っているかを示す。 num_crvs_in_a_loop[0]が外側ループを構成する曲線の数、[1]以降が内側ループ
///   例えば、loop_crvsの num_crvs_in_a_loop[0]番目～num_crvs_in_a_loop[0]+num_crvs_in_a_loop[1]-1番目までが一つ目の内側のループを構成するBezier曲線を示す。
ONGEO_DECL void ONGEO_GetBezierLoops(const ON_BrepFace &face, ON_ClassArray<ON_BezierCurve> &loop_crvs, ON_SimpleArray<int> &num_crvs_in_a_loop);

/// UVカーブを示すベジエ曲線群から、uv点の内外判定を実施する。
/// @param [in] loop_crvs LoopのUVカーブ群
/// @param [in] num_crvs_in_a_loop 各Loopが何個のBezierUVカーブを持っているか
/// @param [in] uv uv点
/// @param [in] tolerance uv空間内での位置トレランス
/// @return true:uv点はループの中、false:uv点はループの外
ONGEO_DECL bool ONGEO_UVPointIsInside(const ON_ClassArray<ON_BezierCurve> &loop_crvs, const ON_SimpleArray<int> &num_crvs_in_a_loop, const ON_2dPoint &uv, double tolerance);

/// UVカーブを示すベジエ曲線群から、uv点の内外判定を実施する。uv点がベジエ曲線上に乗っている場合は、指示したUVベクトルの向きで内外を判定する。
/// @param [in] loop_crvs LoopのUVカーブ群
/// @param [in] num_crvs_in_a_loop 各Loopが何個のBezierUVカーブを持っているか
/// @param [in] uv uv点
/// @param [in] uvdir uvベクトル
/// @param [in] tol_np 最短点計算トレランス
/// @param [in] tol_oncrv 曲線上判定トレランス
/// @return true:uv点はループの中、false:uv点はループの外
ONGEO_DECL bool ONGEO_UVPointIsInside(const ON_ClassArray<ON_BezierCurve> &loop_crvs, const ON_SimpleArray<int> &num_crvs_in_a_loop, const ON_2dPoint &uv, const ON_2dVector &uvdir, double tol_np, double tol_oncrv);

/// UVカーブを示す折れ線群から、uv点の内外判定を実施する。
/// @param [in] loop_pols LoopのUVカーブとなる折れ線群
/// @param [in] num_loop_pols 閉じた折れ線の数
/// @param [in] uv uv点
/// @param [in] tolerance uv空間内での位置トレランス
/// @return true:uv点はループの中、false:uv点はループの外
ONGEO_DECL bool ONGEO_UVPointIsInside(const ON_Polyline *loop_pols, int num_loop_pols, const ON_2dPoint &uv, double tolerance);

ONGEO_DECL ONGEO_SphereTree *ONGEO_NewSphereTree(int num, const ON_NurbsSurface *nbsurfs);
ONGEO_DECL void ONGEO_DeleteSphereTree(ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_CreateTree(ONGEO_SphereTree *, double radius2, double weight_rate, int level);
ONGEO_DECL void ONGEO_SphereTree_RayIntersectTest(const ONGEO_SphereTree *, const ON_3dRay &ray, ON_SimpleArray<ONGEO_SphereTree::Result> &results);
ONGEO_DECL int ONGEO_SphereTree_GetRootNodeIndex(const ONGEO_SphereTree *);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetFirstBezIndexFromBezIndex(const ONGEO_SphereTree *, int bez_index);
ONGEO_DECL int ONGEO_SphereTree_GetNurbsIntervalFromBezIndex(const ONGEO_SphereTree *, const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]);

ONGEO_DECL bool ONGEO_IntersectRayNurbs_Secant(const ON_3dRay &ray, const ON_NurbsSurface &srf, const ON_2dPoint &uv0, const ON_2dPoint &uv1, ON_3dPoint &tuv, ON_3dPoint &ptsrf, ON_3dPoint &ptlin, double tolerance, int max_iter = 50);

/// 多項式関数に値を代入して計算する。
/// @param [in] t 代入する値
/// @param [in] coef 多項式の係数列
/// @param [in] num 係数の個数(=次数+1)
/// @return 計算結果
ONGEO_DECL double ONGEO_Polynomial_Evaluate(double t, const double *coef, int num);

/// 2つの多項式の和を求める。
/// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
/// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
/// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
/// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
/// @param [out] coef_add 多項式1 + 多項式2を示す多項式の係数列
/// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
/// @return true:成功、 false:失敗
ONGEO_DECL bool ONGEO_Polynomial_Add(const double *coef1, int num1, const double *coef2, int num2, double *coef_add, int num_add);

/// 2つの多項式の差を求める。
/// @param [in] coef1 1つめの多項式の係数列 num1 == 0 のときに限り、coef1はヌルポインタでも可
/// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
/// @param [in] coef2 2つめの多項式の係数列 num2 == 0 のときに限り、coef2はヌルポインタでも可
/// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
/// @param [out] coef_sub 多項式1 - 多項式2を示す多項式の係数列
/// @param [in] num 出力多項式の係数の個数 (num1 <= num かつ num2 <= numであること  coef1 または coef2と同じ場所を指していても適切に計算可能)
/// @return true:成功、 false:失敗
ONGEO_DECL bool ONGEO_Polynomial_Subtract(const double *coef1, int num1, const double *coef2, int num2, double *coef_sub, int num_sub);

/// 2つの多項式の積を多項式として展開する。
/// @param [in] coef1 1つめの多項式の係数列
/// @param [in] num1 1つめの多項式の係数の個数(=次数+1)
/// @param [in] coef2 2つめの多項式の係数列
/// @param [in] num2 2つめの多項式の係数の個数(=次数+1)
/// @param [out] coef_mul 2つの多項式の積を展開した多項式の係数列
///	([num1+num2-2]次の多項式となるため、入力側でnum1+num2-1個分の係数の領域を確保する必要がある。)
ONGEO_DECL void ONGEO_Polynomial_Multiply(const double *coef1, int num1, const double *coef2, int num2, double *coef_mul);

/// 多項式の除算をし、商と剰余求める。
/// @param [in] coef1 割られる多項式の係数列
/// @param [in] num1 割られる多項式の係数の個数(=次数+1)
/// @param [in] coef2 割る多項式の係数列
/// @param [in] num2 割る多項式の係数の個数(=次数+1)
/// @param [out] coef_quot 商多項式の係数列 ... 剰余だけが必要な場合はヌルポインタを入れることができる
/// @param [out] num_quot 商多項式の係数の個数 (num_quot >= num1 であること ... 最大次数は coef1の非ゼロ最大次数 - coef2の非ゼロ最大次数 であり、 num1 - num2 とは限らない)
/// @param [out] coef_rem 剰余多項式の係数列
/// @param [out] num_rem 剰余多項式の係数の個数 (num_rem >= num1 かつ num2 - 1であること)
ONGEO_DECL void ONGEO_Polynomial_Divide(const double *coef1, int num1, const double *coef2, int num2, double *coef_quot, int num_quot, double *coef_rem, int num_rem);

/// 多項式を一階微分する。
/// @param [in] coef 多項式の係数列
/// @param [in] num 多項式の係数の個数(=次数+1)
/// @param [out] coef_dif 微分した多項式の係数列
///    (呼び出し元でnum-1個分の領域を確保すること
///     coefと同じ場所に結果を上書きしたい場合は、この引数として&coef[1]を渡し、本関数呼出し後にcoef[0] = 0.0とすること)
ONGEO_DECL void ONGEO_Polynomial_Differential(const double *coef, int num, double *coef_dif);

/// 多項式の係数列からSturm列を生成する。
/// @param coef [in] 多項式の係数列 (次数の大きい順)
/// @param num [in] 多項式の係数の数
/// @param s [out] Sturm列 (呼び出し元でnf*nf個分の領域を確保すること)
ONGEO_DECL void ONGEO_Polynomial_CreateSturmSequence(const double *coef, int num, double *strum);

/// Sturm列を用いて、指定した値での符号反転回数を計算する。
/// @param t [in] 値
/// @param sturm [in] Sturm列
/// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
/// @return 符号反転回数
ONGEO_DECL int ONGEO_Polynomial_NumberOfSignChangesOfSturmSequence(double t, const double *strum, int num);

/// Sturm列を用いて、指定した区間内の多項式の根の数を計算する。
/// @param t1 [in] 区間の下限
/// @param t2 [in] 区間の上限(t1 < t2 であること)
/// @param sturm [in] Sturm列
/// @param num [in] 多項式の係数の数 (Sturm列の配列長はnum*num)
/// @return 解の個数(エラーの場合、負数)
ONGEO_DECL int ONGEO_Polynomial_CalculateNumRoot(double t1, double t2, const double *strum, int num);

/// パスカルの三角形のdim+1行目を計算し、関数側で用意した配列を返す。
/// @param dim [in] 計算対象の次数
/// @return パスカルの三角形のdim+1行目 長さ dim+1の配列を返す。配列は関数が管理するため、delete[]等で解放してはならない。
ONGEO_DECL int *ONGEO_Polynomial_CalcPascalTriangle(int dim);

/// パスカルの三角形を使用して二項係数 mCn を計算する。
/// @param m, n [in] 係数
/// @return 二項係数
ONGEO_DECL int ONGEO_Polynomial_CalcBinomialCoef(int m, int n);

/// Brent法を用いて根を求める。区間 ti1 - ti2 の中で単調増加、または単調減少であること。
/// @param ti1 [in] 狭めたい区間の下限
/// @param ti2 [in] 狭めたい区間の上限 ti1 < ti2であること。
/// @param coef [in] 多項式の係数
/// @param num [in] 多項式の係数の数
/// @param tolerance [in] トレランス
/// @return 根となる値
ONGEO_DECL double ONGEO_Polynomial_FindRootByBrentMethod(double ti1, double ti2, const double *coef, int num, double tolerance);

/// 一変数の秋間補間法により補間値を求めるクラス
struct ONGEO_CLASS ONGEO_Interpolation_Akima_Univariate{
	struct Impl;
	Impl *pimpl;

	ONGEO_Interpolation_Akima_Univariate();
	~ONGEO_Interpolation_Akima_Univariate();

	/// 通過点列をセットする。xは必ず昇順であること。
	bool SetPoints(double *xa, double *ya, int num, int ydim = 1);

	/// xからyの補間値を求める。 y配列は呼び出し元で、ydim分だけ確保すること。
	/// 入力されたyをそのまま戻り値とする。
	double *Evaluate(double x, double y[]);

private:
	ONGEO_Interpolation_Akima_Univariate(const ONGEO_Interpolation_Akima_Univariate &);
	ONGEO_Interpolation_Akima_Univariate &operator =(const ONGEO_Interpolation_Akima_Univariate &);
};

ONGEO_DECL ONGEO_Interpolation_Akima_Univariate *ONGEO_Interpolation_Akima_Univariate_New();
ONGEO_DECL void ONGEO_Interpolation_Akima_Univariate_Delete(ONGEO_Interpolation_Akima_Univariate *);
ONGEO_DECL bool ONGEO_Interpolation_Akima_Univariate_SetPoints(ONGEO_Interpolation_Akima_Univariate *ths, double *xa, double *ya, int num, int ydim = 1);
ONGEO_DECL double *ONGEO_Interpolation_Akima_Univariate_Evaluate(ONGEO_Interpolation_Akima_Univariate *ths, double x, double y[]);

enum ONGEO_NI_ParameterMethod { // 曲線パラメータ [0 - 1] の計算方法
	ONGEO_NI_EquallySpaced, // 等間隔に配置
	ONGEO_NI_ChordLength,   // 制御点列の長さに合わせて配置、最も典型的な方法
	ONGEO_NI_Centripetal    // 制御点列の長さの平方根に合わせて配置、急激な曲がりを含む場合はこちらの方が有利
};
enum ONGEO_NI_Solver{
	ONGEO_NI_Solver_ON_Matrix_Inverse
};

/// 与えられた点列を通るNurbs曲線を作成する。
/// @param [in]  pts     点列
/// @param [in]  pt_cnt  点列の個数
/// @param [in]  order   Nurbs曲線の階数
/// @param [in]  pmethod 曲線パラメータの計算方法
/// @param [put] nc      作成した Nurbs曲線
/// @param [in]  solver  連立一次方程式をどの方法で解くか？
/// @return true:成功、false:失敗
ONGEO_DECL bool ONGEO_FitNurbsCurveToPointArray(const ON_3dPoint *pts, int pt_cnt, int order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsCurve &nc, ONGEO_NI_Solver solver = ONGEO_NI_Solver_ON_Matrix_Inverse);

/// 与えられた点グリッドを通るNurbs曲面を作成する。
/// @param [in]  pts     点グリッド (P_u0v0、 P_u1v0、P_u2v0、...、P_u0v1、P_u1v1、... の順に並べる)
/// @param [in]  u_count 点グリッドのU方向の個数
/// @param [in]  u_count 点グリッドのV方向の個数
/// @param [in]  u_order Nurbs曲面のU方向階数
/// @param [in]  v_order Nurbs曲面のV方向階数
/// @param [in]  pmethod 曲線パラメータの計算方法
/// @param [put] nf      作成した Nurbs曲面
/// @param [in]  solver  連立一次方程式をどの方法で解くか？
/// @return true:成功、false:失敗
ONGEO_DECL bool ONGEO_FitNurbsSurfaceToPointGrid(const ON_3dPoint *pts, int u_count, int v_count, int u_order, int v_order, ONGEO_NI_ParameterMethod pmethod, ON_NurbsSurface &nf, ONGEO_NI_Solver solver = ONGEO_NI_Solver_ON_Matrix_Inverse);

/// IGESデータを扱うクラス
struct ONGEO_CLASS ONGEO_IgesModel{
	ON_String ss;
	struct GlobalSection{
		ON_String parm_delimiter;
		ON_String record_delimiter;
		ON_String product_id_sending;
		ON_String file_name;
		ON_String native_system_id;
		ON_String preprocessor_version;
		int num_of_binbits_for_intrep;
		int max_pow_float;
		int num_of_digits_float;
		int max_pow_double;
		int num_of_digits_double;
		ON_String product_id_receiving;
		double model_scale;
		int unit_flag;
		ON_String unit_name;
		int max_num_lineweight_grad;
		double width_of_max_line_weight;
		ON_String timestamp_filegen;
		double min_resolution_in_unit;
		double approx_max_cod_in_unit;
		ON_String name_author;
		ON_String authors_org;
		int flag_version;
		int flag_draft_std;
		ON_String timestamp_filemod;
		ON_String desc_protocol;
	}gs;
	struct DirectoryEntrySection{
		int entity_type;
		int param_data;
		int structure;
		int line_font;
		int level;
		int view;
		int trans_matrix;
		int label_disp;
		struct{
			int blank_status:4;  // 0:visible, 1:blancked
			int subord_ent_sw:4; // 0:independent, 1:physically dependent, 2:logically dependent, 3: both 1 and 2
			int ent_use_flg:4;   // 0:geometry, 1:annotation, 2:definition, 3:other, 4:logical/positional, 5:2D parametric, 6:construction geometry
			int hierarchy:4;     // 0:global top down, 1:global defer, 2:use hierarchy property
		}stnum;
		int line_weight;
		int color_num;
		int param_line_count;
		int form_num;
		int reserved[2];
		char ent_label[9];
		int ent_subscript;
	};
	ON_SimpleArray<DirectoryEntrySection> des;
	ON_ClassArray<ON_String> ps;

	ONGEO_IgesModel();
	ONGEO_IgesModel(const char *filename);
	ONGEO_IgesModel(const ONGEO_IgesModel &rhs);
	ONGEO_IgesModel &operator =(const ONGEO_IgesModel &rhs);
	bool Load(const char *filename);
	bool Save(const char *filename);
	void Clear();

	/// 指定したインデックスのEntityを空にする(NULLエンティティ化する)。
	/// @param [in] index Directory Entry配列のインデックス
	/// @return true:成功、false:失敗
	bool SetEntityClearOut(int index);

	~ONGEO_IgesModel();
};

ONGEO_DECL ONGEO_IgesModel *ONGEO_NewIgesModel();
ONGEO_DECL ONGEO_IgesModel *ONGEO_NewIgesModel(const char *filename);
ONGEO_DECL void ONGEO_DeleteIgesModel(ONGEO_IgesModel *);
ONGEO_DECL bool ONGEO_IgesModel_Load(ONGEO_IgesModel *, const char *filename);
ONGEO_DECL bool ONGEO_IgesModel_Save(ONGEO_IgesModel *, const char *filename);
ONGEO_DECL void ONGEO_IgesModel_Clear(ONGEO_IgesModel *);

/// 指定したインデックスのEntityを空にする(NULLエンティティ化する)。
/// @param [in] index Directory Entry配列のインデックス
/// @return true:成功、false:失敗
ONGEO_DECL void ONGEO_IgesModel_SetEntityClearOut(ONGEO_IgesModel *, int index);


/// 3dm要素とIges要素との対応情報
struct ONGEO_CLASS ONGEO_IgesTo3dmInfo{
	ONGEO_IgesTo3dmInfo();
	~ONGEO_IgesTo3dmInfo();
	struct Impl;
	Impl *pimpl;

	int DEIndexFromObject(const ON_Object *) const;
private:
	ONGEO_IgesTo3dmInfo(const ONGEO_IgesTo3dmInfo &);
	ONGEO_IgesTo3dmInfo &operator =(const ONGEO_IgesTo3dmInfo &);
};

ONGEO_DECL ONGEO_IgesTo3dmInfo *ONGEO_NewIgesTo3dmInfo();
ONGEO_DECL void ONGEO_DeleteIgesTo3dmInfo(ONGEO_IgesTo3dmInfo *);
ONGEO_DECL int ONGEO_IgesTo3dmInfo_DEIndexFromObject(const ONGEO_IgesTo3dmInfo *, const ON_Object *);

/// IGESデータを読み込み、ONX_Modelに格納する
ONGEO_DECL bool ONGEO_IgesTo3dm(const ONGEO_IgesModel &igs, ONX_Model &onx, ONGEO_IgesTo3dmInfo &info);

/// バイナリSTLを読み込み、ON_Meshに格納する
ONGEO_DECL bool ONGEO_ReadBinarySTL(ON_BinaryArchive &ba, ON_Mesh &mesh, ON_String &header);

/// 一行分のテキスト文字列を読み込む
ONGEO_DECL void ONGEO_ReadLine(ON_BinaryArchive &ba, ON_String &str);

/// テキストSTLを読み込み、ON_Meshに格納する
ONGEO_DECL bool ONGEO_ReadTextSTL(ON_BinaryArchive &ba, ON_Mesh &mesh, ON_String &modelname);

/// 与えられた点列から、DeWall法でドロネー三角形分割を施す
ONGEO_DECL bool ONGEO_DelaunayTriangulation_2D_DeWall(const ON_2dPoint *pts, int num_pts, ON_SimpleArray<int> &simplexes);
ONGEO_DECL bool ONGEO_DelaunayTriangulation_2D_DeWall(const double *pts, int num_pts, ON_SimpleArray<int> &simplexes);
