/*
 * DelaunayTriangulation_DeWall
 * Copylight (C) 2013 mocchi
 * mocchi_2003@yahoo.co.jp
 * License: Boost ver.1
 */

// Reference : P. Cignoni, C. Montani, R. Scopigno,
//   "DeWall: A Fast Divide & Conquer Delaunay Triangulation Algorithm in E^d".
//     Computer Aided Design 1998, Vol 30, Pages 333–341

#define NOMINMAX
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <complex>
#include <limits>
#include <map>

#include "ONGEO.h"
#include "Profile.h"


template <int N> struct traits{
};

template <> struct traits<3>{
	typedef ON_2dPoint ON_Point;
	typedef ON_2dVector ON_Vector;
	static int FII(int j, int i){
		static int fii[3][2] = {{0, 1}, {0, 2}, {1, 2}};
		return fii[j][i];
	}
	static int PII(int i){
		static int pii[3] = {2, 1, 0};
		return pii[i];
	}
	template<typename Ary, typename T>
	static void MakeFirstSimplex(const T &ptA, const T &ptB, ON_SimpleArray<int> &vids_l, ON_SimpleArray<int> &vids_r, const Ary &pts, int idx[3]){
		// 外接円の半径が最小になる点の検出
		int imin, side_min = -1;
		double r_min = std::numeric_limits<double>::max();
		for (int j = 0; j < 2; ++j){
			ON_SimpleArray<int> &vids_i = (j == 0) ? vids_l : vids_r;
			for (int i = 0; i < vids_i.Count(); ++i){
				int index = vids_i[i];
				if (idx[0] == index || idx[1] == index) continue;
				ON_2dPoint pt = pts[index];
				double r = ON_Circle(ptA, ptB, pt).radius;
				if (r_min > r) r_min = r, idx[0] = index, imin = i, side_min = j;
			}
		}
		if (side_min < 0) return;
	}

};

template <> struct traits<4>{
	typedef ON_3dPoint ON_Point;
	typedef ON_3dVector ON_Vector;
	static int FII(int j, int i){
		static int fii[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
		return fii[j][i];
	}
	static int PII(int i){
		static int pii[4] = {3, 2, 1, 0};
		return pii[i];
	}
	template<typename T, typename Ary>
	static void MakeFirstSimplex(const ON_3dPoint &ptA, const ON_3dPoint &ptB, ON_SimpleArray<int> &vids_l, ON_SimpleArray<int> &vids_r, const Ary &pts, int idx[4]){
		traits<3>::MakeFirstSimplex<Ary, ON_Point>(ptA, ptB, vids_l, vids_r, pts, idx);
		// Todo: idx[2]
	}
};

template <int D> struct Region{
	struct Face{
		int fids[D-1];
		void Arrange(){
			for (int j = 0; j < D-1; ++j){
				for (int i = j+1; i < D-1; ++i){
					if (fids[j] > fids[i]) std::swap(fids[j], fids[i]);
				}
			}
		}
		bool operator <(const Face &rhs) const{
			for (int i = 0; i < D-1; ++i){
				if (fids[i] == rhs.fids[i]) continue;
				return (fids[i] < rhs.fids[i]);
			}
			return false;
		}
		bool IsVertexContained(int vid){
			for (int i = 0; i < D-1; ++i) if (fids[i] == vid) return true;
			return false;
		}
	};
	void Arrange(Face &f);
	ON_SimpleArray<int> vids;
	std::map<Face, typename traits<D>::ON_Vector> fids;
	int wall_dim;
	void Swap(Region &rhs){
		ONGEO_Swap(vids, rhs.vids);
		std::swap(fids, rhs.fids);
		std::swap(wall_dim, rhs.wall_dim);
	}
};

// fを無限に伸ばした平面で空間を2分割したとき、
// 球の中心が入った半空間とptが含まれる半空間が同じ時は正、異なる時は負
template <typename Ary>
double DelaunayDistance(Region<3>::Face &f, const ON_2dVector &nrm, ON_2dPoint &pt, Ary &pts){
	ON_2dPoint v[2] = { pts[f.fids[0]], pts[f.fids[1]] };
	ON_2dVector n = nrm * (-1);
	ON_Circle cir(v[0], v[1], pt);
	if (ON_DotProduct(ON_2dPoint(cir.Center()) - v[0], n) < 0) return cir.radius * -1.0;
	return cir.radius;
}

// fを無限に伸ばした平面で空間を2分割したとき、
// 球の中心が入った半空間とptが含まれる半空間が同じ時は正、異なる時は負
template <typename Ary>
double DelaunayDistance(Region<4>::Face &f, const ON_3dVector &nrm, ON_3dPoint &pt, Ary &pts){
	// Todo:
	// ・ptの方向を向くfの法線方向を求める。
	// ・外接球を求める。
	// ・fのいずれかの点から球の中心へ向かうベクトルとの内積を求める。
	// ・外接球の半径に内積の符号をつけてreturnする。
	return -1;
}

template<typename T>
void CalcFaceNormal(const T &pts, const Region<3>::Face &fs, int vidx, ON_2dVector &nrm){
	nrm.PerpendicularTo(pts[fs.fids[1]] - pts[fs.fids[0]]);
	if (ON_DotProduct(nrm, pts[vidx] - pts[fs.fids[0]]) < 0) nrm *= -1;
}

template <int N, typename T>
bool DelaunayTriangulation_DeWall(T &pts, int num_pts, ON_SimpleArray<int> &simplexes){
	int wall_dim = 0;
	double wall_pos;

	ON_ClassArray<Region<N> > stack;
	Region<N> &root = stack.AppendNew();
	root.wall_dim = 0;
	root.vids.SetCapacity(num_pts);
	root.vids.SetCount(num_pts);
	for (int i = 0; i < root.vids.Count(); ++i) root.vids[i] = i;

	while(stack.Count()){
		Region<N> region_c, region_l, region_r;
		region_c.Swap(*stack.Last());
		stack.Remove(stack.Count()-1);

		ON_SimpleArray<int> &vids_c = region_c.vids;
		ON_SimpleArray<int> &vids_l = region_l.vids, &vids_r = region_r.vids;
		std::map<Region<N>::Face, traits<N>::ON_Vector> &fids_l = region_l.fids, &fids_r = region_r.fids;
		int &wall_dim = region_c.wall_dim;

		ON_BoundingBox bb;
		for (int i = 0; i < vids_c.Count(); ++i){
			traits<3>::ON_Point pt = pts[vids_c[i]];
			bb.Set(pt, true);
		}
		wall_pos = (bb.m_max[wall_dim] + bb.m_min[wall_dim]) * 0.5;

		double dist_min = std::numeric_limits<double>::max();
		int side_min = -1;

		std::map<Region<N>::Face, traits<N>::ON_Vector> fids_w;
		{
			traits<N>::ON_Point vtx[N];
			int idx[N];
			// Pointset_Partition
			{
				int imin;
				// 空間の分割とwallに最も近い点の検出
				for (int i = 0; i < vids_c.Count(); ++i){
					int side, ii;
					int index = vids_c[i];
					traits<3>::ON_Point pt = pts[index];
					if (pt[wall_dim] <= wall_pos) vids_l.Append(vids_c[i]), side = 0, ii = vids_l.Count()-1;
					else vids_r.Append(vids_c[i]), side = 1, ii = vids_r.Count()-1;
					double dist = std::abs(wall_pos - pt[wall_dim]);
					if (dist_min > dist) dist_min = dist, side_min = side, idx[0] = index, imin = ii;
				}
				vtx[0] = pts[idx[0]];
			}

			bool add_simplex = false;

			// Make_FirstSimplex
			if (region_c.fids.size() == 0){
				add_simplex = true;
				int imin;
				// 反対側でidx_minに最も近い点の検出
				ON_SimpleArray<int> &vids_opposite = (side_min == 0) ? vids_r : vids_l;
				double dist2_min = std::numeric_limits<double>::max();
				for (int i = 0; i < vids_opposite.Count(); ++i){
					int index = vids_opposite[i];
					double dist2 = (pts[index]-vtx[0]).LengthSquared();
					if (dist2_min > dist2) dist2_min = dist2, idx[1] = index, imin = i;
				}
				vtx[1] = pts[idx[1]];

				traits<N>::MakeFirstSimplex<T, traits<N>::ON_Point>(vtx[0], vtx[1], vids_l, vids_r, pts, idx);
				for (int j = 2; j < N; ++j) vtx[j] = pts[idx[j]];

				// Simplexを追加
				for (int j = 0; j < N; ++j) simplexes.Append(idx[j]);
			}else{
				ON_SimpleArray<int> sii(num_pts);
				sii.SetCount(sii.Capacity());
				// region_cの中にある全ての辺について、新たに設定したwallで分別する。
				for (std::map<Region<N>::Face, traits<N>::ON_Vector>::iterator iter = region_c.fids.begin(); iter != region_c.fids.end(); ++iter){
					const Region<N>::Face &f = iter->first;
					for (int i = 0; i < N-1; ++i){
						int idx = f.fids[i];
						if (sii[idx]) continue;
						double v = pts[idx][wall_dim] - wall_pos;
						sii[idx] = (v < 0) ? 1 : ((v > 0) ? 2 : 3);
					}
				}
				for (std::map<Region<N>::Face, traits<N>::ON_Vector>::iterator iter = region_c.fids.begin(); iter != region_c.fids.end(); ++iter){
					const Region<N>::Face &f = iter->first;
					int side = 0;
					for (int i = 0; i < N-1; ++i) side |= sii[f.fids[i]];

					std::map<Region<N>::Face, traits<N>::ON_Vector> &fids = (side == 1) ? fids_l : ((side == 2) ? fids_r : fids_w);
					fids.insert(std::make_pair(f, iter->second));
				}
			}

			// Simplexを構成する各(D-1)-faceについて、Wallと交わる、または全ての頂点が同一領域に属するかで振り分け
			int fstart = 0;

			// fids_wが空になるまで繰り返し
			for(;;){
				if (add_simplex){

					// vtx、idx には最近作成したsimplexの各頂点情報が格納されている。
					int sii[N] = {0};
					for (int j = 0; j < N; ++j){
						double v = vtx[j][wall_dim] - wall_pos;
						sii[j] = (v < 0) ? 1 : ((v > 0) ? 2 : 3);
					}
					ON_SimpleArray<Region<N>::Face> fs;
					ON_SimpleArray<int> sides;
					ON_ClassArray<traits<N>::ON_Vector> fo; // face の方向(simplex内部へ向いているベクトル)
					// Make_FirstSimplexの場合はSimplexを構成する全てのfが対象。
					// Make_Simplexの場合は、Simplexを構成するfのうち、最初の1つはfids_wから
					//      取り出したものであるため対象外
					for (int j = fstart; j < N; ++j){
						Region<N>::Face &f = fs.AppendNew();
						int &side = sides.AppendNew();
						side = 0;
						int fii;
						for (int i = 0; i < N-1; ++i) fii = traits<N>::FII(j,i), f.fids[i] = idx[fii], side |= sii[fii];
						CalcFaceNormal(pts, f, idx[traits<N>::PII(j)], fo.AppendNew());
						f.Arrange();
					}

					// 追加した simplex を構成する各 face について、対応する fids へ Update
					for (int j = 0; j < fs.Count(); ++j){
						int &side = sides[j];
						Region<N>::Face f = fs[j];
						traits<N>::ON_Vector &nrm = fo[j];
						std::map<Region<N>::Face, traits<N>::ON_Vector> &fids = (side == 1) ? fids_l : ((side == 2) ? fids_r : fids_w);
						std::map<Region<N>::Face, traits<N>::ON_Vector>::iterator iter = fids.find(f);
						// faceが最大2つのSimplexに共有される。よって、既にどこかで使われている
						// faceが見つかった場合は他のSimplexが見つかる可能性がないため、
						// setから削除
						if (iter == fids.end()) fids.insert(std::make_pair(f, nrm));
						else fids.erase(iter);
					}
				}

				if (fids_w.size() == 0) break;

				// fids_wのfaceを1つ取り出し、新たなSimplexを作成する。
				std::map<Region<N>::Face, traits<N>::ON_Vector>::iterator iter_fids_w = fids_w.begin();
				Region<N>::Face f = iter_fids_w->first;
				traits<N>::ON_Vector nrm = iter_fids_w->second;
				traits<N>::ON_Point fv0 = pts[f.fids[0]];
				fids_w.erase(iter_fids_w);

				for (int i = 0; i < N-1; ++i) vtx[i] = pts[f.fids[i]]; 

				// Make_Simplex
				fstart = 1;
				int imin, side_min = -1;
				double r_min = std::numeric_limits<double>::max();
				for (int j = 0; j < 2; ++j){
					ON_SimpleArray<int> &vids_i = (j == 0) ? vids_l : vids_r;
					if (vids_i.Count() == 0) continue;
					for (int i = 0; i < vids_i.Count(); ++i){
						int index = vids_i[i];
						if (f.IsVertexContained(index)) continue;
						traits<N>::ON_Point pt = pts[index];
						if (ON_DotProduct(pt-fv0, nrm) >= 0) continue;

						double r = DelaunayDistance(f, nrm, pt, pts);
						if (r_min > r) r_min = r, idx[N-1] = index, imin = i, side_min = j;
					}
				}

				add_simplex = (side_min != -1);
				if (add_simplex){
					for (int i = 0; i < N-1; ++i) idx[i] = f.fids[i], vtx[i] = pts[idx[i]];

					ON_SimpleArray<int> &vids_i = (side_min == 0) ? vids_l : vids_r;
					vtx[N-1] = pts[idx[N-1]];
					// Simplexを追加
					for (int j = 0; j < N; ++j) simplexes.Append(idx[j]);
				}
			}

			int wall_dim_c = (wall_dim + 1) % 2;
			if (region_l.fids.size()){
				region_l.wall_dim = wall_dim_c;
				Region<N> &r = stack.AppendNew();
				r.Swap(region_l);
			}
			if (region_r.fids.size()){
				region_r.wall_dim = wall_dim_c;
				Region<N> &r = stack.AppendNew();
				r.Swap(region_r);
			}
		}
	}
	return true;
}

bool ONGEO_DelaunayTriangulation_2D_DeWall(const ON_2dPoint *pts, int num_pts, ON_SimpleArray<int> &simplexes){
	PROF("DelaunayTriangulation_2D_DeWall(ON_2dPoint)");
	return DelaunayTriangulation_DeWall<3, const ON_2dPoint *>(pts, num_pts, simplexes);
}

bool ONGEO_DelaunayTriangulation_2D_DeWall(const double *pts, int num_pts, ON_SimpleArray<int> &simplexes){
	PROF("DelaunayTriangulation_2D_DeWall(double array)");
	struct Access{
		const double *pts;
		Access(const double *pts_) : pts(pts_){}
		ON_2dPoint operator [](int idx) const{
			return ON_2dPoint(pts + idx * 2);
		}
	}ac(pts);
	return DelaunayTriangulation_DeWall<3, const Access>(ac, num_pts, simplexes);
}
