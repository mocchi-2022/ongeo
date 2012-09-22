// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#define NOMINMAX
#include "ONGEO.h"
#include <stack>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "Profile.h"

namespace{
	void CreateInitialTree(const ON_SimpleArray<ON_BezierSurface> &bezs, ON_SimpleArray<ONGEO_SphereTree::Node> &nodes){
		nodes.Destroy();
		for (int i = 0; i < bezs.Count(); ++i){
			ONGEO_SphereTree::Node node;
			ONGEO_CalculateRoughBoundingSphere(bezs[i], node.center, node.radius2);
			node.bez_index = i;
			node.direction = bezs[i].IsValid() ? -1 : -3;
			node.childLeft = node.childRight = -1;
			nodes.Append(node);
		}
		if (bezs.Count() == 1){
			return;
		}
		std::vector<int> cbuf(bezs.Count()), tbuf(bezs.Count());
		for (int i = 0; i < bezs.Count(); ++i) cbuf[i] = i;
		int num_errnode = 0;
		for (int i = 0; i < bezs.Count(); ++i){
			if (nodes[i].direction == -3){
				std::swap(cbuf[i], cbuf[num_errnode]);
				++num_errnode;
			}
		}

		std::stack<std::pair<int, std::pair<int, int> > > ranges;
		ranges.push(std::make_pair(-1, std::make_pair(num_errnode, bezs.Count()-1)));
		while(ranges.size()){
			ON_BoundingBox bb;
			int pnode = ranges.top().first;
			std::pair<int, int> range = ranges.top().second;
			ranges.pop();

			for (int i = range.first; i <= range.second; ++i){
				bb.Union(bezs[cbuf[i]].BoundingBox());
			}
			ONGEO_SphereTree::Node node_cur;
			node_cur.center = bb.Center();
			node_cur.radius2 = (bb.m_min-node_cur.center).LengthSquared();
			node_cur.bez_index = -2;
			node_cur.direction = -2;
			node_cur.childLeft = node_cur.childRight = -1;
			if (pnode >= 0){
				if (nodes[pnode].childLeft == -1) nodes[pnode].childLeft = nodes.Count();
				else nodes[pnode].childRight = nodes.Count();
			}
			if (range.second - range.first == 1){
				node_cur.childLeft = cbuf[range.first];
				node_cur.childRight = cbuf[range.second];
				nodes.Append(node_cur);
				continue;
			}

			double Abuf[9] = {0};
			double *A[] = {&Abuf[0], &Abuf[3], &Abuf[6]};
			// 分散共分散行列を求める。
			ON_3dVector mean = ON_3dVector(0,0,0);
			// まずは各成分の平均
			for (int i = range.first; i <= range.second; ++i){
				mean += nodes[cbuf[i]].center;
			}
			double denom = static_cast<double>(range.second - range.first + 1);
			mean /= denom;
			// 分散・共分散行列作成
			for (int i = range.first; i <= range.second; ++i){
				ON_3dVector dif = nodes[cbuf[i]].center - mean;
				A[0][0] += dif.x * dif.x;
				A[0][1] += dif.x * dif.y;
				A[0][2] += dif.x * dif.z;
				A[1][1] += dif.y * dif.y;
				A[1][2] += dif.y * dif.z;
				A[2][2] += dif.z * dif.z;
			}
			A[0][0] /= denom, A[0][1] /= denom, A[0][2] /= denom;
			A[1][1] /= denom, A[1][2] /= denom, A[2][2] /= denom;

			// 相関行列に変換
			double s[3] = {std::sqrt(A[0][0]), std::sqrt(A[1][1]), std::sqrt(A[2][2])};
			A[0][0] = (s[0] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
			A[1][1] = (s[1] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
			A[2][2] = (s[2] < ON_ZERO_TOLERANCE) ? 0.0 : 1.0;
			double ss;
			A[0][1] = ((ss = s[0]*s[1]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][1] /= ss;
			A[0][2] = ((ss = s[0]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[0][2] /= ss;
			A[1][2] = ((ss = s[1]*s[2]) < ON_ZERO_TOLERANCE) ? 0.0 : A[1][2] /= ss;

			A[1][0] = A[0][1], A[2][0] = A[0][2], A[2][1] = A[1][2];

			// 固有値、固有ベクトル計算
			double l[6];
			ONGEO_EigenValue3_Cardano(A, l);
			// 絶対値が最も大きい固有値を探す
			double lmax = l[0];
			double lamax = std::abs(lmax);
			for (int i = 2; i < 6; i += 2){
				double la = std::abs(l[i]);
				if (lamax < la) lmax = l[i], lamax = la;
			}
			ON_3dVector v;
			ONGEO_EigenVector3(A, lmax, v);
			v.Unitize();

			for (int i = range.first; i <= range.second; ++i){
				tbuf[i] = cbuf[i];
			}

			int g1 = range.first, g2 = range.second;

			// 内積の正負で振り分け
			for (int i = range.first; i <= range.second; ++i){
				int index = tbuf[i];
				double t = ON_DotProduct(v, nodes[index].center - node_cur.center);
				if (t >= 0) cbuf[g1++] = index;
				else cbuf[g2--] = index;
			}
			if (g1 == range.second + 1){
				--g1;
				g2 = g1 - 1;
			}else if (g2 == range.first - 1){
				++g2;
				g1 = g2 + 1;
			}
			if (g1 - range.first == 1){
				node_cur.childLeft = cbuf[range.first];
			}else if (g1 > range.first){
				ranges.push(std::make_pair(nodes.Count(), std::make_pair(range.first, g1-1)));
			}
			if (g1 == range.second){
				node_cur.childRight = cbuf[range.second];
			}else if(g2 < range.second){
				ranges.push(std::make_pair(nodes.Count(), std::make_pair(g1, range.second)));
			}
			nodes.Append(node_cur);
		}
	}
}

int ONGEO_CalculateRoughBoundingSphere(const ON_BezierSurface &src, ON_3dPoint &center, double &radius2){
	ON_BoundingBox bb = src.BoundingBox();
	center = bb.Center();
	radius2 = (bb.m_min-center).LengthSquared();
	return 0;
}

void ONGEO_CalculateMinMaxWeight(const ON_BezierSurface &src, double &wmin, double &wmax){
	wmax = wmin = src.Weight(0,0);
	for (int j = 0; j < src.Order(1); ++j){
		for (int i = 0; i < src.Order(0); ++i){
			wmax = std::max(wmax, src.Weight(i,j));
			wmin = std::min(wmin, src.Weight(i,j));
		}
	}
}

int ONGEO_SphereTree::GetRootNodeIndex() const{
	if (bezs.Count() == 0) return -1;
	if (bezs.Count() == 1) return 0;
	else return bezs.Count();
}

ONGEO_SphereTree::ONGEO_SphereTree(int num, const ON_NurbsSurface *nbsurfs){
	nurbs_index_to_bez_first_index.SetCapacity(num);
	nurbs_index_to_bez_first_index.SetCount(num);
	for (int k = 0; k < num; ++k){
		const ON_NurbsSurface &nbsurf = nbsurfs[k];
		nurbs_index_to_bez_first_index[k] = bezs.Count();
		for (int sj = 0; sj <= nbsurf.m_cv_count[1] - nbsurf.m_order[1]; ++sj){
			for (int si = 0; si <= nbsurf.m_cv_count[0] - nbsurf.m_order[0]; ++si){
				ON_BezierSurface bezroot;
				nbsurf.ConvertSpanToBezier(si, sj, bezroot);
				bezs.Append(bezroot);
			}
		}
	}
	num_root_bezsurfs = bezs.Count();
	CreateInitialTree(bezs, nodes);
}
ONGEO_SphereTree::~ONGEO_SphereTree(){
}

int ONGEO_SphereTree::CreateTree(double radius2, double weight_rate, int level){
	if (radius2 <= 0 && weight_rate <= 0 && level < 0) return 1;

	CreateInitialTree(bezs, nodes);
	const static double zero2 = ON_ZERO_TOLERANCE * ON_ZERO_TOLERANCE;

	struct StackItem{
		ON_BezierSurface bez;
		int index;
		int level;
	};
	std::stack<StackItem> stack;
	for (int i = 0; i < bezs.Count(); ++i){
		StackItem siroot;
		siroot.bez = bezs[i];
		siroot.index = i;
		siroot.level = 0;
		if (!siroot.bez.IsValid()) continue;
		stack.push(siroot);
	}

	while(stack.size()){
		StackItem sicur = stack.top();
		stack.pop();
		ON_BezierSurface &bez = sicur.bez;
		Node &cnode = nodes[sicur.index];

		if (cnode.radius2 < zero2) continue;
		bool b_radius2 = false;
		if (radius2 == 0 || cnode.radius2 < radius2) b_radius2 = true;
		bool b_weight_rate = false;
		if (weight_rate > 0){
			double wmin, wmax;
			ONGEO_CalculateMinMaxWeight(bez, wmin, wmax);
			if (weight_rate > wmax / wmin) b_weight_rate = true;
		}else b_weight_rate = true;
		bool b_level = false;
		if (level < 0 || level == sicur.level) b_level = true;
		if (b_radius2 && b_weight_rate && b_level) continue;

		double ulen2 = 0, vlen2 = 0;
		ON_3dPointArray q;
		int qcount = bez.Order(1)*bez.Order(0);
		q.SetCapacity(qcount);
		q.SetCount(qcount);

		for (int j = 0, ji = 0; j < bez.Order(1); ++j, ji += bez.Order(0)){
			for (int i = 0; i < bez.Order(0); ++i){
				bez.GetCV(i,j,q[ji+i]);
			}
		}

		// u方向とv方向の制御点間の長さを計算
		for (int j = 0, ji = 0; j < bez.Order(1); ++j, ji += bez.Order(0)){
			for (int i = 0; i < bez.Order(0)-1; ++i){
				ulen2 += (q[ji+i+1]-q[ji+i]).LengthSquared();
			}
		}
		for (int i = 0; i < bez.Order(0); ++i){
			for (int j = 0, ji = 0; j < bez.Order(1)-1; ++j, ji += bez.Order(0)){
				vlen2 += (q[ji+i+bez.Order(0)]-q[ji+i]).LengthSquared();
			}
		}

		// ulen2とvlen2の大きいほうを半分に分割する。
		StackItem sileft, siright;
		cnode.direction = (ulen2 > vlen2) ? 0 : 1;
		bez.Split(cnode.direction, 0.5, sileft.bez, siright.bez);
		sileft.level = siright.level = sicur.level + 1;

		Node nleft, nright;
		ONGEO_CalculateRoughBoundingSphere(sileft.bez, nleft.center, nleft.radius2);
		ONGEO_CalculateRoughBoundingSphere(siright.bez, nright.center, nright.radius2);
		nleft.bez_index = nright.bez_index = cnode.bez_index;
		nleft.direction = nright.direction = -1;
		nleft.childLeft = nleft.childRight = -1;
		nright.childLeft = nright.childRight = -1;

		cnode.childLeft = sileft.index = nodes.Count();
		cnode.childRight = siright.index = nodes.Count()+1;
		nodes.Append(nleft);
		nodes.Append(nright);

		stack.push(sileft);
		stack.push(siright);
	}

	return 0;
}

int ONGEO_SphereTree::GetNurbsIndexFromBezIndex(int bez_index) const{
	if (bez_index >= num_root_bezsurfs) return -1;
	const int *begin = &nurbs_index_to_bez_first_index[0], *end = begin + nurbs_index_to_bez_first_index.Count();
	const int *p = std::upper_bound(begin, end, bez_index);
	return p - begin - 1;
}
int ONGEO_SphereTree::GetFirstBezIndexFromBezIndex(int bez_index) const{
	int nbindex = GetNurbsIndexFromBezIndex(bez_index);
	return nurbs_index_to_bez_first_index[nbindex];
}

void ONGEO_SphereTree::RayIntersectTest(const ON_3dRay &ray, ON_SimpleArray<Result> &results) const{
	PROF("ONGEO_SphereTree::RayIntersectTest");
	results.Destroy();
	std::stack<Result> si;

#if 0
	for (int i = 0; i < bezs.Count(); ++i){
		Result rroot;
		rroot.bez_index = i;
		rroot.node_index = i;
		rroot.uint.Set(0,1);
		rroot.vint.Set(0,1);
		si.push(rroot);
	}
#endif
	Result rroot;
	rroot.node_index = GetRootNodeIndex();
	rroot.bez_index = rroot.node_index != 0 ? -2 : 0;
	rroot.uint.Set(0,1);
	rroot.vint.Set(0,1);
	if (rroot.node_index < 0) return;
	si.push(rroot);
	while(si.size()){
		Result rcur = si.top();
		si.pop();
		const Node &cnode = nodes[rcur.node_index];
#if 0
		if (cnode.radius2 < 0){
			results.Append(rcur);
			continue;
		}
#endif
		ON_3dVector vp = cnode.center - ray.m_P;
		double tn = ON_DotProduct(ray.m_V, vp);
		double tn2 = tn * tn;
		double dist2 = vp.LengthSquared() - tn2;
		if (dist2 > cnode.radius2) continue;
		if (tn < 0 && tn2 > cnode.radius2) continue;

		if (cnode.direction == -1){
			results.Append(rcur);
		}else{
			Result rleft, rright;
			if (cnode.bez_index == -2){
				rleft.bez_index = (cnode.childLeft >= bezs.Count()) ? -2 : cnode.childLeft;
				rright.bez_index = (cnode.childRight >= bezs.Count()) ? -2 : cnode.childRight;
			}else{
				rleft.bez_index = rright.bez_index = cnode.bez_index;
			}
			rleft.node_index = cnode.childLeft;
			rright.node_index = cnode.childRight;
			if (cnode.direction == 0){
				rleft.uint.m_t[0] = rcur.uint.m_t[0];
				rleft.uint.m_t[1] = rright.uint.m_t[0] = rcur.uint.Mid();
				rright.uint.m_t[1] = rcur.uint.m_t[1];
				rleft.vint = rright.vint = rcur.vint;
			}else if (cnode.direction == 1){
				rleft.vint.m_t[0] = rcur.vint.m_t[0];
				rleft.vint.m_t[1] = rright.vint.m_t[0] = rcur.vint.Mid();
				rright.vint.m_t[1] = rcur.vint.m_t[1];
				rleft.uint = rright.uint = rcur.uint;
			}else{
				rleft.uint = rright.uint = rcur.uint;
				rleft.vint = rright.vint = rcur.vint;
			}
			si.push(rleft);
			si.push(rright);
		}
	}
}

int ONGEO_SphereTree::GetNurbsIntervalFromBezIndex(const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]) const{
	if (!nbsurfs) return -1;
	int nbsurf_index = GetNurbsIndexFromBezIndex(bez_index);
	const ON_NurbsSurface &nbsurf = nbsurfs[nbsurf_index];
	int sstride = nbsurf.m_cv_count[0] - nbsurf.m_order[0] + 1;
	int bez_index_of_parent_nbsurf = bez_index - GetFirstBezIndexFromBezIndex(bez_index);
	int si = bez_index_of_parent_nbsurf % sstride, sj = bez_index_of_parent_nbsurf / sstride;

	range[0].Set(nbsurf.m_knot[0][si+nbsurf.m_order[0]-2], nbsurf.m_knot[0][si+nbsurf.m_order[0]-1]);
	range[1].Set(nbsurf.m_knot[1][sj+nbsurf.m_order[1]-2], nbsurf.m_knot[1][sj+nbsurf.m_order[1]-1]);
	return nbsurf_index;
}

ONGEO_SphereTree *ONGEO_NewSphereTree(int num, const ON_NurbsSurface *nbsurfs){
	return new ONGEO_SphereTree(num, nbsurfs);
}
void ONGEO_DeleteSphereTree(ONGEO_SphereTree *this_){
	delete this_;
}
int ONGEO_SphereTree_CreateTree(ONGEO_SphereTree *this_, double radius2, double weight_rate, int level){
	return this_->CreateTree(radius2, weight_rate, level);
}
void ONGEO_SphereTree_RayIntersectTest(const ONGEO_SphereTree *this_, const ON_3dRay &ray, ON_SimpleArray<ONGEO_SphereTree::Result> &results){
	this_->RayIntersectTest(ray, results);
}
int ONGEO_SphereTree_GetRootNodeIndex(const ONGEO_SphereTree *this_){
	return this_->GetRootNodeIndex();
}

int ONGEO_SphereTree_GetNurbsIndexFromBezIndex(const ONGEO_SphereTree *this_, int bez_index){
	return this_->GetNurbsIndexFromBezIndex(bez_index);
}

int ONGEO_SphereTree_GetFirstBezIndexFromBezIndex(const ONGEO_SphereTree *this_, int bez_index){
	return this_->GetFirstBezIndexFromBezIndex(bez_index);
}

int ONGEO_SphereTree_GetNurbsIntervalFromBezIndex(const ONGEO_SphereTree *this_, const ON_NurbsSurface *nbsurfs, int bez_index, ON_Interval range[2]){
	if (!this_) return -1;
	return this_->GetNurbsIntervalFromBezIndex(nbsurfs, bez_index, range);
}
