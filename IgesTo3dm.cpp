// Copyright mocchi 2013. mocchi_2003@yahoo.co.jp
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#include "ONGEO.h"

#include <cstring>
#include <cmath>
#include <map>

namespace{
	bool Token(const ONGEO_IgesModel &igs, const ON_String &src, int &idx, ON_String &token){
		const char *s = src.Array();
		if (!s) return false;
		const char *sc = s + idx;
		const char *token_end = std::strchr(sc, igs.gs.parm_delimiter[0]);
		if (!token_end) token_end = std::strchr(sc, igs.gs.record_delimiter[0]);
		if (!token_end) token_end = s + src.Length();
		token = src.Mid(idx, static_cast<int>(token_end - sc));
		int new_idx = static_cast<int>(token_end - s);
		if (new_idx < src.Length() - 1) ++new_idx;
		if (new_idx == idx) return false;
		idx = new_idx;
		return true;
	}
	int DEPtr2Index(int deptr){
		return (deptr > 0) ? (deptr - 1) / 2 : -1;
	}

	class DEIndex;
	ON_SimpleArray<DEIndex *> deindices;
	class DEIndex : public ON_UserData{
		ON_OBJECT_DECLARE( DEIndex );
	public:
		int deidx;
		DEIndex() : ON_UserData(){
			m_userdata_uuid = DEIndex::m_DEIndex_class_id.Uuid();
			m_application_uuid = ON_opennurbs5_id;
			m_userdata_copycount = 1;
			deindices.Append(this);
		}
		~DEIndex(){ }
		DEIndex(const DEIndex &rhs);
		DEIndex &operator =(const DEIndex &rhs){
			ON_UserData::operator =(rhs);
			deidx = rhs.deidx;
			return *this;
		}
	};
	ON_OBJECT_IMPLEMENT( DEIndex, ON_UserData, "66AF7BEC-8D83-497d-A293-EB7968049EAA" );

	template <typename T> T *SafeDup(ON_SimpleArray<ON_Object *> &objs, int idx){
		if (idx < 0 || idx >= objs.Count()) return 0;
		T *obj = T::Cast(objs[idx]);
		if (!obj) return 0;
		T *objo = obj->Duplicate();
		return objo;
	}
	template <typename T> T *SafeDup(ON_Object * obj_){
		T *obj = T::Cast(obj_);
		if (!obj) return 0;
		T *objo = obj->Duplicate();
		return objo;
	}
}

// ===== ONGEO_IgesTo3dmInfo =====
struct ONGEO_IgesTo3dmInfo::Impl{
	typedef std::pair<const ON_Object *, int> ObjIdxPair;

	// ** Commitより前の処理 **
	void AddTempDEIndexToObject(ON_Object *obj, int idx){
		DEIndex *index = new DEIndex();
		index->deidx = idx;
		obj->AttachUserData(index);
	}
	// AddTempDeIndexToObjectしてからCommitDEIndicesするまでの間、有効
	int GetTempDEIndexFromObject(const ON_Object *obj){
		ON_UserData *ud = obj->FirstUserData();
		DEIndex *index;
		while(ud && !(index = DEIndex::Cast(ud))) ud = ud->Next();
		if (!ud || !index) return -1;
		return index->deidx;
	}

	// CommitDEIndices実施後は、Commitより前の処理はできなくなる。

	// ** Commitより以降の処理 **
	void CommitDEIndices(){
		for (int i = 0; i < deindices.Count(); ++i){
			DEIndex *index = deindices[i];
			ON_Object *owner = index->Owner();
			obj2DEidx.Append(ObjIdxPair(owner, index->deidx));
			owner->DetachUserData(index);
			delete index;
		}
		deindices.Empty();
		obj2DEidx.QuickSort(Compare);

		// NULLポインタは一番後ろにまとまる。これを削除する。
		for (int i = obj2DEidx.Count() - 1; i >= 0 && !obj2DEidx.Last()->first; --i){
			obj2DEidx.Remove(i);
		}
	}

	void DeleteDEIndexToObject(ON_Object *obj){
		ObjIdxPair key(obj, 0);
		int aindex = obj2DEidx.BinarySearch(&key, Compare);
		if (aindex < 0) return;
		obj2DEidx[aindex].first = 0;
	}

	static int Compare(const ObjIdxPair *lhs, const ObjIdxPair *rhs){
		// ソート後、NULLポインタが配列の最後尾にくるように、降順ソートする。
#if 0
		if ( lhs->first == rhs->first ) return 0;
		if ( !lhs->first ) return 1;
		if ( !rhs->first ) return -1;
#endif
		if( lhs->first > rhs->first ) return -1;
		if( lhs->first < rhs->first ) return  1;
		return 0;
	}

private:
	ON_ClassArray<ObjIdxPair> obj2DEidx;
	friend struct ONGEO_IgesTo3dmInfo;
};

ONGEO_IgesTo3dmInfo::ONGEO_IgesTo3dmInfo(){
	pimpl = new Impl;
}
ONGEO_IgesTo3dmInfo::~ONGEO_IgesTo3dmInfo(){
	delete pimpl;
}

int ONGEO_IgesTo3dmInfo::DEIndexFromObject(const ON_Object *obj) const{
	Impl::ObjIdxPair key(obj, 0);
	int aindex = pimpl->obj2DEidx.BinarySearch(&key, Impl::Compare);
	if (aindex < 0) return -1;
	return pimpl->obj2DEidx[aindex].second;
}

ONGEO_IgesTo3dmInfo *ONGEO_NewIgesTo3dmInfo(){
	return new ONGEO_IgesTo3dmInfo();
}
void ONGEO_DeleteIgesTo3dmInfo(ONGEO_IgesTo3dmInfo *oiinfo){
	delete oiinfo;
}
int ONGEO_IgesTo3dmInfo_DEIndexFromObject(const ONGEO_IgesTo3dmInfo *oiinfo, const ON_Object *obj){
	if (!oiinfo) return -1;
	return oiinfo->DEIndexFromObject(obj);
}

// IgesTo3dm で3dmに変換できるIges Entity
// Type:100 Circular Arc                    済
// Type:102 Composite Curve                 済
// Type:104 Conic Arc                       未
// Type:106 Copious Data                    済 (form 1,2,3,11,12,13,63)
// Type:108 Plane                           済 (form 0,1)
// Type:110 Line                            済 (form 0,1,2 ただし1,2は有限直線として解釈する)
// Type:112 Parametric Spline Curve         未
// Type:116 Point                           済 (Subfigure Definition Entity は無視)
// Type:120 Surface of Revolution           未
// Type:122 Tabulated Cylinder              済
// Type:124 Transformation Matrix           済
// Type:126 Rational B-Spline Curve         済
// Type:128 Rational B-Spline Surface       済
// Type:141 Boundary
// Type:142 Curve on a Parametric Surface   済
// Type:143 Bounded Surface
// Type:144 Trimmed Surface                 済
// Type:158 Sphere                          済
// Type:186 Manifold Solid B-Rep Object
// Type:314 Color Definition                済
// Type:402 Associativity Instance          済 (form 1,7,14,15)
// Type:502 Vertex
// Type:504 Edge
// Type:508 Loop
// Type:510 Face
// Type:514 Shell

bool ONGEO_IgesTo3dm(const ONGEO_IgesModel &igs, ONX_Model &onx, ONGEO_IgesTo3dmInfo &info){
	ON_SimpleArray<ON_Object *> objs(igs.des.Count());
	objs.SetCount(objs.Capacity());

	std::map<int, ON_Xform> xforms;

	ON_ObjectArray<ON_Group> &groups = onx.m_group_table;
	std::map<int, ON_SimpleArray<int> > oidx2gidx; // first:グループに属するオブジェクトIndex、 second:グループIndex列

	ON_ObjectArray<ON_Layer> &layers = onx.m_layer_table;
	std::map<int, int> layerptr2idx; // first: 負;レイヤーをあらわすDEPtr 正;レイヤー番号、 second:レイヤーインデックス+1

	ON_SimpleArray<int> converted(objs.Count());
	for (int i = 0; i < converted.Capacity(); ++i) converted.Append(0);

	ON::unit_system unit_igs2onx[] = {
		ON::no_unit_system, ON::inches, ON::millimeters, ON::no_unit_system, ON::feet, ON::miles,
		ON::meters, ON::kilometers, ON::mils, ON::microns, ON::centimeters, ON::microinches
	};

	ON_Color col[9] = {
		ON_Color(0,0,0), ON_Color(0,0,0), ON_Color(255,0,0), ON_Color(0,255,0), ON_Color(0,0,255),
		ON_Color(255,255,0), ON_Color(255,0,255), ON_Color(0,255,255), ON_Color(255,255,255)
	};

	onx.m_settings.m_ModelUnitsAndTolerances.m_unit_system.m_unit_system = 
		unit_igs2onx[(igs.gs.unit_flag < 0 || igs.gs.unit_flag >= 11) ? 0 : igs.gs.unit_flag];

	// Todo:transformの処理
	for (int i = 0; i < igs.des.Count(); ++i){
		ON_SimpleArray<int> stack;
		stack.Append(i);
		while(stack.Count()){
			int ii = *stack.Last();
			stack.SetCount(stack.Count()-1);
			if (converted[ii]) continue;

			const ONGEO_IgesModel::DirectoryEntrySection &de = igs.des[ii];

			ON_Xform *xform = 0;
			if (de.trans_matrix > 0){
				int ti = (de.trans_matrix - 1) / 2;
				if (!converted[ti]){
					stack.Append(ii); stack.Append(ti);
					continue;
				}
				xform = &xforms[ti];
			}

			ON_3dmObjectAttributes att;
			if (de.stnum.hierarchy != 1){
				att.SetColorSource(de.color_num ? ON::color_from_object : ON::color_from_layer);
				if (de.color_num > 0 && de.color_num < 9) att.m_color = col[de.color_num];
				else if (de.color_num < 0){
					int ci = ((-de.color_num) - 1) / 2;
					if(!converted[ci]){
						stack.Append(ii); stack.Append(ci);
						continue;
					}else{
						ON_3dmObjectAttributes *attc = ON_3dmObjectAttributes::Cast(objs[ci]);
						if (attc) att.m_color = attc->m_color;
					}
				}
			}else att.SetColorSource(ON::color_from_layer);

			int index = 0;
			ON_Object *obj = 0;
			ON_String token;
			const ON_String &prm = igs.ps[ii];
			Token(igs, prm, index, token);

			if (de.entity_type == 100){
				static ON_ArcCurve prt_arc;
				ON_ArcCurve *arc = prt_arc.Duplicate();
				ON_Plane pln(ON_3dPoint(0,0,0), ON_3dPoint(1,0,0), ON_3dPoint(0,1,0));
				ON_2dPoint pt[2];
				double dx[2] = {0}, dy[2] = {0};
				if (Token(igs, prm, index, token)) pln.origin.z = std::strtod(token, 0);
				if (Token(igs, prm, index, token)) pln.origin.x = std::strtod(token, 0);
				if (Token(igs, prm, index, token)) pln.origin.y = std::strtod(token, 0);
				for (int h = 0; h < 2; ++h){
					if (Token(igs, prm, index, token)) pt[h].x = std::strtod(token, 0);
					if (Token(igs, prm, index, token)) pt[h].y = std::strtod(token, 0);
					dx[h] = pt[h].x - pln.origin.x;
					dy[h] = pt[h].y - pln.origin.y;
				}
				double radius = std::sqrt(dx[0] * dx[0] + dy[0] * dy[0]);
				ON_Interval angles(std::atan2(dy[0], dx[0]), std::atan2(dy[1], dx[1]));
				if (angles.IsDecreasing()) angles.m_t[0] -= ON_PI * 2.0;
				arc->m_arc.Create(ON_Circle(pln, radius), angles);
				if (xform) arc->Transform(*xform);
				obj = arc;
			}else if (de.entity_type == 102){
				int num = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) : 0;
				bool defined_forward = false;
				ON_SimpleArray<int> ccrv_i(num); ccrv_i.SetCount(num);
				for (int h = 0; h < num; ++h){
					if (Token(igs, prm, index, token)) ccrv_i[h] = DEPtr2Index(std::strtol(token, 0, 0));
					if (!converted[ccrv_i[h]]){
						if (!defined_forward) stack.Append(ii), defined_forward = true;
						stack.Append(ccrv_i[h]);
					}
				}
				if (defined_forward) continue;
				static ON_PolyCurve prt_pcrv;
				ON_PolyCurve *pcrv = prt_pcrv.Duplicate();
				ON_Curve *crv = 0;
				for (int h = 0; h < num; ++h){
					if (crv = SafeDup<ON_Curve>(objs, ccrv_i[h])) pcrv->Append(crv);
				}
				if (xform) pcrv->Transform(*xform);
				obj = pcrv;
			}else if (de.entity_type == 104){
				// Todo:
			}else if (de.entity_type == 106){
				static ON_PointCloud prt_pc;
				static ON_PolylineCurve prt_plc;
				ON_3dPointArray *pa = 0;
				ON_PointCloud *pc = 0;
				ON_PolylineCurve *plc = 0;
				int pts_type;
				if (de.form_num == 1 || de.form_num == 2 || de.form_num == 3 ||
					de.form_num == 11 || de.form_num == 12 || de.form_num == 13 || de.form_num == 63){
					if (de.form_num < 10) pc = prt_pc.Duplicate(), pa = &pc->m_P;
					else plc = prt_plc.Duplicate(), pa = &plc->m_pline;
					int t = Token(igs, prm, index, token) ? std::strtol(token, 0, 0) : 0;
					if (!t && de.form_num != 63) pts_type = de.form_num % 10;
					else pts_type = t;
				}
				if (pa){
					int num = Token(igs, prm, index, token) ? std::strtol(token, 0, 0) : 0;
					pa->SetCapacity(num + (de.form_num == 63) ? 1 : 0);
					if (plc) plc->m_t.SetCapacity(pa->Capacity());
					if (pts_type == 1){
						double z = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
						for (int h = 0; h < num; ++h){
							double x = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							double y = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							pa->Append(ON_3dPoint(x, y, z));
						}
					}else if (pts_type == 2){
						for (int h = 0; h < num; ++h){
							double x = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							double y = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							double z = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							pa->Append(ON_3dPoint(x, y, z));
						}
					}else if (pts_type == 3){
						for (int h = 0; h < num; ++h){
							double x = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							double y = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							double z = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
							pa->Append(ON_3dPoint(x, y, z));
							for (int hh = 0; hh < 3; ++hh) Token(igs, prm, index, token);
						}
					}
					if (plc){
						if (de.form_num == 63) pa->Append(*pa[0]);
						for (int h = 0; h < pa->Count(); ++h) plc->m_t.Append(h);
						if (xform) plc->Transform(*xform);
						obj = plc;
					}else if (pc){
						if (xform) pc->Transform(*xform);
						obj = pc;
					}
				}
			}else if (de.entity_type == 108){
				double equation[4] = {0};
				for (int h = 0; h < 4; ++h){
					if (Token(igs, prm, index, token)) equation[h] = std::strtod(token, 0);
				}
				int closed_crv = (Token(igs, prm, index, token)) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1;
				if (closed_crv >= 0 && !converted[closed_crv]){
					stack.Append(ii); stack.Append(closed_crv);
					continue;
				}
				static ON_PlaneSurface prt_ps;
				if (de.form_num == 0 || closed_crv < 0){
					ON_PlaneSurface *ps = prt_ps.Duplicate();
					ps->m_plane = ON_Plane(equation);
					for (int h = 0; h < 2; ++h){
						ps->SetExtents(h, ON_Interval(-198000,198000));
						ps->SetDomain(h, -198000, 198000);
					}
					if (xform) ps->Transform(*xform);
					obj = ps;
				}else if (de.form_num == 1 && closed_crv >= 0){
					ON_SimpleArray<ON_Curve *> boundary(1);
					boundary.Append(SafeDup<ON_Curve>(objs, closed_crv));
					if (boundary[0]){
						boundary[0] = boundary[0]->Duplicate();
						ON_Brep *brep = ON_BrepTrimmedPlane(ON_Plane(equation), boundary, true);
						if (brep){
							if (xform) brep->Transform(*xform);
							obj = brep;
						}
					}
				}
			}else if (de.entity_type == 110){
				static ON_LineCurve prt_lin;
				if (de.form_num == 0 || de.form_num == 1 || de.form_num == 2){
					ON_LineCurve *lin = prt_lin.Duplicate();
					ON_3dPoint pt[2];
					for (int h = 0; h < 2; ++h){
						if (Token(igs, prm, index, token)) pt[h].x = std::strtod(token, 0);
						if (Token(igs, prm, index, token)) pt[h].y = std::strtod(token, 0);
						if (Token(igs, prm, index, token)) pt[h].z = std::strtod(token, 0);
					}
					lin->SetStartPoint(pt[0]), lin->SetEndPoint(pt[1]);
					if (xform) lin->Transform(*xform);
					obj = lin;
				}
			}else if (de.entity_type == 116){
				static ON_Point prt_pt;
				ON_Point *pt = prt_pt.Duplicate();
				if (Token(igs, prm, index, token)) pt->point.x = std::strtod(token, 0);
				if (Token(igs, prm, index, token)) pt->point.y = std::strtod(token, 0);
				if (Token(igs, prm, index, token)) pt->point.z = std::strtod(token, 0);
				if (xform) pt->Transform(*xform);
				obj = pt;
			}else if (de.entity_type == 122){
				if (Token(igs, prm, index, token)){
					int dircrv = DEPtr2Index(std::strtol(token, 0, 0));
					if (dircrv >= 0){
						if (!converted[dircrv]){
							stack.Append(ii), stack.Append(dircrv);
							continue;
						}
						ON_3dVector v;
						ON_Curve *crv = ON_Curve::Cast(objs[dircrv]);
						if (crv){
							for (int h = 0; h < 3; ++h) if (Token(igs, prm, index, token)) v[h] = std::strtod(token, 0);
							ON_NurbsCurve nc1, nc2;
							crv->NurbsCurve(&nc1);
							nc2 = nc1;
							nc2.Translate(v-nc1.PointAtStart());
							ON_NurbsSurface *srf = ON_NurbsSurface::New();
							if (srf) srf->CreateRuledSurface(nc1, nc2);
							obj = srf;
						}
					}
				}
			}else if (de.entity_type == 124){
				ON_Xform &xform = xforms[ii];
				for (int ri = 0; ri < 3; ++ri){
					for (int ci = 0; ci < 4; ++ci){
						if (Token(igs, prm, index, token)) xform.m_xform[ri][ci] = std::strtod(token, 0);
					}
				}
				xform.m_xform[3][3] = 1.0;
				obj = 0;
			}else if (de.entity_type == 126){
				int num_ctrlpt = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) + 1 : 0;
				int degree     = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) : 0; // 次数
				Token(igs, prm, index, token);
				Token(igs, prm, index, token);
				bool rational  = (Token(igs, prm, index, token)) ? ((std::strtol(token, 0, 0)) ? false : true) : false;
				Token(igs, prm, index, token);
				ON_NurbsCurve *nc = ON_NurbsCurve::New(de.stnum.ent_use_flg == 5 ? 2 : 3, rational, degree + 1, num_ctrlpt);

				// ノットベクトル
				Token(igs, prm, index, token);
				for (int h = 0; h < nc->KnotCount(); ++h){
					if (Token(igs, prm, index, token)) nc->m_knot[h] = std::strtod(token, 0);
					else if (h > 0) nc->m_knot[h] = nc->m_knot[h-1];
					else nc->m_knot[h] = 0;
				}
				Token(igs, prm, index, token);

				// Weight
				for (int h = 0; h < num_ctrlpt; ++h){
					double w = (Token(igs, prm, index, token)) ? std::strtod(token, 0) : 1.0;
					if (rational) nc->CV(h)[nc->m_dim] = w;
				}

				// XYZ
				for (int h = 0; h < num_ctrlpt; ++h){
					double *cv = nc->CV(h);
					double w = rational ? cv[nc->m_dim] : 1.0;
					for (int hh = 0; hh < 3; ++hh){
						double v = (Token(igs, prm, index, token)) ? std::strtod(token, 0) : 0.0;
						if (hh < nc->m_dim) cv[hh] = v * w;
					}
				}

				// Limits
				ON_Interval lim;
				if (Token(igs, prm, index, token)) lim.m_t[0] = std::strtod(token, 0);
				if (Token(igs, prm, index, token)) lim.m_t[1] = std::strtod(token, 0);
				nc->SetDomain(lim.m_t[0], lim.m_t[1]);

				if (xform) nc->Transform(*xform);
				obj = nc;
			}else if (de.entity_type == 128){
				int num_ctrlpt1 = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) + 1 : 0;
				int num_ctrlpt2 = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) + 1 : 0;
				int degree1     = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) : 0; // 次数
				int degree2     = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) : 0; // 次数
				Token(igs, prm, index, token);
				Token(igs, prm, index, token);
				bool rational  = (Token(igs, prm, index, token)) ? ((std::strtol(token, 0, 0)) ? false : true) : false;
				Token(igs, prm, index, token);
				Token(igs, prm, index, token);

				ON_NurbsSurface *nf = ON_NurbsSurface::New(de.stnum.ent_use_flg == 5 ? 2 : 3, rational, degree1 + 1, degree2 + 1, num_ctrlpt1, num_ctrlpt2);

				// ノットベクトル
				for (int ki = 0; ki < 2; ++ki){
					Token(igs, prm, index, token);
					for (int h = 0; h < nf->KnotCount(ki); ++h){
						if (Token(igs, prm, index, token)) nf->m_knot[ki][h] = std::strtod(token, 0);
						else if (h > 0) nf->m_knot[ki][h] = nf->m_knot[ki][h-1];
						else nf->m_knot[ki][h] = 0;
					}
					Token(igs, prm, index, token);
				}
				
				// Weight
				for (int hv = 0; hv < num_ctrlpt2; ++hv){
					for (int hu = 0; hu < num_ctrlpt1; ++hu){
						double w = (Token(igs, prm, index, token)) ? std::strtod(token, 0) : 1.0;
						if (rational) nf->CV(hu, hv)[nf->m_dim] = w;
					}
				}

				// XYZ
				for (int hv = 0; hv < num_ctrlpt2; ++hv){
					for (int hu = 0; hu < num_ctrlpt1; ++hu){
						double *cv = nf->CV(hu, hv);
						double w = rational ? cv[nf->m_dim] : 1.0;
						for (int hh = 0; hh < 3; ++hh){
							double v = (Token(igs, prm, index, token)) ? std::strtod(token, 0) : 0.0;
							if (hh < nf->m_dim) cv[hh] = v * w;
						}
					}
				}

				// Limits
				for (int h = 0; h < 2; ++h){
					ON_Interval lim;
					if (Token(igs, prm, index, token)) lim.m_t[0] = std::strtod(token, 0);
					if (Token(igs, prm, index, token)) lim.m_t[1] = std::strtod(token, 0);
					nf->SetDomain(h, lim.m_t[0], lim.m_t[1]);
				}

				if (xform) nf->Transform(*xform);
				obj = nf;
			}else if (de.entity_type == 142){
				ON_SimpleArray<int> fwddef;

				Token(igs, prm, index, token);
				Token(igs, prm, index, token);
				int bi = (Token(igs, prm, index, token)) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1;
				if (bi >= 0 && !converted[bi]) fwddef.Append(bi);
				int ci = (Token(igs, prm, index, token)) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1;
				if (ci >= 0 && !converted[ci]) fwddef.Append(ci);
				if (fwddef.Count()){
					stack.Append(ii), stack.Append(fwddef.Count(), fwddef);
					continue;
				}

				static ON_CurveOnSurface cs_prt;
				ON_CurveOnSurface *cs = cs_prt.Duplicate();
				cs->m_c2 = SafeDup<ON_Curve>(objs, bi);
				cs->m_c3 = SafeDup<ON_Curve>(objs, ci);
				if (xform) if (cs->m_c3) cs->m_c3->Transform(*xform);
				obj = cs;
			}else if (de.entity_type == 144){
				ON_SimpleArray<int> fwddef;

				int si = Token(igs, prm, index, token) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1;
				if (si >= 0 && !converted[si]) fwddef.Append(si);

				int num_outer = Token(igs, prm, index, token) ? std::strtol(token, 0, 0) : 0;
				int num_inner = Token(igs, prm, index, token) ? std::strtol(token, 0, 0) : 0;
				int boi = Token(igs, prm, index, token) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1;
				if (boi >= 0 && !converted[boi]) fwddef.Append(boi);

				int num_loops = num_inner + ((boi >= 0) ? 1 : 0);
				ON_SimpleArray<int> bis(num_loops);
				bis.Append(boi);
				for (int h = 0; h < num_inner; ++h){
					bis.Append(Token(igs, prm, index, token) ? DEPtr2Index(std::strtol(token, 0, 0)) : -1);
					if (*bis.Last() >= 0 && !converted[*bis.Last()]) fwddef.Append(*bis.Last());
				}

				if (fwddef.Count()){
					stack.Append(ii), stack.Append(fwddef.Count(), fwddef);
					continue;
				}

				if (si >= 0){
					ON_Surface *srf = SafeDup<ON_Surface>(objs[si]);
					if (srf){
						ON_Brep *brep = ON_Brep::New();
						ON_BrepFace &f = (boi >= 0) ? (brep->m_S.Append(srf), brep->NewFace(0)) : (brep->Create(srf), brep->m_F[0]);
						for (int h = 0; h < bis.Count(); ++h){
							ON_CurveOnSurface *cs = (bis[h] >= 0) ? ON_CurveOnSurface::Cast(objs[bis[h]]) : 0;
							if (!cs) continue;

							ON_BrepLoop &l = brep->NewLoop(h ? l.inner : l.outer, f);
							if (cs->m_c3){
								ON_BrepVertex &v1 = brep->NewVertex(cs->m_c3->PointAtStart());
								ON_BrepVertex &v2 = cs->m_c3->IsClosed() ? v1 : brep->NewVertex(cs->m_c3->PointAtEnd());
								v1.m_tolerance = v2.m_tolerance = 0;
								ON_BrepEdge &e = brep->NewEdge(v1, v2, brep->AddEdgeCurve(cs->m_c3->Duplicate()));
								e.m_tolerance = 0;
								if (cs->m_c2){
									ON_BrepTrim &t = brep->NewTrim(e, false, l, brep->AddTrimCurve(cs->m_c2->Duplicate()));
									t.m_tolerance[0] = t.m_tolerance[1] = 0;
								}
							}
						}
						// ループの向きの調整
						{
#if 0
							int count_loop_crvs, capacity_loop_crvs;
							ON_ClassArray<ON_BezierCurve> loop_crvs;
							ON_SimpleArray<int> num_crvs_in_a_loop;
							ONGEO_GetBezierLoops(f, loop_crvs, num_crvs_in_a_loop);
							count_loop_crvs = loop_crvs.Count();
							capacity_loop_crvs = loop_crvs.Capacity();

							ON_BezierCurve *lckeep = loop_crvs.KeepArray();
							ON_ClassArray<ON_BezierCurve> loop_crvs_s;
							ON_SimpleArray<int> num_crvs_in_a_loop_s(1);
							num_crvs_in_a_loop_s.AppendNew();
							int sum = 0;
							ON_SimpleArray<int> inout(count_loop_crvs);
							ON_2dPoint uv = ON_2dPoint(brep->m_L[0].m_pbox.Min()) + ON_2dVector(-1, -1);
							for (int h = 0; h < num_crvs_in_a_loop.Count(); ++h){
								num_crvs_in_a_loop_s[0] = num_crvs_in_a_loop[h];
								loop_crvs_s.SetArray(lckeep+sum, num_crvs_in_a_loop_s[0], num_crvs_in_a_loop_s[0]);

								inout.Append(ONGEO_UVPointIsInside(loop_crvs_s, num_crvs_in_a_loop_s, uv, ON_ZERO_TOLERANCE) ? 1 : 0);

								loop_crvs_s.KeepArray();
								sum += num_crvs_in_a_loop_s[0];
							}
							loop_crvs.SetArray(lckeep, count_loop_crvs, capacity_loop_crvs);
							for (int hl = 0; hl < brep->m_L.Count(); ++hl){
								ON_BrepLoop &l = brep->m_L[hl];
								if ((l.m_type == l.outer && inout[hl] == 1) || (l.m_type == l.inner && inout[hl] == 0)){
									brep->FlipLoop(l);
								}
							}
#else
							for (int hl = 0; hl < brep->m_L.Count(); ++hl){
								ON_BrepLoop &l = brep->m_L[hl];
								ON_BrepLoop::TYPE type = brep->ComputeLoopType(l);
								if ((l.m_type == l.outer && type == l.inner) || (l.m_type == l.inner && type == l.outer)){
									brep->FlipLoop(l);
								}
							}
#endif
						}
						obj = brep;
					}
				}
			}else if (de.entity_type == 158){
				double radius = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
				double x = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
				double y = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
				double z = Token(igs, prm, index, token) ? std::strtod(token, 0) : 0;
				ON_Brep *brep = ON_BrepSphere(ON_Sphere(ON_3dPoint(x,y,z), radius));
				if (brep && xform) brep->Transform(*xform);
				obj = brep;
			}else if (de.entity_type == 314){
				ON_3dmObjectAttributes *att_ = att.Duplicate();
				double rgb[3] = {0};
				if (Token(igs, prm, index, token)) rgb[0] = std::strtod(token, 0) / 100.0;
				if (Token(igs, prm, index, token)) rgb[1] = std::strtod(token, 0) / 100.0;
				if (Token(igs, prm, index, token)) rgb[2] = std::strtod(token, 0) / 100.0;
				att_->m_color.SetFractionalRGB(rgb[0], rgb[1], rgb[2]);
				obj = att_;
			}else if (de.entity_type == 402){
				if (de.form_num == 1 || de.form_num == 7 || de.form_num == 14 || de.form_num == 15){
					ON_Group &g = groups.AppendNew();
					g.SetGroupName(de.ent_label);
					g.SetGroupIndex(groups.Count() - 1);
					int num = (Token(igs, prm, index, token)) ? std::strtol(token, 0, 0) : 0;
					for (int h = 0; h < num; ++h){
						if (Token(igs, prm, index, token))
							oidx2gidx[DEPtr2Index(std::strtol(token, 0, 0))].Append(g.GroupIndex());
					}
				}
			}else{
			}

			if (obj && (de.stnum.subord_ent_sw & 1) == 0){
				int &layer_idx = layerptr2idx[de.level];
				if (layer_idx == 0){
					ON_Layer &layer = layers.AppendNew();
					layer.SetLayerIndex(layers.Count()-1);
					ON_String layername;
					layername.Format("IGES level %d", de.level);
					layer.SetLayerName(layername);
					layer_idx = layers.Count();
				}

				ONX_Model_Object &onxobj = onx.m_object_table.AppendNew();
				onxobj.m_object = obj;
				onxobj.m_bDeleteObject = true;
				onxobj.m_attributes = att;
				onxobj.m_attributes.m_layer_index = layer_idx;
				const char *el = de.ent_label;
				while(*el == ' ' && *el != '\0') ++el;
				onxobj.m_attributes.m_name = ON_wString(el);
			}
			info.pimpl->AddTempDEIndexToObject(obj, ii);
			objs[ii] = obj;
			converted[ii] = 1;
		}
	}

	for (int i = 0; i < onx.m_object_table.Count(); ++i){
		ONX_Model_Object &onxobj = onx.m_object_table[i];
		int idx = info.pimpl->GetTempDEIndexFromObject(onxobj.m_object);
		if (idx < 0) continue;
		std::map<int, ON_SimpleArray<int> >::iterator iter = oidx2gidx.find(idx);
		if (iter == oidx2gidx.end()) continue;
		for (int h = 0; h < iter->second.Count(); ++h){
			onxobj.m_attributes.AddToGroup(iter->second[h]);
		}
	}
	info.pimpl->CommitDEIndices();

	ONX_Model objdelete;
	for (int i = 0; i < objs.Count(); ++i){
		if (!objs[i] || (igs.des[i].stnum.subord_ent_sw & 1) == 0) continue;
		info.pimpl->DeleteDEIndexToObject(objs[i]);
		ONX_Model_Object &obj = objdelete.m_object_table.AppendNew();
		obj.m_object = objs[i];
		obj.m_bDeleteObject = true;
	}
	info.pimpl->CommitDEIndices();
	return true;
}

