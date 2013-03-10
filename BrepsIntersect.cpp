#include "ONGEO.h"
#include "Profile.h"

ONGEO_BrepsRayIntersect::ONGEO_BrepsRayIntersect(const ON_Brep **breps, int num_breps){
	for (int j = 0; j < num_breps; ++j){
		const ON_Brep *brep = breps[j];
		const ON_BrepFaceArray &fs = brep->m_F;
		for (int i = 0; i < fs.Count(); ++i){
			const ON_BrepFace &f = fs[i];
			ON_NurbsSurface &nbsurf = nbsurfs.AppendNew();
			f.SurfaceOf()->NurbsSurface(&nbsurf);
			nbs2brep.Append(j);

			LoopsInAFace &lif = loops_faces.AppendNew();
			ONGEO_GetBezierLoops(f, lif.loop_crvs, lif.num_crvs_in_a_loop);

			if (lif.loop_crvs.Count() == 0 || lif.num_crvs_in_a_loop.Count() == 0) continue;
			int numcrvs_outerloop = lif.num_crvs_in_a_loop[0];
			if (!numcrvs_outerloop) continue;
			ON_BoundingBox bb;
			ONGEO_CalculateTightBoundingBox(lif.loop_crvs.First(), numcrvs_outerloop, ON_ZERO_TOLERANCE, bb);
			nbsurf.Trim(0, ON_Interval(bb.m_min.x, bb.m_max.x));
			nbsurf.Trim(1, ON_Interval(bb.m_min.y, bb.m_max.y));
		}
	}
	st = ONGEO_NewSphereTree(nbsurfs.Count(), &nbsurfs[0]);
	st->CreateTree(0, 1.02, -1);
}

ONGEO_BrepsRayIntersect::~ONGEO_BrepsRayIntersect(){
	ONGEO_DeleteSphereTree(st);
}

ONGEO_BrepsRayIntersect::ONGEO_BrepsRayIntersect(const ONGEO_BrepsRayIntersect &rhs){
	st = ONGEO_NewSphereTree(0, 0);
	(*this) = rhs;
}

ONGEO_BrepsRayIntersect &ONGEO_BrepsRayIntersect::operator =(const ONGEO_BrepsRayIntersect &rhs){
	PROF("construct ONGEO_BrepsRayIntersect");
	nbsurfs = rhs.nbsurfs;
	loops_faces = rhs.loops_faces;
	nbs2brep = rhs.nbs2brep;
	*st = *rhs.st;
	return *this;
}

void ONGEO_BrepsRayIntersect::Run(const ON_3dRay &ray, ON_SimpleArray<Result> &results, bool (*Test)(void *data, TestStage stage, int nbsurf_index, const Result *candidate), void *data) const{
	PROF("ONGEO_BrepsRayIntersect::Run");
	ON_SimpleArray<ONGEO_SphereTree::Result> stres;
	st->RayIntersectTest(ray, stres);

	for (int j = 0; j < stres.Count(); ++j){
		ONGEO_SphereTree::Result &stre = stres[j];
		ON_BezierSurface bez = st->bezs[stre.bez_index];
		bez.Trim(0, stre.uint);
		bez.Trim(1, stre.vint);

		ON_Interval range[2];
		int nbsurf_index = st->GetNurbsIntervalFromBezIndex(nbsurfs, stre.bez_index, range);
		range[0] = range[0].ParameterAt(stre.uint);
		range[1] = range[1].ParameterAt(stre.vint);

		if (Test && !Test(data, BeforeIntersectRayBezier, nbsurf_index, 0)) continue;

		ON_3dPointArray tuvs, ptsrf, ptlin;
		ONGEO_IntersectRayBezier_QuasiInterpolating(ray, bez, tuvs, ptsrf, ptlin, 1e-4);
		for (int i = 0; i < tuvs.Count(); ++i){
			const LoopsInAFace &lif = loops_faces[nbsurf_index];
			bool inside = true;
			Result &res = results.AppendNew();
			res.uv.Set(range[0].ParameterAt(tuvs[i][1]), range[1].ParameterAt(tuvs[i][2]));
			res.codlin = ptlin[i];
			res.codsrf = ptsrf[i];
			res.t = tuvs[i][0];
			res.nbs_index = nbsurf_index;

			if (Test && !Test(data, BeforeUVPointIsInside, nbsurf_index, &res)) continue;
			if (!ONGEO_UVPointIsInside(lif.loop_crvs, lif.num_crvs_in_a_loop, res.uv, 1e-4)){
				results.Remove();
				continue;
			}
			if (Test && !Test(data, AfterUVPointIsInside, nbsurf_index, &res)){
				results.Remove();
				continue;
			}
		}
	}
}

ONGEO_BrepsRayIntersect *ONGEO_New_BrepsRayIntersect(const ON_Brep **breps, int num_breps){
	return new ONGEO_BrepsRayIntersect(breps, num_breps);
}

void ONGEO_BrepsRayIntersect_Run(const ONGEO_BrepsRayIntersect *bri, const ON_3dRay &ray, ON_SimpleArray<ONGEO_BrepsRayIntersect::Result> &result, bool (*Test)(void *data, ONGEO_BrepsRayIntersect::TestStage stage, int nbsurf_index, const ONGEO_BrepsRayIntersect::Result *candidate), void *data){
	bri->Run(ray, result, Test, data);
}
