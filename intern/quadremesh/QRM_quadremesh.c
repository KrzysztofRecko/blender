/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s): Alexander Pinzon Fernandez
 *                 Krzysztof Recko
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file intern/quadremesh/MOD_quadremesh.c
 *  \ingroup quadremesh
 */

#include "MEM_guardedalloc.h"

#include "DNA_modifier_types.h"
#include "DNA_object_types.h"

#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_linklist.h"
#include "BLI_memarena.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_mesh_mapping.h"

#include "PIL_time.h"

#include "MOD_util.h"
#include "quadremesh_system.h"

#define QR_MAKEPOLYS
#define QR_SHOWORIENTATIONS

#if 0 /* UNUSED ROUTINES */
/**
* Random point, P, uniformly from within triangle ABC, method given by
* Robert Osada, Thomas Funkhouser, Bernard Chazelle, and David Dobkin. 2002. Shape distributions. ACM Trans. Graph. 21,
* 4 (October 2002), 807-832. DOI=10.1145/571647.571648 http://doi.acm.org/10.1145/571647.571648
* a,b,c are the triangle points
* r1, r2, are the randon numbers betwen [0, 1]
* P = (1 − sqrt(r1)) A + sqrt(r1)(1 − r2) B + sqrt(r1) * r2 C
*/
static void uniformRandomPointWithinTriangle(float r[3], float a[3], float b[3], float c[3])
{
	float va, vb, vc;
	float pa[3], pb[3], pc[3];
	float r1, r2;
	r1 = BLI_frand();
	r2 = BLI_frand();
	va = 1.0f - sqrtf(r1);
	vb = sqrtf(r1) * (1.0f - r2);
	vc = sqrtf(r1) * r2;
	mul_v3_v3fl(pa, a, va);
	mul_v3_v3fl(pb, b, vb);
	mul_v3_v3fl(pc, c, vc);
}

static void uniformRandomPointWithinFace(float r[3], LaplacianSystem *sys, int indexf){
	int *vin;
	vin = sys->faces[indexf];
	uniformRandomPointWithinTriangle(r, sys->co[vin[0]], sys->co[vin[1]], sys->co[vin[2]]);
}

/*
* Compute the normal curvature
* k = dot(2*no, (pi - pj)) / (|pi - pj|)^2
* no = normal on vertex pi
* pi - pj is a vector direction on this case the gradient field direction
* the gradient field direction on some vertex is computed how the average of the faces around vertex
*/
static void computeSampleDistanceFunctions(LaplacianSystem *sys, float user_h, float user_alpha) {
	int i, j, *fin, lin;
	float avg1[3], avg2[3], no[3], k1, k2, h1, h2;

	sys->h1 = MEM_mallocN(sizeof(float) * im->num_verts, "QuadRemeshH1");
	sys->h2 = MEM_mallocN(sizeof(float) * im->num_verts, "QuadRemeshH2");

	for (i = 0; i < sys->total_verts; i++) {
		zero_v3(avg1);
		zero_v3(avg2);

		fin = sys->ringf_map[i].indices;
		lin = sys->ringf_map[i].count;

		for (j = 0; j < lin; j++) {
			add_v3_v3(avg1, sys->gf1[j]);
			add_v3_v3(avg2, sys->gf2[j]);
		}

		mul_v3_fl(avg1, 1.0f / ((float)lin));
		mul_v3_fl(avg2, 1.0f / ((float)lin));

		copy_v3_v3(no, sys->no[i]);
		mul_v3_fl(no, 2.0f);

		k1 = dot_v3v3(no, avg1) / dot_v3v3(avg1, avg1);
		k2 = dot_v3v3(no, avg2) / dot_v3v3(avg2, avg2);

		h1 = user_h / (1.0f + user_alpha * (logf(1.0f + k1*k1)));
		h2 = user_h / (1.0f + user_alpha * (logf(1.0f + k2*k2)));

		sys->h1[i] = h1;
		sys->h2[i] = h2;
	}
}
bool nextPointFromVertex(float r_co[3], int *r_face, int *r_edge, GradientFlowSystem *gfsys, int in_v)
{
	int i, f, vf1, vf2;
	float gf[3], a[3], b[3], alfa, beta, gamma, dummy[3];
	QuadRemeshSystem *sys = gfsys->sys;

	for (i = 0; i < sys->ringf_map[in_v].count; i++) {
		f = sys->ringf_map[in_v].indices[i];
		if (!getSecondAndThirdVert(&vf1, &vf2, sys, f, in_v)) continue;

		sub_v3_v3v3(a, sys->co[vf1], sys->co[in_v]);
		sub_v3_v3v3(b, sys->co[vf2], sys->co[in_v]);
		copy_v3_v3(gf, gfsys->gfield[f]);

		alfa = angle_v3v3(a, b);
		beta = angle_v3v3(a, gf);
		gamma = angle_v3v3(gf, b);

		if (beta + gamma >= alfa - 0.001f && beta + gamma <= alfa + 0.001f) {
			/* vertex on this face */
			add_v3_v3(gf, sys->co[in_v]);
			isect_line_line_v3(sys->co[vf1], sys->co[vf2], sys->co[in_v], gf, r_co, dummy);
			*r_edge = getEdgeFromVerts(sys, vf1, vf2);
			*r_face = f;

			return true;
		}
	}

	return false;
}

/*
* alpha is degree of anisotropic curvature sensitivity
* h is the desired distance
* return ve[0] number of vertices
* return ve[1] number of edges
*/
void estimateNumberGFVerticesEdges(int ve[2], QuadRemeshSystem *sys, float h)
{
	int i, totalv, totale;
	float area = 0.0f;
	float sqrtarea;
	for (i = 0; i < sys->total_faces; i++) {
		area += area_tri_v3(sys->co[sys->faces[i][0]], sys->co[sys->faces[i][1]], sys->co[sys->faces[i][2]]);
	}
	sqrtarea = sqrtf(area);
	if (h > 0.0f) {
		totalv = ((sqrtarea / h) + 1.0f);
		totale = totalv * sqrtarea * 2.0f;
		totalv = totalv * totalv;
	}
	else{
		totalv = sqrtarea + 1.0f;
		totale = totalv * sqrtarea * 2.0f;
		totalv = totalv * totalv;
	}
	ve[0] = totalv;
	ve[1] = totale;
}

static bool isPointOnSegment(float in_v[3], float in_l1[3], float in_l2[3])
{
	float a[3], b[3];
	sub_v3_v3v3(a, in_v, in_l1);
	sub_v3_v3v3(b, in_l2, in_v);
	return dot_v3v3(a, b) >= 0.0f;
}


static int getFaceFromTwoEdges(InputMesh *im, int in_e1, int in_e2)
{
	if (im->faces_edge[in_e1][0] == im->faces_edge[in_e2][0])
		return im->faces_edge[in_e2][0];

	return im->faces_edge[in_e2][1];
}

static int getEdgeOppositeToVertex(InputMesh *im, int in_v, int in_f)
{
	int i;

	for (i = 0; i < 3; i++)
		if (im->faces[in_f][i] == in_v)
			break;

	return getEdgeFromVerts(im, im->faces[in_f][(i + 1) % 3], im->faces[in_f][(i + 2) % 3]);
}

static float getSamplingDistanceFunctionOnFace(GradientFlowSystem *gfsys, int in_f, float in_co[3])
{
	float h1, h2, h3, uv[2];
	InputMesh *im = &gfsys->sys->input_mesh;

	resolve_tri_uv_v3(uv, in_co, im->co[im->faces[in_f][0]], im->co[im->faces[in_f][1]], im->co[im->faces[in_f][2]]);

	h1 = gfsys->hfunction[im->faces[in_f][0]];
	h2 = gfsys->hfunction[im->faces[in_f][1]];
	h3 = gfsys->hfunction[im->faces[in_f][2]];

	return uv[0] * h1 + uv[1] * h2 + (1.0f - uv[0] - uv[1]) * h3;
}
#endif // UNUSED ROUTINES

#ifdef QR_GENERATE
/* MESH GENERATION */

static void makeFeatureEdges(OutputMesh *om, InputMesh *im)
{
	int i, j, f1, f2;
	float angle;
	GFPoint p1, p2;

	for (i = 0; i < im->num_edges; i++) {
		f1 = im->faces_edge[i][0];
		f2 = im->faces_edge[i][1];
		if (f1 == QR_NO_FACE || f2 == QR_NO_FACE)
			angle = M_PI * 2.0f;
		else
			angle = angle_normalized_v3v3(im->no[f1], im->no[f2]);

		if (angle >= M_PI_2 * 0.7f) {
			p1.type = eVert;
			p1.v = im->edges[i][0];
			copy_v3_v3(p1.co, im->co[im->edges[i][0]]);
			addGFPoint(im, om, &p1);

			for (j = 0; j < im->ringe_map[p1.v].count; j++)
				insertOnQREdge(om, &om->ringe[im->ringe_map[p1.v].indices[j]], p1.id);

			p2.type = eVert;
			p2.v = im->edges[i][1];
			copy_v3_v3(p2.co, im->co[im->edges[i][1]]);
			addGFPoint(im, om, &p2);

			for (j = 0; j < im->ringe_map[p2.v].count; j++)
				insertOnQREdge(om, &om->ringe[im->ringe_map[p2.v].indices[j]], p2.id);

			linkOnQREdge(om, GFSYSNONE, &om->ringe[i], p1.id, p2.id);
		}
	}
}

static void deleteDegenerateVerts(OutputMesh *om)
{
	int i, n;
	bool *donevert, doverts = true;
	QRVertID *vertmap;
	QRVert *newverts;
	QRDiskLink *it;

	donevert = MEM_callocN(sizeof(QRVertID) * om->num_verts, __func__);

	while (doverts) {
		doverts = false;
		for (i = 0; i < om->num_verts; i++) {
			if (donevert[i])
				continue;

			donevert[i] = true;

			if (om->verts[i].num_links == 1) {
				doverts = true;
				donevert[om->verts[i].link->v] = false;
				unlinkVerts(om, om->verts[i].link);
			}
			else if (om->verts[i].num_links == 2) {
				doverts = true;
				donevert[om->verts[i].link->v] = false;
				donevert[om->verts[i].link->next->v] = false;
				dissolveVert(om, i);
			}
		}
	}

	MEM_SAFE_FREE(donevert);
	vertmap = MEM_mallocN(sizeof(QRVertID) * om->num_verts, __func__);

	for (i = 0, n = 0; i < om->num_verts; i++) {
		if (om->verts[i].num_links == 0)
			vertmap[i] = -1;
		else
			vertmap[i] = n++;
	}

	newverts = MEM_callocN(sizeof(QRVert) * n, "Newlinks");

	for (i = 0, n = 0; i < om->num_verts; i++) {
		if (vertmap[i] == -1) continue;

		memcpy(&newverts[n], &om->verts[i], sizeof(QRVert));

		it = newverts[n].link;
		do {
			it->v = vertmap[it->v];
			it = it->next;
		} while (it != newverts[n].link);
		n++;
	}

	MEM_SAFE_FREE(om->verts);
	MEM_SAFE_FREE(vertmap);
	om->num_verts = n;
	om->verts = newverts;
}

static void hideEdgesOnFaces(OutputMesh *om, InputMesh *im)
{
	int f;
	QREdge *e;
	QREdgeLink *itel;
	LinkNode *it;

	for (f = 0; f < im->num_faces; f++) {
		for (it = om->ringf[f]; it; it = it->next) {
			e = (QREdge*)it->link;
			for (itel = e->v1; itel; itel = itel->next)
				if (itel->elink) {
					itel->elink->e = 0;
					itel->elink->brother->e = 0;
				}
		}
	}
}

static void makeEdges(OutputMesh *om, MEdge *r_edges)
{
	int i, j, e;
	QRDiskLink *it;

	for (i = 0, e = 0; i < om->num_verts; i++) {
		for (j = 0, it = om->verts[i].link;
			 j < om->verts[i].num_links;
			 j++, it = it->next)
		{
			if (it->e < 0) {
				it->brother->e = it->e = e;

				r_edges[e].v1 = i;
				r_edges[e].v2 = it->v;
				r_edges[e].bweight = 0;
				r_edges[e].crease = 0;
				r_edges[e].flag = ME_EDGEDRAW;
				e++;
			}
		}
	}
}

static void makePolys(OutputMesh *om, MPoly *r_polys, MLoop *r_loops)
{
	int i, j, s, num_loops, num_polys, max_num_polys;
	QRVertID v;
	QRDiskLink *it, *lit;

	num_loops = 0;
	num_polys = 0;
	max_num_polys = 2 + om->num_edges - om->num_verts;

	for (i = 0; i < om->num_verts; i++) {
		for (j = 0, it = om->verts[i].link;
			 j < om->verts[i].num_links;
			 j++, it = it->next)
		{
			if (it->poly_on_right) continue;

			s = num_loops;
			v = i;
			lit = it;
			do {
				r_loops[num_loops].v = v;
				r_loops[num_loops].e = lit->e;
				num_loops++;

				v = lit->v;
				lit->poly_on_right = true;
				lit = lit->brother->next;
			} while (v != i);
			
			r_polys[num_polys].loopstart = s;
			r_polys[num_polys].totloop = num_loops - s;
			r_polys[num_polys].mat_nr = 0;
			r_polys[num_polys].flag = 0;
			if (++num_polys == max_num_polys) {
				printf("Num polys: %d / %d  ", num_polys, max_num_polys);
				return;
			}
		}
	}

	printf("Num polys: %d / %d  ", num_polys, max_num_polys);
}

#if 0
static void makeNormals(OutputMesh *om)
{
	int i, tv = om->vert;
	float npos[3];

	for (i = 0; i < tv; i++) {
		normal_short_to_float_v3(npos, om->verts[i].no);
		mul_v3_fl(npos, 0.1f);
		add_v3_v3(npos, om->verts[i].co);
		addEdge(om, i, addVert(om, npos, NULL));
	}
}
#endif // 0

static void initOutputMesh(OutputMesh *om, InputMesh *im)
{
	om->num_verts = om->alloc_verts = om->num_edges = 0;

	om->memarena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Output Mesh");
	om->vonvs = MEM_mallocN(sizeof(QRVertID) * im->num_verts, __func__);
	copy_vn_i(om->vonvs, im->num_verts, -1);
	om->ringe = MEM_callocN(sizeof(QREdge) * im->num_edges, "GFListEdges");
	om->ringf = MEM_callocN(sizeof(LinkNode *) * im->num_faces, "GFListFaces");
}

static void freeOutputMesh(OutputMesh *om)
{
	om->num_edges = om->num_verts = om->alloc_verts = 0;
	MEM_SAFE_FREE(om->verts);
	MEM_SAFE_FREE(om->vonvs);
	MEM_SAFE_FREE(om->ringe);
	MEM_SAFE_FREE(om->ringf);
	if (om->memarena) {
		BLI_memarena_free(om->memarena);
		om->memarena = NULL;
	}
}
#endif

static DerivedMesh *makeOrientationsMesh(QuadRemeshSystem *sys)
{
	int i, j;
	float mid[3], gra1[3], gra2[3];
	DerivedMesh *ret;
	InputMesh *im = &sys->input_mesh;
	MVert *verts;
	MEdge *edges;
#if 0
	ret = CDDM_new(im->num_verts * 5, im->num_verts * 4, 0, 0, 0);

	verts = ret->getVertArray(ret);
	edges = ret->getEdgeArray(ret);

	for (i = 0; i < im->num_verts; i++) {
		copy_v3_v3(mid, im->co[i]);
		copy_v3_v3(verts[i * 5].co, mid);

		for (j = 0; j < 4; j++) {
			madd_v3_v3v3fl(gra1, mid, sys->gfsys[0]->gf[i], 0.1f);
			copy_v3_v3(verts[i * 5 + j + 1].co, gra1);
			cross_v3_v3v3(gra2, sys->gfsys[0]->gf[i], sys->input_mesh.vno[i]);
			copy_v3_v3(sys->gfsys[0]->gf[i], gra2);

			edges[i * 4 + j].v1 = i * 5;
			edges[i * 4 + j].v2 = i * 5 + j + 1;
			edges[i * 4 + j].flag = ME_EDGEDRAW;
		}
	}
#else
	ret = CDDM_new(im->num_faces * 4, im->num_faces * 2, 0, 0, 0);

	verts = ret->getVertArray(ret);
	edges = ret->getEdgeArray(ret);

	for (i = 0; i < im->num_faces; i++) {
		mid_v3_v3v3v3(mid, im->co[im->faces[i][0]], im->co[im->faces[i][1]], im->co[im->faces[i][2]]);

		for (j = 0; j < 4; j++) {
			madd_v3_v3v3fl(gra1, mid, sys->cf[i][j], 0.04f);
			copy_v3_v3(verts[i * 4 + j].co, gra1);
		}

		edges[i * 2].v1 = i * 4;
		edges[i * 2].v2 = i * 4 + 2;
		edges[i * 2].flag = ME_EDGEDRAW;

		edges[i * 2 + 1].v1 = i * 4 + 1;
		edges[i * 2 + 1].v2 = i * 4 + 3;
		edges[i * 2 + 1].flag = ME_EDGEDRAW;
	}
#endif

	return ret;
}

#ifdef QR_GENERATE
static DerivedMesh *remesh(QuadRemeshSystem *sys)
{
	int i, num_loops, num_polys;
	double start_time;
	MVert *verts;
	OutputMesh *om = &sys->output_mesh;
	DerivedMesh *ret;

	start_time = PIL_check_seconds_timer();

	sys->rng = BLI_rng_new(sys->qmd->rng_seed);

	initOutputMesh(om, &sys->input_mesh);

	makeFeatureEdges(om, &sys->input_mesh);
	computeFlowLines(sys);
	//hideEdgesOnFaces(om, &sys->input_mesh);
	//makeNormals(om);

	if (!om->num_verts) {
		freeOutputMesh(om);
		BLI_rng_free(sys->rng);

		return NULL;
	}

#ifdef QR_MAKEPOLYS
	deleteDegenerateVerts(om);

	num_loops = om->num_edges * 2;
	num_polys = 2 + om->num_edges - om->num_verts;
#else
	num_loops = num_polys = 0;
#endif

	ret = CDDM_new(om->num_verts, om->num_edges, 0, num_loops, num_polys);

	verts = ret->getVertArray(ret);

	for (i = 0; i < om->num_verts; i++) {
		copy_v3_v3(verts[i].co, om->verts[i].co);
		normal_float_to_short_v3(verts[i].no, om->verts[i].no);
		verts[i].bweight = 0;
		verts[i].flag = 0;
	}

	makeEdges(om, ret->getEdgeArray(ret));

#ifdef QR_MAKEPOLYS
	makePolys(om, ret->getPolyArray(ret), ret->getLoopArray(ret));

	CDDM_recalc_tessellation(ret);
	//CDDM_calc_edges_tessface(ret);
#endif

	freeOutputMesh(om);
	BLI_rng_free(sys->rng);

	printf("Mesh generation time: %f ms\n", (PIL_check_seconds_timer() - start_time) * 1000.0f);

	return ret;
}
#endif

DerivedMesh *makeResultMesh(QuadRemeshSystem *sys, Object *ob, DerivedMesh *in)
{
	/* Get input mesh */
	if (!sys->input_mesh.is_alloc) {
		sys->qmd->flag |= MOD_QUADREMESH_INPUT_DIRTY;
	}

	if (sys->qmd->flag & MOD_QUADREMESH_INPUT_DIRTY) {
		if (sys->input_mesh.is_alloc) {
			freeInputMesh(&sys->input_mesh);
		}

		getInput(sys, ob, in);
		sys->qmd->flag |= MOD_QUADREMESH_FIELD_DIRTY;
	}

	/* Compute harmonic field */
	if (!sys->is_alloc) {
		sys->qmd->flag |= MOD_QUADREMESH_FIELD_DIRTY;
	}

	if (sys->qmd->flag & MOD_QUADREMESH_FIELD_DIRTY) {
		if (sys->is_alloc) {
			MEM_SAFE_FREE(sys->U_field);
			MEM_SAFE_FREE(sys->cf);

			sys->is_alloc = false;
		}

		//getHarmonicGradients(sys);
		getPrincipalCurvatures(sys);

		if (!sys->has_solution) {
			modifier_setError((ModifierData*)sys->qmd, "No soultion found given these constraints");
			return in;
		}

		sys->qmd->flag |= MOD_QUADREMESH_REMESH;
	}

#ifdef QR_GENERATE
	/* Remesh */
	if (sys->qmd->flag & MOD_QUADREMESH_REMESH) {
		if (sys->cache_mesh) {
			sys->cache_mesh->release(sys->cache_mesh);
		}

		sys->cache_mesh = remesh(sys);
	}
#endif

	sys->qmd->flag = 0;

#ifndef QR_SHOWORIENTATIONS
	if (sys->cache_mesh) {
		return CDDM_copy(sys->cache_mesh);
	}
#else
	return makeOrientationsMesh(sys);
#endif

	return NULL;
}