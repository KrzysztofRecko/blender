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
 * Author: Alexander Pinzon Fernandez
 * All rights reserved.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_quadremesh.c
 *  \ingroup modifiers
 */

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_string.h"
#include "BLI_rand.h"

#include "MEM_guardedalloc.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_particle.h"
#include "BKE_deform.h"

#include "MOD_util.h"
#include "MOD_quadremesh_geom.h"


#ifdef WITH_OPENNL

static void deleteLaplacianSystem(LaplacianSystem *sys)
{
	deleteGradientFlowSystem(sys->gfsys1);
	deleteGradientFlowSystem(sys->gfsys2);
	MEM_SAFE_FREE(sys->faces);
	MEM_SAFE_FREE(sys->edges);
	MEM_SAFE_FREE(sys->faces_edge);
	MEM_SAFE_FREE(sys->co);
	MEM_SAFE_FREE(sys->no);
	MEM_SAFE_FREE(sys->constraints);
	MEM_SAFE_FREE(sys->weights);
	MEM_SAFE_FREE(sys->U_field);
	MEM_SAFE_FREE(sys->h1);
	MEM_SAFE_FREE(sys->h2);
	MEM_SAFE_FREE(sys->gf1);
	MEM_SAFE_FREE(sys->gf2);
	MEM_SAFE_FREE(sys->ringf_indices);
	MEM_SAFE_FREE(sys->ringv_indices);
	MEM_SAFE_FREE(sys->ringe_indices);
	MEM_SAFE_FREE(sys->ringf_map);
	MEM_SAFE_FREE(sys->ringv_map);
	MEM_SAFE_FREE(sys->ringe_map);
	if (sys->context) {
		nlDeleteContext(sys->context);
	}
	MEM_SAFE_FREE(sys->mvert);
	MEM_SAFE_FREE(sys);
}

static void createFaceRingMap(LaplacianSystem *sys)
{
	int i, j, totalr = 0;
	int *index_iter;

	sys->ringf_map = MEM_callocN(sizeof(MeshElemMap) * sys->total_verts, "DeformRingMap");

	for (i = 0; i < sys->total_faces; i++) {
		for (j = 0; j < 3; j++) {
			sys->ringf_map[sys->faces[i][j]].count++;
			totalr++;
		}
	}

	sys->ringf_indices = MEM_callocN(sizeof(int) * totalr, "DeformRingIndex");
	index_iter = sys->ringf_indices;
	for (i = 0; i < sys->total_verts; i++) {
		sys->ringf_map[i].indices = index_iter;
		index_iter += sys->ringf_map[i].count;
		sys->ringf_map[i].count = 0;
	}

	for (i = 0; i < sys->total_faces; i++) {
		for (j = 0; j < 3; j++) {
			MeshElemMap *map = &sys->ringf_map[sys->faces[i][j]];
			map->indices[map->count++] = i;
		}
	}
}

static void createVertRingMap(LaplacianSystem *sys)
{
	int i, totalr = 0;
	int *index_iter;

	sys->ringv_map = MEM_callocN(sizeof(MeshElemMap) * sys->total_verts, "DeformNeighborsMap");

	for (i = 0; i < sys->total_edges; i++) {
		sys->ringv_map[sys->edges[i][0]].count++;
		sys->ringv_map[sys->edges[i][1]].count++;
		totalr += 2;
	}

	sys->ringv_indices = MEM_callocN(sizeof(int) * totalr, "DeformNeighborsIndex");
	index_iter = sys->ringv_indices;
	for (i = 0; i < sys->total_verts; i++) {
		sys->ringv_map[i].indices = index_iter;
		index_iter += sys->ringv_map[i].count;
		sys->ringv_map[i].count = 0;
	}

	for (i = 0; i < sys->total_edges; i++) {
		MeshElemMap *map1 = &sys->ringv_map[sys->edges[i][0]];
		MeshElemMap *map2 = &sys->ringv_map[sys->edges[i][1]];
		map1->indices[map1->count++] = sys->edges[i][1];
		map2->indices[map2->count++] = sys->edges[i][0];
	}
}

static void createEdgeRingMap(LaplacianSystem *sys)
{
	int i, totalr = 0;
	int *index_iter;

	sys->ringe_map = MEM_callocN(sizeof(MeshElemMap) * sys->total_verts, "DeformNeighborsMap");

	for (i = 0; i < sys->total_edges; i++) {
		sys->ringe_map[sys->edges[i][0]].count++;
		sys->ringe_map[sys->edges[i][1]].count++;
		totalr += 2;
	}

	sys->ringe_indices = MEM_callocN(sizeof(int) * totalr, "DeformNeighborsIndex");
	index_iter = sys->ringe_indices;
	for (i = 0; i < sys->total_verts; i++) {
		sys->ringe_map[i].indices = index_iter;
		index_iter += sys->ringe_map[i].count;
		sys->ringe_map[i].count = 0;
	}

	for (i = 0; i < sys->total_edges; i++) {
		MeshElemMap *map1 = &sys->ringe_map[sys->edges[i][0]];
		MeshElemMap *map2 = &sys->ringe_map[sys->edges[i][1]];
		map1->indices[map1->count++] = i;
		map2->indices[map2->count++] = i;
	}
}

static void computeFacesAdjacentToEdge(int fs[2], LaplacianSystem *sys, int indexe)
{
	int i, v1, v2, counter;
	int *fidn, numf;
	int *vin;
	v1 = sys->edges[indexe][0];
	v2 = sys->edges[indexe][1];
	numf = sys->ringf_map[v1].count;
	fidn = sys->ringf_map[v1].indices;
	counter = 0;
	fs[0] = -1;
	fs[1] = -1;

	for (i = 0; i < numf && counter < 2; i++) {
		vin = sys->faces[fidn[i]];
		if (vin[0] == v2 || vin[1] == v2 || vin[2] == v2) {
			fs[counter++] = fidn[i];
		}
	}
}

static void createFacesByEdge(LaplacianSystem *sys){
	int ei;
	for (ei = 0; ei < sys->total_edges; ei++) {
		computeFacesAdjacentToEdge(sys->faces_edge[ei], sys, ei);
	}
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

static void initLaplacianMatrix(LaplacianSystem *sys)
{
	float v1[3], v2[3], v3[3], no[3];
	float w2, w3;
	int j, fi;
	unsigned int idv1, idv2, idv3;

	for (fi = 0; fi < sys->total_faces; fi++) {
		const unsigned int *vidf = sys->faces[fi];

		idv1 = vidf[0];
		idv2 = vidf[1];
		idv3 = vidf[2];

		normal_tri_v3(no, sys->co[idv1], sys->co[idv2], sys->co[idv3]);
		copy_v3_v3(sys->no[fi], no);

		for (j = 0; j < 3; j++) {
			idv1 = vidf[j];
			idv2 = vidf[(j + 1) % 3];
			idv3 = vidf[(j + 2) % 3];

			copy_v3_v3(v1, sys->co[idv1]);
			copy_v3_v3(v2, sys->co[idv2]);
			copy_v3_v3(v3, sys->co[idv3]);

			w2 = cotangent_tri_weight_v3(v3, v1, v2);
			w3 = cotangent_tri_weight_v3(v2, v3, v1);

			if (sys->constraints[idv1] == 1) {
				nlMatrixAdd(idv1, idv1, w2 + w3);
			}
			else  {
				nlMatrixAdd(idv1, idv2, -w2);
				nlMatrixAdd(idv1, idv3, -w3);
				nlMatrixAdd(idv1, idv1, w2 + w3);
			}

		}
	}
	
}

static void computeScalarField(LaplacianSystem *sys)
{
	int vid, i, n;
	n = sys->total_verts;

#ifdef OPENNL_THREADING_HACK
	modifier_opennl_lock();
#endif
	
	nlNewContext();
	sys->context = nlGetCurrent();

	nlSolverParameteri(NL_NB_VARIABLES, n);
	nlSolverParameteri(NL_SYMMETRIC, NL_FALSE);
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_ROWS, n);
	nlSolverParameteri(NL_NB_RIGHT_HAND_SIDES, 1);
	nlBegin(NL_SYSTEM);
	
	for (i = 0; i < n; i++) {
		nlSetVariable(0, i, 0);
	}
		
	nlBegin(NL_MATRIX);
	
	initLaplacianMatrix(sys);

	for (i = 0; i < n; i++) {
		if (sys->constraints[i] == 1) {
			nlRightHandSideSet(0, i, sys->weights[i]);
		}
		else {
			nlRightHandSideSet(0, i, 0);
		}
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	
	if (nlSolveAdvanced(NULL, NL_TRUE)) {
		sys->has_solution = true;

		for (vid = 0; vid < sys->total_verts; vid++) {
			sys->U_field[vid] = nlGetVariable(0, vid);
		}	
	}
	else {
		sys->has_solution = false;
	}
#ifdef OPENNL_THREADING_HACK
	modifier_opennl_unlock();
#endif
}

/**
* Compute the gradient fields
*
* xi, xj, xk, are the vertices of the face
* ui, uj, uk, are the values of scalar fields for every vertex of the face
* n is the normal of the face.
* gf1 is the unknown field gradient 1.
* gf2 is the unknown field gradient 2.
*
* |xj - xi|         |uj - ui|
* |xk - xj| * gf1 = |uk - uj|
* |   nf  |         |   0   |
*
* gf2 = cross(n, gf1)
*/
static void computeGradientFields(LaplacianSystem * sys)
{
	int fi, i, j, k;
	float val, a[3][3], u[3], inv_a[3][3], gf1[3], g[3], w[3];

	for (fi = 0; fi < sys->total_faces; fi++) {
		const unsigned int *vidf = sys->faces[fi];
		i = vidf[0];
		j = vidf[1];
		k = vidf[2];
		sub_v3_v3v3(a[0], sys->co[j], sys->co[i]);
		sub_v3_v3v3(a[1], sys->co[k], sys->co[j]);
		copy_v3_v3(a[2], sys->no[fi]);

		/* Correct way*/
		transpose_m3(a);
		u[0] = sys->U_field[j] - sys->U_field[i];
		u[1] = sys->U_field[k] - sys->U_field[j];
		u[2] = 0;
		invert_m3_m3(inv_a, a);
		//mul_v3_m3v3(sys->gf1[fi], inv_a, u);
		mul_v3_m3v3(gf1, inv_a, u);

		/* Project Gradient fields on face*/
		normalize_v3_v3(g, gf1);
		val = dot_v3v3(g, sys->no[fi]);
		mul_v3_v3fl(u, sys->no[fi], val);
		sub_v3_v3v3(w, g, u);
		normalize_v3_v3(sys->gf1[fi], w);

		cross_v3_v3v3(g, sys->no[fi], sys->gf1[fi]);
		normalize_v3_v3(sys->gf2[fi], g);
		//cross_v3_v3v3(sys->gf2[fi], sys->no[fi], sys->gf1[fi]);
	}
}

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
	vb = sqrtf(r1) * ( 1.0f - r2);
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

static LaplacianSystem* initSystem(QuadRemeshModifierData *qmd, Object *ob, DerivedMesh *dm)
{
	int i, defgrp_index;
	float wpaint;
	MDeformVert *dvert = NULL;
	MDeformVert *dv = NULL, *d;
	LaplacianSystem *sys = MEM_callocN(sizeof(LaplacianSystem), "QuadRemeshCache");

	DM_ensure_tessface(dm);

	/* Get vertices */
	sys->total_verts = dm->getNumVerts(dm);
	sys->co = MEM_mallocN(sizeof(float[3]) * sys->total_verts, "QuadRemeshCoordinates");
	MVert *arrayvect = dm->getVertArray(dm);
	for (i = 0; i < dm->getNumVerts(dm); i++) {
		copy_v3_v3(sys->co[i], arrayvect[i].co);
	}

	/* Get edges */
	sys->edges = MEM_mallocN(sizeof(int[2]) * dm->getNumEdges(dm) * 2, "QuadRemeshEdges");
	MEdge *arrayedge = dm->getEdgeArray(dm);
	for (sys->total_edges = 0, i = 0; i < dm->getNumEdges(dm); i++) {
		sys->edges[sys->total_edges][0] = arrayedge[i].v1;
		sys->edges[sys->total_edges][1] = arrayedge[i].v2;
		sys->total_edges++;
	}

	/* Get faces */
	MFace *tessface = dm->getTessFaceArray(dm);
	sys->faces = MEM_mallocN(sizeof(int[3]) * dm->getNumTessFaces(dm) * 2, "QuadRemeshFaces");
	for (sys->total_faces = 0, i = 0; i < dm->getNumTessFaces(dm); i++) {
		sys->faces[sys->total_faces][0] = tessface[i].v1;
		sys->faces[sys->total_faces][1] = tessface[i].v2;
		sys->faces[sys->total_faces][2] = tessface[i].v3;
		sys->total_faces++;

		if (tessface[i].v4 != 0) {
			sys->faces[sys->total_faces][0] = tessface[i].v1;
			sys->faces[sys->total_faces][1] = tessface[i].v3;
			sys->faces[sys->total_faces][2] = tessface[i].v4;
			sys->total_faces++;

			sys->edges[sys->total_edges][0] = tessface[i].v1;
			sys->edges[sys->total_edges][1] = tessface[i].v3;
			sys->total_edges++;
		}
	}

	/* Get features */
	sys->constraints = MEM_callocN(sizeof(int) * sys->total_verts, __func__);
	sys->weights = MEM_callocN(sizeof(float) * sys->total_verts, __func__);

	modifier_get_vgroup(ob, dm, qmd->anchor_grp_name, &dvert, &defgrp_index);
	BLI_assert(dvert != NULL);

	dv = dvert;
	for (i = 0, sys->total_features = 0; i < sys->total_verts; i++, dv++) {
		d = defvert_find_index(dv, defgrp_index);
		wpaint = defvert_find_weight(dv, defgrp_index);

		if (d && wpaint < 0.19 || wpaint > 0.89) {
			sys->constraints[i] = 1;
			sys->weights[i] = -1.0f + wpaint * 2.0f;
			sys->total_features++;
		}
	}

	/* Allocate rest of memory */
	sys->faces_edge = MEM_mallocN(sizeof(int[2]) * sys->total_edges, "QuadRemeshFacesEdge");
	sys->no = MEM_callocN(sizeof(float[3]) * sys->total_faces, "QuadRemeshNormals");

	sys->gf1 = MEM_mallocN(sizeof(float[3]) * sys->total_faces, "QuadRemeshGradientField1");
	sys->gf2 = MEM_mallocN(sizeof(float[3]) * sys->total_faces, "QuadRemeshGradientField2");
	sys->U_field = MEM_mallocN(sizeof(float) * sys->total_verts, "QuadRemeshUField");
	sys->h1 = MEM_mallocN(sizeof(float) * sys->total_verts, "QuadRemeshH1");
	sys->h2 = MEM_mallocN(sizeof(float) * sys->total_verts, "QuadRemeshH2");
	sys->gfsys1 = sys->gfsys2 = NULL;
	
	createFaceRingMap(sys);
	createVertRingMap(sys);
	createEdgeRingMap(sys);

	createFacesByEdge(sys);

	return sys;
}

static LaplacianSystem *QuadRemeshModifier_do(QuadRemeshModifierData *qmd, Object *ob, DerivedMesh *dm)
{
	int i;
	LaplacianSystem *sys = NULL;
	
	int defgrp_index;
	MDeformVert *dvert = NULL;
	MDeformVert *dv = NULL;
	float mmin = 1000, mmax = 0;
	float y;
	int x;

	qmd->flag |= MOD_QUADREMESH_COMPUTE_FLOW;
	qmd->flag |= MOD_QUADREMESH_REMESH;
	strcpy(qmd->anchor_grp_name, "Group");

	if (qmd->flag & MOD_QUADREMESH_COMPUTE_FLOW) {
		if (strlen(qmd->anchor_grp_name) >= 1) {
			if (qmd->cache_system) {
				sys = qmd->cache_system;
				deleteLaplacianSystem(sys);
			}
			qmd->cache_system = initSystem(qmd, ob, dm);
			sys = qmd->cache_system;

			computeScalarField(sys);

			if (sys->has_solution) {
				computeGradientFields(sys);
				computeSampleDistanceFunctions(sys, 2.0f, 10.0f);

				/* normalization of vgroup weights */
				/*if (!defgroup_find_name(ob, "QuadRemeshFlow")) {
					BKE_defgroup_new(ob, "QuadRemeshFlow");
					modifier_get_vgroup(ob, dm, "QuadRemeshFlow", &dvert, &defgrp_index);
					BLI_assert(dvert != NULL);
					dv = dvert;

					for (i = 0; i < numVerts; i++) {
						mmin = min_ff(mmin, sys->U_field[i]);
						mmax = max_ff(mmax, sys->U_field[i]);
					}

					for (i = 0; i < numVerts; i++) {
						y = (sys->U_field[i] - mmin) / (mmax - mmin);
						x = y * 60;
						y = (x % 2 == 0 ? 0.1 : 0.9);
						defvert_add_index_notest(dv, defgrp_index, y);
						dv++;
					}
				}*/
			}
		}
		
		qmd->flag &= ~MOD_QUADREMESH_COMPUTE_FLOW;
	}

	if (qmd->flag & MOD_QUADREMESH_REMESH && qmd->cache_system) {
		sys = qmd->cache_system;
		if (sys->has_solution) {
			sys->h = 2.0f;
			computeFlowLines(sys);
		}
		qmd->flag &= ~MOD_QUADREMESH_REMESH;
	}

	if (qmd->cache_system) {
		sys = qmd->cache_system;
	}

	return sys;
}


#else  /* WITH_OPENNL */
static void QuadRemeshModifier_do(
        QuadRemeshModifierData *lmd, Object *ob, DerivedMesh *dm,
        float (*vertexCos)[3], int numVerts)
{
	(void)lmd, (void)ob, (void)dm, (void)vertexCos, (void)numVerts;
}
#endif  /* WITH_OPENNL */

static void initData(ModifierData *md)
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;
	lmd->anchor_grp_name[0] = '\0';
	lmd->flag = 0;
	lmd->cache_system = NULL;
}

static void copyData(ModifierData *md, ModifierData *target)
{
	QuadRemeshModifierData *qmd = (QuadRemeshModifierData *)md;
	QuadRemeshModifierData *tqmd = (QuadRemeshModifierData *)target;
	modifier_copyData_generic(md, target);
	tqmd->cache_system = NULL;

}

static bool isDisabled(ModifierData *md, int UNUSED(useRenderParams))
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;
	//if (lmd->anchor_grp_name[0]) return 0;
	//return 1;
	return 0;
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *md)
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;
	CustomDataMask dataMask = 0;
	if (lmd->anchor_grp_name[0]) dataMask |= CD_MASK_MDEFORMVERT;
	return dataMask;
}

#if 0
static void deformVerts(ModifierData *md, Object *ob, DerivedMesh *derivedData,
						float(*vertexCos)[3], int numVerts, ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *dm = get_dm(ob, NULL, derivedData, NULL, false, false);

	QuadRemeshModifier_do((QuadRemeshModifierData *)md, ob, dm);
	if (dm != derivedData) {
		dm->release(dm);
	}
}
#endif // 0


static DerivedMesh *applyModifier(ModifierData *md,
	Object *ob,
	DerivedMesh *dm,
	ModifierApplyFlag UNUSED(flag))
{
	//DerivedMesh *dm2 = get_dm(ob, NULL, dm, NULL, false, false);
	//QuadRemeshModifier_do((QuadRemeshModifierData *)md, ob, dm, (void *)dm->getVertArray(dm), dm->getNumVerts(dm));

	MVert *arrayvect;
	MEdge *arrayedge;
	int i;
	DerivedMesh *result;
	
	LaplacianSystem *sys = QuadRemeshModifier_do((QuadRemeshModifierData *)md, ob, dm);
	
	if (sys) {
		/*result = CDDM_new(gfsys->totalf * 2, gfsys->totalf, 0, 0, 0);
		arrayvect = result->getVertArray(result);
		for (i = 0; i < gfsys->totalf; i++) {
			float cent[3], v[3];
			cent_tri_v3(cent, sys->co[sys->faces[i][0]], sys->co[sys->faces[i][1]], sys->co[sys->faces[i][2]]);
			mul_v3_fl(gfsys->gfield[i], 0.1f);
			add_v3_v3v3(v, cent, gfsys->gfield[i]);
			copy_v3_v3(arrayvect[i * 2].co, v);
			copy_v3_v3(arrayvect[i * 2 + 1].co, cent);
		}
		arrayedge = result->getEdgeArray(result);
		for (i = 0; i < gfsys->totalf; i++) {
			arrayedge[i].v1 = i * 2;
			arrayedge[i].v2 = i * 2 + 1;
			arrayedge[i].flag |= ME_EDGEDRAW;
		}*/

		result = CDDM_new(sys->totvert,
						  sys->gfsys1->totedge + sys->gfsys2->totedge,
						  0, 0, 0);
		arrayvect = result->getVertArray(result);
		for (i = 0; i < sys->totvert; i++) {
			copy_v3_v3(arrayvect[i].co, sys->mvert[i].co);
		}

		arrayedge = result->getEdgeArray(result);
		for (i = 0; i < sys->gfsys1->totedge; i++) {
			arrayedge[i].v1 = sys->gfsys1->medge[i].v1;
			arrayedge[i].v2 = sys->gfsys1->medge[i].v2;
			arrayedge[i].flag |= ME_EDGEDRAW;
		}
		for (i = 0; i < sys->gfsys2->totedge; i++) {
			arrayedge[i + sys->gfsys1->totedge].v1 = sys->gfsys2->medge[i].v1;
			arrayedge[i + sys->gfsys1->totedge].v2 = sys->gfsys2->medge[i].v2;
			arrayedge[i + sys->gfsys1->totedge].flag |= ME_EDGEDRAW;
		}
	}
	else{
		result = dm;
	}
	
	//CDDM_calc_edges_tessface(result);

	return result;
}

static void freeData(ModifierData *md)
{
#ifdef WITH_OPENNL
	QuadRemeshModifierData *qmd = (QuadRemeshModifierData *)md;
	LaplacianSystem *sys = (LaplacianSystem *)qmd->cache_system;
	if (sys) {
		deleteLaplacianSystem(sys);
	}
#endif
}

ModifierTypeInfo modifierType_QuadRemesh = {
	/* name */              "QuadRemesh",
	/* structName */        "QuadRemeshModifierData",
	/* structSize */        sizeof(QuadRemeshModifierData),
	/* type */              eModifierTypeType_Nonconstructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
							eModifierTypeFlag_AcceptsCVs |
							eModifierTypeFlag_SupportsEditmode,
	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          freeData,
	/* isDisabled */        isDisabled,
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
