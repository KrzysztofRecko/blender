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

/** \file blender/modifiers/intern/MOD_quadremesh_geom.c
 *  \ingroup modifiers
 */

#include "MEM_guardedalloc.h"
#include "MOD_quadremesh_geom.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_string.h"
#include "BLI_rand.h"
#include "BLI_linklist.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_particle.h"
#include "BKE_deform.h"

#include "MOD_util.h"

//#define QR_SHOWQUERIES
#define QR_LINELIMIT 40000

void estimateNumberGFVerticesEdges(int ve[2], LaplacianSystem *sys, float h)
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

GradientFlowMesh *newGradientFlowMesh(int totalvert, int totaledge)
{
	GradientFlowMesh *gfmesh = MEM_mallocN(sizeof(GradientFlowMesh), "GradientFlowMesh");
	gfmesh->mvert = MEM_mallocN(sizeof(GFVert) * totalvert, "GFVert");
	gfmesh->medge = MEM_mallocN(sizeof(GFEdge) * totaledge, "GFEdge");
	gfmesh->totvert = 0;
	gfmesh->totedge = 0;
	gfmesh->allocvert = totalvert;
	gfmesh->allocedge = totaledge;

	return gfmesh;
}

void deleteGradientFlowMesh(GradientFlowMesh * gfmesh) 
{
	if (gfmesh) {
		MEM_SAFE_FREE(gfmesh->mvert);
		MEM_SAFE_FREE(gfmesh->medge);
		MEM_SAFE_FREE(gfmesh);
	}
}

GFVertID addVertGFMesh(GradientFlowMesh *gfmesh, float co[3])
{
	if (gfmesh->totvert == gfmesh->allocvert) {
		gfmesh->allocvert = gfmesh->allocvert * 2 + 10;
		gfmesh->mvert = MEM_reallocN(gfmesh->mvert, sizeof(GFVert) * gfmesh->allocvert);
	}

	copy_v3_v3(gfmesh->mvert[gfmesh->totvert].co, co);
	gfmesh->totvert++;

	return gfmesh->totvert - 1;
}

GFVertID addGFVertGFMesh(GradientFlowMesh *gfmesh, GFVert gfvert)
{
	return addVertGFMesh(gfmesh, gfvert.co);
}

GFEdgeID addEdgeGFMesh(GradientFlowMesh *gfmesh, int index_v1, int index_v2, int index_face)
{
	if (gfmesh->totedge == gfmesh->allocedge) {
		gfmesh->allocedge += MOD_QUADREMESH_ALLOC_BLOCK;
		gfmesh->medge = MEM_reallocN(gfmesh->medge, sizeof(GFEdge) * gfmesh->allocedge);
	}

	gfmesh->medge[gfmesh->totedge].face = index_face;
	gfmesh->medge[gfmesh->totedge].v1 = index_v1;
	gfmesh->medge[gfmesh->totedge].v2 = index_v2;
	gfmesh->totedge++;

	return gfmesh->totedge - 1;
}

GFEdgeID addGFEdgeGFMesh(GradientFlowMesh *gfmesh, GFEdge gfedge)
{
	return addEdgeGFMesh(gfmesh, gfedge.v1, gfedge.v2, gfedge.face);
}

/*
* List of vertices from original mesh with special features (edge dihedral angle less that 90) to be preserves
* return the size of array
*/
int *findFeaturesOnMesh(int size[2], LaplacianSystem *sys)
{
	int i, f1, f2, total;
	float angle;
	int *listverts = MEM_callocN(sizeof(int) * sys->total_verts, __func__);
	int *listdest = NULL;
	total = 0;

	for (i = 0; i < sys->total_edges; i++) {
		f1 = sys->faces_edge[i][0];
		f2 = sys->faces_edge[i][1];
		angle = angle_normalized_v3v3(sys->no[f1], sys->no[f2]);
		if (angle >= M_PI_2) {
			listverts[sys->edges[i][0]] = 1;
			listverts[sys->edges[i][1]] = 1;
		}
	}

	for (i = 0; i < sys->total_verts; i++) {
		if (sys->constraints[i] == 1) {
			listverts[i] = 1;
		}
	}

	for (i = 0; i<sys->total_verts; i++) {
		if (listverts[i] == 1) {
			total++;
		}
	}
	if (total > 0) {
		listdest = MEM_mallocN(sizeof(int)* total, __func__);
	}
	total = 0;
	for (i = 0; i<sys->total_verts; i++) {
		if (listverts[i] == 1) {
			listdest[total++] = i;
		}
	}
	MEM_SAFE_FREE(listverts);
	size[0] = total;
	return listdest;
}

void addSeedToQueue(Heap *aheap, float in_co[3], bool is_vertex, int in_val, float weight)
{
	GFSeed *seed = MEM_mallocN(sizeof(GFSeed), __func__);
	copy_v3_v3(seed->co, in_co);
	seed->val = in_val;
	seed->type = is_vertex ? eSeedVert : eSeedFace;

	BLI_heap_insert(aheap, weight, seed);
}

GradientFlowSystem *newGradientFlowSystem(LaplacianSystem *sys, float *mhfunction, float(*mgfield)[3])
{
	int ve[2], i;
	int *lverts, sizeverts[2];
	GradientFlowSystem *gfsys = MEM_callocN(sizeof(GradientFlowSystem), "GradientFlowSystem");
	lverts = NULL;

	estimateNumberGFVerticesEdges(ve, sys, sys->h);
	
	gfsys->sys = sys;
	gfsys->mesh = newGradientFlowMesh(ve[0], ve[1]);
	
	gfsys->ringf_list = MEM_callocN(sizeof(LinkNode *) * sys->total_faces, "GFListFaces");
	gfsys->ringe_list = MEM_callocN(sizeof(LinkNode *) * sys->total_edges, "GFListEdges");
	
	gfsys->heap_seeds = BLI_heap_new();
	
	gfsys->totalf = sys->total_faces;
	gfsys->totale = sys->total_edges;
	
	lverts = findFeaturesOnMesh(sizeverts, sys);
	
	for (i = 0; i < sizeverts[0]; i++) {
		addSeedToQueue(gfsys->heap_seeds, sys->co[lverts[i]], true, lverts[i], 0.0f);
	}
	
	gfsys->hfunction = mhfunction;
	gfsys->gfield = mgfield;
	MEM_SAFE_FREE(lverts);
	
	return gfsys;
}

void deleteGradientFlowSystem(GradientFlowSystem *gfsys) 
{
	int i;
	if (gfsys) {
		for (i = 0; i < gfsys->totalf; i++) {
			BLI_linklist_free(gfsys->ringf_list[i], NULL);
		}
		for (i = 0; i < gfsys->totale; i++) {
			BLI_linklist_free(gfsys->ringe_list[i], NULL);
		}
		deleteGradientFlowMesh(gfsys->mesh);
		MEM_SAFE_FREE(gfsys->ringf_list);
		MEM_SAFE_FREE(gfsys->ringe_list);
		BLI_heap_free(gfsys->heap_seeds, MEM_freeN);
		MEM_SAFE_FREE(gfsys);
	}
}

GFVertID addVertGFSystem(GradientFlowSystem *gfsys, float co[3])
{
	GFVertID ret = addVertGFMesh(gfsys->mesh, co);

	return ret;
}

GFVertID addGFVertGFSystem(GradientFlowSystem *gfsys, GFVert gfvert)
{
	return addVertGFSystem(gfsys, gfvert.co);
}

GFEdgeID addEdgeGFSystem(GradientFlowSystem *gfsys, int index_v1, int index_v2, int in_f)
{
	GFEdgeID ret = addEdgeGFMesh(gfsys->mesh, index_v1, index_v2, in_f);

	BLI_linklist_prepend(&gfsys->ringf_list[in_f], (void*)ret);

	return ret;
}

GFEdgeID addGFEdgeGFSystem(GradientFlowSystem *gfsys, GFEdge gfedge)
{
	return addEdgeGFSystem(gfsys, gfedge.v1, gfedge.v2, gfedge.face);
}

#if 0
int addEdgeTwoFacesGFSystem(GradientFlowSystem *gfsys, int index_v1, int index_v2, int index_face1, int index_face2)
{
	int pos = addEdgeGFMesh(gfsys->mesh, index_v1, index_v2, index_face1);
	if (index_face1 >= 0) {
		if (gfsys->ringf_list[index_face1]) {
			addNodeGFList(gfsys->ringf_list[index_face1], pos);
		}
		else{
			gfsys->ringf_list[index_face1] = newGFList(pos);
		}
	}
	if (index_face2 >= 0) {
		if (gfsys->ringf_list[index_face2]) {
			addNodeGFList(gfsys->ringf_list[index_face2], pos);
		}
		else{
			gfsys->ringf_list[index_face2] = newGFList(pos);
		}
	}
	return pos;

}
#endif // 0


GFSeed *getTopSeedFromQueue(struct Heap *aheap)
{
	GFSeed *vert = BLI_heap_popmin(aheap);

	return vert;
}

bool isOnSegmentLine(float p1[3], float p2[3], float q[3]){
	if (fabsf(len_v3v3(p1, q) + len_v3v3(p2, q) - len_v3v3(p1, p2)) < MOD_QUADREMESH_MIN_LEN) {
		return true;
	}
	return false;
}

/*
* Return 1 if the intersections exist
* Return -1 if the intersections does not exist
*/
bool intersecionLineSegmentWithVector(float r[3], float p1[3], float p2[3], float ori[3], float dir[3])
{
	float v[3], i1[3], i2[3];
	int i;

	add_v3_v3v3(v, ori, dir);
	i = isect_line_line_v3(p1, p2, ori, v, i1, i2);
	if (i == 0) {
		sub_v3_v3v3(i1, p1, ori);
		normalize_v3(i1);
		if (equals_v3v3(i1, dir)) {
			copy_v3_v3(r, p1);
		}
		else {
			copy_v3_v3(r, p2);
		}
	}
	else {
		sub_v3_v3v3(v, i1, ori);
		normalize_v3(v);
		if (equals_v3v3(v, dir)) {
			if (isOnSegmentLine(p1, p2, i1)) {
				copy_v3_v3(r, i1);
			}
			else{
				return false;
			}
		}
		else {
			return false;
		}
	}
	return true;
}

bool intersectionVectorWithTriangle(float r[3], float p1[3], float p2[3], float p3[3], float ori[3], float dir[3])
{
	if (intersecionLineSegmentWithVector(r, p1, p2, ori, dir) == 1) {
		return true;
	}
	else if (intersecionLineSegmentWithVector(r, p2, p3, ori, dir) == 1) {
		return true;
	}
	else if (intersecionLineSegmentWithVector(r, p3, p1, ori, dir) == 1) {
		return true;
	}
	return false;

}

int getEdgeFromVerts(LaplacianSystem *sys, int v1, int v2)
{
	int *eidn, nume, i;
	nume = sys->ringe_map[v1].count;
	eidn = sys->ringe_map[v1].indices;
	for (i = 0; i < nume; i++) {
		if (sys->edges[eidn[i]][0] == v2 || sys->edges[eidn[i]][1] == v2){
			return eidn[i];
		}
	}

	return -1;
}

int getOtherFaceAdjacentToEdge(LaplacianSystem *sys, int oldface, int e)
{
	if (sys->faces_edge[e][0] == oldface) {
		return sys->faces_edge[e][1];
	}

	return sys->faces_edge[e][0];
}

/* Project Gradient fields on face*/
void projectVectorOnFace(float r[3], float no[3], float dir[3])
{
	float g[3], val, u[3], w[3];
	normalize_v3_v3(g, dir);
	val = dot_v3v3(g, no);
	mul_v3_v3fl(u, no, val);
	sub_v3_v3v3(w, g, u);
	normalize_v3_v3(r, w);
}

int getThirdVert(LaplacianSystem *sys, int oldface, int v1, int v2)
{
	int a, b;

	if      (v1 == sys->faces[oldface][0]) a = 0;
	else if (v1 == sys->faces[oldface][1]) a = 1;
	else if (v1 == sys->faces[oldface][2]) a = 2;
	else a = 4;

	if      (v2 == sys->faces[oldface][0]) b = 0;
	else if (v2 == sys->faces[oldface][1]) b = 1;
	else if (v2 == sys->faces[oldface][2]) b = 2;
	else b = 4;

	BLI_assert(a + b < 4 && a != b);
	return sys->faces[oldface][3 - a - b];
}

bool getSecondAndThirdVert(int *r1, int *r2, LaplacianSystem *sys, int face, int v)
{
	int a;

	if      (v == sys->faces[face][0]) a = 0;
	else if (v == sys->faces[face][1]) a = 1;
	else if (v == sys->faces[face][2]) a = 2;
	else return false;

	*r1 = sys->faces[face][(a + 1) % 3];
	*r2 = sys->faces[face][(a + 2) % 3];

	return true;
}

/**
 * /return true - different point reached
 *         false - the same point reached
 */
bool nextPoint(float r_co[3], int *r_edge, GradientFlowSystem *gfsys, int in_f, float in_co[3], float in_dir[3])
{
	int i, pick = -1, v = -1;
	bool is_on_vertex = false;
	float a[3][3], b[3][3], c[2][3], co2[3],
		dummy[3], zero[3] = { 0.0f, 0.0f, 0.0f };
	LaplacianSystem *sys = gfsys->sys;

	/* check if direction is coplanar to triangle */
	/* check if point is inside triangle */
	/* check if triangle is degenerate */

	add_v3_v3v3(co2, in_dir, in_co); /* second point on direction */

	for (i = 0; i < 3; i++)
		sub_v3_v3v3(a[i], in_co, sys->co[sys->faces[in_f][i]]);

	for (i = 0; i < 3; i++)
		if (compare_v3v3(a[i], zero, 0.00001f))
			v = i;

	if (v != -1) {
		is_on_vertex = true;
		normalize_v3_v3(c[0], a[(v + 1) % 3]);
		normalize_v3_v3(c[1], a[(v + 2) % 3]);
		add_v3_v3(c[0], c[1]);
		mul_v3_v3fl(a[v], c[0], -0.5f);
	}

	for (i = 0; i < 3; i++) {
		cross_v3_v3v3(b[i], in_dir, a[i]);
		dummy[i] = dot_v3v3(b[i], sys->no[in_f]);
	}

	for (i = 0; i < 3; i++)
		if (dummy[i] < 0.0f && dummy[(i + 1) % 3] >= 0.0f)
			pick = i;

	if (pick == -1)
		*r_edge = 0;
	else
		*r_edge = getEdgeFromVerts(sys, sys->faces[in_f][pick], sys->faces[in_f][(pick + 1) % 3]);

	isect_line_line_v3(sys->co[sys->faces[in_f][pick]], sys->co[sys->faces[in_f][(pick + 1) % 3]], in_co, co2, r_co, dummy);

	if (is_on_vertex && pick != (v + 1) % 3)
		return false;
	else
		return true;

#if 0
		if (compare_v3v3(r_co, in_co, 0.001f)) { /* result too close to input */
			sub_v3_v3v3(dummy, sys->co[sys->edges[*r_edge][1]], sys->co[sys->edges[*r_edge][0]]);
			if (dot_v3v3(dummy, in_dir) < 0.0f) {
				copy_v3_v3(r_co, sys->co[sys->edges[*r_edge][0]]);
				*r_edge = sys->edges[*r_edge][0];
			}
			else {
				copy_v3_v3(r_co, sys->co[sys->edges[*r_edge][1]]);
				*r_edge = sys->edges[*r_edge][1];
			}
			return false;
		}
		else 


			return true;

	else { /* input on vertex v of face */
		sub_v3_v3v3(b[0], sys->co[sys->faces[in_f][(v + 1) % 3]], in_co);
		sub_v3_v3v3(b[1], sys->co[sys->faces[in_f][(v + 2) % 3]], in_co);

		cross_v3_v3v3(c[0], in_dir, b[0]);
		cross_v3_v3v3(c[1], b[1], in_dir);

		if (dot_v3v3(c[0], sys->no[in_f]) < 0.0f) { /* direction points away from triangle */
			copy_v3_v3(r_co, in_co);
			*r_edge = getEdgeFromVerts(sys, sys->faces[in_f][(v + 1) % 3], sys->faces[in_f][v]);

			return false;
		}
		else if (dot_v3v3(c[1], sys->no[in_f]) < 0.0f) {
			copy_v3_v3(r_co, in_co);
			*r_edge = getEdgeFromVerts(sys, sys->faces[in_f][(v + 2) % 3], sys->faces[in_f][v]);

			return false;
		}
		else { /* direction points towards an edge */
			isect_line_line_v3(sys->co[sys->faces[in_f][(v + 1) % 3]], sys->co[sys->faces[in_f][(v + 2) % 3]], in_co, co2, r_co, dummy);
			*r_edge = getEdgeFromVerts(sys, sys->faces[in_f][(v + 1) % 3], sys->faces[in_f][(v + 2) % 3]);

			return true;
		}
	}
#endif

#if 0

	iv1 = sys->edges[in_e][0];
	iv2 = sys->edges[in_e][1];
	ov1 = getThirdVert(sys, in_f, iv1, iv2);

	sub_v3_v3v3(a, sys->co[ov1], in_co);
	sub_v3_v3v3(e, sys->co[iv1], sys->co[iv2]);

	cross_v3_v3v3(b, a, e);
	if (dot_v3v3(b, sys->no[in_f]) < 0.0f) {
		mul_v3_fl(e, -1.0f);
		iv1 = sys->edges[in_e][1];
		iv2 = sys->edges[in_e][0];
	}

	cross_v3_v3v3(b, a, in_dir);
	cross_v3_v3v3(c, e, in_dir);

	if (dot_v3v3(b, sys->no[in_f]) > 0.0f) {
		*r_edge = getEdgeFromVerts(sys, iv1, ov1);
		if (dot_v3v3(c, sys->no[in_f]) > 0.0f) { /* vertex on iv1 */
			copy_v3_v3(r_co, sys->co[iv1]);
			return false;
		}
		else {                                   /* vertex on edge iv1 - ov1 */
			add_v3_v3v3(c, in_dir, in_co);
			isect_line_line_v3(sys->co[iv1], sys->co[ov1], in_co, c, r_co, dummy);
		}
	}
	else {
		*r_edge = getEdgeFromVerts(sys, iv2, ov1);
		if (dot_v3v3(c, sys->no[in_f]) > 0.0f) { /* vertex on iv2 */
			copy_v3_v3(r_co, sys->co[iv2]);
			return false;
		}
		else {                                   /* vertex on edge iv2 - ov1 */
			add_v3_v3v3(c, in_dir, in_co);
			isect_line_line_v3(sys->co[iv2], sys->co[ov1], in_co, c, r_co, dummy);
}
	}

	return true;
#endif // 0

}

#if 0
bool nextPointFromVertex(float r_co[3], int *r_face, int *r_edge, GradientFlowSystem *gfsys, int in_v)
{
	int i, f, vf1, vf2;
	float gf[3], a[3], b[3], alfa, beta, gamma, dummy[3];
	LaplacianSystem *sys = gfsys->sys;

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
#endif

bool intersectSegmentWithOthersOnFace(GradientFlowSystem *gfsys, float in_a[3], float in_b[3], int in_f)
{
	int e;
	float dummy[3], lambda;
	LinkNode *iter;
	LaplacianSystem *sys = gfsys->sys;
	GFVert *vertices = gfsys->mesh->mvert;
	GFEdge *edges = gfsys->mesh->medge;

	for (iter = gfsys->ringf_list[in_f]; iter; iter = iter->next) {
		e = (int)iter->link;
		if (isect_line_line_strict_v3(in_a, in_b, vertices[edges[e].v1].co, vertices[edges[e].v2].co, dummy, &lambda))
			return true;
	}

	return false;
}

/**
 * 0 - intersection found
 * 1 - intersection not found
 * 2 - wrong direction for this face 
 */
int queryDirection(GradientFlowSystem *gfsys, float in_co[3], int in_f, float in_dir[3], float dist, float maxdist, bool make_seed)
{
	int e, oldf;
	float c[3], len, actlen, oldco[3], newco[3], newco2[3], dir[3];
	LaplacianSystem *sys = gfsys->sys;

	copy_v3_v3(dir, in_dir);
	copy_v3_v3(oldco, in_co);
	oldf = in_f;
	actlen = 0.0f;

	while (1) {
		len = dot_v3v3(dir, sys->no[in_f]);
		mul_v3_v3fl(c, sys->no[in_f], len);
		sub_v3_v3(dir, c);
		if (normalize_v3(dir) < FLT_EPSILON) 
			return 2;
		if (dot_v3v3(sys->no[oldf], sys->no[in_f]) < 0.0f) mul_v3_fl(dir, -1.0f);

		if (!nextPoint(newco, &e, gfsys, in_f, oldco, dir)) return 2;
		oldf = in_f;

		sub_v3_v3v3(c, newco, oldco);
		len = len_v3(c);
		actlen += len;

		copy_v3_v3(newco2, newco);

		if (actlen - len < dist) {
			if (actlen > dist) {
				mul_v3_v3fl(newco2, c, (dist - actlen + len) / len);
				add_v3_v3(newco2, oldco);
			}

#ifdef QR_SHOWQUERIES
			int vf1, vf2;
			vf1 = addVertGFSystem(gfsys, oldco);
			vf2 = addVertGFSystem(gfsys, newco2);
			addEdgeGFSystem(gfsys, vf1, vf2, 0);
#endif

			if (intersectSegmentWithOthersOnFace(gfsys, oldco, newco2, in_f)) return 0;
		}
		else if (!make_seed)
			return 1;

		if (make_seed && actlen > maxdist) {
			mul_v3_v3fl(newco2, c, (maxdist - actlen + len) / len);
			add_v3_v3(newco2, oldco);

#ifdef QR_SHOWQUERIES
			int vf1, vf2;
			vf1 = addVertGFSystem(gfsys, oldco);
			vf2 = addVertGFSystem(gfsys, newco2);
			addEdgeGFSystem(gfsys, vf1, vf2, 0);
#endif

			addSeedToQueue(gfsys->heap_seeds, newco2, false, in_f, 0.0f);
			return 1;
		}
#ifdef QR_SHOWQUERIES
		else {
			int vf1, vf2;
			vf1 = addVertGFSystem(gfsys, oldco);
			vf2 = addVertGFSystem(gfsys, newco2);
			addEdgeGFSystem(gfsys, vf1, vf2, 0);
		}
#endif

		in_f = getOtherFaceAdjacentToEdge(sys, in_f, e);
		copy_v3_v3(oldco, newco);
	}
}

void getOrthogonalDirection(float r_dir[3], GradientFlowSystem *gfsys, float in_dir[3], int in_f)
{
	LaplacianSystem *sys = gfsys->sys;

	/*if (in_ev < 0) {
		in_ev = -in_ev - 1;
		for (i = 0; i < sys->ringf_map[in_ev].count; i++) {
			f = sys->ringf_map[in_ev].indices[i];
			add_v3_v3(no, sys->no[f]);
		}
		mul_v3_fl(no, 1.0f / (float)sys->ringf_map[in_ev].count);
	}
	else {
		add_v3_v3(no, sys->no[sys->faces_edge[in_ev][0]]);
		add_v3_v3(no, sys->no[sys->faces_edge[in_ev][1]]);
		mul_v3_fl(no, 0.5f);
	}*/

	cross_v3_v3v3(r_dir, sys->no[in_f], in_dir);
	normalize_v3(r_dir);
}

bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], float in_newco[3], int in_f, float dist, float maxdist)
{
	int d;
	bool make_seed = BLI_frand() > 0.75f;
	float seg[3], dir[3];
	LaplacianSystem *sys = gfsys->sys;

	sub_v3_v3v3(seg, in_oldco, in_newco);
	getOrthogonalDirection(dir, gfsys, seg, in_f);

	for (d = 0; d < 2; d++) {
		if (!queryDirection(gfsys, in_newco, in_f, dir, dist, maxdist, make_seed))
			return false;
		mul_v3_fl(dir, -1.0f);
	}

	return true;
}

#if 0

float getSamplingDistanceFunctionOnFace(LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexface)
{
	float h1, h2, h3;
	h1 = gfsys->hfunction[sys->faces[indexface][0]];
	h2 = gfsys->hfunction[sys->faces[indexface][1]];
	h3 = gfsys->hfunction[sys->faces[indexface][2]];
	return min_fff(h1, h2, h3);
}

float getMaxSamplingDistanceFunctionOnFace(LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexface)
{
	float h1, h2, h3;
	h1 = gfsys->hfunction[sys->faces[indexface][0]];
	h2 = gfsys->hfunction[sys->faces[indexface][1]];
	h3 = gfsys->hfunction[sys->faces[indexface][2]];
	return max_fff(h1, h2, h3);
}

#endif // 0

void computeGFLine(GradientFlowSystem *gfsys, GFSeed *in_seed)
{
	int i, f, d, e, v, newe;
	bool is_vertex;
	float newco[3], oldco[3], gf[3], dir = 1.0f;
	LaplacianSystem *sys = gfsys->sys;
	GFVertID oldgfv, newgfv, seed;

	seed = -1;

	for (d = 0; d < 2; d++) {
		if (in_seed->type == eSeedVert) {
			is_vertex = true;
			v = in_seed->val;
			f = 0;
		}
		else {
			is_vertex = false;
			f = in_seed->val;
			e = -1;
		}

		oldgfv = -1;
		copy_v3_v3(oldco, in_seed->co);

		while (1) {
			if (is_vertex) {
				for (i = 0; i < sys->ringf_map[v].count; i++) {
					f = sys->ringf_map[v].indices[i];

					mul_v3_v3fl(gf, gfsys->gfield[f], dir);
					if (nextPoint(newco, &e, gfsys, f, oldco, gf)) {
						is_vertex = false;
						break;
					}
				}
				if (is_vertex) break;
			}
			else {
				mul_v3_v3fl(gf, gfsys->gfield[f], dir);
				if (!nextPoint(newco, &newe, gfsys, f, oldco, gf)) {
					is_vertex = true;
					v = compare_v3v3(newco, sys->co[sys->edges[e][0]], 0.000001f) ? sys->edges[e][0] : sys->edges[e][1];
					continue;
				}

				if (newe == e) {
					if (dir * sys->U_field[sys->edges[e][0]] < dir * sys->U_field[sys->edges[e][1]])
						v = sys->edges[e][0];
					else
						v = sys->edges[e][1];

					copy_v3_v3(newco, sys->co[v]);
					is_vertex = true;
				}
				else e = newe;

				if (compare_v3v3(newco, oldco, 0.000001f))
					break;
			}

			if (!checkPoint(gfsys, oldco, newco, f, 0.08f, 0.12f)) break;

			if (seed == -1) seed = addVertGFSystem(gfsys, in_seed->co);
			if (oldgfv == -1) oldgfv = seed;
			newgfv = addVertGFSystem(gfsys, newco);
			addEdgeGFSystem(gfsys, oldgfv, newgfv, f);
			oldgfv = newgfv;
			copy_v3_v3(oldco, gfsys->mesh->mvert[oldgfv].co);

			f = getOtherFaceAdjacentToEdge(sys, f, e);
		}
		dir = -dir;
	}
}

void computeFlowLines(LaplacianSystem *sys) {
	GFSeed *seed;
	int comp = 0;
	
	if (sys->gfsys1) deleteGradientFlowSystem(sys->gfsys1);
	if (sys->gfsys2) deleteGradientFlowSystem(sys->gfsys2);
	sys->gfsys1 = newGradientFlowSystem(sys, sys->h1, sys->gf1);
	sys->gfsys2 = newGradientFlowSystem(sys, sys->h2, sys->gf2);

	while (!BLI_heap_is_empty(sys->gfsys1->heap_seeds)) {
		seed = getTopSeedFromQueue(sys->gfsys1->heap_seeds);
		if (++comp < QR_LINELIMIT)
			computeGFLine(sys->gfsys1, seed);
		MEM_SAFE_FREE(seed);
	}

	comp = 0;
	
	while (!BLI_heap_is_empty(sys->gfsys2->heap_seeds)) {
		seed = getTopSeedFromQueue(sys->gfsys2->heap_seeds);
		if (++comp < QR_LINELIMIT)
			computeGFLine(sys->gfsys2, seed);
		MEM_SAFE_FREE(seed);
	}
	
}