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

GFVertID addVertGFMesh(GradientFlowMesh *gfmesh, float co[3], int index_edge)
{
	if (gfmesh->totvert == gfmesh->allocvert) {
		gfmesh->allocvert += MOD_QUADREMESH_ALLOC_BLOCK;
		gfmesh->mvert = MEM_reallocN(gfmesh->mvert, sizeof(GFVert) * gfmesh->allocvert);
	}

	copy_v3_v3(gfmesh->mvert[gfmesh->totvert].co, co);
	gfmesh->mvert[gfmesh->totvert].idev = index_edge;
	gfmesh->totvert++;

	return gfmesh->totvert - 1;
}

GFVertID addGFVertGFMesh(GradientFlowMesh *gfmesh, GFVert gfvert)
{
	return addVertGFMesh(gfmesh, gfvert.co, gfvert.idev);
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

	/*for (i = 0; i < sys->total_edges; i++) {
		f1 = sys->faces_edge[i][0];
		f2 = sys->faces_edge[i][1];
		angle = angle_normalized_v3v3(sys->no[f1], sys->no[f2]);
		if (angle >= M_PI_2) {
			listverts[sys->edges[i][0]] = 1;
			listverts[sys->edges[i][1]] = 1;
		}
	}*/

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

void addSeedToQueue(Heap *aheap, float value, GFVert *vert)
{
	BLI_heap_insert(aheap, value, vert);
}

GradientFlowSystem *newGradientFlowSystem(LaplacianSystem *sys, float *mhfunction, float(*mgfield)[3])
{
	int ve[2], i;
	int *lverts, sizeverts[2];	
	GFVert *mv;
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
		mv = MEM_mallocN(sizeof(GFVert), __func__);
		copy_v3_v3(mv->co, sys->co[lverts[i]]);
		mv->idev = -lverts[i] - 1;
		addSeedToQueue(gfsys->heap_seeds, 0.0f, mv);
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

GFVertID addVertGFSystem(GradientFlowSystem *gfsys, float co[3], int index_edge)
{
	GFVertID ret = addVertGFMesh(gfsys->mesh, co, index_edge);

	return ret;
}

GFVertID addGFVertGFSystem(GradientFlowSystem *gfsys, GFVert gfvert)
{
	return addVertGFSystem(gfsys, gfvert.co, gfvert.idev);
}

GFEdgeID addEdgeGFSystem(GradientFlowSystem *gfsys, int index_v1, int index_v2, int index)
{
	GFEdgeID ret = addEdgeGFMesh(gfsys->mesh, index_v1, index_v2, index);

	if (index >= 0)
		BLI_linklist_prepend(&gfsys->ringf_list[index], (void*)ret);
	else
		BLI_linklist_prepend(&gfsys->ringe_list[-index - 1], (void*)ret);

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


GFVert *getTopSeedFromQueue(struct Heap *aheap)
{
	GFVert *vert = BLI_heap_popmin(aheap);

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

#if 0
/*
* rind[0] return >0 if solution found
* rind[1] return index or face or edge
* r[3] direction 
*/
void computeGradientDirectionOnVert(int rind[2], float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexvert)
{
	int i, numv, *indv, ind_maxu, op1, op2, index_edge, fs[2];
	float maxu;
	numv = sys->ringv_map[indexvert].count;
	indv = sys->ringv_map[indexvert].indices;
	ind_maxu = indexvert;
	maxu = sys->U_field[indexvert];
	for (i = 0; i < numv; i++) {
		if (sys->U_field[indv[i]] > maxu) {
			ind_maxu = indv[i];
			maxu = sys->U_field[indv[i]];
		}
	}

	if (ind_maxu == indexvert) {
		rind[0] = GRA_DIR_ON_NONE;
		rind[1] = -1;
		zero_v3(r);
		return;
	}

	index_edge = getEdgeFromVerts(sys, indexvert, ind_maxu);
	copy_v2_v2_int(fs, sys->faces_edge[index_edge]);
	op1 = getThirdVert(sys, fs[0], index_edge);
	op2 = getThirdVert(sys, fs[1], index_edge);

	if (sys->U_field[op1] <= sys->U_field[indexvert] && sys->U_field[op2] <= sys->U_field[indexvert]) {
		rind[0] = GRA_DIR_ON_EDGE;
		rind[1] = index_edge;
		sub_v3_v3v3(r, sys->co[ind_maxu], sys->co[indexvert]);
		return;
	}

	if (sys->U_field[op1] > sys->U_field[op2]) {
		rind[0] = GRA_DIR_ON_FACE;
		rind[1] = fs[0];
		copy_v3_v3(r, gfsys->gfield[fs[0]]);
		return;
	}

	if (sys->U_field[op2] > sys->U_field[op1]) {
		rind[0] = GRA_DIR_ON_FACE;
		rind[1] = fs[1];
		copy_v3_v3(r, gfsys->gfield[fs[1]]);
		return;
	}
	rind[0] = GRA_DIR_ON_EDGE;
	rind[1] = index_edge;
	sub_v3_v3v3(r, sys->co[ind_maxu], sys->co[indexvert]);
	return;

}

int getGradientDirectionType(int u1, int u2, int u3, int u4)
{
	if (u1 >= u2 && u1 >= u3 && u1 >= u4) {
		if (u2 >= u3 && u2 >= u4) {
			return GRA_DIR_ON_EDGE;
		}
		if (u3 >= u4) {
			return GRA_DIR_ON_FACE;
		}
		if (u4 >= u3) {
			return GRA_DIR_ON_FACE;
		}
	}
	if (u2 >= u1 && u2 >= u3 && u2 >= u4) {
		if (u1 >= u3 && u1 >= u4) {
			return GRA_DIR_ON_EDGE;
		}
		if (u3 >= u4) {
			return GRA_DIR_ON_FACE;
		}
		if (u4 >= u3) {
			return GRA_DIR_ON_FACE;
		}
	}
	if (u3 >= u4) {
		return GRA_DIR_ON_FACE;
	}
	if (u4 >= u3) {
		return GRA_DIR_ON_FACE;
	}

	return GRA_DIR_ON_NONE;
}

void computeGradientDirectionOnEdge(int rtype, int rind, float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexedge)
{
	if (u1 >= u2 && u1 >= u3 && u1 >= u4) {
		if (u2 >= u3 && u2 >= u4) {
			rtype = GRA_DIR_ON_EDGE;
			rind = indexedge;
			sub_v3_v3v3(r, sys->co[iv2], sys->co[iv1]);
			return;
		}
		if (u3 >= u4) {
			rtype = GRA_DIR_ON_FACE;
			rind = fs[0];
			copy_v3_v3(r, gfsys->gfield[fs[0]]);
			return;
		}
		if (u4 >= u3) {
			rtype = GRA_DIR_ON_FACE;
			rind = fs[1];
			copy_v3_v3(r, gfsys->gfield[fs[1]]);
			return;
		}
	}
	if (u2 >= u1 && u2 >= u3 && u2 >= u4) {
		if (u1 >= u3 && u1 >= u4) {
			rtype = GRA_DIR_ON_EDGE;
			rind = indexedge;
			sub_v3_v3v3(r, sys->co[iv1], sys->co[iv2]);
			return;
		}
		if (u3 >= u4) {
			rtype = GRA_DIR_ON_FACE;
			rind = fs[0];
			copy_v3_v3(r, gfsys->gfield[fs[0]]);
			return;
		}
		if (u4 >= u3) {
			rtype = GRA_DIR_ON_FACE;
			rind = fs[1];
			copy_v3_v3(r, gfsys->gfield[fs[1]]);
			return;
		}
	}
	if (u3 >= u4) {
		rtype = GRA_DIR_ON_FACE;
		rind = fs[0];
		copy_v3_v3(r, gfsys->gfield[fs[0]]);
		return;
	}
	if (u4 >= u3) {
		rtype = GRA_DIR_ON_FACE;
		rind = fs[1];
		copy_v3_v3(r, gfsys->gfield[fs[1]]);
		return;
	}

	rtype = GRA_DIR_ON_NONE;
}

void computeGradientDirectionOnEdgeInverse(int rind[2], float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexedge)
{
	int fs[2], iv1, iv2, ov1, ov2;
	float u1, u2, u3, u4;
	copy_v2_v2_int(fs, sys->faces_edge[indexedge]);
	iv1 = sys->edges[indexedge][0];
	iv2 = sys->edges[indexedge][1];
	ov1 = getThirdVert(sys, fs[0], indexedge);
	ov2 = getThirdVert(sys, fs[1], indexedge);
	u1 = sys->U_field[iv1];
	u2 = sys->U_field[iv2];
	u3 = sys->U_field[ov1];
	u4 = sys->U_field[ov2];

	if (u1 <= u2 && u1 <= u3 && u1 <= u4) {
		if (u2 <= u3 && u2 <= u4) {
			rind[0] = GRA_DIR_ON_EDGE;
			rind[1] = indexedge;
			sub_v3_v3v3(r, sys->co[iv2], sys->co[iv1]);
			return;
		}
		if (u3 <= u4) {
			rind[0] = GRA_DIR_ON_FACE;
			rind[1] = fs[0];
			copy_v3_v3(r, gfsys->gfield[fs[0]]);
			return;
		}
		if (u4 <= u3) {
			rind[0] = GRA_DIR_ON_FACE;
			rind[1] = fs[1];
			copy_v3_v3(r, gfsys->gfield[fs[1]]);
			return;
		}
	}
	if (u2 <= u1 && u2 <= u3 && u2 <= u4) {
		if (u1 <= u3 && u1 <= u4) {
			rind[0] = GRA_DIR_ON_EDGE;
			rind[1] = indexedge;
			sub_v3_v3v3(r, sys->co[iv1], sys->co[iv2]);
			return;
		}
		if (u3 <= u4) {
			rind[0] = GRA_DIR_ON_FACE;
			rind[1] = fs[0];
			copy_v3_v3(r, gfsys->gfield[fs[0]]);
			return;
		}
		if (u4 <= u3) {
			rind[0] = GRA_DIR_ON_FACE;
			rind[1] = fs[1];
			copy_v3_v3(r, gfsys->gfield[fs[1]]);
			return;
		}
	}
	if (u3 <= u4) {
		rind[0] = GRA_DIR_ON_FACE;
		rind[1] = fs[0];
		copy_v3_v3(r, gfsys->gfield[fs[0]]);
		return;
	}
	if (u4 <= u3) {
		rind[0] = GRA_DIR_ON_FACE;
		rind[1] = fs[1];
		copy_v3_v3(r, gfsys->gfield[fs[1]]);
		return;
	}

	rind[0] = GRA_DIR_ON_NONE;

}

/*
* ori coordinate of origin point
* dir direction to made query
* indexface Face in original mesh
* maxradius lenght to made query on dir direction
*/

int nearGFEdgeInGFMeshFromEdge(LaplacianSystem *sys, GradientFlowSystem *gfsys, float ori[3], float dir[3], int indexedge, float maxradius)
{
	int fs[2], res, v1, v2, v3;
	float d1[3], d2[3], r[3];

	if (is_zero_v3(dir)) {
		return -1;
	}

	copy_v2_v2_int(fs, sys->faces_edge[indexedge]);
	projectVectorOnFace(d1, sys->no[fs[0]], dir);
	mul_v3_v3fl(d2, d1, -1.0f);
	v1 = sys->edges[indexedge][0];
	v2 = sys->edges[indexedge][1];
	v3 = getThirdVert(sys, fs[0], indexedge);

	if (intersecionLineSegmentWithVector(r, sys->co[v1], sys->co[v3], ori, d1) 
		|| intersecionLineSegmentWithVector(r, sys->co[v2], sys->co[v3], ori, d1)) {
		res = nearGFEdgeInGFMesh(sys, gfsys, ori, d1, fs[0], maxradius);
		if (res >=0 ) {
			return res;
		}
	} 
	else if (intersecionLineSegmentWithVector(r, sys->co[v1], sys->co[v3], ori, d2) 
		|| intersecionLineSegmentWithVector(r, sys->co[v2], sys->co[v3], ori, d2)) {
		res = nearGFEdgeInGFMesh(sys, gfsys, ori, d2, fs[0], maxradius);
		if (res >=0 ) {
			return res;
		}
	}
	projectVectorOnFace(d1, sys->no[fs[1]], dir);
	mul_v3_v3fl(d2, d1, -1.0f);
	v3 = getThirdVert(sys, fs[1], indexedge);
	if (intersecionLineSegmentWithVector(r, sys->co[v1], sys->co[v3], ori, d1) 
		|| intersecionLineSegmentWithVector(r, sys->co[v2], sys->co[v3], ori, d1)) {
		res = nearGFEdgeInGFMesh(sys, gfsys, ori, d1, fs[1], maxradius);
		if (res >=0 ) {
			return res;
		}
	} 
	else if (intersecionLineSegmentWithVector(r, sys->co[v1], sys->co[v3], ori, d2) 
		|| intersecionLineSegmentWithVector(r, sys->co[v2], sys->co[v3], ori, d2)) {
		res = nearGFEdgeInGFMesh(sys, gfsys, ori, d2, fs[1], maxradius);
		if (res >=0 ) {
			return res;
		}
	}

	return -1;
}

bool nearGFEdgeInGFMesh(LaplacianSystem *sys, GradientFlowSystem *gfsys, float ori[3], float dir[3], int indexface, float maxradius)
{
	GFList *gfl;
	int index_gfedge, iv1, iv2, eid, newface;
	int res;
	float i1[3], v[3], r[3];
	add_v3_v3v3(v, ori, dir);
	/* Query on flow lines inside face[indexface]*/
	if (gfsys->ringf_list[indexface]) {
		gfl = gfsys->ringf_list[indexface];
		while (gfl) {
			index_gfedge = gfl->index;
			if (index_gfedge >= 0 && index_gfedge < gfsys->mesh->totedge) {
				iv1 = gfsys->mesh->medge[index_gfedge].v1;
				iv2 = gfsys->mesh->medge[index_gfedge].v2;


				if (intersecionLineSegmentWithVector(i1, gfsys->mesh->mvert[iv1].co, gfsys->mesh->mvert[iv2].co, ori, v)) {
					if (len_v3v3(i1, ori) < maxradius){
						return index_gfedge;
					}
				}
			}
			gfl = gfl->next;
		}
	}
	else {
		/*Do not flow lines found, then search on adjacent faces*/
		if (intersecionLineSegmentWithVector(i1, sys->co[sys->faces[indexface][0]], sys->co[sys->faces[indexface][1]],
			ori, v)) {
			res = len_v3v3(i1, ori);
			if (res > maxradius) {
				return -1;
			}
			else {
				eid = getEdgeFromVerts(sys, sys->faces[indexface][0], sys->faces[indexface][1]);
				newface = getOtherFaceAdjacentToEdge(sys, indexface, eid);
				//projectVectorOnFace(r, sys->no[indexface], dir);
				cross_v3_v3v3(r, sys->no[indexface], gfsys->gfield[indexface]);
				return nearGFEdgeInGFMesh(sys, gfsys, i1, r, newface, maxradius - res);
			}
		}
		else if (intersecionLineSegmentWithVector(i1, sys->co[sys->faces[indexface][1]], sys->co[sys->faces[indexface][2]],
			ori, v)) {
			res = len_v3v3(i1, ori);
			if (res > maxradius) {
				return -1;
			}
			else {
				eid = getEdgeFromVerts(sys, sys->faces[indexface][1], sys->faces[indexface][2]);
				newface = getOtherFaceAdjacentToEdge(sys, indexface, eid);
				//projectVectorOnFace(r, sys->no[indexface], dir);
				cross_v3_v3v3(r, sys->no[indexface], gfsys->gfield[indexface]);
				return nearGFEdgeInGFMesh(sys, gfsys, i1, r, newface, maxradius - res);
			}
		}
		else if (intersecionLineSegmentWithVector(i1, sys->co[sys->faces[indexface][2]], sys->co[sys->faces[indexface][0]],
			ori, v)) {
			res = len_v3v3(i1, ori);
			if (res > maxradius) {
				return -1;
			}
			else {
				eid = getEdgeFromVerts(sys, sys->faces[indexface][2], sys->faces[indexface][0]);
				newface = getOtherFaceAdjacentToEdge(sys, indexface, eid);
				//projectVectorOnFace(r, sys->no[indexface], dir);
				cross_v3_v3v3(r, sys->no[indexface], gfsys->gfield[indexface]);
				return nearGFEdgeInGFMesh(sys, gfsys, i1, r, newface, maxradius - res);
			}
		}
	}
	return -1;
}
#endif

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

bool nextPointFromEdge(float r_co[3], int *r_edge, GradientFlowSystem *gfsys, int in_f, int in_e, float in_co[3], float in_dir[3])
{
	int iv1, iv2, ov1;
	float a[3], b[3], c[3], 
		alfa, beta, gamma, dummy[3];
	LaplacianSystem *sys = gfsys->sys;

	iv1 = sys->edges[in_e][0];
	iv2 = sys->edges[in_e][1];

	ov1 = getThirdVert(sys, in_f, iv1, iv2);

	sub_v3_v3v3(a, sys->co[iv1], in_co);
	sub_v3_v3v3(b, sys->co[ov1], in_co);

	alfa = angle_v3v3(a, in_dir);
	beta = angle_v3v3(a, b);
	gamma = angle_v3v3(in_dir, b);

	if (beta > alfa) {
		if (gamma > beta) {
			/* vertex on iv1 */
			copy_v3_v3(r_co, sys->co[iv1]);
			*r_edge = -iv1 - 1;
		}
		else {
			/* vertex on edge iv1 - ov1 */
			add_v3_v3v3(c, in_dir, in_co);
			isect_line_line_v3(sys->co[iv1], sys->co[ov1], in_co, c, r_co, dummy);
			*r_edge = getEdgeFromVerts(sys, iv1, ov1);
		}
	}
	else {
		if (gamma > (float)M_PI - beta) {
			/* vertex on iv2 */
			copy_v3_v3(r_co, sys->co[iv2]);
			*r_edge = -iv2 - 1;
		}
		else {
			/* vertex on edge iv2 - ov1 */
			add_v3_v3v3(c, in_dir, in_co);
			isect_line_line_v3(sys->co[iv2], sys->co[ov1], in_co, c, r_co, dummy);
			*r_edge = getEdgeFromVerts(sys, iv2, ov1);
		}
	}

	return true;
}

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
int queryDirection(GradientFlowSystem *gfsys, float in_co[3], int in_e, int in_f, float in_dir[3], float maxdist)
{
	int i, e;
	float c[3], len, co[3], dir[3];
	LaplacianSystem *sys = gfsys->sys;
	copy_v3_v3(dir, in_dir);

	len = dot_v3v3(dir, sys->no[in_f]);
	mul_v3_v3fl(c, sys->no[in_f], len);
	sub_v3_v3(dir, c);
	normalize_v3(dir);

	nextPointFromEdge(co, &e, gfsys, in_f, in_e, in_co, dir);

	if (e < 0)
		return 2;

	sub_v3_v3v3(c, co, in_co);
	len = len_v3(c);

	if (maxdist < len) {
		mul_v3_fl(c, maxdist / len);
		add_v3_v3v3(co, c, in_co);
		maxdist = 0.0f;
	}
	else
		maxdist -= len;

	/*int vf1, vf2;
	vf1 = addVertGFSystem(gfsys, in_co, 0);
	vf2 = addVertGFSystem(gfsys, co, 0);
	addEdgeGFSystem(gfsys, vf1, vf2, 0);*/

	if (intersectSegmentWithOthersOnFace(gfsys, in_co, co, in_f))
		return 0;
		
	if (maxdist > 0.0f) {
		in_f = getOtherFaceAdjacentToEdge(sys, in_f, e);
		i = queryDirection(gfsys, co, e, in_f, dir, maxdist);
		if (i == 2) return 1;
		else return i;
	}
	else
		return 1;
}

void getOrthogonalDirection(float r_dir[3], GradientFlowSystem *gfsys, float in_dir[3], int in_ev)
{
	int i, f;
	float no[3];
	LaplacianSystem *sys = gfsys->sys;

	zero_v3(no);

	if (in_ev < 0) {
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
	}

	cross_v3_v3v3(r_dir, no, in_dir);
	normalize_v3(r_dir);
}

bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], float in_newco[3], int in_ev, float maxdist)
{
	int i, f, d, vf1, vf2;
	float seg[3], dir[3];
	LaplacianSystem *sys = gfsys->sys;

	sub_v3_v3v3(seg, in_oldco, in_newco);
	getOrthogonalDirection(dir, gfsys, seg, in_ev);

	for (d = 0; d < 2; d++) {
		if (in_ev < 0) {
			for (i = 0; i < sys->ringf_map[-in_ev - 1].count; i++) {
				f = sys->ringf_map[-in_ev - 1].indices[i];
				if (!getSecondAndThirdVert(&vf1, &vf2, sys, f, -in_ev - 1)) continue;
#if 0
					float co[3], a[3], b[3], alfa, beta, gamma;
					len = dot_v3v3(dir, sys->no[f]);
					mul_v3_v3fl(b, sys->no[f], len);
					sub_v3_v3v3(c, dir, b);
					normalize_v3(c);

					sub_v3_v3v3(a, sys->co[in_ev], sys->co[vf1]);
					sub_v3_v3v3(b, sys->co[in_ev], sys->co[vf2]);

					alfa = angle_v3v3(a, b);
					beta = angle_v3v3(a, dir);
					gamma = angle_v3v3(dir, b);

					if (beta + gamma >= alfa - 0.001f && beta + gamma <= alfa + 0.001f) {
						if (!queryDirection(gfsys, in_newco, getEdgeFromVerts(sys, in_ev, vf2), f, c, maxdist))
							return false;
						break;
					}
#endif // 0
				if (!queryDirection(gfsys, in_newco, getEdgeFromVerts(sys, -in_ev - 1, vf2), f, dir, maxdist))
					return false;
			}
		}
		else {
			if (!queryDirection(gfsys, in_newco, in_ev, sys->faces_edge[in_ev][0], dir, maxdist))
				return false;
			if (!queryDirection(gfsys, in_newco, in_ev, sys->faces_edge[in_ev][1], dir, maxdist))
				return false;
		}
		mul_v3_fl(dir, -1.0f);
	}
	return true;
}

bool nextPointFlowLine(GFVertID *r, GradientFlowSystem *gfsys, GFVertID in_id, int *face)
{
	float co[3];
	int edge;
	LaplacianSystem *sys = gfsys->sys;
	GFVert *in = &gfsys->mesh->mvert[in_id];

	if (in->idev < 0) {
		if (!nextPointFromVertex(co, face, &edge, gfsys, -in->idev - 1)) return false;
		if (!checkPoint(gfsys, in->co, co, edge, 0.1f)) return false;

		*r = addVertGFSystem(gfsys, co, edge);
		addEdgeGFSystem(gfsys, in_id, *r, *face);
	}
	else {
		*face = getOtherFaceAdjacentToEdge(sys, *face, in->idev);
		nextPointFromEdge(co, &edge, gfsys, *face, in->idev, in->co, gfsys->gfield[*face]);

		if (!checkPoint(gfsys, in->co, co, edge, 0.1f)) return false;

		*r = addVertGFSystem(gfsys, co, edge);
		addEdgeGFSystem(gfsys, in_id, *r, edge < 0 ? edge : *face);
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

void computeGFLine(LaplacianSystem *sys, GradientFlowSystem *gfsys, GFVert *gfvert_seed)
{
	GFVertID oldgfv, newgfv;
	int face = 0;

	oldgfv = addVertGFSystem(gfsys, gfvert_seed->co, gfvert_seed->idev);

	while (nextPointFlowLine(&newgfv, gfsys, oldgfv, &face)) {
		oldgfv = newgfv;
	}
}

int computeNewSeed(float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexf, float ori[3], float dir[3], float mh)
{
	float v1[3], v2[3], v3[3], res[3], d1[3], d2[3], lenr;
	int fs[2], indf1, i1, i2, i3;
	bool intersection;

	if (mh <= 0) {
		return -1;
	}

	i1 = sys->faces[indexf][0];
	i2 = sys->faces[indexf][1];
	i3 = sys->faces[indexf][2];

	copy_v3_v3(v1, sys->co[i1]);
	copy_v3_v3(v2, sys->co[i2]);
	copy_v3_v3(v3, sys->co[i3]);
	normalize_v3_v3(d1, dir);

	//printf("computeNewSeed 0\n");
	intersection = intersecionLineSegmentWithVector(res, v1, v2, ori, d1);
	lenr = intersection ? len_v3v3(res, ori) : 0.0f;
	if (intersection && lenr > MOD_QUADREMESH_MIN_LEN) {
		//printf("computeNewSeed 1\n");
		lenr = len_v3v3(res, ori);
		//printf("computeNewSeed 2\n");
		if (lenr >= mh) {
			copy_v3_v3(r, res);
			//printf("computeNewSeed 3\n");
			return getEdgeFromVerts(sys, i1, i2);
		} 
		else {
			//printf("computeNewSeed 4\n");
			copy_v2_v2_int(fs, sys->faces_edge[getEdgeFromVerts(sys, i1, i2)]);
			indf1 = fs[0] == indexf ? fs[1] : fs[0];
			projectVectorOnFace(d2, sys->no[indf1], d1);
			//printf("computeNewSeed len=%f \t5\n, lenr");
			return computeNewSeed(r, sys, gfsys, indf1, res, d2, mh - lenr);
		}
	} 
	else {
		intersection = intersecionLineSegmentWithVector(res, v2, v3, ori, d1);
		lenr = intersection ? len_v3v3(res, ori) : 0.0f;
		if (intersection && lenr > MOD_QUADREMESH_MIN_LEN) {
			//printf("computeNewSeed 6\n");
			lenr = len_v3v3(res, ori);
			//printf("computeNewSeed 7\n");
			if (lenr >= mh) {
				copy_v3_v3(r, res);
				//printf("computeNewSeed 8\n");
				return getEdgeFromVerts(sys, i2, i3);
			}
			else {
				//printf("computeNewSeed 9\n");
				copy_v2_v2_int(fs, sys->faces_edge[getEdgeFromVerts(sys, i2, i3)]);
				indf1 = fs[0] == indexf ? fs[1] : fs[0];
				projectVectorOnFace(d2, sys->no[indf1], d1);
				//printf("computeNewSeed len=%f \t10\n, lenr");
				return computeNewSeed(r, sys, gfsys, indf1, res, d2, mh - lenr);
			}
		}
		else  {
			intersection = intersecionLineSegmentWithVector(res, v3, v1, ori, d1);
			lenr = intersection ? len_v3v3(res, ori) : 0.0f;
			if (intersection && lenr > MOD_QUADREMESH_MIN_LEN) {
				//printf("computeNewSeed 10\n");
				lenr = len_v3v3(res, ori);
				//printf("computeNewSeed 11\n");
				if (lenr >= mh) {
					//printf("computeNewSeed 12\n");
					copy_v3_v3(r, res);
					return getEdgeFromVerts(sys, i3, i1);
				}
				else {
					//printf("computeNewSeed 13\n");
					copy_v2_v2_int(fs, sys->faces_edge[getEdgeFromVerts(sys, i3, i1)]);
					indf1 = fs[0] == indexf ? fs[1] : fs[0];
					projectVectorOnFace(d2, sys->no[indf1], d1);
					//printf("computeNewSeed 14\n");
					//printf("computeNewSeed len=%f \t15\n, lenr");
					return computeNewSeed(r, sys, gfsys, indf1, res, d2, mh - lenr);
				}
			}
		}
	}
	//printf("computeNewSeed 15\n");
	return -1;
}

void computeFlowLines(LaplacianSystem *sys) {
	GFVert *seed;
	
	if (sys->gfsys1) deleteGradientFlowSystem(sys->gfsys1);
	sys->gfsys1 = newGradientFlowSystem(sys, sys->h1, sys->gf1);
	//sys->gfsys2 = newGradientFlowSystem(sys, sys->h2, sys->gf2);

	while (!BLI_heap_is_empty(sys->gfsys1->heap_seeds)) {
		seed = getTopSeedFromQueue(sys->gfsys1->heap_seeds);
		computeGFLine(sys, sys->gfsys1, seed);
		MEM_SAFE_FREE(seed);
	}
	
	//while (!BLI_heap_is_empty(gfsys2->heap_seeds)) {
	//	seed = getTopSeedFromQueue(gfsys2->heap_seeds);
	//	computeGFLine(sys, gfsys2, seed);
	//	//MEM_SAFE_FREE(seed);
	//}
	
}