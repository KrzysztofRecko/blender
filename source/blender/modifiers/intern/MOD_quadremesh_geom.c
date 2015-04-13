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
#define QR_SAMPLING_RATE 0.03f
#define QR_MINDIST 0.04f
#define QR_SEEDDIST 0.08f

#if 0
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

void deleteGradientFlowMesh(GradientFlowMesh * gfmesh) 
{
	if (gfmesh) {
		MEM_SAFE_FREE(gfmesh->mvert);
		MEM_SAFE_FREE(gfmesh->medge);
		MEM_SAFE_FREE(gfmesh);
	}
}

#endif // 0

GFVertID addVert(LaplacianSystem *sys, float co[3])
{
	if (sys->totvert == sys->allocvert) {
		sys->allocvert = sys->allocvert * 2 + 10;
		sys->mvert = MEM_reallocN(sys->mvert, sizeof(GFVert) * sys->allocvert);
	}

	copy_v3_v3(sys->mvert[sys->totvert].co, co);
	sys->mvert[sys->totvert].e[0] = -1;
	sys->mvert[sys->totvert].e[1] = -1;
	sys->mvert[sys->totvert].e[2] = -1;
	sys->mvert[sys->totvert].e[3] = -1;
	sys->totvert++;

	return sys->totvert - 1;
}

GFEdgeID addEdge(GradientFlowSystem *gfsys, int index_v1, int index_v2, int index_face)
{
	if (gfsys->totedge == gfsys->allocedge) {
		gfsys->allocedge += MOD_QUADREMESH_ALLOC_BLOCK;
		gfsys->medge = MEM_reallocN(gfsys->medge, sizeof(GFEdge) * gfsys->allocedge);
	}

	gfsys->medge[gfsys->totedge].v1 = index_v1;
	gfsys->medge[gfsys->totedge].v2 = index_v2;
	gfsys->totedge++;

	BLI_linklist_prepend(&gfsys->ringf_list[index_face], (void*)(gfsys->totedge - 1));

	return gfsys->totedge - 1;
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

	//estimateNumberGFVerticesEdges(ve, sys, sys->h);
	
	gfsys->sys = sys;
	//gfsys->mesh = newGradientFlowMesh(ve[0], ve[1]);
	
	gfsys->ringf_list = MEM_callocN(sizeof(LinkNode *) * sys->total_faces, "GFListFaces");
	
	gfsys->heap_seeds = BLI_heap_new();
	
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
		for (i = 0; i < gfsys->sys->total_faces; i++) {
			BLI_linklist_free(gfsys->ringf_list[i], NULL);
		}
		MEM_SAFE_FREE(gfsys->ringf_list);
		BLI_heap_free(gfsys->heap_seeds, MEM_freeN);
		MEM_SAFE_FREE(gfsys->medge);
		MEM_SAFE_FREE(gfsys);
	}
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

#if 0
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
#endif

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

#if 0
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
#endif

/**
 * /return 0 - all good
 *         1 - new vertex very close to 1st vertex of result edge
 *         2 -           -||-           2nd vertex of result edge
 *         3 - other error (eg. input outside of triangle)
 */
int nextPoint(float r_co[3], int *r_edge, GradientFlowSystem *gfsys, int in_f, float in_co[3], float in_dir[3])
{
	int i, pick = -1, v = -1;
	bool is_on_vertex = false;
	float a[3][3], b[3][3], c[2][3], co2[3],
		dummy[3];
	LaplacianSystem *sys = gfsys->sys;

	/* check if direction is coplanar to triangle */
	/* check if point is inside triangle */
	/* check if triangle is degenerate */

	add_v3_v3v3(co2, in_dir, in_co); /* second point on direction */

	for (i = 0; i < 3; i++) {
		sub_v3_v3v3(a[i], in_co, sys->co[sys->faces[in_f][i]]);
		if (dot_v3v3(a[i], a[i]) < FLT_EPSILON)
			v = i;
	}

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
		return 3;

	*r_edge = getEdgeFromVerts(sys, sys->faces[in_f][pick], sys->faces[in_f][(pick + 1) % 3]);

	isect_line_line_v3(sys->co[sys->faces[in_f][pick]], sys->co[sys->faces[in_f][(pick + 1) % 3]], in_co, co2, r_co, dummy);

	if (len_squared_v3v3(r_co, sys->co[sys->edges[*r_edge][0]]) < FLT_EPSILON)
		return 1;
	if (len_squared_v3v3(r_co, sys->co[sys->edges[*r_edge][1]]) < FLT_EPSILON)
		return 2;

	return 0;
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
	GFVert *vertices = sys->mvert;
	GFEdge *edges = gfsys->medge;

	for (iter = gfsys->ringf_list[in_f]; iter; iter = iter->next) {
		e = (int)iter->link;
		/* TODO: Maybe grow the segment a little to pick up more stuff? */
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

		if (nextPoint(newco, &e, gfsys, in_f, oldco, dir)) return 2;
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
		
		if (actlen > dist && !make_seed)
			return 1;

#ifdef QR_SHOWQUERIES
		if (actlen > dist && actlen < maxdist) {
			int vf1, vf2;
			vf1 = addVertGFSystem(gfsys, oldco);
			vf2 = addVertGFSystem(gfsys, newco);
			addEdgeGFSystem(gfsys, vf1, vf2, 0);
		}
#endif

		if (actlen > maxdist) {
			mul_v3_v3fl(newco2, c, (maxdist - actlen + len) / len);
			add_v3_v3(newco2, oldco);

#ifdef QR_SHOWQUERIES
			int vf1, vf2;
			vf1 = addVertGFSystem(gfsys, oldco);
			vf2 = addVertGFSystem(gfsys, newco2);
			addEdgeGFSystem(gfsys, vf1, vf2, 0);
#endif

			addSeedToQueue(gfsys->heap_seeds, newco2, false, in_f, -maxdist);
			return 1;
		}

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
	if (dot_v3v3(seg, seg) < FLT_EPSILON) return true;

	getOrthogonalDirection(dir, gfsys, seg, in_f);

	for (d = 0; d < 2; d++) {
		if (!queryDirection(gfsys, in_newco, in_f, dir, dist, maxdist, make_seed))
			return false;
		mul_v3_fl(dir, -1.0f);
	}

	return true;
}

float getSamplingDistanceFunctionOnFace(GradientFlowSystem *gfsys, int in_f, float in_co[3])
{
	float h1, h2, h3, uv[2];
	LaplacianSystem *sys = gfsys->sys;

	resolve_tri_uv_v3(uv, in_co, sys->co[sys->faces[in_f][0]], sys->co[sys->faces[in_f][1]], sys->co[sys->faces[in_f][2]]);

	h1 = gfsys->hfunction[sys->faces[in_f][0]];
	h2 = gfsys->hfunction[sys->faces[in_f][1]];
	h3 = gfsys->hfunction[sys->faces[in_f][2]];

	return uv[0] * h1 + uv[1] * h2 + (1.0f - uv[0] - uv[1]) * h3;
}

void addSegmentToLine(GradientFlowSystem *gfsys, GFLine *line, int in_f, float in_co[3])
{
	int i = gfsys == gfsys->sys->gfsys1 ? 0 : 1;
	GFVertID newv;

	if (line->d == 1) i += 2;

	newv = addVert(gfsys->sys, in_co);

	if (line->d == 0) addEdge(gfsys, line->end, newv, in_f);
	else addEdge(gfsys, newv, line->end, in_f);

	gfsys->sys->mvert[newv].e[i] = line->end;
	gfsys->sys->mvert[line->end].e[(i + 2) % 4] = newv;

	line->end = newv;
}

bool changeLineDirection(GradientFlowSystem *gfsys, GFLine *line)
{
	int i;

	if (line->seed != -1) {
		/* flush queue */
		for (i = 0; i < line->num_q; i++)
			addSegmentToLine(gfsys, line, line->qf[i], line->qco[i]);

		line->num_q = 0;
		copy_v3_v3(line->oldco, gfsys->sys->mvert[line->seed].co);
	}

	/* reset to original state */
	line->end = line->seed;
	copy_v3_v3(line->lastchk, line->oldco);

	return ++line->d != 2;
}

void initGFLine(GradientFlowSystem *gfsys, GFLine *line, GFSeed *in_seed)
{
	line->seed = line->end = -1;
	line->d = 0;
	
	copy_v3_v3(line->oldco, in_seed->co);
}

bool addPointToLine(GradientFlowSystem *gfsys, GFLine *line, int in_f, float in_newco[3])
{
	const float chklen = QR_SAMPLING_RATE; /* sampling rate */

	int i;
	float seg[3], newchk[3];
	float curlen;

	/* qco[0] - first point after last checked
	 * qco[num_q - 1] - last added point */

	sub_v3_v3v3(seg, in_newco, line->oldco);
	curlen = len_v3(seg);

	if (line->seed == -1) {
		if (!checkPoint(gfsys, in_newco, line->oldco, in_f, QR_MINDIST, QR_SEEDDIST)) return false;

		line->end = line->seed = addVert(gfsys->sys, line->oldco);
		line->lastchklen = line->qlen = 0.0f;
		line->num_q = 0;
		copy_v3_v3(line->lastchk, line->oldco);
	}

	while (line->qlen + curlen > line->lastchklen + chklen) {
		mul_v3_v3fl(newchk, seg, (line->lastchklen + chklen - line->qlen) / curlen);
		add_v3_v3(newchk, line->oldco);

		if (!checkPoint(gfsys, line->lastchk, newchk, in_f, QR_MINDIST, QR_SEEDDIST)) {
			if (line->num_q == 0) addSegmentToLine(gfsys, line, in_f, line->lastchk);
			else addSegmentToLine(gfsys, line, line->qf[0], line->lastchk);
			
			line->num_q = 0;
			return false;
		}
		
		/* flush queue */
		for (i = 0; i < line->num_q; i++)
			addSegmentToLine(gfsys, line, line->qf[i], line->qco[i]);

		line->num_q = 0;
		copy_v3_v3(line->lastchk, newchk);
		line->lastchklen += chklen;
	}

	if (line->num_q == 10) return false;

	copy_v3_v3(line->oldco, in_newco);
	copy_v3_v3(line->qco[line->num_q], in_newco);
	line->qf[line->num_q] = in_f;
	line->qlen += curlen;
	line->num_q++;

	return true;
}

void computeGFLine(GradientFlowSystem *gfsys, GFSeed *in_seed)
{
	int i, f, e, v, r, newe;
	bool is_vertex;
	float newco[3], gf[3], dir = 1.0f;
	LaplacianSystem *sys = gfsys->sys;
	GFLine line;

	initGFLine(gfsys, &line, in_seed);

	do {
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

		while (1) {
			if (is_vertex) {
				for (i = 0; i < sys->ringf_map[v].count; i++) {
					f = sys->ringf_map[v].indices[i];

					mul_v3_v3fl(gf, gfsys->gfield[f], dir);
					if (!nextPoint(newco, &e, gfsys, f, line.oldco, gf)) {
						is_vertex = false;
						break;
					}
				}
				if (is_vertex) break;
			}
			else {
				mul_v3_v3fl(gf, gfsys->gfield[f], dir);
				r = nextPoint(newco, &newe, gfsys, f, line.oldco, gf);
				
				if (r == 1) {
					is_vertex = true;
					v = sys->edges[newe][0];
				}
				else if (r == 2) {
					is_vertex = true;
					v = sys->edges[newe][1];
				}
				else if (r == 3) break;

				if (newe == e) {
					if (dir * sys->U_field[sys->edges[e][0]] < dir * sys->U_field[sys->edges[e][1]])
						v = sys->edges[e][0];
					else
						v = sys->edges[e][1];

					copy_v3_v3(newco, sys->co[v]);
					is_vertex = true;
				}

				e = newe;
			}
			
			if (!addPointToLine(gfsys, &line, f, newco))	break;

			f = getOtherFaceAdjacentToEdge(sys, f, e);
		} /* while (1) */

		dir = -dir;
	} while (changeLineDirection(gfsys, &line));
}

void computeFlowLines(LaplacianSystem *sys) {
	GFSeed *seed;
	int comp = 0;
	
	if (sys->gfsys1) deleteGradientFlowSystem(sys->gfsys1);
	if (sys->gfsys2) deleteGradientFlowSystem(sys->gfsys2);
	MEM_SAFE_FREE(sys->mvert);
	sys->allocvert = 0;
	sys->totvert = 0;
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

void insertOnEdge(GradientFlowSystem *gfsys, GFVertID in_vid, GFEdgeID in_e)
{
	int i = gfsys == gfsys->sys->gfsys1 ? 0 : 1;
	GFVert *verts = gfsys->sys->mvert;
	GFVertID va = gfsys->medge[in_e].v1;
	GFVertID vb = gfsys->medge[in_e].v2;
	GFVertID vc = verts[va].e[i + 2];

	while (vc != vb) {
		if (IS_EQF(len_v3v3(verts[va].co, verts[in_vid].co) + len_v3v3(verts[vc].co, verts[in_vid].co), 
			       len_v3v3(verts[va].co, verts[vc].co))) break;

		va = vc;
		vc = verts[vc].e[i + 2];
	}

	verts[va].e[i + 2] = in_vid;
	verts[vc].e[i] = in_vid;
	verts[in_vid].e[i + 2] = vc;
	verts[in_vid].e[i] = va;
}

void generateIntersectionsOnFaces(LaplacianSystem *sys)
{
	int f;
	float isection[3], lambda;
	GFEdgeID e1, e2;
	LinkNode *iter1, *iter2;
	GFVertID newv;

	for (f = 0; f < sys->total_faces; f++) {
		for (iter1 = sys->gfsys1->ringf_list[f]; iter1; iter1 = iter1->next) {
			e1 = (GFEdgeID)iter1->link;

			for (iter2 = sys->gfsys2->ringf_list[f]; iter2; iter2 = iter2->next) {
				e2 = (GFEdgeID)iter2->link;

				if (isect_line_line_strict_v3(sys->mvert[sys->gfsys1->medge[e1].v1].co,
					                          sys->mvert[sys->gfsys1->medge[e1].v2].co,
					                          sys->mvert[sys->gfsys2->medge[e2].v1].co,
					                          sys->mvert[sys->gfsys2->medge[e2].v2].co,
					                          isection, &lambda))
				{
					newv = addVert(sys, isection);
					insertOnEdge(sys->gfsys1, newv, e1);
					insertOnEdge(sys->gfsys2, newv, e2);
				}
			}
		}
	}
}

void deleteDegenerateVerts(LaplacianSystem *sys)
{
	int i, n, m;
	GFVert *verts = sys->mvert;

	for (i = 0; i < sys->totvert; i++) {
		for (n = 0, m = 0; n < 4; n++)
			if (verts[i].e[n] >= 0) m++;

		if (m == 1) {
			for (n = 0; n < 4; n++)
				if (verts[i].e[n] >= 0) break;

			verts[verts[i].e[n]].e[(n + 2) % 4] = -1;
			verts[i].e[n] = -1;
			verts[i].e[0] = -2;
		}
		else if (m == 2) {
			if (verts[i].e[0] >= 0 && verts[i].e[2] >= 0) n = 0;
			else if (verts[i].e[1] >= 0 && verts[i].e[3] >= 0) n = 1;
			else continue;

			verts[verts[i].e[n]].e[n + 2] = verts[i].e[n + 2];
			verts[verts[i].e[n + 2]].e[n] = verts[i].e[n];
			verts[i].e[n] = -1;
			verts[i].e[n + 2] = -1;
			verts[i].e[0] = -2;
		}
	}

	/*for (i = 0, n = 0; i < sys->totvert; i++) {
		if (verts[i].e[0] == -2) continue;

		copy_v3_v3(verts[n].co, verts[i].co);
		copy_v4_v4_int(verts[n].e, verts[i].e);
		n++;
	}

	sys->totvert = n;*/
}

void makeEdges(LaplacianSystem *sys)
{
	int i, j, k;

	sys->totedge = sys->allocedge = 0;
	MEM_SAFE_FREE(sys->medge);

	for (i = 0; i < sys->totvert; i++) {
		for (j = 0; j < 4; j++) {
			if (sys->mvert[i].e[j] >= 0) {
				if (sys->totedge == sys->allocedge) {
					sys->allocedge = sys->allocedge * 2 + 10;
					sys->medge = MEM_reallocN(sys->medge, sys->allocedge * sizeof(GFEdge));
				}
				sys->medge[sys->totedge].v1 = i;
				sys->medge[sys->totedge].v2 = sys->mvert[i].e[j];
				sys->mvert[i].e[j] = -1;
				sys->mvert[sys->mvert[i].e[j]].e[(j + 2) % 4] = -1;
				sys->totedge++;
			}
		}
	}
}

void makePolys(LaplacianSystem *sys)
{
	int i, j, d, s, dd;
	GFVertID v, n;

	sys->totloop = sys->allocloop = 0;
	sys->totpolys = sys->allocpolys = 0;
	MEM_SAFE_FREE(sys->loops);
	MEM_SAFE_FREE(sys->polys);

	for (i = 0; i < sys->totvert; i++) {
		for (j = 0; j < 4; j++) {
			if (sys->mvert[i].e[j] >= 0) {
				d = j;
				v = i;
				s = sys->totloop;
				do {
					if (sys->totloop == sys->allocloop) {
						sys->allocloop = sys->allocloop * 2 + 10;
						sys->loops = MEM_reallocN(sys->loops, sys->allocloop * sizeof(MLoop));
					}
					sys->loops[sys->totloop].v = v;
					sys->loops[sys->totloop].e = 0;
					sys->totloop++;

					n = sys->mvert[v].e[d];
					sys->mvert[v].e[d] = -1;
					v = n;
					dd = d;
					d = (d + 1) % 4;
					while (d != dd && sys->mvert[v].e[d] < 0) d = (d + 3) % 4;
					if (sys->mvert[v].e[d] < 0) break;
				} while (v != i);

				if (sys->totpolys == sys->allocpolys) {
					sys->allocpolys = sys->allocpolys * 2 + 10;
					sys->polys = MEM_reallocN(sys->polys, sys->allocpolys * sizeof(MPoly));
				}
				sys->polys[sys->totpolys].loopstart = s;
				sys->polys[sys->totpolys].totloop = sys->totloop - s;
				sys->polys[sys->totpolys].mat_nr = 0;
				sys->polys[sys->totpolys].flag = 0;
				sys->totpolys++;
			}
		}
	}
}

void generateMesh(LaplacianSystem *sys)
{
	MVert *arrayvect;
	MEdge *arrayedge;
	MPoly *polys;
	MLoop *loops;
	int i;

#if 0
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

	/*sys->resultDM = CDDM_new(sys->totvert,
					sys->gfsys1->totedge + sys->gfsys2->totedge,
					0, 0, 0);
					arrayvect = sys->resultDM->getVertArray(sys->resultDM);
					for (i = 0; i < sys->totvert; i++) {
					copy_v3_v3(arrayvect[i].co, sys->mvert[i].co);
					}

					arrayedge = sys->resultDM->getEdgeArray(sys->resultDM);
					for (i = 0; i < sys->gfsys1->totedge; i++) {
					arrayedge[i].v1 = sys->gfsys1->medge[i].v1;
					arrayedge[i].v2 = sys->gfsys1->medge[i].v2;
					arrayedge[i].flag |= ME_EDGEDRAW;
					}
					for (i = 0; i < sys->gfsys2->totedge; i++) {
					arrayedge[i + sys->gfsys1->totedge].v1 = sys->gfsys2->medge[i].v1;
					arrayedge[i + sys->gfsys1->totedge].v2 = sys->gfsys2->medge[i].v2;
					arrayedge[i + sys->gfsys1->totedge].flag |= ME_EDGEDRAW;
					}*/
#endif // 0

	generateIntersectionsOnFaces(sys);
	deleteDegenerateVerts(sys);
	makePolys(sys);

	sys->resultDM = CDDM_new(sys->totvert, 0, 0, sys->totloop, sys->totpolys);
	arrayvect = sys->resultDM->getVertArray(sys->resultDM);
	for (i = 0; i < sys->totvert; i++) {
		copy_v3_v3(arrayvect[i].co, sys->mvert[i].co);
	}

	/*arrayedge = sys->resultDM->getEdgeArray(sys->resultDM);
	for (i = 0; i < sys->totedge; i++) {
		arrayedge[i].v1 = sys->medge[i].v1;
		arrayedge[i].v2 = sys->medge[i].v2;
		arrayedge[i].flag |= ME_EDGEDRAW;
	}*/

	loops = sys->resultDM->getLoopArray(sys->resultDM);
	memcpy(loops, sys->loops, sys->totloop * sizeof(MLoop));
	polys = sys->resultDM->getPolyArray(sys->resultDM);
	memcpy(polys, sys->polys, sys->totpolys * sizeof(MPoly));
}