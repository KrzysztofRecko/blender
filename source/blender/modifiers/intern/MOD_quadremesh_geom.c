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


#if 0 /* UNUSED ROUTINES */
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

/*
* alpha is degree of anisotropic curvature sensitivity
* h is the desired distance
* return ve[0] number of vertices
* return ve[1] number of edges
*/
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
#endif // UNUSED ROUTINES

static MVertID addVert(OutputMesh *om, float in_co[3], float in_no[3])
{
	if (om->totvert == om->allocvert) {
		om->allocvert = om->allocvert * 2 + 10;
		om->verts = MEM_reallocN(om->verts, sizeof(MVert) * om->allocvert);
		om->vlinks = MEM_reallocN(om->vlinks, sizeof(QREdgeLinkList) * om->allocvert);
	}

	copy_v3_v3(om->verts[om->totvert].co, in_co);
	if (in_no)
		normal_float_to_short_v3(om->verts[om->totvert].no, in_no);

	om->vlinks[om->totvert].link = NULL;
	om->vlinks[om->totvert].num_links = 0;
	om->totvert++;

	return om->totvert - 1;
}

static QREdgeLink *insertLink(OutputMesh *om, MVertID in_a, MVertID in_b)
{
	float vno[3];
	QREdgeLink *it, *l = BLI_memarena_alloc(om->memarena, sizeof(QREdgeLink));

	l->e = -1;
	l->v = in_b;
	l->poly_on_right = false;

	if (om->vlinks[in_a].num_links == 0) {
		om->vlinks[in_a].link = l;
		l->next = l;
		l->prev = l;
		l->ang = 0.0f;
		copy_v3_v3(om->vlinks[in_a].vec, om->verts[in_b].co);
	}
	else {
		normal_short_to_float_v3(vno, om->verts[in_a].no);
		l->ang = angle_signed_on_axis_v3v3v3_v3(om->vlinks[in_a].vec, om->verts[in_a].co, om->verts[in_b].co, vno);

		if (l->ang < om->vlinks[in_a].link->ang) {
			it = om->vlinks[in_a].link->prev;
			om->vlinks[in_a].link = l;
		}
		else
			for (it = om->vlinks[in_a].link; it->next != om->vlinks[in_a].link; it = it->next)
				if (it->next->ang > l->ang && it->ang < l->ang) break;

		it->next->prev = l;
		l->next = it->next;
		l->prev = it;
		it->next = l;
	}
	om->vlinks[in_a].num_links++;

	return l;
}

static QREdgeLink *linkVerts(OutputMesh *om, MVertID in_v1, MVertID in_v2)
{
	QREdgeLink *l1, *l2;

	l1 = insertLink(om, in_v1, in_v2);
	l2 = insertLink(om, in_v2, in_v1);

	l1->brother = l2;
	l2->brother = l1;

	return l1;
}

static void deleteLink(OutputMesh *om, QREdgeLinkList *ll, QREdgeLink *l)
{
	if (l->next == l) {
		ll->link = NULL;
	}
	else {
		if (ll->link == l) ll->link = l->next;
		l->prev->next = l->next;
		l->next->prev = l->prev;
	}
	ll->num_links--;
}

static void unlinkVerts(OutputMesh *om, QREdgeLink *l)
{
	QREdgeLinkList *ll1, *ll2;

	ll2 = &om->vlinks[l->v];
	ll1 = &om->vlinks[l->brother->v];

	deleteLink(om, ll1, l);
	deleteLink(om, ll2, l->brother);
}

static GFEdgeID addGFEdge(GradientFlowSystem *gfsys, MVertID in_v1, MVertID in_v2, int in_f)
{
	if (gfsys->totedge == gfsys->allocedge) {
		gfsys->allocedge = gfsys->allocedge * 2 + 10;
		gfsys->medge = MEM_reallocN(gfsys->medge, sizeof(GFEdge) * gfsys->allocedge);
	}

	gfsys->medge[gfsys->totedge].v1 = in_v1;
	gfsys->medge[gfsys->totedge].v2 = in_v2;

	BLI_linklist_prepend_arena(&gfsys->medge[gfsys->totedge].elinks,
						       (void*)linkVerts(&gfsys->sys->output_mesh, in_v1, in_v2), gfsys->memarena);
	BLI_linklist_prepend_arena(&gfsys->ringf_list[in_f], (void*)gfsys->totedge, gfsys->memarena);
	
	gfsys->totedge++;
	return gfsys->totedge - 1;
}

/* SEED QUEUE */

static void addSeedToQueue(Heap *aheap, float in_co[3], bool is_vertex, int in_val, float weight)
{
	GFSeed *seed = MEM_mallocN(sizeof(GFSeed), __func__);
	copy_v3_v3(seed->co, in_co);
	seed->val = in_val;
	seed->type = is_vertex ? eSeedVert : eSeedFace;

	BLI_heap_insert(aheap, weight, seed);
}

static GFSeed *getTopSeedFromQueue(struct Heap *aheap)
{
	GFSeed *vert = BLI_heap_popmin(aheap);

	return vert;
}

/*
* List of vertices from original mesh with special features (edge dihedral angle less that 90) to be preserves
* return the size of array
*/
static int *findFeaturesOnMesh(InputMesh *im, int size[2])
{
	int i, f1, f2, total;
	float angle;
	int *listverts = MEM_callocN(sizeof(int) * im->num_verts, __func__);
	int *listdest = NULL;
	total = 0;

	for (i = 0; i < im->num_edges; i++) {
		f1 = im->faces_edge[i][0];
		f2 = im->faces_edge[i][1];
		angle = angle_normalized_v3v3(im->no[f1], im->no[f2]);
		if (angle >= M_PI_2) {
			listverts[im->edges[i][0]] = 1;
			listverts[im->edges[i][1]] = 1;
		}
	}

	for (i = 0; i < im->num_verts; i++) {
		if (im->constraints[i] == 1) {
			listverts[i] = 1;
		}
	}

	for (i = 0; i < im->num_verts; i++) {
		if (listverts[i] == 1) {
			total++;
		}
	}
	if (total > 0) {
		listdest = MEM_mallocN(sizeof(int)* total, __func__);
	}
	total = 0;
	for (i = 0; i < im->num_verts; i++) {
		if (listverts[i] == 1) {
			listdest[total++] = i;
		}
	}
	MEM_SAFE_FREE(listverts);
	size[0] = total;
	return listdest;
}

/* GRADIENT FLOW SYSTEM MANAGEMENT */

GradientFlowSystem *newGradientFlowSystem(LaplacianSystem *sys, float *mhfunction, float(*mgfield)[3])
{
	int ve[2], i;
	int *lverts, sizeverts[2];
	GradientFlowSystem *gfsys = MEM_callocN(sizeof(GradientFlowSystem), "GradientFlowSystem");
	gfsys->memarena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Gradient FLow System");
	lverts = NULL;

	//estimateNumberGFVerticesEdges(ve, sys, sys->h);
	
	gfsys->sys = sys;
	//gfsys->mesh = newGradientFlowMesh(ve[0], ve[1]);
	
	gfsys->ringf_list = MEM_callocN(sizeof(LinkNode *) * sys->input_mesh.num_faces, "GFListFaces");
	
	gfsys->heap_seeds = BLI_heap_new();
	
	lverts = findFeaturesOnMesh(&sys->input_mesh, sizeverts);
	
	for (i = 0; i < sizeverts[0]; i++) {
		addSeedToQueue(gfsys->heap_seeds, sys->input_mesh.co[lverts[i]], true, lverts[i], 0.0f);
	}
	
	gfsys->hfunction = mhfunction;
	gfsys->gfield = mgfield;
	MEM_SAFE_FREE(lverts);
	
	return gfsys;
}

void deleteGradientFlowSystem(GradientFlowSystem *gfsys) 
{
	if (gfsys) {
		BLI_memarena_free(gfsys->memarena);
		BLI_heap_free(gfsys->heap_seeds, MEM_freeN);
		MEM_SAFE_FREE(gfsys->medge);
		MEM_SAFE_FREE(gfsys);
	}
}

/* GENERAL PURPOSE FUNCTIONS FOR GETTING AROUND THE FACES */

static int getEdgeFromVerts(InputMesh *im, int v1, int v2)
{
	int *eidn, nume, i;
	nume = im->ringe_map[v1].count;
	eidn = im->ringe_map[v1].indices;
	for (i = 0; i < nume; i++) {
		if (im->edges[eidn[i]][0] == v2 || im->edges[eidn[i]][1] == v2){
			return eidn[i];
		}
	}

	return -1;
}

static int getOtherFaceAdjacentToEdge(InputMesh *im, int in_f, int in_e)
{
	if (im->faces_edge[in_e][0] == in_f) {
		return im->faces_edge[in_e][1];
	}

	return im->faces_edge[in_e][0];
}

/**
 * /return 0 - all good
 *         1 - new vertex very close to 1st vertex of result edge
 *         2 -           -||-           2nd vertex of result edge
 *         3 - other error (eg. input outside of triangle)
 */
static int nextPoint(float r_co[3], int *r_edge, GradientFlowSystem *gfsys, int in_f, float in_co[3], float in_dir[3])
{
	int i, pick = -1, v = -1;
	bool is_on_vertex = false;
	float a[3][3], b[3][3], c[2][3], co2[3], dummy[3];
	InputMesh *im = &gfsys->sys->input_mesh;

	/* check if direction is coplanar to triangle */
	/* check if point is inside triangle */
	/* check if triangle is degenerate */

	add_v3_v3v3(co2, in_dir, in_co); /* second point on direction */

	for (i = 0; i < 3; i++) {
		sub_v3_v3v3(a[i], in_co, im->co[im->faces[in_f][i]]);
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
		dummy[i] = dot_v3v3(b[i], im->no[in_f]);
	}

	for (i = 0; i < 3; i++)
		if (dummy[i] < 0.0f && dummy[(i + 1) % 3] >= 0.0f)
			pick = i;

	if (pick == -1)
		return 3;

	*r_edge = getEdgeFromVerts(im, im->faces[in_f][pick], im->faces[in_f][(pick + 1) % 3]);

	isect_line_line_v3(im->co[im->faces[in_f][pick]], im->co[im->faces[in_f][(pick + 1) % 3]], in_co, co2, r_co, dummy);

	if (len_squared_v3v3(r_co, im->co[im->edges[*r_edge][0]]) < FLT_EPSILON)
		return 1;
	if (len_squared_v3v3(r_co, im->co[im->edges[*r_edge][1]]) < FLT_EPSILON)
		return 2;

	return 0;
}

static bool intersectSegmentWithOthersOnFace(GradientFlowSystem *gfsys, float in_a[3], float in_b[3], int in_f)
{
	int e;
	float dummy[3], lambda;
	LinkNode *iter;
	OutputMesh *om = &gfsys->sys->output_mesh;

	for (iter = gfsys->ringf_list[in_f]; iter; iter = iter->next) {
		e = (int)iter->link;
		/* TODO: Maybe grow the segment a little to pick up more stuff? */
		if (isect_line_line_strict_v3(in_a, in_b, om->verts[gfsys->medge[e].v1].co,
			                          om->verts[gfsys->medge[e].v2].co, dummy, &lambda))
			return true;
	}

	return false;
}

/**
 * 0 - intersection found
 * 1 - intersection not found
 * 2 - wrong direction for this face 
 */
static int queryDirection(GradientFlowSystem *gfsys, float in_co[3], int in_f, float in_dir[3], float dist, float maxdist, bool make_seed)
{
	int e, oldf;
	float c[3], len, actlen, oldco[3], newco[3], newco2[3], dir[3];
	InputMesh *im = &gfsys->sys->input_mesh;

	copy_v3_v3(dir, in_dir);
	copy_v3_v3(oldco, in_co);
	oldf = in_f;
	actlen = 0.0f;

	while (1) {
		len = dot_v3v3(dir, im->no[in_f]);
		mul_v3_v3fl(c, im->no[in_f], len);
		sub_v3_v3(dir, c);
		if (normalize_v3(dir) < FLT_EPSILON) 
			return 2;
		if (dot_v3v3(im->no[oldf], im->no[in_f]) < 0.0f) mul_v3_fl(dir, -1.0f);

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

		in_f = getOtherFaceAdjacentToEdge(im, in_f, e);
		copy_v3_v3(oldco, newco);
	}
}

static bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], float in_newco[3], int in_f, float dist, float maxdist)
{
	int d;
	bool make_seed = BLI_frand() > 0.75f;
	float seg[3], dir[3];
	InputMesh *im = &gfsys->sys->input_mesh;

	sub_v3_v3v3(seg, in_oldco, in_newco);
	if (dot_v3v3(seg, seg) < FLT_EPSILON) return true;

	cross_v3_v3v3(dir, im->no[in_f], seg);
	normalize_v3(dir);

	for (d = 0; d < 2; d++) {
		if (!queryDirection(gfsys, in_newco, in_f, dir, dist, maxdist, make_seed))
			return false;
		mul_v3_fl(dir, -1.0f);
	}

	return true;
}

/* SAMPLING DISTANCE FUNCTION, UNUSED RIGHT NOW */

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

/* GFLINE STUFF */

static void addSegmentToLine(GradientFlowSystem *gfsys, GFLine *line, int in_f, float in_co[3])
{
	MVertID newv;

	newv = addVert(&gfsys->sys->output_mesh, in_co, gfsys->sys->input_mesh.no[in_f]);
	addGFEdge(gfsys, line->end, newv, in_f);

	line->end = newv;
}

static bool changeLineDirection(GradientFlowSystem *gfsys, GFLine *line)
{
	int i;

	if (line->seed != -1) {
		/* flush queue */
		for (i = 0; i < line->num_q; i++)
			addSegmentToLine(gfsys, line, line->qf[i], line->qco[i]);

		line->num_q = 0;
		copy_v3_v3(line->oldco, gfsys->sys->output_mesh.verts[line->seed].co);
	}

	/* reset to original state */
	line->end = line->seed;
	copy_v3_v3(line->lastchk, line->oldco);

	return ++line->d != 2;
}

static void initGFLine(GradientFlowSystem *gfsys, GFLine *line, GFSeed *in_seed)
{
	line->seed = line->end = -1;
	line->d = 0;
	
	copy_v3_v3(line->oldco, in_seed->co);
}

static bool addPointToLine(GradientFlowSystem *gfsys, GFLine *line, int in_f, float in_newco[3])
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

		line->end = line->seed = addVert(&gfsys->sys->output_mesh, line->oldco, gfsys->sys->input_mesh.no[in_f]);
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

static void computeGFLine(GradientFlowSystem *gfsys, GFSeed *in_seed)
{
	int i, f, e, v, r, newe;
	bool is_vertex;
	float newco[3], gf[3], dir = 1.0f;
	InputMesh *im = &gfsys->sys->input_mesh;
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
				for (i = 0; i < im->ringf_map[v].count; i++) {
					f = im->ringf_map[v].indices[i];

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
					v = im->edges[newe][0];
				}
				else if (r == 2) {
					is_vertex = true;
					v = im->edges[newe][1];
				}
				else if (r == 3) break;

				if (newe == e) {
					if (dir * gfsys->sys->U_field[im->edges[e][0]] < dir * gfsys->sys->U_field[im->edges[e][1]])
						v = im->edges[e][0];
					else
						v = im->edges[e][1];

					copy_v3_v3(newco, im->co[v]);
					is_vertex = true;
				}

				e = newe;
			}
			
			if (!addPointToLine(gfsys, &line, f, newco))	break;

			f = getOtherFaceAdjacentToEdge(im, f, e);
		} /* while (1) */

		dir = -dir;
	} while (changeLineDirection(gfsys, &line));
}

static void computeFlowLines(LaplacianSystem *sys) {
	GFSeed *seed;
	int comp = 0;
	
	if (sys->gfsys1) deleteGradientFlowSystem(sys->gfsys1);
	if (sys->gfsys2) deleteGradientFlowSystem(sys->gfsys2);
	sys->gfsys1 = newGradientFlowSystem(sys, sys->h1, sys->gf1);
	sys->gfsys2 = newGradientFlowSystem(sys, sys->h2, sys->gf2);

	while (!BLI_heap_is_empty(sys->gfsys1->heap_seeds)) {
		seed = getTopSeedFromQueue(sys->gfsys1->heap_seeds);
		//if (++comp < QR_LINELIMIT)
			computeGFLine(sys->gfsys1, seed);
		MEM_SAFE_FREE(seed);
	}

	comp = 0;
	
	while (!BLI_heap_is_empty(sys->gfsys2->heap_seeds)) {
		seed = getTopSeedFromQueue(sys->gfsys2->heap_seeds);
		//if (++comp < QR_LINELIMIT)
			computeGFLine(sys->gfsys2, seed);
		MEM_SAFE_FREE(seed);
	}
}

/* MESH GENERATION */

static void insertOnEdge(GradientFlowSystem *gfsys, MVertID in_vid, GFEdgeID in_e)
{
	MVert *verts = gfsys->sys->output_mesh.verts;
	MVertID va = gfsys->medge[in_e].v1;
	MVertID vb = gfsys->medge[in_e].v2;
	LinkNode *it = gfsys->medge[in_e].elinks;
	float len = len_v3v3(verts[va].co, verts[in_vid].co);

	while (((QREdgeLink*)it->link)->v != vb) {
		if (len_v3v3(verts[va].co, verts[((QREdgeLink*)it->link)->v].co) > len) break;

		it = it->next;
	}
	va = ((QREdgeLink*)it->link)->brother->v;
	vb = ((QREdgeLink*)it->link)->v;

	unlinkVerts(&gfsys->sys->output_mesh, (QREdgeLink*)it->link);
	it->link = (void*)linkVerts(&gfsys->sys->output_mesh, va, in_vid);
	BLI_linklist_prepend_arena(&it->next, (void*)linkVerts(&gfsys->sys->output_mesh, in_vid, vb), gfsys->memarena);
}

static void generateIntersectionsOnFaces(LaplacianSystem *sys)
{
	int f;
	float isection[3], lambda;
	GFEdgeID e1, e2;
	LinkNode *iter1, *iter2;
	MVertID newv;

	for (f = 0; f < sys->input_mesh.num_faces; f++) {
		for (iter1 = sys->gfsys1->ringf_list[f]; iter1; iter1 = iter1->next) {
			e1 = (GFEdgeID)iter1->link;

			for (iter2 = sys->gfsys2->ringf_list[f]; iter2; iter2 = iter2->next) {
				e2 = (GFEdgeID)iter2->link;

				if (isect_line_line_strict_v3(sys->output_mesh.verts[sys->gfsys1->medge[e1].v1].co,
					                          sys->output_mesh.verts[sys->gfsys1->medge[e1].v2].co,
					                          sys->output_mesh.verts[sys->gfsys2->medge[e2].v1].co,
					                          sys->output_mesh.verts[sys->gfsys2->medge[e2].v2].co,
					                          isection, &lambda))
				{
					newv = addVert(&sys->output_mesh, isection, sys->input_mesh.no[f]);
					insertOnEdge(sys->gfsys1, newv, e1);
					insertOnEdge(sys->gfsys2, newv, e2);
				}
			}
		}
	}
}

static void deleteDegenerateVerts(OutputMesh *om)
{
	int i, n, m;
	MVertID a, b;
	MVert *verts = om->verts;

	for (i = 0; i < om->totvert; i++) {
		m = om->vlinks[i].num_links;

		if (m == 0) {
			verts[i].flag = ME_HIDE;
		}
		else if (m == 1) {
			verts[i].flag = ME_HIDE;
			unlinkVerts(om, om->vlinks[i].link);
		}
		else if (m == 2) {
			verts[i].flag = ME_HIDE;
			a = om->vlinks[i].link->v;
			b = om->vlinks[i].link->next->v;
			unlinkVerts(om, om->vlinks[i].link);
			unlinkVerts(om, om->vlinks[i].link);
			linkVerts(om, a, b);
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

static int addEdge(OutputMesh *om, unsigned int in_v1, unsigned int in_v2)
{
	if (om->totedge == om->allocedge) {
		om->allocedge = om->allocedge * 2 + 10;
		om->edges = MEM_reallocN(om->edges, om->allocedge * sizeof(MEdge));
	}
	om->edges[om->totedge].v1 = in_v1;
	om->edges[om->totedge].v2 = in_v2;
	om->edges[om->totedge].flag = ME_EDGEDRAW;
	
	return om->totedge++;
}

static void makeEdges(OutputMesh *om)
{
	int i, j;
	QREdgeLink *it;

	for (i = 0; i < om->totvert; i++) {
		if (om->verts[i].flag & ME_HIDE) continue;
		for (j = 0, it = om->vlinks[i].link; j < om->vlinks[i].num_links; j++, it = it->next) {
			if (it->e < 0 && !(om->verts[it->v].flag & ME_HIDE)) {
				it->brother->e = it->e = addEdge(om, i, it->v);
			}
		}
	}
}

static int addLoop(OutputMesh *om, unsigned int in_v, unsigned int in_e)
{
	if (om->totloop == om->allocloop) {
		om->allocloop = om->allocloop * 2 + 10;
		om->loops = MEM_reallocN(om->loops, om->allocloop * sizeof(MLoop));
	}
	om->loops[om->totloop].v = in_v;
	om->loops[om->totloop].e = in_e;
	
	return om->totloop++;
}

static int addPoly(OutputMesh *om, int in_ls, int in_tl)
{
	if (om->totpolys == om->allocpolys) {
		om->allocpolys = om->allocpolys * 2 + 10;
		om->polys = MEM_reallocN(om->polys, om->allocpolys * sizeof(MPoly));
	}
	om->polys[om->totpolys].loopstart = in_ls;
	om->polys[om->totpolys].totloop = in_tl;
	om->polys[om->totpolys].mat_nr = 0;
	om->polys[om->totpolys].flag = 0;

	return om->totpolys++;
}

static void makePolys(OutputMesh *om)
{
	int i, j, d, s, dd;
	MVertID v, n;
	QREdgeLink *it, *lit;

	for (i = 0; i < om->totvert; i++) {
		if (om->verts[i].flag & ME_HIDE) continue;

		for (j = 0, it = om->vlinks[i].link; j < om->vlinks[i].num_links; j++, it = it->next) {
			if (it->poly_on_right) continue;

			s = om->totloop;
			v = i;
			lit = it;
			do {
				addLoop(om, v, lit->e);
				v = lit->v;
				lit->poly_on_right = true;
				lit = lit->brother->next;
			} while (v != i);
			
			addPoly(om, s, om->totloop - s);
		}
	}
}

void makeNormals(OutputMesh *om)
{
	int i, tv = om->totvert;
	float npos[3];

	for (i = 0; i < tv; i++) {
		normal_short_to_float_v3(npos, om->verts[i].no);
		mul_v3_fl(npos, 0.1f);
		add_v3_v3(npos, om->verts[i].co);
		addEdge(om, i, addVert(om, npos, NULL));
	}
}

void freeOutputMesh(OutputMesh *om)
{
	MEM_SAFE_FREE(om->verts);
	MEM_SAFE_FREE(om->edges);
	MEM_SAFE_FREE(om->polys);
	MEM_SAFE_FREE(om->loops);
	if (om->memarena) {
		BLI_memarena_free(om->memarena);
		om->memarena = NULL;
	}
}

void generateMesh(LaplacianSystem *sys)
{
	OutputMesh *om = &sys->output_mesh;

	sys->output_mesh.memarena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Output Mesh");

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
					copy_v3_v3(arrayvect[i].co, sys->verts[i].co);
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

	computeFlowLines(sys);
	generateIntersectionsOnFaces(sys);
	//deleteDegenerateVerts(om);
	makeEdges(om);
	//makePolys(om);
	makeNormals(om);
}

DerivedMesh *getResultMesh(LaplacianSystem *sys)
{
	MVert *verts;
	MEdge *edges;
	MPoly *polys;
	MLoop *loops;
	OutputMesh *om = &sys->output_mesh;
	DerivedMesh *ret;

	if (!om->totvert) return NULL;

	ret = CDDM_new(om->totvert, om->totedge, 0, om->totloop, om->totpolys);

	verts = ret->getVertArray(ret);
	memcpy(verts, om->verts, om->totvert * sizeof(MVert));
	edges = ret->getEdgeArray(ret);
	memcpy(edges, om->edges, om->totedge * sizeof(MEdge));
	loops = ret->getLoopArray(ret);
	memcpy(loops, om->loops, om->totloop * sizeof(MLoop));
	polys = ret->getPolyArray(ret);
	memcpy(polys, om->polys, om->totpolys * sizeof(MPoly));

	CDDM_recalc_tessellation(ret);
	CDDM_calc_edges_tessface(ret);
	//ret->dirty |= DM_DIRTY_NORMALS;

	return ret;
}