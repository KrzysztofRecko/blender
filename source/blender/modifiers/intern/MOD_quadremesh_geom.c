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

#include "PIL_time.h"

#include "MOD_util.h"

//#define QR_SHOWQUERIES
#define QR_LINELIMIT 1000000


#if 0 /* UNUSED ROUTINES */
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

/* MATH STUFF */

static bool isectLines(const float v1[3], const float v2[3],
                       const float v3[3], const float v4[3],
                       float vi[3], float *r_lambda)
{
	float a[3], b[3], c[3], ab[3], cb[3], ca[3];//, dir1[3], dir2[3];
	float d, div;

	sub_v3_v3v3(c, v3, v1);
	sub_v3_v3v3(a, v2, v1);
	sub_v3_v3v3(b, v4, v3);

	//normalize_v3_v3(dir1, a);
	//normalize_v3_v3(dir2, b);
	//d = dot_v3v3(dir1, dir2);
	//if (d == 1.0f || d == -1.0f) {
	//	/* colinear */
	//	return false;
	//}

	cross_v3_v3v3(ab, a, b);
	d = dot_v3v3(c, ab);
	div = dot_v3v3(ab, ab);

	/* test zero length line */
	if (UNLIKELY(div == 0.0f)) {
		return false;
	}
	/* test if the two lines are coplanar */
	else if (d > -0.000001f && d < 0.000001f) {
		float f1, f2;
		cross_v3_v3v3(cb, c, b);
		cross_v3_v3v3(ca, c, a);

		f1 = dot_v3v3(cb, ab) / div;
		f2 = dot_v3v3(ca, ab) / div;

		if (f1 >= 0.0f && f1 <= 1.0f &&
		    f2 >= 0.0f && f2 <= 1.0f)
		{
			if (vi) {
				mul_v3_fl(a, f1);
				add_v3_v3v3(vi, v1, a);
			}

			if (r_lambda) *r_lambda = f1;

			return true; /* intersection found */
		}
	}
	
	return false;
}

static bool isPointOnSegment(float in_v[3], float in_l1[3], float in_l2[3])
{
	float a[3], b[3];
	sub_v3_v3v3(a, in_v, in_l1);
	sub_v3_v3v3(b, in_l2, in_v);
	return dot_v3v3(a, b) >= 0.0f;
}

float absoluteAngleAxis(const float v1[3], const float v2[3], const float axis[3])
{
	float tproj[3];
	float angle, angle2, dot;

	//angle = angle_normalized_v3v3(v1, v2);

	dot = dot_v3v3(v1, v2);
	CLAMP(dot, -1.0f, 1.0f);
	angle = acos(dot);
	//if (!IS_EQF(angle - angle2, 0.0f))
		//printf("%.8f  |  %.8f, %.8f, %.8f\n", angle - angle2, angle, angle2, dot_v3v3(v1, v2_proj));

	/* calculate the sign (reuse 'tproj') */
	cross_v3_v3v3(tproj, v2, v1);
	if (dot_v3v3(tproj, axis) < 0.0f) {
		angle = ((float)(M_PI * 2.0f)) - angle;
	}

	return angle;
}

static void getPerpendicularToNormal(float r[3], float in_a[3], float in_b[3], float in_no[3])
{
	float proj[3];

	sub_v3_v3v3(r, in_b, in_a);
	project_v3_v3v3(proj, r, in_no);
	sub_v3_v3(r, proj);
}

/* QRVert STUFF */

static QRDiskLink *getLink(OutputMesh *om, QRVertID in_v1, QRVertID in_v2)
{
	QRDiskLink *it;

	it = om->verts[in_v1].link;
	while(it) {
		if (it->v == in_v2)
			return it;
		
		it = it->next;
		if (it == om->verts[in_v1].link)
			break;
	}

	return NULL;
}

static QRDiskLink *insertLink(OutputMesh *om, QRVertID in_a, QRVertID in_b)
{
	float vec[3];
	QRDiskLink *it, *l;

	l = BLI_memarena_alloc(om->memarena, sizeof(QRDiskLink));

	l->e = -1;
	l->v = in_b;
	l->poly_on_right = false;

	getPerpendicularToNormal(vec, om->verts[in_a].co, om->verts[in_b].co, om->verts[in_a].no);
	normalize_v3(vec);

	if (om->verts[in_a].num_links == 0) {
		om->verts[in_a].link = l;
		copy_v3_v3(om->verts[in_a].vec, vec);

		l->next = l;
		l->prev = l;
		l->ang = 0.0f;
	}
	else {
		l->ang = absoluteAngleAxis(om->verts[in_a].vec, vec, om->verts[in_a].no);

		if (l->ang <= om->verts[in_a].link->ang) {
			it = om->verts[in_a].link->prev;
			om->verts[in_a].link = l;
		}
		else {
			for (it = om->verts[in_a].link; it->next != om->verts[in_a].link; it = it->next) {
				if (it->next->ang > l->ang)
					break;
			}
		}

		it->next->prev = l;
		l->next = it->next;
		l->prev = it;
		it->next = l;
	}
	om->verts[in_a].num_links++;

	return l;
}

static QRDiskLink *linkVerts(OutputMesh *om, QRVertID in_v1, QRVertID in_v2)
{
	QRDiskLink *l1, *l2;

	l1 = getLink(om, in_v1, in_v2);
	if (l1)
		return l1;

	l1 = insertLink(om, in_v1, in_v2);
	l2 = insertLink(om, in_v2, in_v1);

	l1->brother = l2;
	l2->brother = l1;

	om->num_edges++;

	return l1;
}

static void deleteLink(OutputMesh *om, QRVert *ll, QRDiskLink *l)
{
	BLI_assert(ll->num_links != 0);

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

static void unlinkVerts(OutputMesh *om, QRDiskLink *l)
{
	QRVert *ll1, *ll2;

	ll2 = &om->verts[l->v];
	ll1 = &om->verts[l->brother->v];

	deleteLink(om, ll1, l);
	deleteLink(om, ll2, l->brother);

	om->num_edges--;
}

/* ADDING STUFF TO OUTPUT MESH */

static QRVertID addVert(OutputMesh *om, float in_co[3], float in_no[3])
{
	if (om->num_verts == om->alloc_verts) {
		om->alloc_verts = om->alloc_verts * 2 + 10;
		om->verts = MEM_reallocN(om->verts, sizeof(QRVert) * om->alloc_verts);
	}

	copy_v3_v3(om->verts[om->num_verts].co, in_co);
	if (in_no)
		copy_v3_v3(om->verts[om->num_verts].no, in_no);

	om->verts[om->num_verts].link = NULL;
	om->verts[om->num_verts].num_links = 0;
	
	return om->num_verts++;
}

/* SEED QUEUE */

static void addSeedToQueue(Heap *aheap, float in_co[3], GFPointType in_type, int in_val, float weight)
{
	GFPoint *seed = MEM_mallocN(sizeof(GFPoint), __func__);
	copy_v3_v3(seed->co, in_co);
	seed->type = in_type;
	seed->e = seed->f = seed->v = seed->id = -1;
	if (in_type == eVert)
		seed->v = in_val;
	else
		seed->f = in_val;
	

	BLI_heap_insert(aheap, weight, seed);
}

static GFPoint *getTopSeedFromQueue(struct Heap *aheap)
{
	GFPoint *p = BLI_heap_popmin(aheap);

	return p;
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

static void getInitialSeeds(GradientFlowSystem *gfsys)
{
	int i, *lverts, sizeverts[2];

	lverts = findFeaturesOnMesh(&gfsys->sys->input_mesh, sizeverts);
	
	for (i = 0; i < sizeverts[0]; i++)
		addSeedToQueue(gfsys->seeds, gfsys->sys->input_mesh.co[lverts[i]], eVert, lverts[i], 0.0f);
	
	MEM_SAFE_FREE(lverts);
}

/* QREdge ROUTINES */

static void appendOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl = BLI_memarena_alloc(om->memarena, sizeof(QREdgeLink));

	newl->v = in_v;
	newl->dist = in_dist;
	newl->elink = NULL;
	newl->next = NULL;

	if (in_e->v2 == NULL) {
		in_e->v1 = in_e->v2 = newl;
	}
	else {
		in_e->v2->next = newl;
		in_e->v2 = newl;
	}
}

static void prependOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl = BLI_memarena_alloc(om->memarena, sizeof(QREdgeLink));

	newl->v = in_v;
	newl->dist = in_dist;
	newl->next = in_e->v1;
	newl->elink = NULL;

	in_e->v1 = newl;
	if (in_e->v2 == NULL)
		in_e->v2 = newl;
}

static void insertAfterOnQREdge(OutputMesh *om, QREdge *in_e, QREdgeLink *in_l, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl = BLI_memarena_alloc(om->memarena, sizeof(QREdgeLink));

	newl->v = in_v;
	newl->dist = in_dist;
	newl->next = in_l->next;

	if (in_l->elink) {
		unlinkVerts(om, in_l->elink);
		in_l->elink = linkVerts(om, in_l->v, in_v);
		newl->elink = linkVerts(om, in_v, in_l->next->v);
	}
	else newl->elink = NULL;

	in_l->next = newl;
}

static void insertOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_vid)
{
	float tmp, vec[3];
	QRVert *verts = om->verts;
	QREdgeLink *it;

	if (in_e->num_links == 0) {
		appendOnQREdge(om, in_e, in_vid, 0.0f);
		copy_v3_v3(in_e->orig, verts[in_vid].co);
	}
	else if (in_e->num_links == 1) {
		sub_v3_v3v3(in_e->dir, verts[in_vid].co, in_e->orig);
		tmp = len_v3(in_e->dir);

		if (tmp < FLT_EPSILON)
			return;

		appendOnQREdge(om, in_e, in_vid, tmp);
		mul_v3_fl(in_e->dir, 1.0f / tmp);
	}
	else {
		sub_v3_v3v3(vec, verts[in_vid].co, in_e->orig);
		tmp = dot_v3v3(vec, in_e->dir);

		//if (IS_EQF(tmp, 0.0f))
			//return;

		if (tmp > in_e->v2->dist)
			appendOnQREdge(om, in_e, in_vid, tmp);
		else if (tmp < in_e->v1->dist)
			prependOnQREdge(om, in_e, in_vid, tmp);
		else {
			for (it = in_e->v1; it->next && it->next->dist < tmp; it = it->next);
			if (in_vid == it->v || in_vid == it->next->v) return;
			insertAfterOnQREdge(om, in_e, it, in_vid, tmp);
		}
	}
	in_e->num_links++;
}

static void linkOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_v1, QRVertID in_v2)
{
	QREdgeLink *it;

	for (it = in_e->v1; it && it->next && it->v != in_v1; it = it->next)
		if (it->v == in_v2) {
			in_v2 = in_v1;
			break;
		}

	for (; it && it->next && it->v != in_v2; it = it->next)
		if (!it->elink)
			it->elink = linkVerts(om, it->v, it->next->v);
}

static QREdge *addQREdgeToList(OutputMesh *om, LinkNode **list, QRVertID in_v1, QRVertID in_v2)
{
	QREdge *newe = BLI_memarena_calloc(om->memarena, sizeof(QREdge));

	insertOnQREdge(om, newe, in_v1);
	insertOnQREdge(om, newe, in_v2);
	linkOnQREdge(om, newe, in_v1, in_v2);
	BLI_linklist_prepend_arena(list, (void*)newe, om->memarena);
	
	return newe;
}

static bool isectSegmentWithQREdge(OutputMesh *om, float in_a[3], float in_b[3], QREdge *in_e)
{
	return isectLines(in_a, in_b, om->verts[in_e->v1->v].co, om->verts[in_e->v2->v].co, NULL, NULL);
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
static int nextPoint(float r_co[3], int *r_edge, InputMesh *im, int in_f, float in_co[3], float in_dir[3])
{
	int i, pick = -1, v = -1;
	bool is_on_vertex = false;
	float a[3][3], b[3][3], c[2][3], co2[3], dummy[3];

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
	QREdge *e;
	LinkNode *iter;

	for (iter = gfsys->ringf[in_f]; iter; iter = iter->next) {
		e = (QREdge*)iter->link;
		
		if (isectSegmentWithQREdge(&gfsys->sys->output_mesh, in_a, in_b, e))
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

	while (in_f != QR_NO_FACE) {
		len = dot_v3v3(dir, im->no[in_f]);
		mul_v3_v3fl(c, im->no[in_f], len);
		sub_v3_v3(dir, c);
		if (normalize_v3(dir) < FLT_EPSILON) 
			return 2;
		if (dot_v3v3(im->no[oldf], im->no[in_f]) < 0.0f) mul_v3_fl(dir, -1.0f);

		if (nextPoint(newco, &e, im, in_f, oldco, dir)) return 2;
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

			addSeedToQueue(gfsys->seeds, newco2, eFace, in_f, -maxdist);
			return 1;
		}

		in_f = getOtherFaceAdjacentToEdge(im, in_f, e);
		copy_v3_v3(oldco, newco);
	}
}

static bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], float in_newco[3], int in_f, float dist, float maxdist)
{
	int d;
	bool make_seed = BLI_frand() > gfsys->sys->qmd->seeding_probability;
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

#if 0
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
#endif // 0

/* GFLINE STUFF */

static void addGFPoint(InputMesh *im, OutputMesh *om, GFPoint *in_p)
{
	int i;
	float no[3];
	QRVertID newv;

	if (in_p->type == eVert) {
		if (om->vonvs[in_p->v] == -1) {
			zero_v3(no);
			for (i = 0; i < im->ringf_map[in_p->v].count; i++)
				add_v3_v3(no, im->no[im->ringf_map[in_p->v].indices[i]]);

			mul_v3_fl(no, 1.0f / (float)im->ringf_map[in_p->v].count);

			newv = addVert(om, in_p->co, no);
			om->vonvs[in_p->v] = newv;
		}
		else
			newv = om->vonvs[in_p->v];
	}
	else if (in_p->type == eEdge) {
		if (im->faces_edge[in_p->e][0] == -1)
			copy_v3_v3(no, im->no[im->faces_edge[in_p->e][1]]);
		else if (im->faces_edge[in_p->e][1] == -1)
			copy_v3_v3(no, im->no[im->faces_edge[in_p->e][0]]);
		else {
			add_v3_v3v3(no, im->no[im->faces_edge[in_p->e][0]],
						    im->no[im->faces_edge[in_p->e][1]]);
			mul_v3_fl(no, 0.5f);
		}

		newv = addVert(om, in_p->co, no);
	}
	else {
		copy_v3_v3(no, im->no[in_p->f]);
		newv = addVert(om, in_p->co, no);
	}

	in_p->id = newv;
}

static void addSegmentToLine(GFLine *line, GFPoint *in_p)
{
	int i;
	OutputMesh *om = &line->gfsys->sys->output_mesh;
	InputMesh *im = &line->gfsys->sys->input_mesh;

	addGFPoint(im, om, in_p);

	if (in_p->type == eVert) {
		for (i = 0; i < im->ringe_map[in_p->v].count; i++)
			insertOnQREdge(om, &om->ringe[im->ringe_map[in_p->v].indices[i]], in_p->id);

		if (line->end.type == eEdge)
			linkOnQREdge(om, &om->ringe[line->end.e], line->end.id, in_p->id);

		//addQREdge(gfsys, line->end, newv, in_p->f);
	}
	else if (in_p->type == eEdge) {
		insertOnQREdge(om, &om->ringe[in_p->e], in_p->id);

		if (line->end.type == eEdge && line->end.e == in_p->e) {
			linkOnQREdge(om, &om->ringe[in_p->e], line->end.id, in_p->id);
		}
		else {
			addQREdgeToList(om, &line->gfsys->ringf[in_p->f], line->end.id, in_p->id);
		}
	}
	else {
		addQREdgeToList(om, &line->gfsys->ringf[in_p->f], line->end.id, in_p->id);
	}

	memcpy(&line->end, in_p, sizeof(GFPoint));
}

static void flushGFLineQueue(GFLine *line)
{
	int i;

	for (i = 0; i < line->num_q; i++)
		addSegmentToLine(line, &line->q[i]);

	line->num_q = 0;
}

static bool enqueuePoint(GFLine *line, GFPoint *in_p)
{
	if (line->num_q == 10)
		return false;

	memcpy(&line->q[line->num_q], in_p, sizeof(GFPoint));
	line->num_q++;

	return true;
}

static bool changeLineDirection(GFLine *line)
{
	flushGFLineQueue(line);

	/* reset to original state */
	memcpy(&line->end, &line->seed, sizeof(GFPoint));
	memcpy(&line->lastchkp, &line->seed, sizeof(GFPoint));
	line->oldp = &line->seed;
	line->lastchklen = line->qlen = 0.0f;

	return ++line->d != 2;
}

static bool initGFLine(GradientFlowSystem *gfsys, GFLine *line, GFPoint *in_seed)
{
	const float maxdist = gfsys->sys->qmd->max_line_dist;
	const float seeddist = maxdist * 2.0f;

	int i;
	float old[3];
	OutputMesh *om = &gfsys->sys->output_mesh;
	InputMesh *im = &gfsys->sys->input_mesh;

	if (in_seed->type != eVert) {
		add_v3_v3v3(old, in_seed->co, gfsys->gf[in_seed->f]);
		if (!checkPoint(gfsys, old, in_seed->co, in_seed->f, maxdist, seeddist))
			return false;
	}
		
	addGFPoint(im, om, in_seed);

	memcpy(&line->seed, in_seed, sizeof(GFPoint));
	memcpy(&line->end, in_seed, sizeof(GFPoint));
	memcpy(&line->lastchkp, in_seed, sizeof(GFPoint));
	line->oldp = &line->seed;
	line->gfsys = gfsys;

	line->d = line->num_q = 0;
	line->lastchklen = line->qlen = 0.0f;

	if (line->seed.type == eVert) {
		for (i = 0; i < im->ringe_map[line->seed.v].count; i++)
			insertOnQREdge(om, &om->ringe[im->ringe_map[line->seed.v].indices[i]], line->seed.id);
	}

	return true;
}

static bool addPointToLine(GFLine *line, GFPoint *in_p)
{
	const float chklen = line->gfsys->sys->qmd->sampling_interval;
	const float maxdist = line->gfsys->sys->qmd->max_line_dist;
	const float seeddist = maxdist * 2.0f;

	float seg[3];
	float curlen;
	GFPoint newchkp;

	/* qco[0] - first point after last checked
	 * qco[num_q - 1] - last added point */

	sub_v3_v3v3(seg, in_p->co, line->oldp->co);
	curlen = len_v3(seg);

	while (line->qlen + curlen > line->lastchklen + chklen) {
		mul_v3_v3fl(newchkp.co, seg, (line->lastchklen + chklen - line->qlen) / curlen);
		add_v3_v3(newchkp.co, line->oldp->co);
		if (line->oldp->type == eEdge && line->oldp->e == in_p->e) {
			newchkp.type = eEdge;
			newchkp.e = in_p->e;
		}
		else
			newchkp.type = eFace;
		newchkp.f = in_p->f;

		if (!checkPoint(line->gfsys, line->lastchkp.co, newchkp.co, newchkp.f, maxdist, seeddist)) {
			addSegmentToLine(line, &line->lastchkp);
			line->num_q = 0;

			return false;
		}
		
		flushGFLineQueue(line);
		memcpy(&line->lastchkp, &newchkp, sizeof(GFPoint));
		line->lastchklen += chklen;
	}

	if(!enqueuePoint(line, in_p))
		return false;

	line->qlen += curlen;
	line->oldp = &line->q[line->num_q - 1];

	return true;
}

static void computeGFLine(GFLine *line)
{
	int i, r, newe;
	float gf[3], dir = 1.0f;
	LaplacianSystem *sys = line->gfsys->sys;
	InputMesh *im = &sys->input_mesh;
	GFPoint p;

	do {
		memcpy(&p, &line->seed, sizeof(GFPoint));

		do {
			if (p.type == eVert) {
				for (i = 0; i < im->ringf_map[p.v].count; i++) {
					p.f = im->ringf_map[p.v].indices[i];

					mul_v3_v3fl(gf, line->gfsys->gf[p.f], dir);
					if (!nextPoint(p.co, &p.e, im, p.f, line->oldp->co, gf)) {
						p.type = eEdge;
						break;
					}
				}
				if (p.type == eVert) break;
			}
			else {
				mul_v3_v3fl(gf, line->gfsys->gf[p.f], dir);
				r = nextPoint(p.co, &newe, im, p.f, line->oldp->co, gf);

				if (r == 1) {
					p.type = eVert;
					p.v = im->edges[newe][0];
				}
				else if (r == 2) {
					p.type = eVert;
					p.v = im->edges[newe][1];
				}
				else if (r == 3)
					break;
				else
					p.type = eEdge;

				if (newe == p.e) {
					if (dir * sys->U_field[im->edges[p.e][0]] < dir * sys->U_field[im->edges[p.e][1]])
						p.v = im->edges[p.e][0];
					else
						p.v = im->edges[p.e][1];

					copy_v3_v3(p.co, im->co[p.v]);
					p.type = eVert;
				}

				p.e = newe;
			}
			
			if (!addPointToLine(line, &p)) break;
			p.f = getOtherFaceAdjacentToEdge(im, p.f, p.e);
		} while (p.f != QR_NO_FACE);

		dir = -dir;
	} while (changeLineDirection(line));
}

static void computeFlowLines(LaplacianSystem *sys) {
	GFPoint *seed;
	GFLine line;
	int s, comp;
	
	for (s = 0; s < 2; s++) {
		getInitialSeeds(sys->gfsys[s]);

		comp = 0;
		while (!BLI_heap_is_empty(sys->gfsys[s]->seeds)) {
			seed = getTopSeedFromQueue(sys->gfsys[s]->seeds);

			if (++comp < QR_LINELIMIT && initGFLine(sys->gfsys[s], &line, seed))
				computeGFLine(&line);

			MEM_SAFE_FREE(seed);
		}
	}
}

/* MESH GENERATION */

static void generateIntersectionsOnFaces(OutputMesh *om, InputMesh *im)
{
	int f;
	float isection[3], lambda;
	QREdge *e1, *e2;
	LinkNode *iter1, *iter2;
	QRVertID newv;

	for (f = 0; f < im->num_faces; f++) {
		for (iter1 = om->ringf[0][f]; iter1; iter1 = iter1->next) {
			e1 = (QREdge*)iter1->link;

			for (iter2 = om->ringf[1][f]; iter2; iter2 = iter2->next) {
				e2 = (QREdge*)iter2->link;

				if (isectLines(om->verts[e1->v1->v].co, om->verts[e1->v2->v].co,
					           om->verts[e2->v1->v].co, om->verts[e2->v2->v].co,
					           isection, &lambda))
				{
					newv = addVert(om, isection, im->no[f]);
					insertOnQREdge(om, e1, newv);
					insertOnQREdge(om, e2, newv);
				}
			}
		}
	}
}

static void deleteDegenerateVerts(OutputMesh *om)
{
	int i, m, n;
	QRVertID a, b, *vertmap;
	QRVert *newverts;
	QRDiskLink *it;

	for (i = 0; i < om->num_verts; i++) {
		m = om->verts[i].num_links;

		if (m == 1)
			unlinkVerts(om, om->verts[i].link);
		else if (m == 2) {
			a = om->verts[i].link->v;
			b = om->verts[i].link->next->v;
			unlinkVerts(om, om->verts[i].link);
			unlinkVerts(om, om->verts[i].link);
			linkVerts(om, a, b);
		}
	}

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
		for (it = om->ringf[0][f]; it; it = it->next) {
			e = (QREdge*)it->link;
			for (itel = e->v1; itel; itel = itel->next)
				if (itel->elink) {
					itel->elink->e = 0;
					itel->elink->brother->e = 0;
				}
		}

		for (it = om->ringf[1][f]; it; it = it->next) {
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
	int i, j, s, num_loops, num_polys;
	QRVertID v;
	QRDiskLink *it, *lit;

	num_loops = 0;
	num_polys = 0;

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
			num_polys++;
		}
	}
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
	fill_vn_i(om->vonvs, im->num_verts, -1);
	om->ringe = MEM_callocN(sizeof(QREdge) * im->num_edges, "GFListEdges");
	om->ringf[0] = MEM_callocN(sizeof(LinkNode *) * im->num_faces, "GFListFaces");
	om->ringf[1] = MEM_callocN(sizeof(LinkNode *) * im->num_faces, "GFListFaces");
}

void freeOutputMesh(OutputMesh *om)
{
	om->num_edges = om->num_verts = om->alloc_verts = 0;
	MEM_SAFE_FREE(om->verts);
	MEM_SAFE_FREE(om->vonvs);
	MEM_SAFE_FREE(om->ringe);
	MEM_SAFE_FREE(om->ringf[0]);
	MEM_SAFE_FREE(om->ringf[1]);
	if (om->memarena) {
		BLI_memarena_free(om->memarena);
		om->memarena = NULL;
	}
}

DerivedMesh *makeResultMesh(LaplacianSystem *sys)
{
	int i;
	double start_time;
	MVert *verts;
	OutputMesh *om = &sys->output_mesh;
	DerivedMesh *ret;

	start_time = PIL_check_seconds_timer();

	initOutputMesh(om, &sys->input_mesh);
	sys->gfsys[0]->ringf = om->ringf[0];
	sys->gfsys[1]->ringf = om->ringf[1];

	computeFlowLines(sys);
	generateIntersectionsOnFaces(om, &sys->input_mesh);
	deleteDegenerateVerts(om);
	//hideEdgesOnFaces(om, &sys->input_mesh);
	//makeNormals(om);

	if (!om->num_verts)
		return NULL;

	ret = CDDM_new(om->num_verts, om->num_edges, 0,
				   om->num_edges * 2, 2 + om->num_edges - om->num_verts);

	verts = ret->getVertArray(ret);

	for (i = 0; i < om->num_verts; i++) {
		copy_v3_v3(verts[i].co, om->verts[i].co);
		normal_float_to_short_v3(verts[i].no, om->verts[i].no);
		verts[i].bweight = 0;
		verts[i].flag = 0;
	}

	makeEdges(om, ret->getEdgeArray(ret));
	makePolys(om, ret->getPolyArray(ret), ret->getLoopArray(ret));

	freeOutputMesh(om);
	CDDM_recalc_tessellation(ret);
	//CDDM_calc_edges_tessface(ret);
	//ret->dirty |= DM_DIRTY_NORMALS;

	printf("Mesh generation time: %f ms", (PIL_check_seconds_timer() - start_time) * 1000.0f);

	return ret;
}