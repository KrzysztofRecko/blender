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
#define QR_LINELIMIT 100000


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

static bool isPointOnSegment(float in_v[3], float in_l1[3], float in_l2[3])
{
	float a[3], b[3];
	sub_v3_v3v3(a, in_v, in_l1);
	sub_v3_v3v3(b, in_l2, in_v);
	return dot_v3v3(a, b) >= 0.0f;
}
#endif // UNUSED ROUTINES

/* MATH STUFF */

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

static void getNormalAtEdge(float r_no[3], InputMesh *im, int in_e)
{
	if (im->faces_edge[in_e][0] == QR_NO_FACE)
		copy_v3_v3(r_no, im->no[im->faces_edge[in_e][1]]);
	else if (im->faces_edge[in_e][1] == QR_NO_FACE)
		copy_v3_v3(r_no, im->no[im->faces_edge[in_e][0]]);
	else {
		add_v3_v3v3(r_no, im->no[im->faces_edge[in_e][0]],
						  im->no[im->faces_edge[in_e][1]]);
		mul_v3_fl(r_no, 0.5f);
	}
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

	sub_v3_v3v3(vec, om->verts[in_b].co, om->verts[in_a].co);
	project_plane_v3_v3v3(vec, vec, om->verts[in_a].no);
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
			break;
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
	newl->gfsysid = 0;

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
	newl->gfsysid = 0;

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
	newl->gfsysid = in_l->gfsysid;

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

static void linkOnQREdge(OutputMesh *om, GFSysID sys_id, QREdge *in_e, QRVertID in_v1, QRVertID in_v2)
{
	QREdgeLink *it;

	for (it = in_e->v1; it && it->next && it->v != in_v1; it = it->next)
		if (it->v == in_v2) {
			in_v2 = in_v1;
			break;
		}

	for (; it && it->next && it->v != in_v2; it = it->next) {
		if (!it->elink)
			it->elink = linkVerts(om, it->v, it->next->v);

		it->gfsysid |= sys_id;
	}
}

static QREdge *addQREdgeToFace(OutputMesh *om, InputMesh *im, GFSysID sys_id, int in_f, QRVertID in_v1, QRVertID in_v2)
{
	float isection[3];
	QREdge *e;
	QREdge *newe = BLI_memarena_calloc(om->memarena, sizeof(QREdge));
	QRVertID newv;
	LinkNode *it;

	insertOnQREdge(om, newe, in_v1);
	insertOnQREdge(om, newe, in_v2);
	linkOnQREdge(om, sys_id, newe, in_v1, in_v2);
	BLI_linklist_prepend_arena(&om->ringf[in_f], (void*)newe, om->memarena);

	for (it = om->ringf[in_f]->next; it; it = it->next) {
		e = (QREdge*)it->link;

		if ((e->v1->gfsysid & sys_id) != 0)
			continue;

		if (isect_seg_seg_unsafe_v3(isection,
			                        om->verts[e->v1->v].co, om->verts[e->v2->v].co,
					                om->verts[in_v1].co, om->verts[in_v2].co))
		{
			newv = addVert(om, isection, im->no[in_f]);
			insertOnQREdge(om, e, newv);
			insertOnQREdge(om, newe, newv);
		}
	}
	
	return newe;
}

static bool isectSegmentWithOthersOnFace(OutputMesh *om, GFSysID sys_id, float in_a[3], float in_b[3], int in_f)
{
	QREdge *e;
	LinkNode *iter;

	for (iter = om->ringf[in_f]; iter; iter = iter->next) {
		e = (QREdge*)iter->link;
		
		if ((e->v1->gfsysid & sys_id) == 0)
			continue;

		if (isect_seg_seg_unsafe_v3(NULL, in_a, in_b, om->verts[e->v1->v].co, om->verts[e->v2->v].co))
			return true;
	}

	return false;
}

static bool isectPointWithQREdge(OutputMesh *om, GFSysID sys_id, float in_co[3], int in_e)
{
	float vec[3], tmp;
	QREdge *e = &om->ringe[in_e];
	QREdgeLink *it;

	if (!e || e->num_links < 2)
		return false;

	sub_v3_v3v3(vec, in_co, e->orig);
	tmp = dot_v3v3(vec, e->dir);

	if (tmp <= e->v2->dist && tmp >= e->v1->dist) {
		for (it = e->v1; it->next && it->next->dist < tmp; it = it->next);
		
		if (it->elink != NULL && (it->gfsysid & sys_id) != 0)
			return true;
	}

	return false;
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

static void makePoint(GFPoint *r_p, GFPointType p_type, int p_vef, float p_co[3])
{
	r_p->type = p_type;
	r_p->v = p_vef;
	r_p->id = -1;
	copy_v3_v3(r_p->co, p_co);
}

static void makeSegment(GFSegment *r, InputMesh *im, GFPoint *prev, int in_f)
{
	if (prev->type == eEdge && r->p.type == eEdge) {
		if (prev->e == r->p.e) { /* eEdgeSegment */
			r->type = eEdgeSegment;
			r->e = r->p.e;
		}
		else { /* eFaceSegment */
			r->type = eFaceSegment;
			r->f = in_f; //getFaceFromTwoEdges(im, prev->e, p_vef);
		}
	}
	else if (prev->type == eVert && r->p.type == eEdge) { /* likely eFaceSegment */
		if (im->edges[r->p.e][0] == prev->v || im->edges[r->p.e][1] == prev->v) {
			r->type = eEdgeSegment;
			r->e = r->p.e;
		}
		else {
			r->type = eFaceSegment;
			r->f = in_f;
		}
	}
	else if (prev->type == eEdge && r->p.type == eVert) { /* dunno */
		if (im->edges[prev->e][0] == r->p.v || im->edges[prev->e][1] == r->p.v) {
			r->type = eEdgeSegment;
			r->e = prev->e;
		}
		else {
			r->type = eFaceSegment;
			r->f = in_f;
		}
	}
	else if (prev->type == eVert && r->p.type == eVert) { /* eEdgeSegment */
		r->type = eEdgeSegment;
		r->e = getEdgeFromVerts(im, prev->v, r->p.v);
	}
	else { /* prev on a face or this on a face */
		r->type = eFaceSegment;
		r->f = in_f; //p_type == eFace ? p_vef : prev->f;
	}
}

/**
 * /return true - all ok
 *         false - error
 */
static bool nextPointOnFace(GFPoint *r_p, InputMesh *im, GFPoint *in_p, int in_f, float in_dir[3])
{
	int i, pick = -1, v = -1, e;
	bool is_on_vertex = false;
	float a[3][3], b[3][3], c[2][3], co2[3], result[3], dummy[3];

	/* check if direction is coplanar to triangle */
	/* check if point is inside triangle */
	/* check if triangle is degenerate */

	add_v3_v3v3(co2, in_dir, in_p->co); /* second point on direction */

	for (i = 0; i < 3; i++) {
		sub_v3_v3v3(a[i], in_p->co, im->co[im->faces[in_f][i]]);
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
		return false;

	e = getEdgeFromVerts(im, im->faces[in_f][pick], im->faces[in_f][(pick + 1) % 3]);

	isect_line_line_v3(im->co[im->faces[in_f][pick]],
					   im->co[im->faces[in_f][(pick + 1) % 3]],
					   in_p->co, co2, result, dummy);

	if (len_squared_v3v3(result, im->co[im->edges[e][0]]) < FLT_EPSILON)
		makePoint(r_p, eVert, im->edges[e][0], im->co[im->edges[e][0]]);
	else if (len_squared_v3v3(result, im->co[im->edges[e][1]]) < FLT_EPSILON)
		makePoint(r_p, eVert, im->edges[e][1], im->co[im->edges[e][1]]);
	else
		makePoint(r_p, eEdge, e, result);

	return true;
}

/**
 * 0 - intersection found
 * 1 - intersection not found
 * 2 - wrong direction for this face
 */
static int queryDirection(GradientFlowSystem *gfsys, GFPoint *in_p, int in_f, float in_dir[3], float dist, float maxdist, bool make_seed)
{
	int oldf;
	float c[3], len, actlen, prevlen, chkco[3], seedco[3], dir[3];
	InputMesh *im = &gfsys->sys->input_mesh;
	OutputMesh *om = &gfsys->sys->output_mesh;
	GFPoint oldp, newp;

	copy_v3_v3(dir, in_dir);
	memcpy(&oldp, in_p, sizeof(GFPoint));
	oldf = in_f;
	actlen = prevlen = 0.0f;

	while (in_f != QR_NO_FACE && oldp.type != eVert) {
		project_plane_v3_v3v3(dir, dir, im->no[in_f]);

		if (normalize_v3(dir) < FLT_EPSILON) 
			break;

		if (dot_v3v3(im->no[oldf], im->no[in_f]) < 0.0f)
			mul_v3_fl(dir, -1.0f);

		if (!nextPointOnFace(&newp, im, &oldp, in_f, dir))
			break;
		oldf = in_f;

		sub_v3_v3v3(c, newp.co, oldp.co);
		len = len_v3(c);
		actlen = prevlen + len;

		if (actlen > maxdist && make_seed) {
			madd_v3_v3v3fl(seedco, oldp.co, c, (maxdist - prevlen) / len);

			addSeedToQueue(gfsys->seeds, seedco, eFace, in_f, -maxdist);
		}
		else
			copy_v3_v3(seedco, newp.co);

		if (actlen > dist && prevlen < dist)
			madd_v3_v3v3fl(chkco, oldp.co, c, (dist - prevlen) / len);
		else
			copy_v3_v3(chkco, newp.co);

#ifdef QR_SHOWQUERIES
		int vf1, vf2;
		if (make_seed) {
			vf1 = addVert(&gfsys->sys->output_mesh, oldp.co, NULL);
			vf2 = addVert(&gfsys->sys->output_mesh, seedco, NULL);
		}
		else {
			vf1 = addVert(&gfsys->sys->output_mesh, oldp.co, NULL);
			vf2 = addVert(&gfsys->sys->output_mesh, chkco, NULL);
		}
		linkVerts(&gfsys->sys->output_mesh, vf1, vf2);
#endif

		if (prevlen < dist) {
			if (isectSegmentWithOthersOnFace(om, gfsys->id, oldp.co, chkco, in_f))
				return 0;

			if (newp.type == eEdge && actlen <= dist) {
				if (isectPointWithQREdge(om, gfsys->id, chkco, newp.e))
					return 0;
			}
		}
		
		if ((actlen > dist && !make_seed) || actlen > maxdist)
			return 1;

		if (newp.type == eVert)
			break;

		in_f = getOtherFaceAdjacentToEdge(im, in_f, newp.e);
		memcpy(&oldp, &newp, sizeof(GFPoint));
		prevlen = actlen;
	}

	return 2;
}

static bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], GFPoint *in_p, float dist, float maxdist)
{
	bool make_seed = BLI_rng_get_float(gfsys->sys->rng) < gfsys->sys->qmd->seeding_probability;

	int d, f;
	float seg[3], dir[3], no[3];
	InputMesh *im = &gfsys->sys->input_mesh;

	sub_v3_v3v3(seg, in_oldco, in_p->co);

	if (dot_v3v3(seg, seg) < FLT_EPSILON)
		return true;

	if (in_p->type == eEdge) {
		getNormalAtEdge(no, im, in_p->e);
		f = im->faces_edge[in_p->e][0];
		if (f == QR_NO_FACE)
			f = im->faces_edge[in_p->e][1];
		if (f == QR_NO_FACE)
			return true;
	}
	else {
		copy_v3_v3(no, im->no[in_p->f]);
		f = in_p->f;
	}

	cross_v3_v3v3(dir, no, seg);
	normalize_v3(dir);

	for (d = 0; d < 2; d++) {
		if (!queryDirection(gfsys, in_p, f, dir, dist, maxdist, make_seed))
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
		getNormalAtEdge(no, im, in_p->e);
		newv = addVert(om, in_p->co, no);
	}
	else { /* eFace */
		copy_v3_v3(no, im->no[in_p->f]);
		newv = addVert(om, in_p->co, no);
	}

	in_p->id = newv;
}

static void addSegmentToLine(GFLine *line, GFSegment *in_s)
{
	int i;
	OutputMesh *om = &line->gfsys->sys->output_mesh;
	InputMesh *im = &line->gfsys->sys->input_mesh;

	addGFPoint(im, om, &in_s->p);

	if (in_s->p.type == eVert) {
		for (i = 0; i < im->ringe_map[in_s->p.v].count; i++)
			insertOnQREdge(om, &om->ringe[im->ringe_map[in_s->p.v].indices[i]], in_s->p.id);
	}
	else if (in_s->p.type == eEdge)
		insertOnQREdge(om, &om->ringe[in_s->p.e], in_s->p.id);

	if (in_s->type == eEdgeSegment)
		linkOnQREdge(om, line->gfsys->id, &om->ringe[in_s->e], line->end.id, in_s->p.id);
	else
		addQREdgeToFace(om, im, line->gfsys->id, in_s->f, line->end.id, in_s->p.id);

	memcpy(&line->end, &in_s->p, sizeof(GFPoint));
}

static void flushGFLineQueue(GFLine *line)
{
	int i;

	for (i = 0; i < line->num_q; i++)
		addSegmentToLine(line, &line->q[i]);

	line->num_q = 0;
}

static bool enqueueSegment(GFLine *line, GFSegment *in_s)
{
	if (line->num_q == QR_GFLINE_QSIZE)
		return false;

	memcpy(&line->q[line->num_q], in_s, sizeof(GFSegment));
	line->num_q++;

	return true;
}

static bool changeLineDirection(GFLine *line)
{
	flushGFLineQueue(line);

	/* reset to original state */
	memcpy(&line->end, &line->seed, sizeof(GFPoint));
	line->lastchkp = &line->seed;
	line->oldp = &line->seed;
	line->lastchklen = line->qlen = 0.0f;
	line->lastchks.type = eNoSegment;

	return ++line->d != 2;
}

static bool initGFLine(GradientFlowSystem *gfsys, GFLine *line, GFPoint *in_seed)
{
	const float maxdist = gfsys->sys->qmd->max_line_dist;
	const float seeddist = maxdist * QR_MAXDIST_TO_SEEDDIST;

	int i;
	float old[3];
	OutputMesh *om = &gfsys->sys->output_mesh;
	InputMesh *im = &gfsys->sys->input_mesh;

	if (in_seed->type == eFace) {
		add_v3_v3v3(old, in_seed->co, gfsys->gf[in_seed->f]);
		if (!checkPoint(gfsys, old, in_seed, maxdist, seeddist))
			return false;
	}
		
	addGFPoint(im, om, in_seed);

	memcpy(&line->seed, in_seed, sizeof(GFPoint));
	memcpy(&line->end, in_seed, sizeof(GFPoint));
	line->lastchkp = &line->seed;
	line->oldp = &line->seed;
	line->gfsys = gfsys;
	line->lastchks.type = eNoSegment;

	line->d = line->num_q = 0;
	line->lastchklen = line->qlen = 0.0f;

	if (line->seed.type == eVert) {
		for (i = 0; i < im->ringe_map[line->seed.v].count; i++)
			insertOnQREdge(om, &om->ringe[im->ringe_map[line->seed.v].indices[i]], line->seed.id);
	}

	return true;
}

static bool nextLineSegment(GFLine *line, GFSegment *in_s)
{
	const float chklen = line->gfsys->sys->qmd->sampling_interval;
	const float maxdist = line->gfsys->sys->qmd->max_line_dist;
	const float seeddist = maxdist * QR_MAXDIST_TO_SEEDDIST;

	float seg[3], co[3];
	float curlen;
	GFSegment newchks;
	InputMesh *im = &line->gfsys->sys->input_mesh;

	/* qco[0] - first point after last checked
	 * qco[num_q - 1] - last added point */

	sub_v3_v3v3(seg, in_s->p.co, line->oldp->co);
	curlen = len_v3(seg);

	while (line->qlen + curlen > line->lastchklen + chklen) {
		mul_v3_v3fl(co, seg, (line->lastchklen + chklen - line->qlen) / curlen);
		add_v3_v3(co, line->oldp->co);

		if (in_s->type == eEdgeSegment)
			makePoint(&newchks.p, eEdge, in_s->e, co);
		else
			makePoint(&newchks.p, eFace, in_s->f, co);

		makeSegment(&newchks, im, line->oldp, in_s->f);

		if (!checkPoint(line->gfsys, line->lastchkp->co, &newchks.p, maxdist, seeddist)) {
			if (line->lastchks.type != eNoSegment)
				addSegmentToLine(line, &line->lastchks);
			line->num_q = 0;

			return false;
		}
		
		flushGFLineQueue(line);
		memcpy(&line->lastchks, &newchks, sizeof(GFSegment));
		line->lastchkp = &line->lastchks.p;
		line->lastchklen += chklen;
	}

	if(!enqueueSegment(line, in_s))
		return false;

	line->qlen += curlen;
	line->oldp = &line->q[line->num_q - 1].p;

	return true;
}

static void computeGFLine(GFLine *line)
{
	int i, f = QR_NO_FACE;
	float gf[3], dir = 1.0f;
	LaplacianSystem *sys = line->gfsys->sys;
	InputMesh *im = &sys->input_mesh;
	GFSegment news;

	do {
		if (line->oldp->type == eFace)
			f = line->oldp->f;

		do {
			if (line->oldp->type == eVert) {
				for (i = 0; i < im->ringf_map[line->oldp->v].count; i++) {
					f = im->ringf_map[line->oldp->v].indices[i];

					mul_v3_v3fl(gf, line->gfsys->gf[f], dir);
					if (!nextPointOnFace(&news.p, im, line->oldp, f, gf))
						continue;

					if (news.p.type == eEdge)
						break;
				}
				if (news.p.type != eEdge)
					break;
			}
			else {
				mul_v3_v3fl(gf, line->gfsys->gf[f], dir);
				if (!nextPointOnFace(&news.p, im, line->oldp, f, gf))
					break;
				
				if (news.p.type == eEdge && line->oldp->type == eEdge && line->oldp->e == news.p.e) {
					if (dir * sys->U_field[im->edges[news.p.e][0]] < dir * sys->U_field[im->edges[news.p.e][1]])
						makePoint(&news.p, eVert, im->edges[news.p.e][0], im->co[im->edges[news.p.e][0]]);
					else
						makePoint(&news.p, eVert, im->edges[news.p.e][1], im->co[im->edges[news.p.e][1]]);
				}
			}

			makeSegment(&news, im, line->oldp, f);
			if (!nextLineSegment(line, &news))
				break;
				
			if (line->oldp->type == eEdge)
				f = getOtherFaceAdjacentToEdge(im, f, line->oldp->e);
		} while (f != QR_NO_FACE);

		dir = -dir;
	} while (changeLineDirection(line));
}

/* MESH GENERATION */

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

static void makeFeatureEdges(OutputMesh *om, InputMesh *im)
{
	int i, j, f1, f2;
	float angle;
	GFPoint p1, p2;

	for (i = 0; i < im->num_edges; i++) {
		f1 = im->faces_edge[i][0];
		f2 = im->faces_edge[i][1];
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
	om->ringf = MEM_callocN(sizeof(LinkNode *) * im->num_faces, "GFListFaces");
}

void freeOutputMesh(OutputMesh *om)
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

DerivedMesh *makeResultMesh(LaplacianSystem *sys)
{
	int i;
	double start_time;
	MVert *verts;
	OutputMesh *om = &sys->output_mesh;
	DerivedMesh *ret;

	start_time = PIL_check_seconds_timer();

	sys->rng = BLI_rng_new(PIL_check_seconds_timer_i());

	initOutputMesh(om, &sys->input_mesh);
	sys->gfsys[0]->id = GFSYS1;
	sys->gfsys[1]->id = GFSYS2;

	makeFeatureEdges(om, &sys->input_mesh);
	computeFlowLines(sys);
	deleteDegenerateVerts(om);
	//hideEdgesOnFaces(om, &sys->input_mesh);
	//makeNormals(om);

	if (!om->num_verts)
		return NULL;

	ret = CDDM_new(om->num_verts, om->num_edges, 0, //0, 0);
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
	BLI_rng_free(sys->rng);

	CDDM_recalc_tessellation(ret);
	//CDDM_calc_edges_tessface(ret);
	//ret->dirty |= DM_DIRTY_NORMALS;

	printf("Mesh generation time: %f ms\n", (PIL_check_seconds_timer() - start_time) * 1000.0f);

	return ret;
}