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
#ifdef QR_GENERATE
#include "MEM_guardedalloc.h"

#include "DNA_modifier_types.h"

#include "BLI_utildefines.h"
#include "BLI_math.h"
#include "BLI_heap.h"
#include "BLI_rand.h"

#include "BKE_mesh_mapping.h"

#include "quadremesh_util.h"

#define QR_GFLINE_QSIZE 64
#define QR_LINELIMIT 100000
//#define QR_SHOWQUERIES

typedef struct GFLine {
	struct GradientFlowSystem *gfsys;
	GFPoint end, seed;
	GFPoint *oldp, *lastchkp;
	GFSegment lastchks;
	int d;

	float lastchklen;

	float qlen;
	int num_q;
	GFSegment q[QR_GFLINE_QSIZE];
} GFLine;


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
static int queryDirection(GradientFlowSystem *gfsys, GFPoint *in_p, int in_f, float in_dir[3],
						  float dist, float maxdist, bool make_seed)
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

		oldf = in_f;
		in_f = getOtherFaceAdjacentToEdge(im, in_f, newp.e);
		memcpy(&oldp, &newp, sizeof(GFPoint));
		prevlen = actlen;
	}

	return 2;
}

static bool checkPoint(GradientFlowSystem *gfsys, float in_oldco[3], GFPoint *in_p, float dist, float maxdist)
{
	bool make_seed = BLI_rng_get_float(gfsys->sys->rng) < QR_SEEDPROB;

	int d, f;
	float seg[3], dir[3], no[3];
	InputMesh *im = &gfsys->sys->input_mesh;

	sub_v3_v3v3(seg, in_oldco, in_p->co);

	if (len_squared_v3(seg) < FLT_EPSILON)
		return true;

	if (in_p->type == eEdge) {
		getNormalAtEdge(no, im, in_p->e);

		if ((f = im->faces_edge[in_p->e][0]) == QR_NO_FACE)	
			if ((f = im->faces_edge[in_p->e][1]) == QR_NO_FACE)
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

	if (line->d == 0 && line->end.id > line->seed.id + 3 &&
		line->end.type == eFace && line->seed.type == eFace &&
		line->end.f == line->seed.f)
	{
		addQREdgeToFace(&line->gfsys->sys->output_mesh,
						&line->gfsys->sys->input_mesh,
						line->gfsys->id, line->end.f,
						line->end.id, line->seed.id);
		return false;
	}

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
	const float chklen = QR_SAMPLING_RATE;
	const float maxdist = line->gfsys->sys->qmd->max_line_dist;
	const float seeddist = maxdist * QR_MAXDIST_TO_SEEDDIST;//(1.0f + BLI_rng_get_float(line->gfsys->sys->rng) * (QR_MAXDIST_TO_SEEDDIST - 1.0f));

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
	QuadRemeshSystem *sys = line->gfsys->sys;
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

void computeFlowLines(QuadRemeshSystem *sys) {
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

#endif
