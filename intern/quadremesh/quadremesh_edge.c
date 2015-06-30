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
#include <float.h>

#include "BLI_utildefines.h"
#include "BLI_memarena.h"
#include "BLI_linklist.h"
#include "BLI_math.h"

#include "quadremesh_util.h"

static QREdgeLink *newQREdgeLink(OutputMesh *om, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl;
	newl = BLI_memarena_alloc(om->memarena, sizeof(QREdgeLink));

	newl->v = in_v;
	newl->dist = in_dist;
	newl->elink = NULL;
	newl->next = NULL;
	newl->gfsysid = 0;

	return newl;
}

static void appendOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl;
	newl = newQREdgeLink(om, in_v, in_dist);

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
	QREdgeLink *newl;
	newl = newQREdgeLink(om, in_v, in_dist);

	newl->next = in_e->v1;
	in_e->v1 = newl;

	if (in_e->v2 == NULL)
		in_e->v2 = newl;
}

static void insertAfterOnQREdge(OutputMesh *om, QREdge *in_e, QREdgeLink *in_l, QRVertID in_v, float in_dist)
{
	QREdgeLink *newl;
	newl = newQREdgeLink(om, in_v, in_dist);

	newl->next = in_l->next;
	newl->gfsysid = in_l->gfsysid;

	if (in_l->elink) {
		unlinkVerts(om, in_l->elink);
		in_l->elink = linkVerts(om, in_l->v, in_v);
		newl->elink = linkVerts(om, in_v, in_l->next->v);
	}

	in_l->next = newl;
}

void insertOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_vid)
{
	float tmp, vec[3];
	QREdgeLink *it;

	if (in_e->num_links == 0) {
		appendOnQREdge(om, in_e, in_vid, 0.0f);
		copy_v3_v3(in_e->orig, om->verts[in_vid].co);
	}
	else if (in_e->num_links == 1) {
		sub_v3_v3v3(in_e->dir, om->verts[in_vid].co, in_e->orig);
		tmp = len_v3(in_e->dir);

		if (tmp < FLT_EPSILON)
			return;

		appendOnQREdge(om, in_e, in_vid, tmp);
		mul_v3_fl(in_e->dir, 1.0f / tmp);
	}
	else {
		sub_v3_v3v3(vec, om->verts[in_vid].co, in_e->orig);
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
			//if (IS_EQF(it->dist, tmp) || IS_EQF(it->next->dist, tmp)) return;
			insertAfterOnQREdge(om, in_e, it, in_vid, tmp);
		}
	}
	in_e->num_links++;
}

void linkOnQREdge(OutputMesh *om, GFSysID sys_id, QREdge *in_e, QRVertID in_v1, QRVertID in_v2)
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

QREdge *addQREdgeToFace(OutputMesh *om, InputMesh *im, GFSysID sys_id, int in_f, QRVertID in_v1, QRVertID in_v2)
{
	float isection[3];
	QREdge *e, *newe;
	QRVertID newv;
	LinkNode *it;

	newe = BLI_memarena_calloc(om->memarena, sizeof(QREdge));

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

#endif
