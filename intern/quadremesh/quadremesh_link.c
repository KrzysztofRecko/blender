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
#include "BLI_utildefines.h"
#include "BLI_memarena.h"
#include "BLI_math.h"

#include "quadremesh_system.h"

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
		l->ang = angle_signed_on_axis_normalized_fast_v3v3_v3(om->verts[in_a].vec, vec, om->verts[in_a].no);

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

QRDiskLink *linkVerts(OutputMesh *om, QRVertID in_v1, QRVertID in_v2)
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

void deleteLink(OutputMesh *om, QRVert *ll, QRDiskLink *l)
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

void unlinkVerts(OutputMesh *om, QRDiskLink *l)
{
	QRVert *ll1, *ll2;

	ll2 = &om->verts[l->v];
	ll1 = &om->verts[l->brother->v];

	deleteLink(om, ll1, l);
	deleteLink(om, ll2, l->brother);

	om->num_edges--;
}

void dissolveVert(OutputMesh *om, QRVertID in_v)
{
	QRVertID a, b;
	QRDiskLink *al, *bl;

	a = om->verts[in_v].link->v;
	al = om->verts[in_v].link->brother;
	b = om->verts[in_v].link->next->v;
	bl = om->verts[in_v].link->next->brother;

	if (a != b) {
		deleteLink(om, &om->verts[in_v], om->verts[in_v].link);
		deleteLink(om, &om->verts[in_v], om->verts[in_v].link);

		if (getLink(om, a, b)) {
			deleteLink(om, &om->verts[a], al);
			deleteLink(om, &om->verts[b], bl);

			om->num_edges -= 2;
		}
		else {
			al->v = b;
			al->brother = bl;
			bl->v = a;
			bl->brother = al;
			
			om->num_edges--;
		}
	}
	else {
		deleteLink(om, &om->verts[in_v], om->verts[in_v].link);
		deleteLink(om, &om->verts[in_v], om->verts[in_v].link);
		deleteLink(om, &om->verts[a], al);
		deleteLink(om, &om->verts[b], bl);

		om->num_edges -= 2;
	}
}

#endif
