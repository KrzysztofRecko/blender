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

#include "MEM_guardedalloc.h"

#include "DNA_modifier_types.h"

#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_heap.h"
#include "BLI_linklist.h"

#include "BKE_mesh_mapping.h"
#include "BKE_DerivedMesh.h"

#include "quadremesh_system.h"

#ifdef QR_GENERATE
/* ADDING STUFF TO OUTPUT MESH */

QRVertID addVert(OutputMesh *om, float in_co[3], float in_no[3])
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


void addGFPoint(InputMesh *im, OutputMesh *om, GFPoint *in_p)
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
#endif

/*
* List of vertices from original mesh with special features (edge dihedral angle less that 90) to be preserves
* return the size of array
*/
static int *findFeaturesOnMesh(InputMesh *im, int size[2])
{
	int i, f1, f2, total, *listverts, *listdest;
	float angle;

	listverts = MEM_callocN(sizeof(int) * im->num_verts, __func__);
	listdest = NULL;
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
		listdest = MEM_mallocN(sizeof(int) * total, __func__);
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

#ifdef QR_GENERATE
void addSeedToQueue(Heap *aheap, float in_co[3], GFPointType in_type, int in_val, float weight)
{
	GFPoint *seed;

	seed = MEM_mallocN(sizeof(GFPoint), __func__);

	copy_v3_v3(seed->co, in_co);
	seed->type = in_type;
	seed->e = seed->f = seed->v = seed->id = -1;
	if (in_type == eVert)
		seed->v = in_val;
	else
		seed->f = in_val;

	BLI_heap_insert(aheap, weight, seed);
}

void getInitialSeeds(QuadRemeshSystem *gfsys)
{
	int i, *lverts, sizeverts[2];

	lverts = findFeaturesOnMesh(&gfsys->sys->input_mesh, sizeverts);
	
	for (i = 0; i < sizeverts[0]; i++)
		addSeedToQueue(gfsys->seeds, gfsys->sys->input_mesh.co[lverts[i]], eVert, lverts[i], 0.0f);
	
	MEM_SAFE_FREE(lverts);
}

/* SEED QUEUE */

GFPoint *getTopSeedFromQueue(struct Heap *aheap)
{
	return (GFPoint*)BLI_heap_popmin(aheap);
}

bool isectSegmentWithOthersOnFace(OutputMesh *om, GFSysID sys_id, float in_a[3], float in_b[3], int in_f)
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

bool isectPointWithQREdge(OutputMesh *om, GFSysID sys_id, float in_co[3], int in_e)
{
	float vec[3], tmp;
	QREdge *e;
	QREdgeLink *it;

	e = &om->ringe[in_e];

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
#endif

void freeQuadRemeshSystem(QuadRemeshSystem *sys)
{
	if (sys->input_mesh.is_alloc) {
		freeInputMesh(&sys->input_mesh);
	}

	if (sys->is_alloc) {
		MEM_SAFE_FREE(sys->U_field);
		MEM_SAFE_FREE(sys->cf);
		MEM_SAFE_FREE(sys->singularities);
	}

	if (sys->cache_mesh) {
		sys->cache_mesh->release(sys->cache_mesh);
	}
	
	MEM_SAFE_FREE(sys);
}

void initQuadRemeshSystem(QuadRemeshModifierData *qmd)
{
	QuadRemeshSystem *sys;

	sys = MEM_callocN(sizeof(QuadRemeshSystem), "QuadRemeshSystem - Cache");

	sys->qmd = qmd;
	qmd->cache_system = sys;
}