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

#include "DNA_object_types.h"

#include "BLI_math.h"

#include "BKE_mesh_mapping.h"
#include "BKE_DerivedMesh.h"
#include "BKE_deform.h"

#include "MOD_util.h"

#include "quadremesh_system.h"

static void createFaceRingMap(InputMesh *im)
{
	int i, j, totalr = 0;
	int *index_iter;

	im->ringf_map = MEM_callocN(sizeof(MeshElemMap) * im->num_verts, "DeformRingMap");

	for (i = 0; i < im->num_faces; i++) {
		for (j = 0; j < 3; j++) {
			im->ringf_map[im->faces[i][j]].count++;
			totalr++;
		}
	}

	im->ringf_indices = MEM_callocN(sizeof(int) * totalr, "DeformRingIndex");
	index_iter = im->ringf_indices;
	for (i = 0; i < im->num_verts; i++) {
		im->ringf_map[i].indices = index_iter;
		index_iter += im->ringf_map[i].count;
		im->ringf_map[i].count = 0;
	}

	for (i = 0; i < im->num_faces; i++) {
		for (j = 0; j < 3; j++) {
			MeshElemMap *map = &im->ringf_map[im->faces[i][j]];
			map->indices[map->count++] = i;
		}
	}
}

static void createEdgeRingMap(InputMesh *im)
{
	int i, totalr = 0;
	int *index_iter;

	im->ringe_map = MEM_callocN(sizeof(MeshElemMap) * im->num_verts, "DeformNeighborsMap");

	for (i = 0; i < im->num_edges; i++) {
		im->ringe_map[im->edges[i][0]].count++;
		im->ringe_map[im->edges[i][1]].count++;
		totalr += 2;
	}

	im->ringe_indices = MEM_callocN(sizeof(int) * totalr, "DeformNeighborsIndex");
	index_iter = im->ringe_indices;
	for (i = 0; i < im->num_verts; i++) {
		im->ringe_map[i].indices = index_iter;
		index_iter += im->ringe_map[i].count;
		im->ringe_map[i].count = 0;
	}

	for (i = 0; i < im->num_edges; i++) {
		MeshElemMap *map1 = &im->ringe_map[im->edges[i][0]];
		MeshElemMap *map2 = &im->ringe_map[im->edges[i][1]];
		map1->indices[map1->count++] = i;
		map2->indices[map2->count++] = i;
	}
}

static void createVertRingMap(InputMesh *im)
{
	int i, totalr = 0;
	int *index_iter;

	im->ringv_map = MEM_callocN(sizeof(MeshElemMap) * im->num_verts, "DeformNeighborsMap");

	for (i = 0; i < im->num_edges; i++) {
		im->ringv_map[im->edges[i][0]].count++;
		im->ringv_map[im->edges[i][1]].count++;
		totalr += 2;
	}

	im->ringv_indices = MEM_callocN(sizeof(int) * totalr, "DeformNeighborsIndex");
	index_iter = im->ringv_indices;
	for (i = 0; i < im->num_verts; i++) {
		im->ringv_map[i].indices = index_iter;
		index_iter += im->ringv_map[i].count;
		im->ringv_map[i].count = 0;
	}

	for (i = 0; i < im->num_edges; i++) {
		MeshElemMap *map1 = &im->ringv_map[im->edges[i][0]];
		MeshElemMap *map2 = &im->ringv_map[im->edges[i][1]];
		map1->indices[map1->count++] = im->edges[i][1];
		map2->indices[map2->count++] = im->edges[i][0];
	}
}

static void createFacesByEdge(InputMesh *im){
	int e, i, j, v1, v2;
	unsigned int *vin;

	im->faces_edge = MEM_mallocN(sizeof(int[2]) * im->num_edges, "QuadRemeshFacesEdge");

	for (e = 0; e < im->num_edges; e++) {
		v1 = im->edges[e][0];
		v2 = im->edges[e][1];

		im->faces_edge[e][0] = QR_NO_FACE;
		im->faces_edge[e][1] = QR_NO_FACE;

		for (i = 0, j = 0; i < im->ringf_map[v1].count && j < 2; i++) {
			vin = im->faces[im->ringf_map[v1].indices[i]];
			if (vin[0] == v2 || vin[1] == v2 || vin[2] == v2) {
				im->faces_edge[e][j++] = im->ringf_map[v1].indices[i];
			}
		}
	}
}

void freeInputMesh(InputMesh *im)
{
	MEM_SAFE_FREE(im->faces);
	MEM_SAFE_FREE(im->edges);
	MEM_SAFE_FREE(im->faces_edge);
	MEM_SAFE_FREE(im->co);
	MEM_SAFE_FREE(im->no);
	MEM_SAFE_FREE(im->vno);
	MEM_SAFE_FREE(im->ev);
	MEM_SAFE_FREE(im->constraints);
	MEM_SAFE_FREE(im->weights);

	MEM_SAFE_FREE(im->ringf_indices);
	MEM_SAFE_FREE(im->ringe_indices);
	MEM_SAFE_FREE(im->ringv_indices);
	MEM_SAFE_FREE(im->ringf_map);
	MEM_SAFE_FREE(im->ringe_map);
	MEM_SAFE_FREE(im->ringv_map);

	im->is_alloc = false;
}

static void getInputMeshData(InputMesh *im, DerivedMesh *dm)
{
	int i;
	MVert *arrayvect;
	MEdge *arrayedge;
	MFace *tessface;

	DM_ensure_tessface(dm);

	/* Get vertices */
	im->num_verts = dm->getNumVerts(dm);
	im->co = MEM_mallocN(sizeof(float[3]) * im->num_verts, "QuadRemeshCoordinates");
	arrayvect = dm->getVertArray(dm);
	for (i = 0; i < dm->getNumVerts(dm); i++) {
		copy_v3_v3(im->co[i], arrayvect[i].co);
	}

	/* Get edges */
	im->edges = MEM_mallocN(sizeof(int[2]) * dm->getNumEdges(dm) * 2, "QuadRemeshEdges");
	arrayedge = dm->getEdgeArray(dm);
	for (im->num_edges = 0, i = 0; i < dm->getNumEdges(dm); i++) {
		im->edges[im->num_edges][0] = arrayedge[i].v1;
		im->edges[im->num_edges][1] = arrayedge[i].v2;
		im->num_edges++;
	}

	/* Edge vectors */
	im->ev = MEM_mallocN(sizeof(float[3]) * im->num_edges, "QuadRemeshEdgeVectors");
	for (i = 0; i < im->num_edges; i++) {
		sub_v3_v3v3(im->ev[i], im->co[im->edges[i][1]], im->co[im->edges[i][0]]);
		normalize_v3(im->ev[i]);
	}

	/* Get faces */
	tessface = dm->getTessFaceArray(dm);
	im->faces = MEM_mallocN(sizeof(int[3]) * dm->getNumTessFaces(dm) * 2, "QuadRemeshFaces");
	for (im->num_faces = 0, i = 0; i < dm->getNumTessFaces(dm); i++) {
		im->faces[im->num_faces][0] = tessface[i].v1;
		im->faces[im->num_faces][1] = tessface[i].v2;
		im->faces[im->num_faces][2] = tessface[i].v3;
		im->num_faces++;

		if (tessface[i].v4 != 0) {
			im->faces[im->num_faces][0] = tessface[i].v1;
			im->faces[im->num_faces][1] = tessface[i].v3;
			im->faces[im->num_faces][2] = tessface[i].v4;
			im->num_faces++;

			im->edges[im->num_edges][0] = tessface[i].v1;
			im->edges[im->num_edges][1] = tessface[i].v3;
			im->num_edges++;
		}
	}

	/* Compute face and vertex normals */
	im->no = MEM_callocN(sizeof(float[3]) * im->num_faces, "QuadRemeshNormals");
	for (i = 0; i < im->num_faces; i++) {
		normal_tri_v3(im->no[i], im->co[im->faces[i][0]],
					             im->co[im->faces[i][1]],
					             im->co[im->faces[i][2]]);
	}

	im->vno = MEM_callocN(sizeof(float[3]) * im->num_verts, "QuadRemeshNormals");
	for (i = 0; i < im->num_verts; i++) {
		zero_v3(im->vno[i]);
	}
	for (i = 0; i < im->num_faces; i++) {
		add_v3_v3(im->vno[im->faces[i][0]], im->no[i]);
		add_v3_v3(im->vno[im->faces[i][1]], im->no[i]);
		add_v3_v3(im->vno[im->faces[i][2]], im->no[i]);
	}
	for (i = 0; i < im->num_verts; i++) {
		normalize_v3(im->vno[i]);
	}
}

static void getFeatures(InputMesh *im, MDeformVert *dvert, int defgrp_index)
{
	int i;
	float wpaint;
	MDeformVert *dv = NULL;
	MDeformWeight *d;

	BLI_assert(dvert != NULL);

	/* Get features */
	im->constraints = MEM_callocN(sizeof(int) * im->num_verts, __func__);
	im->weights = MEM_callocN(sizeof(float) * im->num_verts, __func__);

	dv = dvert;
	for (i = 0, im->num_features = 0; i < im->num_verts; i++, dv++) {
		d = defvert_find_index(dv, defgrp_index);
		wpaint = defvert_find_weight(dv, defgrp_index);

		if (d && (wpaint < 0.19 || wpaint > 0.89)) {
			im->constraints[i] = 1;
			im->weights[i] = -1.0f + wpaint * 2.0f;
			im->num_features++;
		}
	}
}

void getInput(QuadRemeshSystem *sys, Object *ob, DerivedMesh *dm)
{
	int defgrp_index;
	MDeformVert *dvert;

	//modifier_get_vgroup(ob, dm, sys->qmd->anchor_grp_name, &dvert, &defgrp_index);

	getInputMeshData(&sys->input_mesh, dm);
	//getFeatures(&sys->input_mesh, dvert, defgrp_index);

	createFaceRingMap(&sys->input_mesh);
	createEdgeRingMap(&sys->input_mesh);
	createVertRingMap(&sys->input_mesh);
	createFacesByEdge(&sys->input_mesh);

	sys->input_mesh.is_alloc = true;
}

void getNormalAt(float r_no[3], InputMesh *im, int in_vef, GFPointType in_type)
{
	int i;

	if (in_type == eVert) {
		BLI_assert(in_vef <= im->num_verts);

		zero_v3(r_no);
		for (i = 0; i < im->ringf_map[in_vef].count; i++)
			add_v3_v3(r_no, im->no[im->ringf_map[in_vef].indices[i]]);

		mul_v3_fl(r_no, 1.0f / (float)im->ringf_map[in_vef].count);
	}
	else if (in_type == eEdge) {
		BLI_assert(in_vef <= im->num_edges);
		BLI_assert(!(im->faces_edge[in_vef][0] == QR_NO_FACE &&
				     im->faces_edge[in_vef][1] == QR_NO_FACE));

		if (im->faces_edge[in_vef][0] == QR_NO_FACE)
			copy_v3_v3(r_no, im->no[im->faces_edge[in_vef][1]]);
		else if (im->faces_edge[in_vef][1] == QR_NO_FACE)
			copy_v3_v3(r_no, im->no[im->faces_edge[in_vef][0]]);
		else {
			add_v3_v3v3(r_no, im->no[im->faces_edge[in_vef][0]],
							  im->no[im->faces_edge[in_vef][1]]);
			mul_v3_fl(r_no, 0.5f);
		}
	}
	else { /* eFace */
		BLI_assert(in_vef <= im->num_faces);

		copy_v3_v3(r_no, im->no[in_vef]);
	}
}

int getEdgeFromVerts(InputMesh *im, int v1, int v2)
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

int getOtherFaceAdjacentToEdge(InputMesh *im, int in_f, int in_e)
{
	if (im->faces_edge[in_e][0] == in_f) {
		return im->faces_edge[in_e][1];
	}

	return im->faces_edge[in_e][0];
}

int getVertexOppositeToEdge(InputMesh *im, int in_e, int in_f)
{
	int i;
	for (i = 0; i < 3; i++) {
		if (im->faces[in_f][i] != im->edges[in_e][0] &&
			im->faces[in_f][i] != im->edges[in_e][1])
		{
			return im->faces[in_f][i];
		}
	}
	return 0;
}