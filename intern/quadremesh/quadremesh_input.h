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

#ifndef __QUADREMESH_INPUT_H__
#define __QUADREMESH_INPUT_H__

#include "quadremesh_util.h"

typedef struct QuadRemeshSystem QuadRemeshSystem;
typedef struct Object Object;
typedef struct DerivedMesh DerivedMesh;

typedef struct InputMesh {
	bool is_alloc;

	int num_verts, num_edges, num_faces, num_features;
	float(*co)[3];					/* Original vertex coordinates */
	float(*no)[3];					/* Original face normal */
	float(*vno)[3];
	float(*ev)[3]; /* Normalized edge vectors */

	unsigned int(*faces)[3];		/* Copy of MFace (tessface) v1-v3, v2-v4 */
	unsigned int(*edges)[2];		/* Copy of edges v1-v2 */
	unsigned int(*faces_edge)[2];	/* Faces by edges  */

	int *ringf_indices;				/* Indices of faces per vertex */
	int *ringe_indices;				/* Indices of edges per vertex */
	int *ringv_indices;				/* Indices of verts per vertex */
	struct MeshElemMap *ringf_map;  /* Map of faces per vertex */
	struct MeshElemMap *ringe_map;	/* Map of edges per vertex */
	struct MeshElemMap *ringv_map;	/* Map of verts per vertex */

	int *constraints;				/* Feature points constraints*/
	float *weights;					/* Feature points weights*/
} InputMesh;

void getInput(QuadRemeshSystem *sys, Object *ob, DerivedMesh *dm);
void freeInputMesh(InputMesh *im);

void getNormalAt(float r_no[3], InputMesh *im, int in_vef, GFPointType in_type);
int getEdgeFromVerts(InputMesh *im, int v1, int v2);
int getOtherFaceAdjacentToEdge(InputMesh *im, int in_f, int in_e);
int getVertexOppositeToEdge(InputMesh *im, int in_e, int in_f);

#endif