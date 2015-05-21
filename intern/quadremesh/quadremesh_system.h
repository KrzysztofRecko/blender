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

#ifndef __QUADREMESH_SYSTEM_H__
#define __QUADREMESH_SYSTEM_H__

#include "quadremesh_input.h"

typedef struct OutputMesh {
	struct MemArena *memarena;
	struct FlowLine **lines;
	QRVert **verts;
	QRVert **vonvs;
	struct QREdge *ringe;                  /* QREdges per original edges */
	struct LinkNode **ringf;        /* Lists of QREdge per original faces */
	
	int num_lines, alloc_lines;
	int num_verts, alloc_verts;
	int num_edges;
} OutputMesh;

typedef struct QuadRemeshSystem {
	struct QuadRemeshModifierData *qmd;

	InputMesh input_mesh;
	
	bool has_solution, is_alloc;
	float(*cf)[4][3];               /* Cross field on mesh faces */
	int *singularities;
	float *U_field;					/* Initial scalar field */

	OutputMesh output_mesh;
	struct RNG *rng;

	struct DerivedMesh *cache_mesh;
} QuadRemeshSystem;

QRVert *addVert(OutputMesh *om, InputMesh *im, QRVert *in_v);

#endif