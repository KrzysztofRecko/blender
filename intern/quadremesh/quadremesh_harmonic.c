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
#include "DNA_meshdata_types.h"

#include "BLI_utildefines.h"
#include "BLI_math.h"

#include "BKE_deform.h"

#include "MOD_util.h"

#include "quadremesh_system.h"

#if 0
static void initLaplacianMatrix(InputMesh *im)
{
	float v1[3], v2[3], v3[3];
	float w2, w3;
	int j, fi;
	unsigned int idv1, idv2, idv3;

	for (fi = 0; fi < im->num_faces; fi++) {
		const unsigned int *vidf = im->faces[fi];

		idv1 = vidf[0];
		idv2 = vidf[1];
		idv3 = vidf[2];

		for (j = 0; j < 3; j++) {
			idv1 = vidf[j];
			idv2 = vidf[(j + 1) % 3];
			idv3 = vidf[(j + 2) % 3];

			copy_v3_v3(v1, im->co[idv1]);
			copy_v3_v3(v2, im->co[idv2]);
			copy_v3_v3(v3, im->co[idv3]);

			w2 = cotangent_tri_weight_v3(v3, v1, v2);
			w3 = cotangent_tri_weight_v3(v2, v3, v1);

			if (im->constraints[idv1] == 1) {
				nlMatrixAdd(idv1, idv1, w2 + w3);
			}
			else  {
				nlMatrixAdd(idv1, idv2, -w2);
				nlMatrixAdd(idv1, idv3, -w3);
				nlMatrixAdd(idv1, idv1, w2 + w3);
			}

		}
	}
	
}

static bool computeScalarField(QuadRemeshSystem *sys)
{
	int i;
	InputMesh *im = &sys->input_mesh;
	
	nlNewContext();

	nlSolverParameteri(NL_NB_VARIABLES, im->num_verts);
	nlSolverParameteri(NL_SYMMETRIC, NL_FALSE);
	nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
	nlSolverParameteri(NL_NB_ROWS, im->num_verts);
	nlSolverParameteri(NL_NB_RIGHT_HAND_SIDES, 1);
	nlBegin(NL_SYSTEM);
	
	for (i = 0; i < im->num_verts; i++) {
		nlSetVariable(0, i, 0);
	}
		
	nlBegin(NL_MATRIX);
	
	initLaplacianMatrix(im);

	for (i = 0; i < im->num_verts; i++) {
		if (im->constraints[i] == 1) {
			nlRightHandSideSet(0, i, im->weights[i]);
		}
		else {
			nlRightHandSideSet(0, i, 0);
		}
	}
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);
	
	return nlSolveAdvanced(NULL, NL_TRUE);
}

/**
* Compute the gradient fields
*
* xi, xj, xk, are the vertices of the face
* ui, uj, uk, are the values of scalar fields for every vertex of the face
* n is the normal of the face.
* gf1 is the unknown field gradient 1.
* gf2 is the unknown field gradient 2.
*
* |xj - xi|         |uj - ui|
* |xk - xj| * gf1 = |uk - uj|
* |   nf  |         |   0   |
*
* gf2 = cross(n, gf1)
*/
static void computeGradientFields(QuadRemeshSystem * sys)
{
	int fi, i, j, k;
	float val, a[3][3], u[3], inv_a[3][3], gf1[3], g[3], w[3];
	InputMesh *im = &sys->input_mesh;

	for (fi = 0; fi < im->num_faces; fi++) {
		const unsigned int *vidf = im->faces[fi];
		i = vidf[0];
		j = vidf[1];
		k = vidf[2];
		sub_v3_v3v3(a[0], im->co[j], im->co[i]);
		sub_v3_v3v3(a[1], im->co[k], im->co[j]);
		copy_v3_v3(a[2], im->no[fi]);

		/* Correct way*/
		transpose_m3(a);
		u[0] = sys->U_field[j] - sys->U_field[i];
		u[1] = sys->U_field[k] - sys->U_field[j];
		u[2] = 0;
		invert_m3_m3(inv_a, a);
		//mul_v3_m3v3(sys->gf1[fi], inv_a, u);
		mul_v3_m3v3(gf1, inv_a, u);

		/* Project Gradient fields on face*/
		normalize_v3_v3(g, gf1);
		val = dot_v3v3(g, im->no[fi]);
		mul_v3_v3fl(u, im->no[fi], val);
		sub_v3_v3v3(w, g, u);
		normalize_v3_v3(sys->cf[fi][0], w);

		cross_v3_v3v3(g, im->no[fi], sys->cf[fi][0]);
		copy_v3_v3(sys->cf[fi][1], g);
		
		cross_v3_v3v3(g, im->no[fi], sys->cf[fi][1]);
		copy_v3_v3(sys->cf[fi][2], g);

		cross_v3_v3v3(g, im->no[fi], sys->cf[fi][2]);
		copy_v3_v3(sys->cf[fi][3], g);
	}
}

void getHarmonicGradients(QuadRemeshSystem *sys)
{
	int i;
	InputMesh *im = &sys->input_mesh;

#ifdef OPENNL_THREADING_HACK
	modifier_opennl_lock();
#endif

	sys->has_solution = computeScalarField(sys);

	if (sys->has_solution) {
		sys->U_field = MEM_mallocN(sizeof(float) * im->num_verts, "QuadRemeshUField");

		for (i = 0; i < im->num_verts; i++) {
			sys->U_field[i] = nlGetVariable(0, i);
		}

		sys->cf = MEM_mallocN(sizeof(float[4][3]) * im->num_faces, "QuadRemeshGradientField1");

		computeGradientFields(sys);

		sys->is_alloc = true;
	}

#ifdef OPENNL_THREADING_HACK
	modifier_opennl_unlock();
#endif
}

void getUField(QuadRemeshSystem *sys, MDeformVert *dvert, int defgrp_index)
{
	int i;
	float y, mmin = 1000.0f, mmax = -1000.0f;
	MDeformVert *dv = NULL;

	BLI_assert(dvert != NULL);

	if (sys->has_solution) {
		dv = dvert;

		for (i = 0; i < sys->input_mesh.num_verts; i++) {
			mmin = min_ff(mmin, sys->U_field[i]);
			mmax = max_ff(mmax, sys->U_field[i]);
		}

		for (i = 0; i < sys->input_mesh.num_verts; i++) {
			y = (sys->U_field[i] - mmin) / (mmax - mmin);
			//x = y * 60;
			//y = (x % 2 == 0 ? 0.1 : 0.9);
			defvert_add_index_notest(dv, defgrp_index, y);
			dv++;
		}
	}
}

#endif