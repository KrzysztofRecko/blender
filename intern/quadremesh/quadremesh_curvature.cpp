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
 * Contributor(s): Krzysztof Recko
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

extern "C" {

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "BKE_mesh_mapping.h"

#include "quadremesh_system.h"

}

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace Eigen;

static void alignAtEdge(QuadRemeshSystem *sys, int in_e)
{
	int f1, f2;
	float ev[3], vec[3], tmp[3], best[3], best_dot = -1.0f;
	float ang1, ang2, ang_diff, power_sum, ang_mid;
	InputMesh *im = &sys->input_mesh;

	f1 = im->faces_edge[in_e][0];
	f2 = im->faces_edge[in_e][1];

	if (f1 == QR_NO_FACE || f2 == QR_NO_FACE)
		return;

	//if (sys->gfsys[0]->h[f1] < FLT_EPSILON && sys->gfsys[0]->h[f2] < FLT_EPSILON)
	//	return;

	sub_v3_v3v3(ev, im->co[im->edges[in_e][0]],
					im->co[im->edges[in_e][1]]);
	normalize_v3(ev);

	ang1 = angle_signed_on_axis_normalized_fast_v3v3_v3(ev, sys->cf[f1][0], im->no[f1]);
	while (ang1 > M_PI_4)
		ang1 -= M_PI_2;

	ang2 = angle_signed_on_axis_normalized_fast_v3v3_v3(ev, sys->cf[f2][0], im->no[f2]);
	while (ang2 > M_PI_4)
		ang2 -= M_PI_2;

	ang_diff = ang1 - ang2;

	/*power_sum = sys->gfsys[0]->h[f1] + sys->gfsys[0]->h[f2];
	ang_mid = ang_diff * (sys->gfsys[0]->h[f2] / power_sum);

	rotate_normalized_v3_v3v3fl(tmp, sys->gfsys[0]->gf[f1], im->no[f1], ang_mid);
	copy_v3_v3(sys->gfsys[0]->gf[f1], tmp);
	rotate_normalized_v3_v3v3fl(tmp, sys->gfsys[0]->gf[f2], im->no[f2], -(ang_diff - ang_mid));
	copy_v3_v3(sys->gfsys[0]->gf[f2], tmp);


	if (sys->gfsys[0]->h[f2] > sys->gfsys[0]->h[f1]) {
		sys->gfsys[0]->h[f1] = power_sum * 0.05f;
	}
	else {
		sys->gfsys[0]->h[f2] = power_sum * 0.05f;
	}*/
}

static void getCurvatures(QuadRemeshSystem *sys)
{
	int f, e;
	float ev[3], angle, vec[3];
	InputMesh *im = &sys->input_mesh;
	Vector3f edge_vec;
	Matrix3f edge_c, *face_curvatures = new Matrix3f[im->num_faces];

	for (f = 0; f < im->num_faces; f++) {
		face_curvatures[f].setZero();
	}

	for (e = 0; e < im->num_edges; e++) {
		edge_c.setZero();

		sub_v3_v3v3(ev, im->co[im->edges[e][0]],
						im->co[im->edges[e][1]]);
		
		normalize_v3(ev);
		angle = angle_normalized_v3v3(im->no[im->faces_edge[e][0]],
									  im->no[im->faces_edge[e][1]]);
		sub_v3_v3v3(vec, im->co[getVertexOppositeToEdge(&sys->input_mesh, e, im->faces_edge[e][0])], im->co[im->edges[e][0]]);
		if (dot_v3v3(vec, im->no[im->faces_edge[e][1]]) > 0.0f) {
			angle = -angle;
		}

		edge_vec << ev[0], ev[1], ev[2];
		edge_c = angle * (edge_vec * edge_vec.transpose());

		if (im->faces_edge[e][0] != QR_NO_FACE)
			face_curvatures[im->faces_edge[e][0]] += edge_c;
		if (im->faces_edge[e][1] != QR_NO_FACE)
			face_curvatures[im->faces_edge[e][1]] += edge_c;
	}

	for (f = 0; f < im->num_faces; f++) {
		//std::cout << "vert_c: " << std::endl << vert_c << std::endl;
		EigenSolver<Matrix3f> es(face_curvatures[f]);
		Matrix3f eigenv = es.pseudoEigenvectors();
		Vector3cf eigenvalues = es.eigenvalues();
		//std::cout << "Eigenvalues: " << std::endl << es.eigenvalues() << std::endl;
		//std::cout << "Pseudo-eigenvectors: " << std::endl << es.pseudoEigenvectors() << std::endl;

		int max = MAX3_PAIR(eigenvalues(0).real(), eigenvalues(1).real(), eigenvalues(2).real(), 0, 1, 2);
		if (eigenvalues(max).real() > FLT_EPSILON) {
			copy_v3_v3(sys->cf[f][0], eigenv.col(max).data());
			normalize_v3(sys->cf[f][0]);
			//sys->gfsys[0]->h[f] = SQUARE(eigenvalues(max).real());
		}
		else {
			sub_v3_v3v3(sys->cf[f][0], im->co[im->faces[f][0]], im->co[im->faces[f][1]]);
			normalize_v3(sys->cf[f][0]);
			//sys->gfsys[0]->h[f] = 0.0f;
		}

		copy_v3_v3(ev, sys->cf[f][0]);
		for (int i = 1; i < 4; i++) {
			cross_v3_v3v3(sys->cf[f][i], ev, sys->input_mesh.no[f]);
			copy_v3_v3(ev, sys->cf[f][i]);
		}
	}

	delete[] face_curvatures;

	//for (int i = 0; i < 20; i++)
		//for (int e = 0; e < im->num_edges; e++) {
			//alignAtEdge(sys, e);
		//}
}

void getPrincipalCurvatures(QuadRemeshSystem *sys)
{
	int i;
	InputMesh *im = &sys->input_mesh;

	sys->has_solution = true;

	sys->U_field = static_cast<float*>(MEM_mallocN(sizeof(float) * im->num_verts, "QuadRemeshUField"));

	for (i = 0; i < im->num_verts; i++) {
		sys->U_field[i] = 0.0f;
	}

	sys->cf = static_cast<float(*)[4][3]>(MEM_mallocN(sizeof(float[4][3]) * im->num_faces, __func__));
	//sys->gfsys[0]->h = static_cast<float*>(MEM_mallocN(sizeof(float) * im->num_faces, __func__));
	//sys->gfsys[1]->h = static_cast<float*>(MEM_mallocN(sizeof(float) * im->num_faces, __func__));

	getCurvatures(sys);

	sys->is_alloc = true;
}