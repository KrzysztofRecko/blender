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

/** \file blender/modifiers/intern/MOD_quadremesh_geom.h
 *  \ingroup modifiers
 */

#ifndef __MOD_QUADREMESH_GEOM_H__
#define __MOD_QUADREMESH_GEOM_H__

#ifdef WITH_OPENNL

#include "ONL_opennl.h"

#include "BKE_mesh_mapping.h"
#include "BLI_heap.h"
#include "BLI_linklist.h"
#include "BKE_cdderivedmesh.h"
#include "BLI_memarena.h"

#define QR_SAMPLING_RATE 0.03f
#define QR_MINDIST 0.04f
#define QR_SEEDPROB 0.75f
#define QR_SEEDDIST 0.08f

typedef int MEdgeID;
typedef int QRVertID;

typedef struct QRDiskLink {
	struct QRDiskLink *next, *prev, *brother;
	MEdgeID e;
	QRVertID v;
	float ang;
	bool poly_on_right;
} QRDiskLink;

typedef struct QRVert {
	QRDiskLink *link;
	float co[3], no[3], vec[3];
	int num_links;
} QRVert;

typedef struct QREdgeLink {
	struct QREdgeLink *next;
	QRVertID v;
	QRDiskLink *elink;
	float dist;
} QREdgeLink;

typedef struct QREdge {
	QREdgeLink *v1, *v2;
	float dir[3], orig[3];
	int num_links;
} QREdge;

typedef enum {
	eVert,
	eEdge,
	eFace
} GFPointType;

typedef struct GFPoint {
	float co[3];
	int f, v, e;
	QRVertID id;
	GFPointType type;
} GFPoint;

/* GradientFlowSysten, one gfsys for every gradient field */
typedef struct GradientFlowSystem {
	struct Heap *seeds;

	LinkNode **ringf;
	//float *h;
	float(*gf)[3];				/* Gradient Field */

	struct LaplacianSystem *sys;
} GradientFlowSystem;

typedef struct GFLine {
	GradientFlowSystem *gfsys;
	GFPoint end, seed, lastchkp;
	GFPoint *oldp;
	int d;

	float lastchklen;

	float qlen;
	int num_q;
	GFPoint q[10];
} GFLine;

typedef struct InputMesh {
	int num_verts, num_edges, num_faces, num_features;
	float(*co)[3];					/* Original vertex coordinates */
	float(*no)[3];					/* Original face normal */

	unsigned int(*faces)[3];		/* Copy of MFace (tessface) v1-v3, v2-v4 */
	unsigned int(*edges)[2];		/* Copy of edges v1-v2 */
	unsigned int(*faces_edge)[2];	/* Faces by edges  */

	int *ringf_indices;				/* Indices of faces per vertex */
	int *ringe_indices;				/* Indices of edges per vertex */
	MeshElemMap *ringf_map;			/* Map of faces per vertex */
	MeshElemMap *ringe_map;			/* Map of edges per vertex */

	int *constraints;				/* Feature points constraints*/
	float *weights;					/* Feature points weights*/
} InputMesh;

typedef struct OutputMesh {
	MemArena *memarena;
	QRVert *verts;
	QRVertID *vonvs;
	QREdge *ringe;                  /* QREdges per original edges */
	LinkNode **ringf[2];               /* Lists of QREdge per original faces */
	
	int num_verts, alloc_verts;
	int num_edges;
} OutputMesh;

typedef struct LaplacianSystem {
	bool has_solution;

	struct QuadRemeshModifierData *qmd;
	InputMesh input_mesh;
	OutputMesh output_mesh;
	GradientFlowSystem *gfsys[2];
	
	float *U_field;					/* Initial scalar field*/

	NLContext *context;				/* System for solve general implicit rotations */
} LaplacianSystem;

void freeOutputMesh(OutputMesh *om);
DerivedMesh *makeResultMesh(LaplacianSystem *sys);

#endif /*openNl*/
#endif /*__MOD_QUADREMESH_GEOM_H__*/