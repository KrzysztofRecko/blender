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

typedef int GFEdgeID;
typedef int MEdgeID;
typedef int MVertID;

typedef struct QREdgeLink {
	struct QREdgeLink *next, *prev, *brother;
	MEdgeID e;
	MVertID v;
	float ang;
	bool poly_on_right;
} QREdgeLink;

typedef struct QREdgeLinkList {
	QREdgeLink *link;
	float vec[3];
	int num_links;
} QREdgeLinkList;

typedef enum {
	eSeedVert,
	eSeedFace
} GFSeedType;

typedef struct GFSeed {
	float co[3];
	int val;
	GFSeedType type;
} GFSeed;

typedef struct GFEdge {
	LinkNode *elinks;
	MVertID v1, v2;	    /* End vert indices */
} GFEdge;

typedef struct GFLine {
	MVertID end, seed;
	int d;
	float oldco[3];

	float lastchk[3];
	float lastchklen;

	float qco[10][3];
	int qf[10];
	float qlen;
	int num_q;
} GFLine;

typedef struct InputMesh {
	int num_verts, num_edges, num_faces, num_features;
	float(*co)[3];					/* Original vertex coordinates */
	float(*no)[3];					/* Original face normal */

	unsigned int(*faces)[3];		/* Copy of MFace (tessface) v1-v3, v2-v4 */

	unsigned int(*edges)[2];		/* Copy of edges v1-v2 */
	unsigned int(*faces_edge)[2];	/* Faces by edges  */

	int *ringf_indices;				/* Indices of faces per vertex */
	int *ringv_indices;				/* Indices of neighbors(vertex) per vertex */
	int *ringe_indices;				/* Indices of edges per vertex */
	MeshElemMap *ringf_map;			/* Map of faces per vertex */
	//MeshElemMap *ringv_map;			/* Map of vertex per vertex  - unused right now */
	MeshElemMap *ringe_map;			/* Map of edges per vertex */

	int *constraints;				/* Feature points constraints*/
	float *weights;					/* Feature points weights*/
} InputMesh;

typedef struct OutputMesh {
	int totvert, allocvert;
	MVert *verts;
	QREdgeLinkList *vlinks;
	MemArena *memarena;

	int totedge, allocedge;
	MEdge *edges;
	int totloop, allocloop;
	MLoop *loops;
	int totpolys, allocpolys;
	MPoly *polys;
} OutputMesh;

typedef struct LaplacianSystem LaplacianSystem;

#define QR_SAMPLING_RATE 0.03f
#define QR_MINDIST 0.04f
#define QR_SEEDPROB 0.75f
#define QR_SEEDDIST 0.08f

/* GradientFlowSysten, one gfsys for every gradient field */
typedef struct GradientFlowSystem {
	LinkNode **ringf_list;			/* Array list of of GFEdge per original face */
	MemArena *memarena;

	struct Heap *heap_seeds;

	float *hfunction;
	float(*gfield)[3];				/* Gradient Field g1 */

	LaplacianSystem *sys;
} GradientFlowSystem;

typedef struct LaplacianSystem {
	bool command_compute_flow;
	bool has_solution;
	bool command_remesh;

	struct QuadRemeshModifierData *qmd;
	InputMesh input_mesh;
	OutputMesh output_mesh;
	
	float *U_field;					/* Initial scalar field*/

	float(*gf1)[3], (*gf2)[3];		/* Gradient Fields g1 and g2 per face*/
	float *h1, *h2;					/* Sampling distance functions h1 and h2 */
	float h;

	NLContext *context;				/* System for solve general implicit rotations */
	GradientFlowSystem *gfsys1, *gfsys2;
} LaplacianSystem;

GradientFlowSystem *newGradientFlowSystem(LaplacianSystem *sys, float *mhfunction, float(*mgfield)[3]);
void deleteGradientFlowSystem(GradientFlowSystem *gfsys);
void generateMesh(LaplacianSystem *sys);
void freeOutputMesh(OutputMesh *om);
DerivedMesh *getResultMesh(LaplacianSystem *sys);

#endif /*openNl*/
#endif /*__MOD_QUADREMESH_GEOM_H__*/