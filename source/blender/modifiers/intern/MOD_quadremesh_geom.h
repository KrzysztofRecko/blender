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

#define MOD_QUADREMESH_ALLOC_BLOCK 256
#define MOD_QUADREMESH_MIN_LEN 0.000001f

#ifdef WITH_OPENNL

#include "ONL_opennl.h"

#include "BKE_mesh_mapping.h"
#include "BLI_heap.h"
#include "BLI_linklist.h"
#include "BKE_cdderivedmesh.h"

typedef int GFEdgeID;
typedef int GFVertID;

typedef struct GFVert {
	float co[3];		/* Vert position */
	GFVertID e[4];      /* 4 neighboring vertices: 2 for each gfsystem */
} GFVert;

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
	GFVertID v1, v2;	    /* Vert indices */
} GFEdge;

typedef struct GFLine {
	GFVertID end, seed;
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
	//MeshElemMap *ringv_map;			/* Map of vertex per vertex */
	MeshElemMap *ringe_map;			/* Map of edges per vertex */

	int *constraints;				/* Feature points constraints*/
	float *weights;					/* Feature points weights*/
} InputMesh;

typedef struct OutputMesh {
	int totvert, allocvert;
	GFVert *mvert;	/* array of verts */
	int totedge, allocedge;
	GFEdge *medge;

	int totloop, allocloop;
	MLoop *loops;
	int totpolys, allocpolys;
	MPoly *polys;
} OutputMesh;

typedef struct LaplacianSystem LaplacianSystem;

/* GradientFlowSysten, one gfsys for every gradient field */
typedef struct GradientFlowSystem {
	int totedge, allocedge;
	GFEdge *medge;	/* array of edges */
	LinkNode **ringf_list;			/* Array list of of GFEdge per original face */

	struct Heap *heap_seeds;

	float *hfunction;
	float(*gfield)[3];				/* Gradient Field g1 */

	LaplacianSystem *sys;
} GradientFlowSystem;

typedef struct LaplacianSystem {
	bool command_compute_flow;
	bool has_solution;
	bool command_remesh;

	InputMesh input_mesh;
	OutputMesh output_mesh;
	
	float *U_field;					/* Initial scalar field*/

	float(*gf1)[3], (*gf2)[3];		/* Gradient Fields g1 and g2 per face*/
	float *h1, *h2;					/* Sampling distance functions h1 and h2 */
	float h;

	NLContext *context;				/* System for solve general implicit rotations */
	GradientFlowSystem *gfsys1, *gfsys2;
	DerivedMesh *resultDM;
} LaplacianSystem;

/*
* alpha is degree of anisotropic curvature sensitivity
* h is the desired distance
* return ve[0] number of vertices
* return ve[1] number of edges
*/
//void estimateNumberGFVerticesEdges(int ve[2], LaplacianSystem *sys, float h);

//void deleteGradientFlowMesh(GradientFlowMesh * gfmesh);
//int addGFVertGFMesh(GradientFlowMesh *gfmesh, GFVert gfvert);
//int addVertGFMesh(GradientFlowMesh *gfmesh, float co[3], int index_edge);
//int addGFEdgeGFMesh(GradientFlowMesh *gfmesh, GFEdge gfedge);
//int addEdgeGFMesh(GradientFlowMesh *gfmesh, int index_v1, int index_v2, int index_face);
//
//GradientFlowSystem *newGradientFlowSystem(LaplacianSystem *sys, float *mhfunction, float(*mgfield)[3]);
void deleteGradientFlowSystem(GradientFlowSystem *gfsys);
//int addGFVertGFSystem(GradientFlowSystem *gfsys, GFVert gfvert);
//int addVertGFSystem(GradientFlowSystem *gfsys, float co[3], int index_edge);
//int addGFEdgeGFSystem(GradientFlowSystem *gfsys, GFEdge gfedge);
//int addEdgeGFSystem(GradientFlowSystem *gfsys, int index_v1, int index_v2, int index_face);
//int addEdgeTwoFacesGFSystem(GradientFlowSystem *gfsys, int index_v1, int index_v2, int index_face1, int index_face2);
//
//int *findFeaturesOnMesh(int size[2], LaplacianSystem *sys);
//void addSeedToQueue(struct Heap *aheap, float value, GFVert *vert);
//GFVert *getTopSeedFromQueue(struct Heap *aheap);
//
//bool isOnSegmentLine(float p1[3], float p2[3], float q[3]);
//bool intersecionLineSegmentWithVector(float r[3], float p1[3], float p2[3], float ori[3], float dir[3]);
//int getEdgeFromVerts(LaplacianSystem *sys, int v1, int v2);
//int getOtherFaceAdjacentToEdge(LaplacianSystem *sys, int oldface, int inde);
//void projectVectorOnFace(float r[3], float no[3], float dir[3]);
//int getDifferentVertexFaceEdge(LaplacianSystem *sys, int oldface, int inde);
//#define GRA_DIR_ON_NONE 0
//#define GRA_DIR_ON_FACE 1
//#define GRA_DIR_ON_EDGE 2
//void computeGradientDirectionOnVert(int rind[2], float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexvert);
//void computeGradientDirectionOnEdgeInverse(int rind[2], float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexedge);
//int nearGFEdgeInGFMesh(LaplacianSystem *sys, GradientFlowSystem *gfsys, float ori[3], float dir[3], int indexface, float maxradius);
//int nearGFEdgeInGFMeshFromEdge(LaplacianSystem *sys, GradientFlowSystem *gfsys, float ori[3], float dir[3], int indexedge, float maxradius);
//int nextPointFlowLine(float r[3], LaplacianSystem *sys, float q[3], int oldface, int inde);
//int nextPointFlowLineInverse(float r[3], LaplacianSystem *sys, float q[3], int oldface, int inde);
//float getSamplingDistanceFunctionOnFace(LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexface);
//float getMaxSamplingDistanceFunctionOnFace(LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexface);
//void computeGFLine(LaplacianSystem *sys, GradientFlowSystem *gfsys, GFVert *gfvert_seed);
//
//int computeNewSeed(float r[3], LaplacianSystem *sys, GradientFlowSystem *gfsys, int indexf, float ori[3], float dir[3], float mh);

void computeFlowLines(LaplacianSystem *sys);
void generateMesh(LaplacianSystem *sys);

#endif /*openNl*/
#endif /*__MOD_QUADREMESH_GEOM_H__*/