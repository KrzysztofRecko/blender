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

#ifndef __QUADREMESH_UTIL_H__
#define __QUADREMESH_UTIL_H__

#define QR_SEEDPROB 0.2f
#define QR_SAMPLING_RATE 0.04f
#define QR_MAXDIST_TO_SEEDDIST 2.0f
#define QR_NO_FACE 0xffffffff

typedef int MEdgeID;
typedef int QRVertID;

typedef enum {
	GFSYS1 = 1 << 0,
	GFSYS2 = 1 << 1,
	GFSYSNONE = 1 << 2
} GFSysID;

typedef struct QRDiskLink {
	struct QRDiskLink *next, *prev, *brother;
	MEdgeID e;
	QRVertID v;
	float ang;
	bool poly_on_right;
} QRDiskLink;

typedef struct QREdgeLink {
	QRVertID v;
	float dist;

	struct QREdgeLink *next;
	struct QRDiskLink *elink;

	GFSysID gfsysid;
} QREdgeLink;

typedef struct QREdge {
	QREdgeLink *v1, *v2;
	float dir[3], orig[3];
	int num_links;
} QREdge;

typedef struct QRVert {
	struct QRDiskLink *link;
	float co[3], no[3], vec[3];
	int num_links;
} QRVert;

typedef enum {
	eVert,
	eEdge,
	eFace
} GFPointType;

typedef struct GFPoint {
	float co[3];
	union {
		int v, e, f;
	};
	QRVertID id;
	GFPointType type;
} GFPoint;

typedef enum {
	eEdgeSegment,
	eFaceSegment,
	eNoSegment
} GFSegmentType;

typedef struct GFSegment {
	GFPoint p;
	union {
		int e, f;
	};
	GFSegmentType type;
} GFSegment;

/* GradientFlowSysten, one gfsys for every gradient field */
typedef struct GradientFlowSystem {
	GFSysID id;
	struct Heap *seeds;
	//float *h;
	float(*gf)[3];				/* Gradient Field */

	struct QuadRemeshSystem *sys;
} GradientFlowSystem;

typedef struct InputMesh {
	bool is_alloc;

	int num_verts, num_edges, num_faces, num_features;
	float(*co)[3];					/* Original vertex coordinates */
	float(*no)[3];					/* Original face normal */

	unsigned int(*faces)[3];		/* Copy of MFace (tessface) v1-v3, v2-v4 */
	unsigned int(*edges)[2];		/* Copy of edges v1-v2 */
	unsigned int(*faces_edge)[2];	/* Faces by edges  */

	int *ringf_indices;				/* Indices of faces per vertex */
	int *ringe_indices;				/* Indices of edges per vertex */
	struct MeshElemMap *ringf_map;  /* Map of faces per vertex */
	struct MeshElemMap *ringe_map;	/* Map of edges per vertex */

	int *constraints;				/* Feature points constraints*/
	float *weights;					/* Feature points weights*/
} InputMesh;

typedef struct OutputMesh {
	struct MemArena *memarena;
	QRVert *verts;
	QRVertID *vonvs;
	struct QREdge *ringe;                  /* QREdges per original edges */
	struct LinkNode **ringf;        /* Lists of QREdge per original faces */
	
	int num_verts, alloc_verts;
	int num_edges;
} OutputMesh;

typedef struct QuadRemeshSystem {
	struct QuadRemeshModifierData *qmd;

	InputMesh input_mesh;
	
	bool has_solution, is_alloc;
	GradientFlowSystem *gfsys[2];
	float *U_field;					/* Initial scalar field*/

	OutputMesh output_mesh;
	struct RNG *rng;

	struct DerivedMesh *cache_mesh;
} QuadRemeshSystem;

void getNormalAtEdge(float r_no[3], InputMesh *im, int in_e);

QRDiskLink *linkVerts(OutputMesh *om, QRVertID in_v1, QRVertID in_v2);
void deleteLink(OutputMesh *om, QRVert *ll, QRDiskLink *l);
void unlinkVerts(OutputMesh *om, QRDiskLink *l);
void dissolveVert(OutputMesh *om, QRVertID in_v);

void insertOnQREdge(OutputMesh *om, QREdge *in_e, QRVertID in_vid);
void linkOnQREdge(OutputMesh *om, GFSysID sys_id, QREdge *in_e, QRVertID in_v1, QRVertID in_v2);
QREdge *addQREdgeToFace(OutputMesh *om, InputMesh *im, GFSysID sys_id, int in_f, QRVertID in_v1, QRVertID in_v2);

QRVertID addVert(OutputMesh *om, float in_co[3], float in_no[3]);
void addGFPoint(InputMesh *im, OutputMesh *om, GFPoint *in_p);

void addSeedToQueue(struct Heap *aheap, float in_co[3], GFPointType in_type, int in_val, float weight);
GFPoint *getTopSeedFromQueue(struct Heap *aheap);

bool isectSegmentWithOthersOnFace(OutputMesh *om, GFSysID sys_id, float in_a[3], float in_b[3], int in_f);
bool isectPointWithQREdge(OutputMesh *om, GFSysID sys_id, float in_co[3], int in_e);

void getInitialSeeds(GradientFlowSystem *gfsys);
void computeFlowLines(QuadRemeshSystem *sys);

void freeInputMesh(InputMesh *im);
void getInput(QuadRemeshSystem *sys, struct Object *ob, struct DerivedMesh *dm);
void getHarmonicGradients(QuadRemeshSystem *sys);

GradientFlowSystem *newGradientFlowSystem(QuadRemeshSystem *sys);
void freeGradientFlowSystem(GradientFlowSystem *gfsys);

#endif