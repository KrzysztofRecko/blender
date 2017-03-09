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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * Contributor(s): Jiri Hnidek <jiri.hnidek@vslib.cz>.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/mball_tessellate.c
 *  \ingroup bke
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>

#include "MEM_guardedalloc.h"

#include "DNA_object_types.h"
#include "DNA_meta_types.h"
#include "DNA_scene_types.h"

#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_string_utils.h"
#include "BLI_utildefines.h"
#include "BLI_memarena.h"
#include "BLI_task.h"

#include "BKE_global.h"

#include "BKE_depsgraph.h"
#include "BKE_scene.h"
#include "BKE_displist.h"
#include "BKE_mball_tessellate.h"  /* own include */

#include "BLI_strict_flags.h"

/* experimental (faster) normal calculation */
// #define USE_ACCUM_NORMAL
// #define MB_LINEAR_CONVERGE

/* Data types */

typedef struct corner {         /* corner of a cube */
	int lat[3];                 /* index within lattice */
	float co[3], value;         /* location and function value */
	struct corner *next;
} CORNER;

typedef struct cube {           /* partitioning cell (cube) */
	int lat[3];                 /* lattice location of cube */
	CORNER *corners[8];         /* eight corners */
} CUBE;

typedef struct cubes {          /* linked list of cubes acting as stack */
	CUBE cube;                  /* a single cube */
	struct cubes *next;         /* remaining elements */
} CUBES;

typedef struct centerlist {     /* list of cube locations */
	int lat[3];                 /* cube location */
	struct centerlist *next;    /* remaining elements */
} CENTERLIST;

typedef struct edgelist {       /* list of edges */
	int a[3], b[3];             /* edge corner ids */
	int vid;                    /* vertex id */
	struct edgelist *next;      /* remaining elements */
} EDGELIST;

typedef struct intlist {        /* list of integers */
	int i;                      /* an integer */
	struct intlist *next;       /* remaining elements */
} INTLIST;

typedef struct intlists {       /* list of list of integers */
	INTLIST *list;              /* a list of integers */
	struct intlists *next;      /* remaining elements */
} INTLISTS;

typedef struct Box {			/* an AABB with pointer to metalelem */
	float min[3], max[3];
	struct MLSmall *ml;
} Box;

typedef struct MLSmall {        /* a compact structure like this improves cache usage */
	float imat[4][4];           /* this is only stuff needed for computing density and BVH build*/
	short type;
	float expx, expy, expz;
	float rad2, s;
	Box *bb;
} MLSmall;

typedef struct MetaballBVHNode {	/* BVH node */
	Box bb[2];						/* AABB of children */
	struct MetaballBVHNode *child[2];
} MetaballBVHNode;

typedef struct chunk CHUNK;

typedef struct process {        /* parameters, storage */
	float thresh, size;			/* mball threshold, single cube size */
	float delta;				/* small delta for calculating normals */
	unsigned int converge_res;	/* converge procedure resolution (more = slower) */
	unsigned int chunk_res, num_chunks;

	CHUNK **chunks;
	bool cubes_left_to_process;

	MLSmall **mainb;			/* array of all metaelems */
	unsigned int totelem, mem;	/* number of metaelems */

	Box allbb;                  /* Bounding box of all metaelems */

	EDGELIST **common_edges;    /* edge and vertex id hash table, common for all chunks */
	ThreadMutex edge_lock;
	float(*co)[3], (*no)[3];    /* surface vertices - positions and normals */
	unsigned int totvertex;		/* memory size */
	unsigned int curvertex;		/* currently added vertices */

	/* memory allocation from common pool */
	MemArena *pgn_elements;
	MemArena *metaballs;        /* metaballs allocated from their own memarena to improve data locality */
} PROCESS;

typedef struct chunk {
	PROCESS *process;

	Box bb;
	int index, lat_pos[3], min_lat[3], max_lat[3];

	MLSmall **mainb;            /* each chunk gets its own list of mballs which intersect chunk's bb */
	unsigned int elem;

	MetaballBVHNode **bvh_queue;
	MetaballBVHNode bvh;        /* also, each chunk makes its own BVH */
	unsigned int bvh_queue_size;

	CUBES *cubes;               /* cubes, centers and corners are local for a chunk */
	CENTERLIST **centers;
	CORNER **corners;
	EDGELIST **edges;           /* edge and vertex id hash table, common for all chunks */

	int(*indices)[4];           /* output indices */
	unsigned int totindex;		/* size of memory allocated for indices */
	unsigned int curindex;		/* number of currently added indices */

	float(*co)[3], (*no)[3];    /* surface vertices - positions and normals */
	unsigned int totvertex;		/* memory size */
	unsigned int curvertex;		/* currently added vertices */
	unsigned int vertex_offset;

	CHUNK *neighbors[6];
	CENTERLIST *cube_queue;
	MemArena *cube_queue_mem;
	ThreadMutex cube_queue_lock;

	MemArena *mem;
} CHUNK;

/* Forward declarations */
static int vertid(CHUNK *chunk, const CORNER *c1, const CORNER *c2);
static void add_cube(CHUNK *chunk, int i, int j, int k);
static void make_face(CHUNK *chunk, int i1, int i2, int i3, int i4);
static void converge(CHUNK *chunk, const CORNER *c1, const CORNER *c2, float r_p[3]);

/* ******************* SIMPLE BVH ********************* */

static void make_union(const Box *a, const Box *b, Box *r_out)
{
	r_out->min[0] = min_ff(a->min[0], b->min[0]);
	r_out->min[1] = min_ff(a->min[1], b->min[1]);
	r_out->min[2] = min_ff(a->min[2], b->min[2]);

	r_out->max[0] = max_ff(a->max[0], b->max[0]);
	r_out->max[1] = max_ff(a->max[1], b->max[1]);
	r_out->max[2] = max_ff(a->max[2], b->max[2]);
}

static void make_box_from_ml(Box *r, MLSmall *ml)
{
	copy_v3_v3(r->max, ml->bb->max);
	copy_v3_v3(r->min, ml->bb->min);
	r->ml = ml;
}

/**
 * Partitions part of mainb array [start, end) along axis s. Returns i,
 * where centroids of elements in the [start, i) segment lie "on the right side" of div,
 * and elements in the [i, end) segment lie "on the left"
 */
static unsigned int partition_mainb(MLSmall **mainb, unsigned int start, unsigned int end, unsigned int s, float div)
{
	unsigned int i = start, j = end - 1;
	div *= 2.0f;

	while (1) {
		while (i < j && div > (mainb[i]->bb->max[s] + mainb[i]->bb->min[s])) i++;
		while (j > i && div < (mainb[j]->bb->max[s] + mainb[j]->bb->min[s])) j--;

		if (i >= j)	break;

		SWAP(MLSmall *, mainb[i], mainb[j]);
		i++;
		j--;
	}

	if (i == start) i++;

	return i;
}

/**
 * Recursively builds a BVH, dividing elements along the middle of the longest axis of allbox.
 */
static void build_bvh_spatial(
        CHUNK *chunk, MetaballBVHNode *node,
        unsigned int start, unsigned int end, const Box *allbox)
{
	unsigned int part, j, s;
	float dim[3], div;

	/* Maximum bvh queue size is number of nodes which are made, equals calls to this function. */
	chunk->bvh_queue_size++;

	VECSUB(dim, allbox->max, allbox->min);

	s = 0;
	if (dim[1] > dim[0] && dim[1] > dim[2]) s = 1;
	else if (dim[2] > dim[1] && dim[2] > dim[0]) s = 2;

	div = allbox->min[s] + (dim[s] / 2.0f);

	part = partition_mainb(chunk->mainb, start, end, s, div);

	make_box_from_ml(&node->bb[0], chunk->mainb[start]);
	node->child[0] = NULL;

	if (part > start + 1) {
		for (j = start; j < part; j++) {
			make_union(chunk->mainb[j]->bb, &node->bb[0], &node->bb[0]);
		}

		node->child[0] = BLI_memarena_alloc(chunk->mem, sizeof(MetaballBVHNode));
		build_bvh_spatial(chunk, node->child[0], start, part, &node->bb[0]);
	}

	node->child[1] = NULL;
	if (part < end) {
		make_box_from_ml(&node->bb[1], chunk->mainb[part]);

		if (part < end - 1) {
			for (j = part; j < end; j++) {
				make_union(chunk->mainb[j]->bb, &node->bb[1], &node->bb[1]);
			}

			node->child[1] = BLI_memarena_alloc(chunk->mem, sizeof(MetaballBVHNode));
			build_bvh_spatial(chunk, node->child[1], part, end, &node->bb[1]);
		}
	}
	else {
		INIT_MINMAX(node->bb[1].min, node->bb[1].max);
	}
}

/* ******************** ARITH ************************* */

/**
 * BASED AT CODE (but mostly rewritten) :
 * C code from the article
 * "An Implicit Surface Polygonizer"
 * by Jules Bloomenthal, jbloom@beauty.gmu.edu
 * in "Graphics Gems IV", Academic Press, 1994
 *
 * Authored by Jules Bloomenthal, Xerox PARC.
 * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
 * Permission is granted to reproduce, use and distribute this code for
 * any and all purposes, provided that this notice appears in all copies.
 */

#define L   0  /* left direction:	-x, -i */
#define R   1  /* right direction:	+x, +i */
#define B   2  /* bottom direction: -y, -j */
#define T   3  /* top direction:	+y, +j */
#define N   4  /* near direction:	-z, -k */
#define F   5  /* far direction:	+z, +k */
#define LBN 0  /* left bottom near corner  */
#define LBF 1  /* left bottom far corner   */
#define LTN 2  /* left top near corner     */
#define LTF 3  /* left top far corner      */
#define RBN 4  /* right bottom near corner */
#define RBF 5  /* right bottom far corner  */
#define RTN 6  /* right top near corner    */
#define RTF 7  /* right top far corner     */

/**
 * the LBN corner of cube (i, j, k), corresponds with location
 * (i-0.5)*size, (j-0.5)*size, (k-0.5)*size)
 */

#define HASHBIT     (5)
#define HASHSIZE    (size_t)(1 << (3 * HASHBIT))   /*! < hash table size (32768) */

#define HASH(i, j, k) ((((( (i) & 31) << 5) | ( (j) & 31)) << 5) | ( (k) & 31) )

#define MB_BIT(i, bit) (((i) >> (bit)) & 1)
// #define FLIP(i, bit) ((i) ^ 1 << (bit)) /* flip the given bit of i */

/* ******************** DENSITY COPMPUTATION ********************* */

/**
 * Computes density from given metaball at given position.
 * Metaball equation is: ``(1 - r^2 / R^2)^3 * s``
 *
 * r = distance from center
 * R = metaball radius
 * s - metaball stiffness
 */
static float densfunc(MLSmall *ball, float x, float y, float z)
{
	float dist2;
	float dvec[3] = {x, y, z};

	mul_m4_v3(ball->imat, dvec);

	switch (ball->type) {
		case MB_BALL:
			/* do nothing */
			break;
		case MB_CUBE:
			if      (dvec[2] > ball->expz)  dvec[2] -= ball->expz;
			else if (dvec[2] < -ball->expz) dvec[2] += ball->expz;
			else                            dvec[2] = 0.0;
			/* fall through */
		case MB_PLANE:
			if      (dvec[1] >  ball->expy) dvec[1] -= ball->expy;
			else if (dvec[1] < -ball->expy) dvec[1] += ball->expy;
			else                            dvec[1] = 0.0;
			/* fall through */
		case MB_TUBE:
			if      (dvec[0] >  ball->expx) dvec[0] -= ball->expx;
			else if (dvec[0] < -ball->expx) dvec[0] += ball->expx;
			else                            dvec[0] = 0.0;
			break;
		case MB_ELIPSOID:
			dvec[0] /= ball->expx;
			dvec[1] /= ball->expy;
			dvec[2] /= ball->expz;
			break;
#if 0
		/* *** deprecated, could be removed?, do-versioned at least *** */
		case MB_TUBEX:
			if      (dvec[0] >  ball->len) dvec[0] -= ball->len;
			else if (dvec[0] < -ball->len) dvec[0] += ball->len;
			else                           dvec[0] = 0.0;
			break;
		case MB_TUBEY:
			if      (dvec[1] >  ball->len) dvec[1] -= ball->len;
			else if (dvec[1] < -ball->len) dvec[1] += ball->len;
			else                           dvec[1] = 0.0;
			break;
		case MB_TUBEZ:
			if      (dvec[2] >  ball->len) dvec[2] -= ball->len;
			else if (dvec[2] < -ball->len) dvec[2] += ball->len;
			else                           dvec[2] = 0.0;
			break;
			/* *** end deprecated *** */
#endif
	}

	/* ball->rad2 is inverse of squared rad */
	dist2 = 1.0f - (len_squared_v3(dvec) * ball->rad2);

	/* ball->s is negative if metaball is negative */
	return (dist2 < 0.0f) ? 0.0f : (ball->s * dist2 * dist2 * dist2);
}

/**
 * Computes density at given position form all metaballs which contain this point in their box.
 * Traverses BVH using a queue.
 */
static float metaball(CHUNK *chunk, float x, float y, float z)
{
	int i;
	float dens = 0.0f;
	unsigned int front = 0, back = 0;
	MetaballBVHNode *node;

	chunk->bvh_queue[front++] = &chunk->bvh;

	while (front != back) {
		node = chunk->bvh_queue[back++];

		for (i = 0; i < 2; i++) {
			if ((node->bb[i].min[0] <= x) && (node->bb[i].max[0] >= x) &&
			    (node->bb[i].min[1] <= y) && (node->bb[i].max[1] >= y) &&
			    (node->bb[i].min[2] <= z) && (node->bb[i].max[2] >= z))
			{
				if (node->child[i])	chunk->bvh_queue[front++] = node->child[i];
				else dens += densfunc(node->bb[i].ml, x, y, z);
			}
		}
	}

	return chunk->process->thresh - dens;
}

/**
 * Adds face to indices, expands memory if needed.
 */
static void make_face(CHUNK *chunk, int i1, int i2, int i3, int i4)
{
	int *cur;

	if (UNLIKELY(chunk->totindex == chunk->curindex)) {
		chunk->totindex += 4096;
		chunk->indices = MEM_reallocN(chunk->indices, sizeof(int[4]) * chunk->totindex);
	}

	cur = chunk->indices[chunk->curindex++];

	/* displists now support array drawing, we treat tri's as fake quad */

	cur[0] = i1;
	cur[1] = i2;
	cur[2] = i3;
	if (i4 == -1) cur[3] = i3;
	else          cur[3] = i4;
}

static void freechunk(CHUNK *chunk)
{
	if (chunk->corners)        MEM_freeN(chunk->corners);
	if (chunk->centers)        MEM_freeN(chunk->centers);
	if (chunk->edges)          MEM_freeN(chunk->edges);
	if (chunk->bvh_queue)      MEM_freeN(chunk->bvh_queue);
	if (chunk->mainb)          MEM_freeN(chunk->mainb);
	if (chunk->co)             MEM_freeN(chunk->co);
	if (chunk->no)             MEM_freeN(chunk->no);
	if (chunk->indices)        MEM_freeN(chunk->indices);
	if (chunk->mem)            BLI_memarena_free(chunk->mem);
	if (chunk->cube_queue_mem) BLI_memarena_free(chunk->cube_queue_mem);

	BLI_mutex_end(&chunk->cube_queue_lock);

	MEM_freeN(chunk);
}

/* Frees allocated memory */
static void freeprocess(PROCESS *process)
{
	unsigned int i;

	if (process->chunks) {
		for (i = 0; i < process->num_chunks; i++)
			if (process->chunks[i]) freechunk(process->chunks[i]);
		MEM_freeN(process->chunks);
	}

	if (process->num_chunks > 1) { /* These are shared if there is only 1 chunk. */
		if (process->mainb)        MEM_freeN(process->mainb);
		if (process->pgn_elements) BLI_memarena_free(process->pgn_elements);
	}
	if (process->common_edges) MEM_freeN(process->common_edges);
	if (process->co)           MEM_freeN(process->co);
	if (process->no)           MEM_freeN(process->no);
	if (process->metaballs)    BLI_memarena_free(process->metaballs);

	BLI_mutex_end(&process->edge_lock);
}

/* **************** POLYGONIZATION ************************ */

/**** Cubical Polygonization (optional) ****/

#define LB  0  /* left bottom edge	*/
#define LT  1  /* left top edge	*/
#define LN  2  /* left near edge	*/
#define LF  3  /* left far edge	*/
#define RB  4  /* right bottom edge */
#define RT  5  /* right top edge	*/
#define RN  6  /* right near edge	*/
#define RF  7  /* right far edge	*/
#define BN  8  /* bottom near edge	*/
#define BF  9  /* bottom far edge	*/
#define TN  10 /* top near edge	*/
#define TF  11 /* top far edge	*/

static INTLISTS *cubetable[256];
static char faces[256];

/* edge: LB, LT, LN, LF, RB, RT, RN, RF, BN, BF, TN, TF */
static int corner1[12] = {
	LBN, LTN, LBN, LBF, RBN, RTN, RBN, RBF, LBN, LBF, LTN, LTF
};
static int corner2[12] = {
	LBF, LTF, LTN, LTF, RBF, RTF, RTN, RTF, RBN, RBF, RTN, RTF
};
static int leftface[12] = {
	B, L, L, F, R, T, N, R, N, B, T, F
};
/* face on left when going corner1 to corner2 */
static int rightface[12] = {
	L, T, N, L, B, R, R, F, B, F, N, T
};
/* face on right when going corner1 to corner2 */

/**
 * triangulate the cube directly, without decomposition
 */
static void docube(CHUNK *chunk, CUBE *cube)
{
	INTLISTS *polys;
	CORNER *c1, *c2;
	int i, index = 0, count, indexar[8];

	/* Determine which case cube falls into. */
	for (i = 0; i < 8; i++) {
		if (cube->corners[i]->value > 0.0f)	index += (1 << i);
	}

	/* Using faces[] table, adds neighbouring cube if surface intersects face in this direction. */
	if (MB_BIT(faces[index], 0)) add_cube(chunk, cube->lat[0] - 1, cube->lat[1], cube->lat[2]);
	if (MB_BIT(faces[index], 1)) add_cube(chunk, cube->lat[0] + 1, cube->lat[1], cube->lat[2]);
	if (MB_BIT(faces[index], 2)) add_cube(chunk, cube->lat[0], cube->lat[1] - 1, cube->lat[2]);
	if (MB_BIT(faces[index], 3)) add_cube(chunk, cube->lat[0], cube->lat[1] + 1, cube->lat[2]);
	if (MB_BIT(faces[index], 4)) add_cube(chunk, cube->lat[0], cube->lat[1], cube->lat[2] - 1);
	if (MB_BIT(faces[index], 5)) add_cube(chunk, cube->lat[0], cube->lat[1], cube->lat[2] + 1);

	/* Using cubetable[], determines polygons for output. */
	for (polys = cubetable[index]; polys; polys = polys->next) {
		INTLIST *edges;

		count = 0;
		/* Sets needed vertex id's lying on the edges. */
		for (edges = polys->list; edges; edges = edges->next) {
			c1 = cube->corners[corner1[edges->i]];
			c2 = cube->corners[corner2[edges->i]];

			indexar[count] = vertid(chunk, c1, c2);
			count++;
		}

		/* Adds faces to output. */
		switch (count) {
			case 3:
				make_face(chunk, indexar[2], indexar[1], indexar[0], -1);
				break;
			case 4:
				make_face(chunk, indexar[3], indexar[2], indexar[1], indexar[0]);
				break;
			case 5:
				make_face(chunk, indexar[3], indexar[2], indexar[1], indexar[0]);
				make_face(chunk, indexar[4], indexar[3], indexar[0], -1);
				break;
			case 6:
				make_face(chunk, indexar[3], indexar[2], indexar[1], indexar[0]);
				make_face(chunk, indexar[5], indexar[4], indexar[3], indexar[0]);
				break;
			case 7:
				make_face(chunk, indexar[3], indexar[2], indexar[1], indexar[0]);
				make_face(chunk, indexar[5], indexar[4], indexar[3], indexar[0]);
				make_face(chunk, indexar[6], indexar[5], indexar[0], -1);
				break;
		}
	}
}
/**
 * return corner with the given lattice location
 * set (and cache) its function value
 */
static CORNER *setcorner(CHUNK *chunk, const int i, const int j, const int k)
{
	/* for speed, do corner value caching here */
	CORNER *c;
	int index;

	/* does corner exist? */
	index = HASH(i, j, k);
	c = chunk->corners[index];

	for (; c != NULL; c = c->next) {
		if (c->lat[0] == i && c->lat[1] == j && c->lat[2] == k) return c;
	}

	c = BLI_memarena_alloc(chunk->mem, sizeof(CORNER));

	c->lat[0] = i;
	c->lat[1] = j;
	c->lat[2] = k;

	c->co[0] = ((float)i - 0.5f) * chunk->process->size;
	c->co[1] = ((float)j - 0.5f) * chunk->process->size;
	c->co[2] = ((float)k - 0.5f) * chunk->process->size;

#ifdef MB_PRECOMPUTE
	c->value = chunk->process->thresh;
#else
	c->value = metaball(chunk, c->co[0], c->co[1], c->co[2]);
#endif

	c->next = chunk->corners[index];
	chunk->corners[index] = c;
 
	return c;
}
/**
 * return next clockwise edge from given edge around given face
 */
static int nextcwedge(int edge, int face)
{
	switch (edge) {
		case LB:
			return (face == L) ? LF : BN;
		case LT:
			return (face == L) ? LN : TF;
		case LN:
			return (face == L) ? LB : TN;
		case LF:
			return (face == L) ? LT : BF;
		case RB:
			return (face == R) ? RN : BF;
		case RT:
			return (face == R) ? RF : TN;
		case RN:
			return (face == R) ? RT : BN;
		case RF:
			return (face == R) ? RB : TF;
		case BN:
			return (face == B) ? RB : LN;
		case BF:
			return (face == B) ? LB : RF;
		case TN:
			return (face == T) ? LT : RN;
		case TF:
			return (face == T) ? RT : LF;
	}
	return 0;
}

/**
 * \return the face adjoining edge that is not the given face
 */
static int otherface(int edge, int face)
{
	int other = leftface[edge];
	return face == other ? rightface[edge] : other;
}

/**
 * create the 256 entry table for cubical polygonization
 */
static void makecubetable(void)
{
	static bool is_done = false;
	int i, e, c, done[12], pos[8];

	if (is_done) return;
	is_done = true;

	for (i = 0; i < 256; i++) {
		for (e = 0; e < 12; e++) done[e] = 0;
		for (c = 0; c < 8; c++) pos[c] = MB_BIT(i, c);
		for (e = 0; e < 12; e++)
			if (!done[e] && (pos[corner1[e]] != pos[corner2[e]])) {
				INTLIST *ints = NULL;
				INTLISTS *lists = MEM_callocN(sizeof(INTLISTS), "mball_intlist");
				int start = e, edge = e;

				/* get face that is to right of edge from pos to neg corner: */
				int face = pos[corner1[e]] ? rightface[e] : leftface[e];

				while (1) {
					edge = nextcwedge(edge, face);
					done[edge] = 1;
					if (pos[corner1[edge]] != pos[corner2[edge]]) {
						INTLIST *tmp = ints;

						ints = MEM_callocN(sizeof(INTLIST), "mball_intlist");
						ints->i = edge;
						ints->next = tmp; /* add edge to head of list */

						if (edge == start) break;
						face = otherface(edge, face);
					}
				}
				lists->list = ints; /* add ints to head of table entry */
				lists->next = cubetable[i];
				cubetable[i] = lists;
			}
	}

	for (i = 0; i < 256; i++) {
		INTLISTS *polys;
		faces[i] = 0;
		for (polys = cubetable[i]; polys; polys = polys->next) {
			INTLIST *edges;

			for (edges = polys->list; edges; edges = edges->next) {
				if (edges->i == LB || edges->i == LT || edges->i == LN || edges->i == LF) faces[i] |= 1 << L;
				if (edges->i == RB || edges->i == RT || edges->i == RN || edges->i == RF) faces[i] |= 1 << R;
				if (edges->i == LB || edges->i == RB || edges->i == BN || edges->i == BF) faces[i] |= 1 << B;
				if (edges->i == LT || edges->i == RT || edges->i == TN || edges->i == TF) faces[i] |= 1 << T;
				if (edges->i == LN || edges->i == RN || edges->i == BN || edges->i == TN) faces[i] |= 1 << N;
				if (edges->i == LF || edges->i == RF || edges->i == BF || edges->i == TF) faces[i] |= 1 << F;
			}
		}
	}
}

void BKE_mball_cubeTable_free(void)
{
	int i;
	INTLISTS *lists, *nlists;
	INTLIST *ints, *nints;

	for (i = 0; i < 256; i++) {
		lists = cubetable[i];
		while (lists) {
			nlists = lists->next;

			ints = lists->list;
			while (ints) {
				nints = ints->next;
				MEM_freeN(ints);
				ints = nints;
			}

			MEM_freeN(lists);
			lists = nlists;
		}
		cubetable[i] = NULL;
	}
}

/**** Storage ****/

/**
 * Inserts cube at lattice i, j, k into hash table, marking it as "done"
 */
static int setcenter(CHUNK *chunk, CENTERLIST *table[], const int i, const int j, const int k)
{
	int index;
	CENTERLIST *newc, *l, *q;

	index = HASH(i, j, k);
	q = table[index];

	for (l = q; l != NULL; l = l->next) {
		if (l->lat[0] == i && l->lat[1] == j && l->lat[2] == k) return 1;
	}

	newc = BLI_memarena_alloc(chunk->mem, sizeof(CENTERLIST));
	newc->lat[0] = i;
	newc->lat[1] = j;
	newc->lat[2] = k;
	newc->next = q;
	table[index] = newc;

	return 0;
}

/**
 * Adds a vertex, expands memory if needed.
 */
static unsigned int addtovertices(CHUNK *chunk, bool is_common)
{
	if (is_common) {
		if (UNLIKELY(chunk->process->curvertex == chunk->process->totvertex)) {
			chunk->process->totvertex += 4096;
			chunk->process->co = MEM_reallocN(chunk->process->co, chunk->process->totvertex * sizeof(float[3]));
			chunk->process->no = MEM_reallocN(chunk->process->no, chunk->process->totvertex * sizeof(float[3]));

		}

		return chunk->process->curvertex++;
	}
	else {
		if (UNLIKELY(chunk->curvertex == chunk->totvertex)) {
			chunk->totvertex += 4096;
			chunk->co = MEM_reallocN(chunk->co, chunk->totvertex * sizeof(float[3]));
			chunk->no = MEM_reallocN(chunk->no, chunk->totvertex * sizeof(float[3]));
		}

		return chunk->curvertex++;
	}
}

static void copytovertices(CHUNK *chunk, const float co[3], const float no[3], const unsigned int vid, bool is_common)
{
	if (is_common) {
		copy_v3_v3(chunk->process->co[vid], co);
		copy_v3_v3(chunk->process->no[vid], no);
	}
	else {
		copy_v3_v3(chunk->co[vid], co);
		copy_v3_v3(chunk->no[vid], no);
	}
}

/**
 * Computes normal from density field at given point.
 *
 * \note Doesn't do normalization!
 */
static void vnormal(CHUNK *chunk, const float point[3], float r_no[3])
{
#ifdef MB_ACCUM_NORMAL
	r_no[0] = r_no[1] = r_no[2] = 0.0f;
#else
	const float delta = chunk->process->delta;
	const float f = metaball(chunk, point[0], point[1], point[2]);

	r_no[0] = metaball(chunk, point[0] + delta, point[1], point[2]) - f;
	r_no[1] = metaball(chunk, point[0], point[1] + delta, point[2]) - f;
	r_no[2] = metaball(chunk, point[0], point[1], point[2] + delta) - f;
#endif

#if 0
	f = normalize_v3(r_no);

	if (0) {
		float tvec[3];

		delta *= 2.0f;

		f = process->function(process, point[0], point[1], point[2]);

		tvec[0] = process->function(process, point[0] + delta, point[1], point[2]) - f;
		tvec[1] = process->function(process, point[0], point[1] + delta, point[2]) - f;
		tvec[2] = process->function(process, point[0], point[1], point[2] + delta) - f;

		if (normalize_v3(tvec) != 0.0f) {
			add_v3_v3(r_no, tvec);
			normalize_v3(r_no);
		}
	}
#endif
}

static void orient_edge(int a[3], int b[3])
{
	if (a[0] > b[0] || a[1] > b[1] || a[2] > b[2]) {
		SWAP(int, a[0], b[0]);
		SWAP(int, a[1], b[1]);
		SWAP(int, a[2], b[2]);
	}
}
static bool equal_v3_v3_int(const int a[3], const int b[3])
{
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
static bool is_boundary_edge(CHUNK *chunk, const int one[3], const int two[3])
{
	if (one[0] == chunk->min_lat[0]     && two[0] == chunk->min_lat[0])     return true;
	if (one[0] == chunk->max_lat[0] + 1 && two[0] == chunk->max_lat[0] + 1) return true;
	if (one[1] == chunk->min_lat[1]     && two[1] == chunk->min_lat[1])     return true;
	if (one[1] == chunk->max_lat[1] + 1 && two[1] == chunk->max_lat[1] + 1) return true;
	if (one[2] == chunk->min_lat[2]     && two[2] == chunk->min_lat[2])     return true;
	if (one[2] == chunk->max_lat[2] + 1 && two[2] == chunk->max_lat[2] + 1) return true;
	return false;
}
static bool setedge(CHUNK *chunk, const CORNER *c1, const CORNER *c2, unsigned int *r_vid, bool *r_is_common)
{
	int one[3], two[3], index;
	EDGELIST *q, **table;
	copy_v3_v3_int(one, c1->lat);
	copy_v3_v3_int(two, c2->lat);
	orient_edge(one, two);

	index = HASH(one[0], one[1], one[2]) + HASH(two[0], two[1], two[2]);

	*r_is_common = is_boundary_edge(chunk, one, two);

	if (*r_is_common) {
		table = chunk->process->common_edges;
		BLI_mutex_lock(&chunk->process->edge_lock);
	}
	else table = chunk->edges;

	q = table[index];

	for (; q != NULL; q = q->next) {
		if (equal_v3_v3_int(one, q->a) && equal_v3_v3_int(two, q->b)) {
			*r_vid = (unsigned int)q->vid;

			if (*r_is_common) BLI_mutex_unlock(&chunk->process->edge_lock);

			return true;
		}
	}

	*r_vid = addtovertices(chunk, *r_is_common);

	if (*r_is_common) q = BLI_memarena_alloc(chunk->process->pgn_elements, sizeof(EDGELIST));
	else              q = BLI_memarena_alloc(chunk->mem, sizeof(EDGELIST));

	copy_v3_v3_int(q->a, one);
	copy_v3_v3_int(q->b, two);
	q->vid = (int)*r_vid;
	q->next = table[index];
	table[index] = q;

	if (*r_is_common) BLI_mutex_unlock(&chunk->process->edge_lock);

	return false;
}

/**
 * \return the id of vertex between two corners.
 *
 * If it wasn't previously computed, does #converge() and adds vertex to process.
 */
static int vertid(CHUNK *chunk, const CORNER *c1, const CORNER *c2)
{
	float v[3], no[3];
	unsigned int vid;
	bool is_common;
	
	if (setedge(chunk, c1, c2, &vid, &is_common)) {
		if (is_common) return -(int)vid - 2; /* for common vertices */
		else return (int)vid;
	}

	converge(chunk, c1, c2, v);  /* position */
	vnormal(chunk, v, no); /* normal */

	copytovertices(chunk, v, no, vid, is_common);

	if (is_common) return -(int)vid - 2;
	else return (int)vid;
}

/**
 * Given two corners, computes approximation of surface intersection point between them.
 * In case of small threshold, do bisection.
 */
static void converge(CHUNK *chunk, const CORNER *c1, const CORNER *c2, float r_p[3])
{
	float tmp, dens;
	unsigned int i;
	float c1_value, c1_co[3];
	float c2_value, c2_co[3];

	if (c1->value < c2->value) {
		c1_value = c2->value;
		copy_v3_v3(c1_co, c2->co);
		c2_value = c1->value;
		copy_v3_v3(c2_co, c1->co);
	}
	else {
		c1_value = c1->value;
		copy_v3_v3(c1_co, c1->co);
		c2_value = c2->value;
		copy_v3_v3(c2_co, c2->co);
	}

#ifndef MB_LINEAR_CONVERGE
	for (i = 0; i < chunk->process->converge_res; i++) {
		interp_v3_v3v3(r_p, c1_co, c2_co, 0.5f);
		dens = metaball(chunk, r_p[0], r_p[1], r_p[2]);

		if (dens > 0.0f) {
			c1_value = dens;
			copy_v3_v3(c1_co, r_p);
		}
		else {
			c2_value = dens;
			copy_v3_v3(c2_co, r_p);
		}
	}
#endif

	tmp = -c1_value / (c2_value - c1_value);
	interp_v3_v3v3(r_p, c1_co, c2_co, tmp);
}

static void next_lattice(int r[3], const float pos[3], const float size)
{
	r[0] = (int)ceil((pos[0] / size) + 0.5f);
	r[1] = (int)ceil((pos[1] / size) + 0.5f);
	r[2] = (int)ceil((pos[2] / size) + 0.5f);
}
static void prev_lattice(int r[3], const float pos[3], const float size)
{
	next_lattice(r, pos, size);
	r[0]--; r[1]--; r[2]--;
}
static void closest_latice(int r[3], const float pos[3], const float size)
{
	r[0] = (int)floorf(pos[0] / size + 1.0f);
	r[1] = (int)floorf(pos[1] / size + 1.0f);
	r[2] = (int)floorf(pos[2] / size + 1.0f);
}

static void draw_cube(CHUNK *chunk, int i, int j, int k)
{
	PROCESS *process = chunk->process;

	float n[3] = { 0.0f, 0.0f, 0.0f };
	float v[8][3];
	int vids[8], l;
	v[0][0] = ((float)i - 0.5f) * process->size;
	v[0][1] = ((float)j - 0.5f) * process->size;
	v[0][2] = ((float)k - 0.5f) * process->size;

	v[1][0] = ((float)(i + 1) - 0.5f) * process->size;
	v[1][1] = ((float)j - 0.5f) * process->size;
	v[1][2] = ((float)k - 0.5f) * process->size;

	v[2][0] = ((float)i - 0.5f) * process->size;
	v[2][1] = ((float)(j + 1) - 0.5f) * process->size;
	v[2][2] = ((float)k - 0.5f) * process->size;
	v[3][0] = ((float)i - 0.5f) * process->size;
	v[3][1] = ((float)j - 0.5f) * process->size;
	v[3][2] = ((float)(k + 1) - 0.5f) * process->size;

	v[4][0] = ((float)(i + 1) - 0.5f) * process->size;
	v[4][1] = ((float)(j + 1) - 0.5f) * process->size;
	v[4][2] = ((float)k - 0.5f) * process->size;

	v[5][0] = ((float)i - 0.5f) * process->size;
	v[5][1] = ((float)(j + 1) - 0.5f) * process->size;
	v[5][2] = ((float)(k + 1) - 0.5f) * process->size;

	v[6][0] = ((float)(i + 1) - 0.5f) * process->size;
	v[6][1] = ((float)j - 0.5f) * process->size;
	v[6][2] = ((float)(k + 1) - 0.5f) * process->size;

	v[7][0] = ((float)(i + 1) - 0.5f) * process->size;
	v[7][1] = ((float)(j + 1) - 0.5f) * process->size;
	v[7][2] = ((float)(k + 1) - 0.5f) * process->size;

	for (l = 0; l < 8; l++) {
		vids[l] = (int)addtovertices(chunk, false);
		copytovertices(chunk, v[l], n, (unsigned int)vids[l], false);
	}

	make_face(chunk, vids[0], vids[2], vids[4], vids[1]);
	make_face(chunk, vids[0], vids[3], vids[5], vids[2]);
	make_face(chunk, vids[0], vids[1], vids[6], vids[3]);
	make_face(chunk, vids[1], vids[4], vids[7], vids[6]);
	make_face(chunk, vids[2], vids[5], vids[7], vids[4]);
	make_face(chunk, vids[3], vids[6], vids[7], vids[5]);
}
#if 0
static void draw_box(PROCESS *process, Box *box)
{
	if (!box) return;

	float n[3] = { 0.0f, 0.0f, 0.0f };
	float v[8][3];
	int vids[8];
	v[0][0] = box->min[0];
	v[0][1] = box->min[1];
	v[0][2] = box->min[2];
	v[1][0] = box->max[0];
	v[1][1] = box->min[1];
	v[1][2] = box->min[2];
	v[2][0] = box->min[0];
	v[2][1] = box->max[1];
	v[2][2] = box->min[2];
	v[3][0] = box->min[0];
	v[3][1] = box->min[1];
	v[3][2] = box->max[2];

	v[4][0] = box->max[0];
	v[4][1] = box->max[1];
	v[4][2] = box->min[2];

	v[5][0] = box->min[0];
	v[5][1] = box->max[1];
	v[5][2] = box->max[2];

	v[6][0] = box->max[0];
	v[6][1] = box->min[1];
	v[6][2] = box->max[2];

	v[7][0] = box->max[0];
	v[7][1] = box->max[1];
	v[7][2] = box->max[2];

	for (int i = 0; i < 8; i++) {
		vids[i] = addtovertices(process);
		copytovertices(process, v[i], n, vids[i]);
	}

	make_face(process, vids[0], vids[2], vids[4], vids[1]);
	make_face(process, vids[0], vids[3], vids[5], vids[2]);
	make_face(process, vids[0], vids[1], vids[6], vids[3]);
	make_face(process, vids[1], vids[4], vids[7], vids[6]);
	make_face(process, vids[2], vids[5], vids[7], vids[4]);
	make_face(process, vids[3], vids[6], vids[7], vids[5]);
}

#endif // 0

static void enqueue_cube(CHUNK *chunk, int i, int j, int k)
{
	CENTERLIST *newc;

	if (chunk != NULL) {
		//draw_cube(chunk, i, j, k);
		chunk->process->cubes_left_to_process = true;

		BLI_mutex_lock(&chunk->cube_queue_lock);
		newc = BLI_memarena_alloc(chunk->cube_queue_mem, sizeof(CENTERLIST));
		newc->lat[0] = i;
		newc->lat[1] = j;
		newc->lat[2] = k;
		newc->next = chunk->cube_queue;
		chunk->cube_queue = newc;

		BLI_mutex_unlock(&chunk->cube_queue_lock);
	}
}
/**
 * Adds cube at given lattice position to cube stack of process.
 */
static void add_cube(CHUNK *chunk, int i, int j, int k)
{
	CUBES *ncube;
	int n;

	if      (i == chunk->min_lat[0] - 1) enqueue_cube(chunk->neighbors[0], i, j, k);
	else if (i == chunk->max_lat[0] + 1) enqueue_cube(chunk->neighbors[1], i, j, k);
	else if (j == chunk->min_lat[1] - 1) enqueue_cube(chunk->neighbors[2], i, j, k);
	else if (j == chunk->max_lat[1] + 1) enqueue_cube(chunk->neighbors[3], i, j, k);
	else if (k == chunk->min_lat[2] - 1) enqueue_cube(chunk->neighbors[4], i, j, k);
	else if (k == chunk->max_lat[2] + 1) enqueue_cube(chunk->neighbors[5], i, j, k);
	else {
		/* test if cube has been found before */
		if (setcenter(chunk, chunk->centers, i, j, k) == 0) {
			/* push cube on stack: */
			ncube = BLI_memarena_alloc(chunk->mem, sizeof(CUBES));
			ncube->next = chunk->cubes;
			chunk->cubes = ncube;

			ncube->cube.lat[0] = i;
			ncube->cube.lat[1] = j;
			ncube->cube.lat[2] = k;

			/* set corners of initial cube: */
			for (n = 0; n < 8; n++) {
				ncube->cube.corners[n] = setcorner(chunk, i + MB_BIT(n, 2), j + MB_BIT(n, 1), k + MB_BIT(n, 0));
			}
		}
	}
}

static void do_queued_cubes(TaskPool *UNUSED(pool), void *data, int UNUSED(threadid))
{
	CHUNK *chunk = (CHUNK*)data;
	CUBE c;

	while (chunk->cube_queue != NULL) {
		add_cube(chunk, UNPACK3(chunk->cube_queue->lat));
		chunk->cube_queue = chunk->cube_queue->next;
	}

	while (chunk->cubes != NULL) {
		c = chunk->cubes->cube;
		chunk->cubes = chunk->cubes->next;

		docube(chunk, &c);
	}
}

#define INSIDE(a, lbn, rtf)	(a[0] >= lbn[0] && a[1] >= lbn[1] && a[2] >= lbn[2] && \
                             a[0] <= rtf[0] && a[1] <= rtf[1] && a[2] <= rtf[2])

/**
 * Find at most 26 cubes to start polygonization from.
 */
static void find_first_points(CHUNK *chunk, const unsigned int em)
{
	const MLSmall *ml;
	int center[3], lbn[3], rtf[3], it[3], dir[3], add[3];
	float tmp[3], a, b;

	ml = chunk->mainb[em];
	mid_v3_v3v3(tmp, ml->bb->min, ml->bb->max);
	closest_latice(center, tmp, chunk->process->size);
	next_lattice(lbn, ml->bb->min, chunk->process->size);
	prev_lattice(rtf, ml->bb->max, chunk->process->size);

	if (!INSIDE(center, chunk->min_lat, chunk->max_lat)) return;
	DO_MIN(chunk->max_lat, rtf);
	DO_MAX(chunk->min_lat, lbn);
	for (dir[0] = -1; dir[0] <= 1; dir[0]++)
		for (dir[1] = -1; dir[1] <= 1; dir[1]++)
			for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
				if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) continue;
				copy_v3_v3_int(it, center);

				b = setcorner(chunk, it[0], it[1], it[2])->value;
				do {
					VECADD(it, it, dir);
					a = b;
					b = setcorner(chunk, it[0], it[1], it[2])->value;
					if (a * b < 0.0f) {
						VECSUB(add, it, dir);
						DO_MIN(it, add);
						enqueue_cube(chunk, add[0], add[1], add[2]);
						break;
					}
				} while (INSIDE(it, lbn, rtf));
			}
}
static void precompute(CHUNK *chunk, const unsigned int em)
{
	CORNER *c;
	int index, it[3], lbn[3], rtf[3];
	MLSmall *ml = chunk->mainb[em];
	bool make;

	next_lattice(lbn, ml->bb->min, chunk->process->size);
	next_lattice(rtf, ml->bb->max, chunk->process->size);

	DO_MIN(chunk->max_lat, rtf);
	DO_MAX(chunk->min_lat, lbn);

	for (it[0] = lbn[0]; it[0] < rtf[0]; it[0]++)
		for (it[1] = lbn[1]; it[1] < rtf[1]; it[1]++)
			for (it[2] = lbn[2]; it[2] < rtf[2]; it[2]++) {
				make = true;
				index = HASH(it[0], it[1], it[2]);
				c = chunk->corners[index];

				for (; c != NULL; c = c->next) {
					if (equal_v3_v3_int(it, c->lat)) {
						c->value -= densfunc(ml, UNPACK3(c->co));
						make = false; 
						break;
					}
				}
				if (make) {
					c = BLI_memarena_alloc(chunk->mem, sizeof(CORNER));

					VECCOPY(c->lat, it);
					c->co[0] = ((float)it[0] - 0.5f) * chunk->process->size;
					c->co[1] = ((float)it[1] - 0.5f) * chunk->process->size;
					c->co[2] = ((float)it[2] - 0.5f) * chunk->process->size;
					c->value = chunk->process->thresh;// - 1.0f;
					c->value -= densfunc(ml, UNPACK3(c->co));

					c->next = chunk->corners[index];
					chunk->corners[index] = c;
				}
			}
}
static void init_chunk_bb(CHUNK *chunk)
{
	unsigned int i;
	float step[3];
	PROCESS *process = chunk->process;

	chunk->lat_pos[0] = chunk->index % (int)process->chunk_res;
	chunk->lat_pos[1] = (chunk->index / (int)process->chunk_res) % (int)process->chunk_res;
	chunk->lat_pos[2] = chunk->index / (int)process->chunk_res / (int)process->chunk_res;

	for (i = 0; i < 3; i++) {
		step[i] = (process->allbb.max[i] - process->allbb.min[i]) / (float)process->chunk_res;
		chunk->bb.min[i] = process->allbb.min[i] + step[i] * (float)chunk->lat_pos[i];
		chunk->bb.max[i] = process->allbb.min[i] + step[i] * (float)(chunk->lat_pos[i] + 1);
	}
	prev_lattice(chunk->max_lat, chunk->bb.max, process->size);
	next_lattice(chunk->min_lat, chunk->bb.min, process->size);
	for (i = 0; i < 3; i++) {
		if (chunk->lat_pos[i] == 0) chunk->min_lat[i]--;
		chunk->bb.max[i] = ((float)chunk->max_lat[i] + 0.5f) * process->size + process->delta * 2.0f;
		chunk->bb.min[i] = ((float)chunk->min_lat[i] - 0.5f) * process->size - process->delta * 2.0f;
	}
}

static void set_chunk_neighbors(CHUNK *chunk)
{
	int i, index, res, ne_lat[6][3];
	res = (int)chunk->process->chunk_res;

	for (i = 0; i < 6; i++) {
		VECCOPY(ne_lat[i], chunk->lat_pos);
	}
	ne_lat[0][0]--;
	ne_lat[1][0]++;
	ne_lat[2][1]--;
	ne_lat[3][1]++;
	ne_lat[4][2]--;
	ne_lat[5][2]++;

	for (i = 0; i < 6; i++) {
		if (ne_lat[i][0] < 0 || ne_lat[i][0] >= res ||
			ne_lat[i][1] < 0 || ne_lat[i][1] >= res ||
			ne_lat[i][2] < 0 || ne_lat[i][2] >= res)
		{
			chunk->neighbors[i] = NULL;
		}
		else {
			index = ne_lat[i][0] + ne_lat[i][1] * res + ne_lat[i][2] * res * res;
			chunk->neighbors[i] = chunk->process->chunks[index];
		}
	}
}

/**
 * Initialize chunk.
 */
static void init_chunk(TaskPool *UNUSED(pool), void *data, int UNUSED(threadid))
{
	unsigned int i;
	Box allbb;
	CHUNK *chunk = (CHUNK *)data;
	PROCESS *process = chunk->process;
	init_chunk_bb(chunk);
	set_chunk_neighbors(chunk);

	if (process->num_chunks > 1) {
		chunk->mainb = MEM_mallocN(sizeof(MLSmall *) * process->totelem, "chunk's metaballs");
		chunk->mem = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Metaball chunk memarena");

		for (i = 0; i < process->totelem; i++) {
			if (process->mainb[i]->bb->min[0] < chunk->bb.max[0] &&
				process->mainb[i]->bb->min[1] < chunk->bb.max[1] &&
				process->mainb[i]->bb->min[2] < chunk->bb.max[2] &&
				process->mainb[i]->bb->max[0] > chunk->bb.min[0] &&
				process->mainb[i]->bb->max[1] > chunk->bb.min[1] &&
				process->mainb[i]->bb->max[2] > chunk->bb.min[2])
			{
				chunk->mainb[chunk->elem++] = process->mainb[i];
			}
		}
	}
	else {
		chunk->mainb = process->mainb;
		chunk->mem = process->pgn_elements;
		chunk->elem = process->totelem;
	}
	if (chunk->elem > 0) {
		copy_v3_v3(allbb.min, chunk->mainb[0]->bb->min);
		copy_v3_v3(allbb.max, chunk->mainb[0]->bb->max);

		for (i = 1; i < chunk->elem; i++) {
			make_union(chunk->mainb[i]->bb, &allbb, &allbb);
		}

		build_bvh_spatial(chunk, &chunk->bvh, 0, chunk->elem, &allbb);

		chunk->centers = MEM_callocN(HASHSIZE * sizeof(CENTERLIST *), "mbproc->centers");
		chunk->corners = MEM_callocN(HASHSIZE * sizeof(CORNER *), "mbproc->corners");
		chunk->edges = MEM_callocN(2 * HASHSIZE * sizeof(EDGELIST *), "mbproc->edges");
		chunk->bvh_queue = MEM_callocN(sizeof(MetaballBVHNode *) * chunk->bvh_queue_size, "Metaball BVH Queue");
		chunk->cube_queue_mem = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Chunk's memory for additional cubes");
		chunk->cubes = NULL;

#ifdef MB_PRECOMPUTE
		for (i = 0; i < chunk->elem; i++) {
			precompute(chunk, i);
		}
#endif
		for (i = 0; i < chunk->elem; i++) {
			find_first_points(chunk, i);
		}
	}
}
/**
 * The main polygonization proc.
 * Initializes process, mutexes, makes cubetable
 * and queues num_chunks tasks.
 */
static void polygonize(PROCESS *process)
{
	unsigned int i;
	TaskScheduler *task_scheduler = BLI_task_scheduler_get();
	TaskPool *task_pool;

	process->chunks = MEM_callocN(sizeof(CHUNK *) * process->num_chunks, "mbproc->chunks");
	process->common_edges = MEM_callocN(2 * HASHSIZE * sizeof(EDGELIST *), "mbproc->edges");
	BLI_mutex_init(&process->edge_lock);
	makecubetable();

	for (i = 0; i < process->num_chunks; i++) {
		process->chunks[i] = MEM_callocN(sizeof(CHUNK), "Chunk");
		process->chunks[i]->index = (int)i;
		process->chunks[i]->process = process;
		BLI_mutex_init(&process->chunks[i]->cube_queue_lock);
	}
	
	task_pool = BLI_task_pool_create(task_scheduler, process);

	for (i = 0; i < process->num_chunks; i++) {
		BLI_task_pool_push(task_pool, init_chunk, (void*)process->chunks[i], false, TASK_PRIORITY_HIGH);
	}
	BLI_task_pool_work_and_wait(task_pool);

	while (process->cubes_left_to_process) {
		process->cubes_left_to_process = false;
		for (i = 0; i < process->num_chunks; i++) {
			BLI_task_pool_push(task_pool, do_queued_cubes, (void*)process->chunks[i], false, TASK_PRIORITY_HIGH);
		}
		BLI_task_pool_work_and_wait(task_pool);
	}

	BLI_task_pool_free(task_pool);
}

static DispList* make_displist(PROCESS *process)
{
	unsigned int i, j, k;
	DispList *dl;
	dl = MEM_callocN(sizeof(DispList), "mballdisp");
	dl->type = DL_INDEX4;
	dl->nr = (int)process->curvertex;

	/* count vertices from all chunks */
	for (i = 0; i < process->num_chunks; i++) {
		process->chunks[i]->vertex_offset = (unsigned int)dl->nr;
		dl->nr += (int)process->chunks[i]->curvertex;
		dl->parts += (int)process->chunks[i]->curindex;
	}

	/* copy vertex coordinates - first common, and then from chunks */
	dl->verts = MEM_callocN(sizeof(float[3]) * (size_t)dl->nr, "verts");
	for (i = 0; i < process->curvertex; i++) {
		copy_v3_v3(&dl->verts[i * 3], process->co[i]);
	}
	for (i = 0; i < process->num_chunks; i++) {
		for (j = 0, k = process->chunks[i]->vertex_offset; j < process->chunks[i]->curvertex; j++, k++) {
			copy_v3_v3(&dl->verts[k * 3], process->chunks[i]->co[j]);
		}
	}

	/* copy indices */
	dl->index = MEM_callocN(sizeof(int[4]) * (size_t)dl->parts, "indices");
	for (i = 0, k = 0; i < process->num_chunks; i++) {
		for (j = 0; j < process->chunks[i]->curindex; j++) {
			const int *ind = process->chunks[i]->indices[j];
			dl->index[k++] = ind[0] < 0 ? -(ind[0]) - 2 : ind[0] + (int)process->chunks[i]->vertex_offset;
			dl->index[k++] = ind[1] < 0 ? -(ind[1]) - 2 : ind[1] + (int)process->chunks[i]->vertex_offset;
			dl->index[k++] = ind[2] < 0 ? -(ind[2]) - 2 : ind[2] + (int)process->chunks[i]->vertex_offset;
			dl->index[k++] = ind[3] < 0 ? -(ind[3]) - 2 : ind[3] + (int)process->chunks[i]->vertex_offset;
		}
	}

	/* copy or make normals */
	dl->nors = MEM_callocN(sizeof(float[3]) * (size_t)dl->nr, "nors");
#ifdef MB_ACCUM_NORMAL
	for (i = 0; i < (unsigned int)dl->parts; i++) {
		float n[3];
		normal_quad_v3(n,
					   &dl->verts[dl->index[i * 4] * 3],
					   &dl->verts[dl->index[i * 4 + 1] * 3],
					   &dl->verts[dl->index[i * 4 + 2] * 3],
					   &dl->verts[dl->index[i * 4 + 3] * 3]);
		add_v3_v3(&dl->nors[dl->index[i * 4] * 3], n);
		add_v3_v3(&dl->nors[dl->index[i * 4 + 1] * 3], n);
		add_v3_v3(&dl->nors[dl->index[i * 4 + 2] * 3], n);
		add_v3_v3(&dl->nors[dl->index[i * 4 + 3] * 3], n);
	}
	for (i = 0; i < (unsigned int)dl->nr; i++) {
		normalize_v3(&dl->nors[i * 3]);
	}
#else
	for (i = 0; i < process->curvertex; i++) {
		normalize_v3(process->no[i]);
		copy_v3_v3(&dl->nors[i * 3], process->no[i]);
					}
	for (i = 0; i < process->num_chunks; i++) {
		for (j = 0, k = process->chunks[i]->vertex_offset; j < process->chunks[i]->curvertex; j++, k++) {
			normalize_v3(process->chunks[i]->no[j]);
			copy_v3_v3(&dl->nors[k * 3], process->chunks[i]->no[j]);
		}
	}
#endif

	return dl;
}
/**
 * Iterates over ALL objects in the scene and all of its sets, including
 * making all duplis(not only metas). Copies metas to mainb array.
 * Computes bounding boxes for building BVH. */
static void init_meta(EvaluationContext *eval_ctx, PROCESS *process, Scene *scene, Object *ob)
{
	Scene *sce_iter = scene;
	Base *base;
	Object *bob;
	MetaBall *mb;
	const MetaElem *ml;
	float obinv[4][4], obmat[4][4];
	unsigned int i;
	int obnr, zero_size = 0;
	char obname[MAX_ID_NAME];
	SceneBaseIter iter;

	copy_m4_m4(obmat, ob->obmat);   /* to cope with duplicators from BKE_scene_base_iter_next */
	invert_m4_m4(obinv, ob->obmat);

	BLI_split_name_num(obname, &obnr, ob->id.name + 2, '.');

	/* make main array */
	BKE_scene_base_iter_next(eval_ctx, &iter, &sce_iter, 0, NULL, NULL);
	while (BKE_scene_base_iter_next(eval_ctx, &iter, &sce_iter, 1, &base, &bob)) {
		if (bob->type == OB_MBALL) {
			zero_size = 0;
			ml = NULL;

			if (bob == ob && (base->flag & OB_FROMDUPLI) == 0) {
				mb = ob->data;

				if (mb->editelems) ml = mb->editelems->first;
				else ml = mb->elems.first;
			}
			else {
				char name[MAX_ID_NAME];
				int nr;

				BLI_split_name_num(name, &nr, bob->id.name + 2, '.');
				if (STREQ(obname, name)) {
					mb = bob->data;

					if (mb->editelems) ml = mb->editelems->first;
					else ml = mb->elems.first;
				}
			}

			/* when metaball object has zero scale, then MetaElem to this MetaBall
			 * will not be put to mainb array */
			if (has_zero_axis_m4(bob->obmat)) {
				zero_size = 1;
			}
			else if (bob->parent) {
				struct Object *pob = bob->parent;
				while (pob) {
					if (has_zero_axis_m4(pob->obmat)) {
						zero_size = 1;
						break;
					}
					pob = pob->parent;
				}
			}

			if (zero_size) {
				while (ml) {
					ml = ml->next;
				}
			}
			else {
				while (ml) {
					if (!(ml->flag & MB_HIDE)) {
						float pos[4][4], rot[4][4], mat[4][4];
						float expx, expy, expz;
						float tempmin[3], tempmax[3];

						BoundBox bb;
						MLSmall *new_ml;

						/* make a copy because of duplicates */
						new_ml = BLI_memarena_alloc(process->metaballs, sizeof(MLSmall));
						new_ml->bb = BLI_memarena_alloc(process->pgn_elements, sizeof(Box));
						new_ml->expx = ml->expx;
						new_ml->expy = ml->expy;
						new_ml->expz = ml->expz;
						new_ml->type = ml->type;

						/* too big stiffness seems only ugly due to linear interpolation
						* no need to have possibility for too big stiffness */
						if (ml->s > 10.0f) new_ml->s = 10.0f;
						else new_ml->s = ml->s;
						/* if metaball is negative, set stiffness negative */
						if (ml->flag & MB_NEGATIVE) new_ml->s = -new_ml->s;
						/* Translation of MetaElem */
						unit_m4(pos);
						pos[3][0] = ml->x;
						pos[3][1] = ml->y;
						pos[3][2] = ml->z;

						/* Rotation of MetaElem is stored in quat */
						quat_to_mat4(rot, ml->quat);

						/* basis object space -> world -> ml object space -> position -> rotation -> ml local space */
						mul_m4_series(mat, obinv, bob->obmat, pos, rot);
						/* ml local space -> basis object space */
						invert_m4_m4(new_ml->imat, mat);
						/* rad2 is inverse of squared radius */
						new_ml->rad2 = 1 / (ml->rad * ml->rad);
						/* initial dimensions = radius */
						expx = ml->rad;
						expy = ml->rad;
						expz = ml->rad;
						switch (ml->type) {
							case MB_BALL:
								break;
							case MB_CUBE: /* cube is "expanded" by expz, expy and expx */
								expz += ml->expz;
								/* fall through */
							case MB_PLANE: /* plane is "expanded" by expy and expx */
								expy += ml->expy;
								/* fall through */
							case MB_TUBE: /* tube is "expanded" by expx */
								expx += ml->expx;
								break;
							case MB_ELIPSOID: /* ellipsoid is "stretched" by exp* */
								expx *= ml->expx;
								expy *= ml->expy;
								expz *= ml->expz;
								break;
						}
						/* untransformed Bounding Box of MetaElem */
						/* TODO, its possible the elem type has been changed and the exp* values can use a fallback */
						copy_v3_fl3(bb.vec[0], -expx, -expy, -expz);  /* 0 */
						copy_v3_fl3(bb.vec[1], +expx, -expy, -expz);  /* 1 */
						copy_v3_fl3(bb.vec[2], +expx, +expy, -expz);  /* 2 */
						copy_v3_fl3(bb.vec[3], -expx, +expy, -expz);  /* 3 */
						copy_v3_fl3(bb.vec[4], -expx, -expy, +expz);  /* 4 */
						copy_v3_fl3(bb.vec[5], +expx, -expy, +expz);  /* 5 */
						copy_v3_fl3(bb.vec[6], +expx, +expy, +expz);  /* 6 */
						copy_v3_fl3(bb.vec[7], -expx, +expy, +expz);  /* 7 */

						/* transformation of Metalem bb */
						for (i = 0; i < 8; i++)
							mul_m4_v3(mat, bb.vec[i]);

						/* find max and min of transformed bb */
						INIT_MINMAX(tempmin, tempmax);
						for (i = 0; i < 8; i++) {
							DO_MINMAX(bb.vec[i], tempmin, tempmax);
						}

						/* set only point 0 and 6 - AABB of Metaelem */
						copy_v3_v3(new_ml->bb->min, tempmin);
						copy_v3_v3(new_ml->bb->max, tempmax);

						/* add new_ml to mainb[] */
						if (process->totelem == process->mem) {
							process->mem = process->mem * 2 + 10;
							process->mainb = MEM_reallocN(process->mainb, sizeof(MLSmall *) * process->mem);
						}
						process->mainb[process->totelem++] = new_ml;
					}
					ml = ml->next;
				}
			}
		}
	}

	/* compute AABB of all Metaelems */
	if (process->totelem > 0) {
		copy_v3_v3(process->allbb.min, process->mainb[0]->bb->min);
		copy_v3_v3(process->allbb.max, process->mainb[0]->bb->max);
		for (i = 1; i < process->totelem; i++)
			make_union(process->mainb[i]->bb, &process->allbb, &process->allbb);
	}
}

void BKE_mball_polygonize(EvaluationContext *eval_ctx, Scene *scene, Object *ob, ListBase *dispbase)
{
	MetaBall *mb = ob->data;
	PROCESS process = {0};
	int num_threads = BLI_system_thread_count();

	/* TODO: Test these numbers experimentally. 
	 * See what happens on a 64-core processor, maybe? */
	if (num_threads < 2) process.chunk_res = 1;
	else if (num_threads < 9) process.chunk_res = 2;
	else process.chunk_res = 3;

	process.num_chunks = (process.chunk_res * process.chunk_res * process.chunk_res);

	process.thresh = mb->thresh;

	if      (process.thresh < 0.001f) process.converge_res = 16;
	else if (process.thresh < 0.01f)  process.converge_res = 8;
	else if (process.thresh < 0.1f)   process.converge_res = 4;
	else                              process.converge_res = 2;

	if ((eval_ctx->mode != DAG_EVAL_RENDER) && (mb->flag == MB_UPDATE_NEVER)) return;
	if ((G.moving & (G_TRANSFORM_OBJ | G_TRANSFORM_EDIT)) && mb->flag == MB_UPDATE_FAST) return;

	if (eval_ctx->mode == DAG_EVAL_RENDER) {
		process.size = mb->rendersize;
	}
	else {
		process.size = mb->wiresize;
		if ((G.moving & (G_TRANSFORM_OBJ | G_TRANSFORM_EDIT)) && mb->flag == MB_UPDATE_HALFRES) {
			process.size *= 2.0f;
		}
	}

	process.delta = process.size * 0.001f;
	process.pgn_elements = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Metaball memarena");
	process.metaballs = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "Metaballs memarena");

	/* initialize all mainb (MetaElems) */
	init_meta(eval_ctx, &process, scene, ob);

	if (process.totelem > 0) {
		/* don't polygonize metaballs with too high resolution (base mball to small)
		* note: Eps was 0.0001f but this was giving problems for blood animation for durian, using 0.00001f */
		if (ob->size[0] > 0.00001f * (process.allbb.max[0] - process.allbb.min[0]) ||
		    ob->size[1] > 0.00001f * (process.allbb.max[1] - process.allbb.min[1]) ||
		    ob->size[2] > 0.00001f * (process.allbb.max[2] - process.allbb.min[2]))
		{
			polygonize(&process);

			BLI_addtail(dispbase, make_displist(&process));
		}
	}
	freeprocess(&process);
}
