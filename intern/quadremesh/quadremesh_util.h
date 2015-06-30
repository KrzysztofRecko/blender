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



void addSeedToQueue(struct Heap *aheap, float in_co[3], GFPointType in_type, int in_val, float weight);
GFPoint *getTopSeedFromQueue(struct Heap *aheap);

#endif