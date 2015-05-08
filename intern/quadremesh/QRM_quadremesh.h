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

/**
 * Known issues:
 *  - caching the u-field
 *  - semi automatic extrema placement
 *  - mesh anisotropy
 */

/** \file intern/quadremesh/MOD_quadremesh.h
 *  \ingroup quadremesh
 */

#ifndef __MOD_QUADREMESH_H__
#define __MOD_QUADREMESH_H__

/* Default values for the modifier */
#define QR_MINDIST 0.04f

typedef struct QuadRemeshSystem QuadRemeshSystem;

void initQuadRemeshSystem(QuadRemeshModifierData *qmd);
struct DerivedMesh *makeResultMesh(QuadRemeshSystem *sys, Object *ob, DerivedMesh *dm);
void freeQuadRemeshSystem(QuadRemeshSystem *sys);

void getUField(QuadRemeshSystem *sys, MDeformVert *dvert, int defgrp_index);

#endif /*__MOD_QUADREMESH_H__*/