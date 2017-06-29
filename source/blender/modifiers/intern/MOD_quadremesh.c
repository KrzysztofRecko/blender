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

/** \file blender/modifiers/intern/MOD_quadremesh.c
 *  \ingroup modifiers
 */

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_string.h"
#include "BLI_rand.h"

#include "MEM_guardedalloc.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_particle.h"
#include "BKE_deform.h"

#include "MOD_util.h"

#ifdef WITH_MOD_QUADREMESH
#include "QRM_quadremesh.h"

static void initData(ModifierData *md)
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;
	lmd->anchor_grp_name[0] = '\0';
	lmd->flag = 0;

	lmd->cache_system = NULL;
	lmd->max_line_dist = QR_MINDIST;
	lmd->rng_seed = QR_SEED;
	lmd->auto_updates = QR_AUTOUPDATES;
}

static void copyData(ModifierData *md, ModifierData *target)
{
	QuadRemeshModifierData *qmd = (QuadRemeshModifierData *)md;
	
	modifier_copyData_generic(md, target);
	qmd->cache_system = NULL;
}

static bool isDisabled(ModifierData *md, int UNUSED(useRenderParams))
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;

	//if (!lmd->anchor_grp_name[0])
		//return true;

	return false;
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *md)
{
	QuadRemeshModifierData *lmd = (QuadRemeshModifierData *)md;
	CustomDataMask dataMask = 0;
	if (lmd->anchor_grp_name[0]) dataMask |= CD_MASK_MDEFORMVERT;
	return dataMask;
}

//#define QR_SHOW_UFIELD

static DerivedMesh *applyModifier(ModifierData *md,
	Object *ob,
	DerivedMesh *dm,
	ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result = dm;
	QuadRemeshModifierData *qmd = (QuadRemeshModifierData *)md;
	QuadRemeshSystem *sys;

	if (qmd->auto_updates) {
		qmd->flag |= MOD_QUADREMESH_INPUT_DIRTY;
	}

	if (!qmd->cache_system) {
		initQuadRemeshSystem(qmd);
	}

	sys = qmd->cache_system;

	result = makeResultMesh(sys, ob, dm);

#ifdef QR_SHOW_UFIELD
	{
		int defgrp_index;
		MDeformVert *dvert = NULL;

		if (!defgroup_find_name(ob, "QuadRemeshFlow")) {
			BKE_defgroup_new(ob, "QuadRemeshFlow");
		}

		modifier_get_vgroup(ob, dm, "QuadRemeshFlow", &dvert, &defgrp_index);
		getUField(sys, dvert, defgrp_index);
	}
#endif

	return result;
}

static void freeData(ModifierData *md)
{
	QuadRemeshModifierData *qmd = (QuadRemeshModifierData *)md;
	QuadRemeshSystem *sys = (QuadRemeshSystem *)qmd->cache_system;

	if (sys) {
		freeQuadRemeshSystem(sys);
	}
}

ModifierTypeInfo modifierType_QuadRemesh = {
	/* name */              "QuadRemesh",
	/* structName */        "QuadRemeshModifierData",
	/* structSize */        sizeof(QuadRemeshModifierData),
	/* type */              eModifierTypeType_Nonconstructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
							eModifierTypeFlag_AcceptsCVs,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          freeData,
	/* isDisabled */        isDisabled,
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};

#else

ModifierTypeInfo modifierType_QuadRemesh = {
	/* name */              "QuadRemesh",
	/* structName */        "QuadRemeshModifierData",
	/* structSize */        sizeof(QuadRemeshModifierData),
	/* type */              eModifierTypeType_Nonconstructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
							eModifierTypeFlag_AcceptsCVs,

	/* copyData */          NULL,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     NULL,
	/* applyModifierEM */   NULL,
	/* initData */          NULL,
	/* requiredDataMask */  NULL,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};

#endif