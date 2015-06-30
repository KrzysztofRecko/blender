# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****

set(INC
	.
	../../source/blender/blenkernel
	../../source/blender/blenlib
	../../source/blender/makesdna
	../../source/blender/modifiers/intern
	../../source/blender/modifiers
	../../source/blender/depsgraph
	../opennl/extern
	../guardedalloc
)

set(INC_SYS
	${EIGEN3_INCLUDE_DIRS}
)

set(SRC
	QRM_quadremesh.c
	quadremesh_util.c
	quadremesh_line.c
	quadremesh_link.c
	quadremesh_edge.c
	quadremesh_input.c
	quadremesh_harmonic.c
	quadremesh_curvature.cpp
	
	QRM_quadremesh.h
	quadremesh_util.h
	quadremesh_input.h
	quadremesh_system.h
	quadremesh_line.h
)

blender_add_lib(bf_intern_quadremesh "${SRC}" "${INC}" "${INC_SYS}")

