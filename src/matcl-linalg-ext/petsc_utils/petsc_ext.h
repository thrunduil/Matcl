/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#pragma once

#if 0
TODO

#include "matcl-linalg/general/config_linalg.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)

#include "petsc.h"
#include "petscmat.h"

#pragma warning(pop)

PetscErrorCode PCKSPGetKSPRef(PC, KSP**);
PetscErrorCode PCBJacobiGetSubKSPRef(PC,PetscInt*,PetscInt*,KSP**[]);
PetscErrorCode PCASMGetSubKSPRef(PC,PetscInt*,PetscInt*,KSP**[]);
PetscErrorCode PCFieldSplitGetSubKSPRef(PC,PetscInt*,KSP**[]);
PetscErrorCode PCGalerkinGetKSPRef(PC,KSP **);

PetscErrorCode PCMGGetSolverDownRef(PC, PetscInt lev, KSP**);
PetscErrorCode PCMGGetSolverUpRef(PC, PetscInt lev,KSP**);
PetscErrorCode PCMGDestroySolverDownKSP(PC, PetscInt lev, KSP* ksp);
PetscErrorCode PCMGDestroySolverUpKSP(PC, PetscInt lev, KSP* ksp);

#endif