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

#if 0
TODO

#include "petsc_ext.h"
#include <../src/ksp/pc/impls/bjacobi/bjacobi.h>
#include <../src/ksp/pc/impls/asm/asm.h>
#include <../src/ksp/pc/impls/fieldsplit/fieldsplit.h>
#include <../src/ksp/pc/impls/galerkin/galerkin.h>
#include <../src/ksp/pc/impls/ksp/pcksp.h>
#include <petsc/private/pcmgimpl.h>

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant

PetscErrorCode PCKSPGetKSPRef(PC pc, KSP** ksp_ref)
{
    PC_KSP* jac = (PC_KSP*)pc->data;

    if (!jac->ksp)
        PCKSPCreateKSP_KSP(pc);

    *ksp_ref = &jac->ksp;
    return 0;
};

PetscErrorCode PCBJacobiGetSubKSPRef(PC pc,PetscInt *n_local,PetscInt *first_local,KSP ***ksp_ref)
{
    PC_BJacobi *jac = (PC_BJacobi*)pc->data;;

    if (!pc->setupcalled)
        SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_WRONGSTATE,"Must call KSPSetUp() or PCSetUp() first");

    if (n_local) 
        *n_local = jac->n_local;
    if (first_local) 
        *first_local = jac->first_local;

    *ksp_ref  = &jac->ksp;
    jac->same_local_solves = PETSC_FALSE;        /* Assume that local solves are now different;
                                                  not necessarily true though!  This flag is
                                                  used only for PCView_BJacobi() */
    PetscFunctionReturn(0);
};

PetscErrorCode PCASMGetSubKSPRef(PC pc,PetscInt *n_local,PetscInt *first_local,KSP ***ksp_ref)
{
    PC_ASM *osm     = (PC_ASM*)pc->data;

    PetscErrorCode ierr;
    
    if (osm->n_local_true < 1)
    {
        SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ORDER,
              "Need to call PCSetUP() on PC (or KSPSetUp() on the outer KSP object) before calling here");
    }

    if (n_local) 
        *n_local = osm->n_local_true;
  
    if (first_local) 
    {
        ierr          = MPI_Scan(&osm->n_local_true,first_local,1,MPIU_INT,MPI_SUM,PetscObjectComm((PetscObject)pc));
        CHKERRQ(ierr);
    
        *first_local -= osm->n_local_true;
    }

    if (ksp_ref) 
    {
        /* Assume that local solves are now different; not necessarily
           true though!  This flag is used only for PCView_ASM() */
        *ksp_ref    = &osm->ksp;
        osm->same_local_solves = PETSC_FALSE;
    }
    PetscFunctionReturn(0);
}

PetscErrorCode PCFieldSplitGetSubKSPRef(PC pc,PetscInt *n,KSP *** subksp_ref)
{
  PC_FieldSplit     *jac = (PC_FieldSplit*)pc->data;
  PetscErrorCode    ierr;
  PetscInt          cnt   = 0;
  PC_FieldSplitLink ilink = jac->head;

  ierr = PetscMalloc1(jac->nsplits,subksp_ref);
  CHKERRQ(ierr);

  while (ilink) 
  {
    (*subksp_ref)[cnt++]    = &ilink->ksp;
    ilink                   = ilink->next;
  }
  if (cnt != jac->nsplits) 
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Corrupt PCFIELDSPLIT object: number of splits in"
               " linked list %D does not match number in object %D",cnt,jac->nsplits);

  if (n)
      *n = jac->nsplits;

  PetscFunctionReturn(0);
}

PetscErrorCode PCMGGetSolverUpRef(PC pc, PetscInt lev, KSP ** ksp_ref)
{
    PC_MG* mg               = (PC_MG*)pc->data;
    PC_MG_Levels**mglevels  = mg->levels;

    *ksp_ref =  &mglevels[lev]->smoothu;
    PetscFunctionReturn(0);
}
PetscErrorCode PCMGGetSolverDownRef(PC pc, PetscInt lev, KSP ** ksp_ref)
{
    PC_MG* mg               = (PC_MG*)pc->data;
    PC_MG_Levels**mglevels  = mg->levels;

    *ksp_ref =  &mglevels[lev]->smoothd;
    PetscFunctionReturn(0);
}

PetscErrorCode PCMGDestroySolverUpKSP(PC pc, PetscInt lev, KSP* ksp)
{
    PC_MG* mg               = (PC_MG*)pc->data;
    PC_MG_Levels**mglevels  = mg->levels;

    if (mglevels[lev]->smoothd == *ksp)
        PetscFunctionReturn(0);

    KSPDestroy(ksp);
    PetscFunctionReturn(0);
}
PetscErrorCode PCMGDestroySolverDownKSP(PC pc, PetscInt lev, KSP* ksp)
{
    PC_MG* mg               = (PC_MG*)pc->data;
    PC_MG_Levels**mglevels  = mg->levels;

    if (mglevels[lev]->smoothu == *ksp)
        PetscFunctionReturn(0);

    KSPDestroy(ksp);
    PetscFunctionReturn(0);
}


PetscErrorCode PCGalerkinGetKSPRef(PC pc,KSP ** ksp_ref)
{
    PC_Galerkin *jac    = (PC_Galerkin*)pc->data;
    *ksp_ref            = &jac->ksp;
    PetscFunctionReturn(0);
}

#pragma warning(pop)

#endif