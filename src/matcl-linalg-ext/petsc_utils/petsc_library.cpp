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

#include "petsc_library.h"

#include "matcl-linalg/general/mmlib_petsc_exception.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)

#include <petsc.h>
#include "petscerror.h"
#include "petscsys.h"

#pragma warning(pop)

#include <iostream>
#include <sstream>
#include <mutex>

namespace matcl { namespace petsc
{

//--------------------------------------------------------------
//                  petsc_initializer
//--------------------------------------------------------------

static const petsc_initializer g_petsc_init;

/// callback function to use with Petsc error handler framework - throws InternalPetscError excpetions
PetscErrorCode petsc_error_handler(MPI_Comm, int line, const char *func,
                    const char *file, PetscErrorCode n, PetscErrorType p, const char *mess,void *ctx)
{
    (void)ctx;
    (void)p;
    (void)n;
    (void)file;
    (void)line;
    (void)func;

    std::ostringstream os;
    os  << "petsc error: " << mess;

    throw error::internal_petsc_error(os.str());
}


petsc_initializer::petsc_initializer()
{
    char help_message;
    PetscInitialize(0,0,0, &help_message);
    PetscPopSignalHandler();
    PetscPushErrorHandler(&petsc_error_handler, nullptr);
};

petsc_initializer::~petsc_initializer()
{
    PetscFinalize();
};

//--------------------------------------------------------------
//                  petsc_lock_helper
//--------------------------------------------------------------
class petsc_lock_helper::petsc_lock_helper_impl
{
    private:
        using try_lock  = std::unique_lock<std::recursive_mutex>;

    public:
        petsc_lock_helper_impl(std::recursive_mutex& mutex_in)
            : m_lock(mutex_in, std::try_to_lock)
        {}

        bool owns_lock() const { return m_lock.owns_lock();}

    private:
        try_lock    m_lock;
};

std::recursive_mutex petsc_mutex;

petsc_lock_helper::petsc_lock_helper()
    : m_lock(new petsc_lock_helper_impl(petsc_mutex))
{
    if (!m_lock->owns_lock())
    {
        // attempted to lock the Petsc world but another thread already locked it
        throw error::multithreading_not_allowed();
    }
}
petsc_lock_helper::~petsc_lock_helper()
{}

}};

#endif