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

#include "matcl-linalg/general/config_linalg.h"

#include <memory>

namespace matcl { namespace petsc
{

/// class which ensures that PetscInitialize and PetscFinalize and other 
/// global settings are setting during the lifetime of morfa_petsc library.
class petsc_initializer
{
    public:
        petsc_initializer();
        ~petsc_initializer();
};

/// class which locks the petsc lock; only one thread may access the petsc world,
/// but it may do so multiple times - hence recursive mutex is used
class petsc_lock_helper
{
    private:
        class petsc_lock_helper_impl;

        using ptr_type  = std::shared_ptr<petsc_lock_helper>;
        using impl_type = std::shared_ptr<petsc_lock_helper_impl>;

    public:
        petsc_lock_helper();
        ~petsc_lock_helper();

    private:
        impl_type   m_lock;

        petsc_lock_helper(const petsc_lock_helper&) = delete;
        petsc_lock_helper& operator=(const petsc_lock_helper&) = delete;
};

}};
