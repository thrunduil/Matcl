/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-scalar/config.h"
#include "matcl-core/matrix/scalar_types.h"

#include <memory>

namespace matcl { namespace details
{
    class local_rand_state_raii;
    class rand_state;
};};

namespace matcl
{
    MATCL_SCALAR_EXPORT 
    details::local_rand_state_raii local_rand_state(unsigned long s);
};

namespace matcl { namespace details
{

#pragma warning(push)
#pragma warning(disable: 4251)  //needs to have dll-interface to be used by clients of class

class MATCL_SCALAR_EXPORT local_rand_state_raii
{
    private:
        std::shared_ptr<class local_rand_state_raii_impl>   m_impl;

    private:
        local_rand_state_raii(unsigned long s);

        local_rand_state_raii(const local_rand_state_raii&) = delete;
        local_rand_state_raii& operator=(const local_rand_state_raii&&) = delete;

        friend details::local_rand_state_raii matcl::local_rand_state(unsigned long s);

    public:
        ~local_rand_state_raii();

        local_rand_state_raii(local_rand_state_raii&&);
};

#pragma warning(pop)

};};