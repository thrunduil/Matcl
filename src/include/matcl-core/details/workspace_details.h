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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace details
{

#pragma warning(push)
#pragma warning(disable :4324) //structure was padded due to alignment specifier

class MATCL_CORE_EXPORT alignas(MATCL_SIMD_ALIGNMENT) workspace_base
{
    private:
        static const size_t MAX_ON_STACK    = 10 * sizeof(double);

    private:
        union data
        {
            char    m_data[MAX_ON_STACK + sizeof(Integer)];
            void*   m_ptr;
        };

        data        m_data;
        size_t      m_size_of;
        size_t      m_total_size_of;

    public:
        workspace_base();
        workspace_base(size_t size_of);

        workspace_base(workspace_base&&);

        ~workspace_base();

        size_t      size_of() const;
        void*       get_ptr();
        void        resize(size_t new_size_of);

    private:
        void        set_magic_number(void* ptr, size_t size_of);
        void        check_magic_number(void* ptr, size_t size_of);
};

#pragma warning(pop)

}}