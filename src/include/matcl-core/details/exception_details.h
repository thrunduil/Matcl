/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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
#include "matcl-core/general/fwd_decls.h"

namespace matcl { namespace details
{

#pragma warning(push)
#pragma warning(disable: 4251)  //needs to have dll-interface to be used by clients of class

class MATCL_CORE_EXPORT enable_warnings_raii
{
    private:
        int         m_old_state;
        bool        m_global;

    public:
        enable_warnings_raii(bool val, bool global);
        ~enable_warnings_raii();

        enable_warnings_raii(enable_warnings_raii&& other);

    private:
        enable_warnings_raii(const enable_warnings_raii&) = delete;
        enable_warnings_raii& operator=(const enable_warnings_raii&&) = delete;
};

#pragma warning(pop)

};};