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
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/details/exception_details.h"
#include <string>
#include <memory>

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl { namespace error
{

// Thread safe wrapper of std::string to use with matcl::exception message handling
// Insures exclusive writing to std::string
class MATCL_CORE_EXPORT safe_string_message
{
    private:
        class safe_string_message_impl;
        std::shared_ptr<safe_string_message_impl> impl;

    public:
        // create empty string
        safe_string_message();
        // copy content from other
        safe_string_message(const safe_string_message& other);
    
        // copy content from other
        safe_string_message&    operator=(const safe_string_message& other);
        // copy content from rhs
        safe_string_message&    operator=(const std::string& rhs);
        // append current string with rhs
        safe_string_message&    operator+=(const std::string& rhs);
        // get char pointer
        const char*             c_str();
};

}}

#pragma warning(pop)