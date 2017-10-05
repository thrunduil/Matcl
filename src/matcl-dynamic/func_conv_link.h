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

#include "matcl-dynamic/function.h"

#include <cassert>

namespace matcl { namespace dynamic { namespace details
{

class fun_conv_link : public evaler
{
    private:	
        std::vector<function>       m_converters;

        fun_conv_link(const fun_conv_link&) = delete;
        fun_conv_link& operator=(const fun_conv_link&) = delete;

    public:
        fun_conv_link(const std::vector<function>& convs);
        ~fun_conv_link() override;

       bool         make_eval(const object** _args, object& ret) const override;
       void         make_eval(const object** _args) const;
       function     make_converter(int n_deduced, const Type deduced[], Type ded_ret,
                        const std::vector<function>&) const override;
};

};};};