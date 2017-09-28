/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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
#include "matcl-core/general/exception.h"

#include <cassert>

namespace matcl { namespace dynamic { namespace details
{

#pragma warning(push)
#pragma warning(disable:4275)	//used as base for dll-interface class
#pragma warning(disable:4251)	//needs to have dll-interface to be used by clients

class MATCL_DYN_EXPORT fun_evaler_conv : public evaler
{
    private:
        using func_vec              = std::vector<function>;

    private:
        const evaler*               m_fun;		
        std::vector<function>       m_converters;
        std::vector<object>         m_deduced;
        Type                        m_deduced_ret;

    public:
        fun_evaler_conv(const evaler* fun, int n_deduced, const Type deduced[],
                        Type deduced_ret, const std::vector<function>& convs);

        bool make_eval(const object** _args, object& ret) const override;

        function make_converter(int, const Type[], Type, const func_vec&) const override
        {
            matcl_assert(false, "Should never be here!");
            return function();
        };
};

#pragma warning(pop)

};};};