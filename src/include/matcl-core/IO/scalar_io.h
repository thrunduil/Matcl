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

#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-core/options/options_disp.h"
#include "matcl-core/IO/base_io.h"

namespace matcl { namespace details
{

struct MATCL_CORE_EXPORT saveload_scalar_helper
{
    static std::ostream&    eval_save(std::ostream& os, Integer A);
    static std::ostream&    eval_save(std::ostream& os, Float A);
    static std::ostream&    eval_save(std::ostream& os, Real A);
    static std::ostream&    eval_save(std::ostream& os, const Complex& A);
    static std::ostream&    eval_save(std::ostream& os, const Float_complex& A);

    static std::istream&    eval_load(std::istream& is, Integer& A);
    static std::istream&    eval_load(std::istream& is, Float& A);
    static std::istream&    eval_load(std::istream& is, Real& A);
    static std::istream&    eval_load(std::istream& is, Complex& A);
    static std::istream&    eval_load(std::istream& is, Float_complex& A);
};

struct MATCL_CORE_EXPORT to_string_scalar_helper
{
    static std::string eval(Integer v);
    static std::string eval(Float v);
    static std::string eval(Real v);
    static std::string eval(const Complex& v);
    static std::string eval(const Float_complex& v);
    static std::string eval(const Object& v);

    template<class Ty>
    static std::string eval(const dynamic::object_type<Ty>& v)   
    {
        return eval(Object(v));
    };
};

}}
