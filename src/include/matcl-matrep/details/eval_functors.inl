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

#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/details/val_struct_codes.h"


namespace matcl
{

template<class derived>
inline bool test_function_templ<derived>::eval(Integer val) const
{
    return this->template eval_impl<Integer>(val);
};

template<class derived>
inline bool test_function_templ<derived>::eval(Real val) const
{
    return this->template eval_impl<Real>(val);
};

template<class derived>
inline bool test_function_templ<derived>::eval(Float val) const
{
    return this->template eval_impl<Float>(val);
};

template<class derived>
inline bool test_function_templ<derived>::eval(const Complex& val) const
{
    return this->template eval_impl<Complex>(val);
};

template<class derived>
inline bool test_function_templ<derived>::eval(const Float_complex& val) const
{
    return this->template eval_impl<Float_complex>(val);
};

template<class derived>
inline bool test_function_templ<derived>::eval(const Object& val) const
{
    return this->template eval_impl<Object>(val);
};

};