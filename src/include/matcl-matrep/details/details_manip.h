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

#include "matcl-matrep/lib_functions/manip.h"

namespace matcl { namespace details
{

template<class T, class S>
struct convert_scalar_imp
{
    static T eval(const S& val);
};

template<class T,class S> 
typename details::enable_if_scalar2<T,S,typename details::promote_scalar<T>::type>::type
convert_scalar(const S& val)
{
    using TP = typename details::promote_scalar<T>::type;
    using SP = typename details::promote_scalar<S>::type;
    return details::convert_scalar_imp<TP,SP>::eval(val);
};

}}