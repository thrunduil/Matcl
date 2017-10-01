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
#include "matcl-core/IO/archive.h"
#include "matcl-core/IO/scalar_io.h"

namespace matcl
{

template<class S>
inline typename md::enable_if_scalar<S, void>::type
disp(const S& A, const disp_stream_ptr& os, const options& opts)
{
    using SP = typename details::promote_scalar<S>::type;
    return details::disp_impl(SP(A),os,opts);
};

template<class S>
inline typename md::enable_if_scalar<S, void>::type
disp_header(const S& A, const disp_stream_ptr& os, const options& opts)
{
    options opts2 = opts;
    opts2.set(matcl::opt::disp::header_only(true));
    return disp(A, os, opts2);
};

template<class T> inline 
typename details::enable_if_scalar<T,void>::type
load(iarchive& ar,T& A)
{
    Matrix tmp;
    load(ar,tmp);

    using TP    = typename details::promote_scalar<T>::type;
    A           = tmp.get_scalar<TP>();
};

template<class T> inline 
typename details::enable_if_scalar<T,void>::type
save(oarchive& ar,const T& A)
{
    return save(ar,Matrix(A));
};

template<class T> inline
typename details::enable_if_matcl_scalar<T,std::ostream&>::type
operator<<(std::ostream& os, const T& A)
{
    return details::saveload_scalar_helper::eval_save(os, A);
};

template<class T> inline
typename details::enable_if_matcl_scalar<T,std::istream&>::type
operator>>(std::istream& is, T& A)
{
    return details::saveload_scalar_helper::eval_load(is, A);
};

template<class S, class Enable> 
inline std::string matcl::to_string(const S& A)
{
    return details::to_string_scalar_helper::eval(A);
};

};