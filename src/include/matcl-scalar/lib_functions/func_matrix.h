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

#include "matcl-scalar/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

// matrix multiplication x * y
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            mmul(const S1& x, const S2& y);

// multiplication x * y; equivalent to mmul function
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            operator*(const S1& x, const S2& y);

// for scalars the kron function is equivalent to mul
template<class S1, class S2, class Enable = typename md::enable_if_scalar2_ntobj<S1,S2,void>::type>
typename md::unify_types_promote<S1,S2>::type
                            kron(const S1& x, const S2& y);

};

#include "matcl-scalar/details/func_matrix.inl"