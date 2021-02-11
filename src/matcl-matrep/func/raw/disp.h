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

#include <iostream>
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/IO/disp_data_provider.h"

#include "matcl-core/details/IO/disp_impl.h"

namespace matcl { namespace raw
{

template<class value_type>
void disp(const disp_stream_ptr& os, const Matrix<value_type,struct_dense>& m);

template<class value_type>
void disp(const disp_stream_ptr& os, const Matrix<value_type,struct_sparse>& m);

template<class value_type>
void disp(const disp_stream_ptr& os, const Matrix<value_type,struct_banded>& m);

template<class value_type>
void disp(const Matrix<value_type,struct_dense>& v);

template<class value_type>
void disp(const Matrix<value_type,struct_sparse>& v);

template<class value_type>
void disp(const Matrix<value_type,struct_banded>& v);

template<class value_type>
void disp(const Matrix<value_type,struct_dense>& v)
{
    return raw::disp(default_disp_stream(),v); 
};

template<class value_type>
void disp(const Matrix<value_type,struct_sparse>& v)
{
    return raw::disp(default_disp_stream(),v); 
};

template<class value_type>
void disp(const Matrix<value_type,struct_banded>& v)
{
    return raw::disp(default_disp_stream(),v); 
};

}};