/*
 * This file is a part of Matrix Computation Library (MATCL)
 *
 * Copyright (c) Pawe³ Kowal 2017 - 2018
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#define INLINE_TYPE
#include "blas.inl"

#include "matcl-blas-lapack/blas_loader/blas_loader.h"

namespace matcl 
{

std::string lapack::loaded_blas_plugin_name()
{
    return raw_blas_lapack::loaded_blas_plugin_name();
}

bool lapack::load_blas_plugin(const std::string& path)
{
    return raw_blas_lapack::load_blas_plugin(path);
}

void lapack::initialize_plugin()
{
    return raw_blas_lapack::initialize_plugin();
};

};