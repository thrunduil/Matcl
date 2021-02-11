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
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-core/details/IO/io_impl.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

template<class val_type,class struct_type>
std::istream& load(std::istream&, Matrix<val_type,struct_type>&);

template<class val_type,class struct_type>
std::ostream& save(std::ostream&, const Matrix<val_type,struct_type>&);

std::ostream& save(std::ostream&, ti::ti_object ti);
std::istream& load(std::istream&, ti::ti_object& ti);

template<class Mat>
struct mm_helper_save
{    
    static void eval(std::ostream& os, const Mat& A, const std::string& comments);
};

struct mm_helper_load
{
    static void eval(std::istream& is, matcl::Matrix& ret, std::string& comments);
};

};};
