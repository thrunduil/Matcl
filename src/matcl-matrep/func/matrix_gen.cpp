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

#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/func/raw/mvgen.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/details/details_manip.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl
{

namespace mrd = matcl::raw::details;
namespace md = matcl::details;

//-----------------------------------------------------------------
//                      sequences
//-----------------------------------------------------------------
Matrix matcl::range(Real s, Real i, Real e)
{
    Matrix ret = matcl::Matrix(raw::range(s,i,e),false);
    return ret;
};

Matrix matcl::range(Real s, Real e)
{
    Matrix ret = matcl::Matrix(raw::range(s,e),false);
    return ret;
};

Matrix matcl::irange(Integer s, Integer i, Integer e)
{
    Matrix ret = matcl::Matrix(raw::irange(s,i,e),false);
    return ret;
};

Matrix matcl::irange(Integer s, Integer e)
{
    Matrix ret = matcl::Matrix(raw::irange(s,e),false);
    return ret;
};

Matrix matcl::frange(Float s, Float i, Float e)
{
    Matrix ret = matcl::Matrix(raw::frange(s,i,e),false);
    return ret;
};

Matrix matcl::frange(Float s, Float e)
{
    Matrix ret = matcl::Matrix(raw::frange(s,e),false);
    return ret;
};

Matrix matcl::linspace(Real s, Real e, Integer n)
{
    Matrix ret = matcl::Matrix(raw::linspace(s,e,n),false);
    return ret;
};

Matrix matcl::flinspace(Float s, Float e, Integer n)
{
    Matrix ret = matcl::Matrix(raw::flinspace(s,e,n),false);
    return ret;
};

Matrix matcl::logspace(Real s, Real e, Integer n)
{
    Matrix ret = matcl::Matrix(raw::logspace(s,e,n),false);
    return ret;
};

Matrix matcl::flogspace(Float s, Float e, Integer n)
{
    Matrix ret = matcl::Matrix(raw::flogspace(s,e,n),false);
    return ret;
};

//-----------------------------------------------------------------
//                      zeroes
//-----------------------------------------------------------------
Matrix matcl::zeros(Integer r, Integer c)
{
    return Matrix(raw::zeros(r,c),false);
};

Matrix matcl::zeros(ti::ti_object ti,Integer r, Integer c)
{
    return Matrix(raw::zeros(ti,r,c),false);
};

Matrix matcl::spzeros(Integer r, Integer c, Integer nnz)
{
    raw::real_sparse out(ti::ti_empty(),r,c,nnz);
    return Matrix(out,false);
};

Matrix matcl::spzeros(ti::ti_object ti,Integer r, Integer c,Integer nnz)
{
    raw::object_sparse out(ti,r,c,nnz);
    return Matrix(out,false);
};

Matrix matcl::bzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    raw::real_band out(ti::ti_empty(),0.0, r,c,fd,ld);
    
    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::bzeros(ti::ti_object ti,Integer r, Integer c,Integer fd, Integer ld)
{
    raw::object_band out(ti, md::default_value<Object>(ti), r,c,fd,ld);

    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);
    
    return Matrix(out,false);
};

Matrix matcl::izeros(Integer r, Integer c)
{
    return Matrix(raw::izeros(r,c),false);
};

Matrix matcl::ispzeros(Integer r, Integer c, Integer nnz)
{
    raw::integer_sparse out(ti::ti_empty(),r,c,nnz);
    return Matrix(out,false);
};

Matrix matcl::ibzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    raw::integer_band out(ti::ti_empty(),0, r,c,fd,ld);
    
    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::czeros(Integer r, Integer c)
{
    return Matrix(raw::czeros(r,c),false);
};

Matrix matcl::cspzeros(Integer r, Integer c, Integer nnz)
{
    raw::complex_sparse out(ti::ti_empty(),r,c,nnz);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::cbzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    raw::complex_band out(ti::ti_empty(), 0.0, r,c,fd,ld);
    
    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::fzeros(Integer r, Integer c)
{
    return Matrix(raw::fzeros(r,c),false);
};

Matrix matcl::fspzeros(Integer r, Integer c, Integer nnz)
{
    raw::float_sparse out(ti::ti_empty(),r,c,nnz);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::fbzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    raw::float_band out(ti::ti_empty(), 0.0f, r,c,fd,ld);
    
    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::fczeros(Integer r, Integer c)
{
    return Matrix(raw::fczeros(r,c),false);
};

Matrix matcl::fcspzeros(Integer r, Integer c, Integer nnz)
{
    raw::float_complex_sparse out(ti::ti_empty(),r,c,nnz);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::fcbzeros(Integer r, Integer c,Integer fd, Integer ld)
{
    raw::float_complex_band out(ti::ti_empty(), 0.0f, r,c,fd,ld);
    
    if (fd <= 0 && ld >= 0)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::zeros(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return izeros(r,c);
        case value_code::v_float:
            return fzeros(r,c);
        case value_code::v_real:
            return zeros(r,c);
        case value_code::v_float_complex:
            return fczeros(r,c);
        case value_code::v_complex:
            return czeros(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("zeros");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::spzeros(Integer r, Integer c,Integer nnz, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ispzeros(r,c,nnz);
        case value_code::v_float:
            return fspzeros(r,c,nnz);
        case value_code::v_real:
            return spzeros(r,c,nnz);
        case value_code::v_float_complex:
            return fcspzeros(r,c,nnz);
        case value_code::v_complex:
            return cspzeros(r,c,nnz);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("spzeros");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
}

Matrix matcl::bzeros(Integer r, Integer c,Integer fd, Integer ld, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ibzeros(r,c,fd,ld);
        case value_code::v_float:
            return fbzeros(r,c,fd,ld);
        case value_code::v_real:
            return bzeros(r,c,fd,ld);
        case value_code::v_float_complex:
            return fcbzeros(r,c,fd,ld);
        case value_code::v_complex:
            return cbzeros(r,c,fd,ld);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("bzeros");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

//-----------------------------------------------------------------
//                      ones
//-----------------------------------------------------------------
Matrix matcl::ones(Integer r, Integer c)
{
    return Matrix(raw::ones(r,c),false);
};

Matrix matcl::ones(ti::ti_object ti,Integer r, Integer c)
{
    return Matrix(raw::ones(ti,r,c),false);
};

Matrix matcl::iones(Integer r, Integer c)
{
    return Matrix(raw::iones(r,c),false);
};

Matrix matcl::spones(Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::ones(r,c)),false);
};

Matrix matcl::spones(ti::ti_object ti,Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::ones(ti,r,c)),false);
};

Matrix matcl::bones(Integer r, Integer c)
{
    return band(ones(r,c));
};

Matrix matcl::bones(ti::ti_object ti,Integer r, Integer c)
{
    return band(ones(ti,r,c));
};

Matrix matcl::ispones(Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::iones(r,c)),false);
};

Matrix matcl::ibones(Integer r, Integer c)
{
    return band(iones(r,c));
};

Matrix matcl::cones(Integer r, Integer c)
{
    return Matrix(raw::cones(r,c),false);
};

Matrix matcl::cspones(Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::cones(r,c)),false);
};

Matrix matcl::cbones(Integer r, Integer c)
{
    return band(cones(r,c));
};

Matrix matcl::fones(Integer r, Integer c)
{
    return Matrix(raw::fones(r,c),false);
};

Matrix matcl::fspones(Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::fones(r,c)),false);
};

Matrix matcl::fbones(Integer r, Integer c)
{
    return band(fones(r,c));
};

Matrix matcl::fcones(Integer r, Integer c)
{
    return Matrix(raw::fcones(r,c),false);
};

Matrix matcl::fcspones(Integer r, Integer c)
{
    return Matrix(raw::sparse(raw::fcones(r,c)),false);
};

Matrix matcl::fcbones(Integer r, Integer c)
{
    return band(fcones(r,c));
};

Matrix matcl::ones(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return iones(r,c);
        case value_code::v_float:
            return fones(r,c);
        case value_code::v_real:
            return ones(r,c);
        case value_code::v_float_complex:
            return fcones(r,c);
        case value_code::v_complex:
            return cones(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("ones");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::spones(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ispones(r,c);
        case value_code::v_float:
            return fspones(r,c);
        case value_code::v_real:
            return spones(r,c);
        case value_code::v_float_complex:
            return fcspones(r,c);
        case value_code::v_complex:
            return cspones(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("spones");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::bones(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ibones(r,c);
        case value_code::v_float:
            return fbones(r,c);
        case value_code::v_real:
            return bones(r,c);
        case value_code::v_float_complex:
            return fcbones(r,c);
        case value_code::v_complex:
            return cbones(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("bones");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

//-----------------------------------------------------------------
//                      eye
//-----------------------------------------------------------------
Matrix matcl::eye(Integer m, Integer n)
{
    return Matrix(raw::eye(m,n),false);
};

Matrix matcl::eye(Integer m)
{
    return Matrix(raw::eye(m),false);
};

Matrix matcl::eye(ti::ti_object ti,Integer m, Integer n)
{
    return Matrix(raw::eye(ti,m,n),false);
};

Matrix matcl::eye(ti::ti_object ti,Integer m)
{
    return Matrix(raw::eye(ti,m),false);
};

Matrix matcl::speye(Integer m, Integer n)
{
    return Matrix(raw::speye(m,n),false);
};

Matrix matcl::speye(Integer m)
{
    return Matrix(raw::speye(m),false);
};

Matrix matcl::speye(ti::ti_object ti,Integer m, Integer n)
{
    return Matrix(raw::speye(ti,m,n),false);
};

Matrix matcl::speye(ti::ti_object ti,Integer m)
{
    return Matrix(raw::speye(ti,m),false);
};

Matrix matcl::beye(Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::beye(m,n,fd,ld),false);
};

Matrix matcl::beye(Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::beye(m,fd,ld),false);
};

Matrix matcl::beye(ti::ti_object ti,Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::beye(ti,m,n,fd,ld),false);
};

Matrix matcl::beye(ti::ti_object ti,Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::beye(ti,m,fd,ld),false);
};

Matrix matcl::ieye(Integer m, Integer n)
{
    return Matrix(raw::ieye(m,n),false);
};

Matrix matcl::ieye(Integer m)
{
    return Matrix(raw::ieye(m),false);
};

Matrix matcl::ispeye(Integer m, Integer n)
{
    return Matrix(raw::ispeye(m,n),false);
};

Matrix matcl::ispeye(Integer m)
{
    return Matrix(raw::ispeye(m),false);
};

Matrix matcl::ibeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::ibeye(m,n,fd,ld),false);
};

Matrix matcl::ibeye(Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::ibeye(m,fd,ld),false);
};

Matrix matcl::ceye(Integer m, Integer n)
{
    return Matrix(raw::ceye(m,n),false);
};

Matrix matcl::ceye(Integer m)
{
    return Matrix(raw::ceye(m),false);
};

Matrix matcl::cspeye(Integer m, Integer n)
{
    return Matrix(raw::cspeye(m,n),false);
};

Matrix matcl::cspeye(Integer m)
{
    return Matrix(raw::cspeye(m),false);
};

Matrix matcl::cbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::cbeye(m,n,fd,ld),false);
};

Matrix matcl::cbeye(Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::cbeye(m,fd,ld),false);
};

Matrix matcl::feye(Integer m, Integer n)
{
    return Matrix(raw::feye(m,n),false);
};

Matrix matcl::feye(Integer m)
{
    return Matrix(raw::feye(m),false);
};

Matrix matcl::fspeye(Integer m, Integer n)
{
    return Matrix(raw::fspeye(m,n),false);
};

Matrix matcl::fspeye(Integer m)
{
    return Matrix(raw::fspeye(m),false);
};

Matrix matcl::fbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::fbeye(m,n,fd,ld),false);
};

Matrix matcl::fbeye(Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::fbeye(m,fd,ld),false);
};

Matrix matcl::fceye(Integer m, Integer n)
{
    return Matrix(raw::fceye(m,n),false);
};

Matrix matcl::fceye(Integer m)
{
    return Matrix(raw::fceye(m),false);
};

Matrix matcl::fcspeye(Integer m, Integer n)
{
    return Matrix(raw::fcspeye(m,n),false);
};

Matrix matcl::fcspeye(Integer m)
{
    return Matrix(raw::fcspeye(m),false);
};

Matrix matcl::fcbeye(Integer m, Integer n, Integer fd, Integer ld)
{
    return Matrix(raw::fcbeye(m,n,fd,ld),false);
};

Matrix matcl::fcbeye(Integer m, Integer fd, Integer ld)
{
    return Matrix(raw::fcbeye(m,fd,ld),false);
};

Matrix matcl::eye(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ieye(r,c);
        case value_code::v_float:
            return feye(r,c);
        case value_code::v_real:
            return eye(r,c);
        case value_code::v_float_complex:
            return fceye(r,c);
        case value_code::v_complex:
            return ceye(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("eye");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::eye(Integer m, matcl::value_code vt)
{
    return eye(m,m,vt);
};

Matrix matcl::speye(Integer r, Integer c, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ispeye(r,c);
        case value_code::v_float:
            return fspeye(r,c);
        case value_code::v_real:
            return speye(r,c);
        case value_code::v_float_complex:
            return fcspeye(r,c);
        case value_code::v_complex:
            return cspeye(r,c);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("speye");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::speye(Integer m, matcl::value_code vt)
{
    return speye(m,m,vt);
};

Matrix matcl::beye(Integer m, Integer n, Integer fd, Integer ld, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:
            return ibeye(m,n,fd,ld);
        case value_code::v_float:
            return fbeye(m,n,fd,ld);
        case value_code::v_real:
            return beye(m,n,fd,ld);
        case value_code::v_float_complex:
            return fcbeye(m,n,fd,ld);
        case value_code::v_complex:
            return cbeye(m,n,fd,ld);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("beye");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::beye(Integer m, Integer fd, Integer ld, matcl::value_code vt)
{
    return beye(m,m,fd,ld,vt);
};

//-----------------------------------------------------------------
//                      diag
//-----------------------------------------------------------------
Matrix matcl::diag(const Matrix& mat, Integer d)
{
    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::integer_dense>(),d),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::float_dense>(),d),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::real_dense>(),d),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::float_complex_dense>(),d),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::complex_dense>(),d),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diag(mat2.get_impl<raw::object_dense>(),d),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

Matrix matcl::bdiag(const Matrix &mat, Integer d)
{
    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::integer_dense>(),d),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::float_dense>(),d),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::real_dense>(),d),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::float_complex_dense>(),d),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::complex_dense>(),d),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiag(mat2.get_impl<raw::object_dense>(),d),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

Matrix matcl::spdiag(const Matrix &mat, Integer d)
{
    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::integer_dense>(),d),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::float_dense>(),d),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::real_dense>(),d),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::float_complex_dense>(),d),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::complex_dense>(),d),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiag(mat2.get_impl<raw::object_dense>(),d),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

Matrix matcl::diags(const Matrix &mat, const Matrix &d, Integer m, Integer n)
{
    const raw::integer_dense& mat_d = d.impl<raw::integer_dense>();

    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::integer_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::float_dense>(),mat_d,m,n),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::real_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::float_complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::diags(mat2.get_impl<raw::object_dense>(),mat_d,m,n),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

Matrix matcl::bdiags(const Matrix &mat, const Matrix &d, Integer m, Integer n)
{
    const raw::integer_dense& mat_d = d.impl<raw::integer_dense>();

    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::integer_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::float_dense>(),mat_d,m,n),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::real_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::float_complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::bdiags(mat2.get_impl<raw::object_dense>(),mat_d,m,n),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

Matrix matcl::spdiags(const Matrix &mat, const Matrix &d, Integer m, Integer n)
{
    const raw::integer_dense& mat_d = d.impl<raw::integer_dense>();

    switch (mat.get_value_code())
    {
        case value_code::v_integer:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::integer_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::float_dense>(),mat_d,m,n),false);
        }
        case value_code::v_real:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::real_dense>(),mat_d,m,n),false);
        }
        case value_code::v_float_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::float_complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_complex:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::complex_dense>(),mat_d,m,n),false);
        }
        case value_code::v_object:
        {
            Matrix mat2 = matcl::full(mat);
            return Matrix(raw::spdiags(mat2.get_impl<raw::object_dense>(),mat_d,m,n),false);
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

//-----------------------------------------------------------------
//                      rand
//-----------------------------------------------------------------

Matrix matcl::rand(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::rand(rand_ptr, m,n),false);
};

Matrix matcl::irand(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::irand(rand_ptr, m,n),false);
};

Matrix matcl::crand(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::crand(rand_ptr, m,n),false);
};

Matrix matcl::frand(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::frand(rand_ptr, m,n),false);
};

Matrix matcl::fcrand(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::fcrand(rand_ptr, m,n),false);
};

Matrix matcl::rand(Integer m, Integer n,matcl::value_code vt, const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            return irand(m,n, rand_ptr);
        case value_code::v_float:
            return frand(m,n, rand_ptr);
        case value_code::v_real:
            return rand(m,n, rand_ptr);
        case value_code::v_float_complex:
            return fcrand(m,n, rand_ptr);
        case value_code::v_complex:
            return crand(m,n, rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("rand");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::randn(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::randn(rand_ptr,m,n),false);
};

Matrix matcl::crandn(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::crandn(rand_ptr,m,n),false);
};

Matrix matcl::frandn(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::frandn(rand_ptr,m,n),false);
};

Matrix matcl::fcrandn(Integer m, Integer n, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::fcrandn(rand_ptr,m,n),false);
};

Matrix matcl::randn(Integer m, Integer n,matcl::value_code vt, const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            throw error::integer_value_type_not_allowed("randn");
        case value_code::v_float:
            return frandn(m,n,rand_ptr);
        case value_code::v_real:
            return randn(m,n,rand_ptr);
        case value_code::v_float_complex:
            return fcrandn(m,n,rand_ptr);
        case value_code::v_complex:
            return crandn(m,n,rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("randn");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::sprand(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::sprand(rand_ptr,m,n,d),false);
};

Matrix matcl::isprand(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::isprand(rand_ptr,m,n,d),false);
};

Matrix matcl::csprand(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::csprand(rand_ptr,m,n,d),false);
};

Matrix matcl::fsprand(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::fsprand(rand_ptr,m,n,d),false);
};

Matrix matcl::fcsprand(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::fcsprand(rand_ptr,m,n,d),false);
};

Matrix matcl::sprand(Integer m, Integer n, Real d,matcl::value_code vt, const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            return isprand(m,n,d,rand_ptr);
        case value_code::v_float:
            return fsprand(m,n,d,rand_ptr);
        case value_code::v_real:
            return sprand(m,n,d,rand_ptr);
        case value_code::v_float_complex:
            return fcsprand(m,n,d,rand_ptr);
        case value_code::v_complex:
            return csprand(m,n,d,rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("sprand");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::sprandn(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return Matrix(raw::sprandn(rand_ptr,m,n,d),false);
};

Matrix matcl::csprandn(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return matcl::Matrix(raw::csprandn(rand_ptr,m,n,d),false);
};

Matrix matcl::fsprandn(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return matcl::Matrix(raw::fsprandn(rand_ptr,m,n,d),false);
};

Matrix matcl::fcsprandn(Integer m, Integer n, Real d, const matcl::rand_state& rand_ptr)
{
    return matcl::Matrix(raw::fcsprandn(rand_ptr,m,n,d),false);
};

Matrix matcl::sprandn(Integer m, Integer n, Real d,matcl::value_code vt, const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            throw error::integer_value_type_not_allowed("sprandn");
        case value_code::v_float:
            return fsprandn(m,n,d,rand_ptr);
        case value_code::v_real:
            return sprandn(m,n,d,rand_ptr);
        case value_code::v_float_complex:
            return fcsprandn(m,n,d,rand_ptr);
        case value_code::v_complex:
            return csprandn(m,n,d,rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("sprandn");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

template<class Value, class Derived>
struct rand_band_helper
{
    using BM = matcl::raw::Matrix<Value,struct_banded>;

    static Matrix eval(Integer m, Integer n, Integer fd, Integer ld, 
                       const matcl::rand_state& rand_ptr)
    {
        BM res(ti::ti_empty(),m,n,fd,ld);

        Value* ptr_res = res.rep_ptr();

        Integer res_size = res.impl_size();
        for (Integer j = 0; j < res_size; ++j)
            ptr_res[j]      = Derived::rand(rand_ptr);

        if (fd == 0)
        {
            if (ld == 0)    res.get_struct().set(predefined_struct_type::diag);
            else            res.get_struct().set(predefined_struct_type::triu);
        }
        else if (ld == 0)
        {
            res.get_struct().set(predefined_struct_type::tril);
        };

        return Matrix(res,false);
    };
};

Matrix matcl::rand_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Real, rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::rand(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::irand_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Integer, rand_impl>
    {
        static Integer rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::irand(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::crand_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Complex, rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::crand(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::frand_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Float, rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::frand(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld, rand_ptr);
};

Matrix matcl::fcrand_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Float_complex, rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::fcrand(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::rand_band(Integer m, Integer n, Integer fd, Integer ld, matcl::value_code vt, 
                 const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            return irand_band(m,n,fd,ld,rand_ptr);
        case value_code::v_float:
            return frand_band(m,n,fd,ld,rand_ptr);
        case value_code::v_real:
            return rand_band(m,n,fd,ld,rand_ptr);
        case value_code::v_float_complex:
            return fcrand_band(m,n,fd,ld,rand_ptr);
        case value_code::v_complex:
            return crand_band(m,n,fd,ld,rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("rand_band");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::randn_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Real, rand_impl>
    {
        static Real rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::randn(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::crandn_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Complex, rand_impl>
    {
        static Complex rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::crandn(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::frandn_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Float, rand_impl>
    {
        static Float rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::frandn(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::fcrandn_band(Integer m, Integer n, Integer fd, Integer ld, const matcl::rand_state& rand_ptr)
{
    struct rand_impl : rand_band_helper<Float_complex, rand_impl>
    {
        static Float_complex rand(const matcl::rand_state& rand_ptr)
        {
            return matcl::fcrandn(rand_ptr);
        };
    };

    return rand_impl::eval(m,n,fd,ld,rand_ptr);
};

Matrix matcl::randn_band(Integer m, Integer n, Integer fd, Integer ld,matcl::value_code vt,
                  const matcl::rand_state& rand_ptr)
{
    switch (vt)
    {
        case value_code::v_integer:
            throw error::integer_value_type_not_allowed("randn_band");
        case value_code::v_float:
            return frandn_band(m,n,fd,ld,rand_ptr);
        case value_code::v_real:
            return randn_band(m,n,fd,ld,rand_ptr);
        case value_code::v_float_complex:
            return fcrandn_band(m,n,fd,ld,rand_ptr);
        case value_code::v_complex:
            return crandn_band(m,n,fd,ld,rand_ptr);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("randn_band");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

permvec matcl::randperm(Integer n, const matcl::rand_state& rand_ptr)
{
    Matrix p(raw::randperm(rand_ptr, n),false);
    return details::pv_constructor::make(p);
};

//-----------------------------------------------------------------
//                      make_matrix
//-----------------------------------------------------------------
Matrix matcl::make_integer_dense(Integer rows,Integer cols)
{
    raw::integer_dense out(ti::ti_empty(),Integer(),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_real_dense(Integer rows,Integer cols)
{
    raw::real_dense out(ti::ti_empty(),Real(),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_float_dense(Integer rows,Integer cols)
{
    raw::float_dense out(ti::ti_empty(),Float(),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_complex_dense(Integer rows,Integer cols)
{
    raw::complex_dense out(ti::ti_empty(),Complex(),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_float_complex_dense(Integer rows,Integer cols)
{
    raw::float_complex_dense out(ti::ti_empty(),Float_complex(),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_object_dense(ti::ti_object ti,Integer rows,Integer cols)
{
    raw::object_dense out(ti,Object(ti),rows,cols);
    out.get_struct().set(predefined_struct_type::diag);
    return Matrix(out,false);
};

Matrix matcl::make_integer_dense(Integer rows,Integer cols,const Integer *arr)
{
    return Matrix(raw::integer_dense(ti::ti_empty(),arr,rows,cols),false);
};

Matrix matcl::make_integer_dense(Integer rows,Integer cols,const Integer *arr, Integer ld)
{
    return Matrix(raw::integer_dense(ti::ti_empty(),arr,rows,cols,ld),false);
};

Matrix matcl::make_real_dense(Integer rows,Integer cols,const Real *arr)
{
    return Matrix(raw::real_dense(ti::ti_empty(),arr,rows,cols),false);
};

Matrix matcl::make_real_dense(Integer rows,Integer cols,const Real *arr, Integer ld)
{
    return Matrix(raw::real_dense(ti::ti_empty(),arr,rows,cols,ld),false);
};

Matrix matcl::make_float_dense(Integer rows,Integer cols,const Float *arr)
{
    return Matrix(raw::float_dense(ti::ti_empty(),arr,rows,cols),false);
};

Matrix matcl::make_float_dense(Integer rows,Integer cols,const Float *arr, Integer ld)
{
    return Matrix(raw::float_dense(ti::ti_empty(),arr,rows,cols,ld),false);
};

Matrix matcl::make_complex_dense(Integer rows,Integer cols,const Complex *arr)
{
    return Matrix(raw::complex_dense(ti::ti_empty(),arr,rows,cols),false);
};

Matrix matcl::make_complex_dense(Integer rows,Integer cols,const Complex *arr, Integer ld)
{
    return Matrix(raw::complex_dense(ti::ti_empty(),arr,rows,cols,ld),false);
};

Matrix matcl::make_float_complex_dense(Integer rows,Integer cols,const Float_complex *arr)
{
    return Matrix(raw::float_complex_dense(ti::ti_empty(),arr,rows,cols),false);
};

Matrix matcl::make_float_complex_dense(Integer rows,Integer cols,const Float_complex *arr, Integer ld)
{
    return Matrix(raw::float_complex_dense(ti::ti_empty(),arr,rows,cols, ld),false);
};

Matrix matcl::make_complex_dense(Integer rows,Integer cols,const Real *ar_r,const Real *ar_i)
{
    return Matrix(raw::complex_dense(ti::ti_empty(),ar_r,ar_i,rows,cols),false);
};

Matrix matcl::make_float_complex_dense(Integer rows,Integer cols,const Float *ar_r,const Float *ar_i)
{
    return Matrix(raw::float_complex_dense(ti::ti_empty(),ar_r,ar_i,rows,cols),false);
};

Matrix matcl::make_object_dense(ti::ti_object ti, Integer rows,Integer cols,const Object *arr)
{
    return Matrix(raw::object_dense(ti,arr,rows,cols),false);
};

Matrix matcl::make_object_dense(ti::ti_object ti, Integer rows,Integer cols,const Object *arr, Integer ld)
{
    return Matrix(raw::object_dense(ti,arr,rows,cols,ld),false);
};

template<class T, class Enable>
Matrix matcl::make_dense_foreign(Integer rows,Integer cols, T* arr, Integer ld)
{
    using Mat       = raw::Matrix<T,struct_dense>;
    using foreign   = typename Mat::foreign;
    return Matrix(Mat(ti::ti_empty(),arr,rows,cols,ld, foreign()),false);
};

Matrix matcl::make_dense_foreign(ti::ti_object ti, Integer rows,Integer cols, Object* arr, 
                                Integer ld)
{
    using Mat       = raw::Matrix<Object,struct_dense>;
    using foreign   = Mat::foreign;
    return Matrix(Mat(ti, arr, rows, cols, ld,foreign()),false);
};

Matrix matcl::make_integer_dense(Integer val,Integer rows,Integer cols)
{
    raw::integer_dense out(ti::ti_empty(),val,rows,cols);
    if (val == 0)
    {
        out.get_struct().set(predefined_struct_type::diag);
    };
    return Matrix(out,false);
};

Matrix matcl::make_real_dense(Real val,Integer rows,Integer cols)
{
    raw::real_dense out(ti::ti_empty(),val,rows,cols);
    if (val == 0.)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::make_float_dense(Float val,Integer rows,Integer cols)
{
    raw::float_dense out(ti::ti_empty(),val,rows,cols);
    if (val == 0.)
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::make_complex_dense(Complex val,Integer rows,Integer cols)
{
    raw::complex_dense out(ti::ti_empty(),val,rows,cols);
    if (mrd::is_zero(val))
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::make_float_complex_dense(Float_complex val,Integer rows,Integer cols)
{
    raw::float_complex_dense out(ti::ti_empty(),val,rows,cols);
    if (mrd::is_zero(val))
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::make_object_dense(const Object& val,Integer rows,Integer cols)
{
    raw::object_dense out(val.get_type(),val,rows,cols);
    if (mrd::is_zero(val))
        out.get_struct().set(predefined_struct_type::diag);

    return Matrix(out,false);
};

Matrix matcl::make_integer_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::integer_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Integer* ptr_res = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_res[i] = 0;

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_real_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::real_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Real* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = 0;

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_float_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::float_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Float* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = 0.0f;

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_complex_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::complex_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Complex* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = 0;

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_float_complex_band(Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::float_complex_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Float_complex* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = 0.0f;

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_object_band(ti::ti_object ti, Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::object_band tmp(ti,rows,cols,fd,ld);
    Object* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        mrd::reset_helper(ptr_tmp[i],Object(ti));

    if (fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_integer_band(Integer val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::integer_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Integer* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = val;

    if (val == 0 && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_real_band(Real val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::real_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Real* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = val;

    if (val == 0 && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_float_band(Float val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::float_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Float* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = val;
    
    if (val == 0.0f && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_complex_band(Complex val,Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::complex_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Complex* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = val;

    if (mrd::is_zero(val) && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_float_complex_band(Float_complex val,Integer rows,Integer cols, 
                                       Integer fd, Integer ld)
{
    raw::float_complex_band tmp(ti::ti_empty(),rows,cols,fd,ld);
    Float_complex* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        ptr_tmp[i] = val;

    if (mrd::is_zero(val) && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_object_band(const Object& val, Integer rows,Integer cols, Integer fd, Integer ld)
{
    raw::object_band tmp(val.get_type(),rows,cols,fd,ld);
    Object* ptr_tmp = tmp.rep_ptr();

    Integer tmp_size = tmp.impl_size();
    for (Integer i = 0; i < tmp_size; ++i)
        mrd::reset_helper(ptr_tmp[i],val);

    if (mrd::is_zero(val) && fd <= 0 && ld >= 0)
        tmp.get_struct().set(predefined_struct_type::diag);

    return Matrix(tmp,false);
};

Matrix matcl::make_integer_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::integer_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_real_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::real_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_float_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::float_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_complex_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::complex_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_float_complex_sparse(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::float_complex_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_object_sparse(ti::ti_object ti, Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::object_sparse(ti,rows,cols,nzmax),false);
};

Matrix matcl::make_integer_sparse(const Integer *trip_r, const Integer *trip_c, const Integer *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::integer_sparse(ti::ti_empty(),trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_real_sparse(const Integer *trip_r, const Integer *trip_c, const Real *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::real_sparse(ti::ti_empty(),trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_float_sparse(const Integer *trip_r, const Integer *trip_c, const Float *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::float_sparse(ti::ti_empty(),trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_complex_sparse(const Integer *trip_r, const Integer *trip_c, const Complex *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::complex_sparse(ti::ti_empty(),trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_complex_sparse(const Integer *trip_r, const Integer *trip_c, 
                                        const Real *trip_re, const Real *trip_im,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::complex_sparse(ti::ti_empty(),trip_r,trip_c,trip_re,trip_im,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_float_complex_sparse(const Integer *trip_r, const Integer *trip_c, 
                                        const Float_complex *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::float_complex_sparse(ti::ti_empty(),trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_float_complex_sparse(const Integer *trip_r, const Integer *trip_c, 
                                        const Float *trip_re, const Float *trip_im,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::float_complex_sparse(ti::ti_empty(),trip_r,trip_c,trip_re,trip_im,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_object_sparse(ti::ti_object ti, const Integer *trip_r, const Integer *trip_c, 
                                        const Object *trip_x,
                                        Integer r, Integer c, Integer nnz0, Integer nzmax)
{
    nzmax = std::max(nnz0, nzmax);
    return Matrix(raw::object_sparse(ti,trip_r,trip_c,trip_x,r,c,nnz0, nzmax),false);
};

Matrix matcl::make_sparse_matrix(const Matrix& trip_r, const Matrix& trip_c, const Matrix& trip_x,
                                        Integer r, Integer c, Integer nzmax0)
{
    Matrix rm = matcl::full(trip_r);
    Matrix cm = matcl::full(trip_c);
    Matrix xm = matcl::full(trip_x);

    const raw::integer_dense& ri = rm.impl<raw::integer_dense>().make_explicit();
    const raw::integer_dense& ci = cm.impl<raw::integer_dense>().make_explicit();

    if (ri.rows() != 1 && ri.cols() != 1)
        throw error::vector_required(ri.rows(),ri.cols());

    if (ci.rows() != 1 && ci.cols() != 1)
        throw error::vector_required(ci.rows(),ci.cols());

    if (xm.rows() != 1 && xm.cols() != 1)
        throw error::vector_required(xm.rows(),xm.cols());

    Integer size    = ri.size();
    Integer xsize   = imult(xm.rows(),xm.cols());
    Integer nz      = size;
    Integer nzmax   = std::max(nz,nzmax0);

    if (ci.size() != size || xsize != size)
        throw error::invalid_vectors_spmat(ri.size(),ci.size(),xsize);		

    switch(trip_x.get_value_code())
    {
        case value_code::v_integer:
        {
            raw::integer_dense x = xm.get_impl<raw::integer_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::integer_sparse(ti::ti_empty(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),false);
        }
        case value_code::v_float:
        {
            raw::float_dense x = xm.get_impl<raw::float_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::float_sparse(ti::ti_empty(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),false);
        }
        case value_code::v_real:
        {
            raw::real_dense x = xm.get_impl<raw::real_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::real_sparse(ti::ti_empty(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),false);
        }
        case value_code::v_float_complex:
        {
            raw::float_complex_dense x = xm.get_impl<raw::float_complex_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::float_complex_sparse(ti::ti_empty(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),
                          false);
        }
        case value_code::v_complex:
        {
            raw::complex_dense x = xm.get_impl<raw::complex_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::complex_sparse(ti::ti_empty(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),false);
        }
        case value_code::v_object:
        {
            raw::object_dense x = xm.get_impl<raw::object_dense>();
            x.assign_to_fresh(x.make_explicit());
            return Matrix(raw::object_sparse(xm.get_type(),ri.ptr(),ci.ptr(),x.ptr(),r,c,nz,nzmax),false);
        }
        default:
        {
            matcl_assert(0,"unknown value type");
            throw;
        }
    };
};

Matrix matcl::make_dense_matrix(Integer rows,Integer cols, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:     return make_integer_dense(rows,cols);
        case value_code::v_float:       return make_float_dense(rows,cols);
        case value_code::v_real:        return make_real_dense(rows,cols);
        case value_code::v_float_complex:   return make_float_complex_dense(rows,cols);
        case value_code::v_complex:     return make_complex_dense(rows,cols);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_dense_matrix");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::make_band_matrix(Integer rows,Integer cols, Integer fd, Integer ld, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:     return make_integer_band(rows,cols,fd,ld);
        case value_code::v_float:       return make_float_band(rows,cols,fd,ld);
        case value_code::v_real:        return make_real_band(rows,cols,fd,ld);
        case value_code::v_float_complex:   return make_float_complex_band(rows,cols,fd,ld);
        case value_code::v_complex:     return make_complex_band(rows,cols,fd,ld);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_band_matrix");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::make_sparse_matrix(Integer rows,Integer cols, Integer nzmax, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:     return make_integer_sparse(rows,cols,nzmax);
        case value_code::v_float:       return make_float_sparse(rows,cols,nzmax);
        case value_code::v_real:        return make_real_sparse(rows,cols,nzmax);
        case value_code::v_float_complex:   return make_float_complex_sparse(rows,cols,nzmax);
        case value_code::v_complex:     return make_complex_sparse(rows,cols,nzmax);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_sparse_matrix");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::make_integer_dense_noinit(Integer rows,Integer cols, Integer*& ptr_data)
{
    raw::integer_dense mat(ti::ti_empty(),rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_real_dense_noinit(Integer rows,Integer cols, Real*& ptr_data)
{
    raw::real_dense mat(ti::ti_empty(),rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_float_dense_noinit(Integer rows,Integer cols, Float*& ptr_data)
{
    raw::float_dense mat(ti::ti_empty(),rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_complex_dense_noinit(Integer rows,Integer cols, Complex*& ptr_data)
{
    raw::complex_dense mat(ti::ti_empty(),rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_float_complex_dense_noinit(Integer rows,Integer cols, Float_complex*& ptr_data)
{
    raw::float_complex_dense mat(ti::ti_empty(),rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_object_dense_noinit(ti::ti_object ti, Integer rows,Integer cols, Object*& ptr_data)
{
    raw::object_dense mat(ti,rows,cols);
    ptr_data = mat.ptr();
    return Matrix(std::move(mat),false);
};

Matrix matcl::make_integer_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::integer_band(ti::ti_empty(),rows,cols,fd, ld),false);
};

Matrix matcl::make_real_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::real_band(ti::ti_empty(),rows,cols,fd, ld),false);
};

Matrix matcl::make_float_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::float_band(ti::ti_empty(),rows,cols,fd,ld),false);
};

Matrix matcl::make_complex_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::complex_band(ti::ti_empty(),rows,cols,fd,ld),false);
};

Matrix matcl::make_float_complex_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::float_complex_band(ti::ti_empty(),rows,cols,fd,ld),false);
};

Matrix matcl::make_object_band_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer fd, Integer ld)
{
    return Matrix(raw::object_band(ti,rows,cols,fd,ld),false);
};

Matrix matcl::make_integer_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::integer_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_real_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::real_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_float_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::float_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_complex_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::complex_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_float_complex_sparse_noinit(Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::float_complex_sparse(ti::ti_empty(),rows,cols,nzmax),false);
};

Matrix matcl::make_object_sparse_noinit(ti::ti_object ti, Integer rows,Integer cols, Integer nzmax)
{
    return Matrix(raw::object_sparse(ti,rows,cols,nzmax),false);
};

Matrix matcl::make_dense_noinit(Integer rows,Integer cols, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:  
        {
            Integer* ptr;
            return make_integer_dense_noinit(rows,cols,ptr);
        }
        case value_code::v_float:
        {
            Float* ptr;
            return make_float_dense_noinit(rows,cols,ptr);
        }
        case value_code::v_real:
        {
            Real* ptr;
            return make_real_dense_noinit(rows,cols, ptr);
        }
        case value_code::v_float_complex:
        {
            Float_complex* ptr;
            return make_float_complex_dense_noinit(rows,cols, ptr);
        }
        case value_code::v_complex:  
        {
            Complex* ptr;
            return make_complex_dense_noinit(rows,cols, ptr);
        }
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_dense_noinit");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::make_band_noinit(Integer rows,Integer cols, Integer fd, Integer ld, 
                                     matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:  return make_integer_band_noinit(rows,cols,fd,ld);
        case value_code::v_float:    return make_float_band_noinit(rows,cols,fd,ld);
        case value_code::v_real:     return make_real_band_noinit(rows,cols,fd,ld);
        case value_code::v_float_complex: return make_float_complex_band_noinit(rows,cols,fd,ld);
        case value_code::v_complex:  return make_complex_band_noinit(rows,cols,fd,ld);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_band_noinit");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

Matrix matcl::make_sparse_noinit(Integer rows,Integer cols, Integer nzmax, matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:  return make_integer_sparse_noinit(rows,cols,nzmax);
        case value_code::v_float:    return make_float_sparse_noinit(rows,cols,nzmax);
        case value_code::v_real:     return make_real_sparse_noinit(rows,cols,nzmax);
        case value_code::v_float_complex:   return make_float_complex_sparse_noinit(rows,cols,nzmax);
        case value_code::v_complex:  return make_complex_sparse_noinit(rows,cols,nzmax);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("make_sparse_noinit");
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

template MATCL_MATREP_EXPORT
Matrix matcl::make_dense_foreign<Integer,void>(Integer rows,Integer cols, Integer* arr, Integer ld);

template MATCL_MATREP_EXPORT
Matrix matcl::make_dense_foreign<Float,void>(Integer rows,Integer cols, Float* arr, Integer ld);

template MATCL_MATREP_EXPORT
Matrix matcl::make_dense_foreign<Real,void>(Integer rows,Integer cols, Real* arr, Integer ld);

template MATCL_MATREP_EXPORT
Matrix matcl::make_dense_foreign<Complex,void>(Integer rows,Integer cols, Complex* arr, Integer ld);

template MATCL_MATREP_EXPORT
Matrix matcl::make_dense_foreign<Float_complex,void>(Integer rows,Integer cols, Float_complex* arr, Integer ld);

};
