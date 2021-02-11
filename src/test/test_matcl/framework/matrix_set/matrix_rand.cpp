/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matrix_rand.h"
#include "rand_sp.h"
#include "rand_str.h"
#include "matcl-matrep/lib_functions/manip.h"


namespace matcl { namespace test
{

static Integer max_int = 1000;

//======================================================================
//                      rand_matrix_1
//======================================================================
Integer rand_matrix_1::rand_scalar_int()
{
    return irand() % max_int;
};

Real rand_matrix_1::rand_scalar_real()
{
    return randn();
};

Float rand_matrix_1::rand_scalar_float()
{
    return frandn();
};

Complex rand_matrix_1::rand_scalar_compl()
{
    return Complex(randn(),randn());
};

Float_complex rand_matrix_1::rand_scalar_fcompl()
{
    return Float_complex(frandn(), frandn());
};

matcl::Matrix rand_matrix_1::rand_dense_int(Integer m,Integer n)
{
    return iround(max_int * randn(m,n));
};

matcl::Matrix rand_matrix_1::rand_dense_real(Integer m,Integer n)
{
    return randn(m,n);
};
matcl::Matrix rand_matrix_1::rand_dense_float(Integer m,Integer n)
{
    return frandn(m,n);
};
matcl::Matrix rand_matrix_1::rand_dense_compl(Integer m,Integer n)
{
    return crandn(m,n);
};
matcl::Matrix rand_matrix_1::rand_dense_fcompl(Integer m,Integer n)
{
    return fcrandn(m,n);
};

matcl::Matrix rand_matrix_1::rand_sparse_int(Integer m,Integer n, Real d)
{
    Matrix tmp = sprandn(m,n,d);
    return iround(max_int*tmp);
};
matcl::Matrix rand_matrix_1::rand_sparse_real(Integer m,Integer n, Real d)
{
    return sprandn(m,n,d);
};
matcl::Matrix rand_matrix_1::rand_sparse_float(Integer m,Integer n, Real d)
{
    return fsprandn(m,n,d);
};
matcl::Matrix rand_matrix_1::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return csprandn(m,n,d);
};
matcl::Matrix rand_matrix_1::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return fcsprandn(m,n,d);
};

matcl::Matrix rand_matrix_1::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return iround(max_int*randn_band(m,n,fd,ld));
};
matcl::Matrix rand_matrix_1::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return randn_band(m,n,fd,ld);
};
matcl::Matrix rand_matrix_1::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return frandn_band(m,n,fd,ld);
};

matcl::Matrix rand_matrix_1::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return crandn_band(m,n,fd,ld);
};
matcl::Matrix rand_matrix_1::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return fcrandn_band(m,n,fd,ld);
};

//======================================================================
//                      rand_matrix_sparse
//======================================================================
rand_matrix_sparse::rand_matrix_sparse(const rand_matrix_ptr& base_gen)
    :m_base_generator(base_gen)
{};

Integer rand_matrix_sparse::rand_scalar_int()
{
    return m_base_generator->rand_scalar_int();
};
Real rand_matrix_sparse::rand_scalar_real()
{
    return m_base_generator->rand_scalar_real();
};
Float rand_matrix_sparse::rand_scalar_float()
{
    return m_base_generator->rand_scalar_float();
};
Complex rand_matrix_sparse::rand_scalar_compl()
{
    return m_base_generator->rand_scalar_compl();
};
Float_complex rand_matrix_sparse::rand_scalar_fcompl()
{
    return m_base_generator->rand_scalar_fcompl();
};
matcl::Matrix rand_matrix_sparse::rand_dense_int(Integer m,Integer n)
{
    return m_base_generator->rand_dense_int(m,n);
};
matcl::Matrix rand_matrix_sparse::rand_dense_real(Integer m,Integer n)
{
    return m_base_generator->rand_dense_real(m,n);
};
matcl::Matrix rand_matrix_sparse::rand_dense_float(Integer m,Integer n)
{
    return m_base_generator->rand_dense_float(m,n);
};
matcl::Matrix rand_matrix_sparse::rand_dense_compl(Integer m,Integer n)
{
    return m_base_generator->rand_dense_compl(m,n);
};
matcl::Matrix rand_matrix_sparse::rand_dense_fcompl(Integer m,Integer n)
{
    return m_base_generator->rand_dense_fcompl(m,n);
};

matcl::Matrix rand_matrix_sparse::rand_sparse_int(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_int(m,n,d);
};
matcl::Matrix rand_matrix_sparse::rand_sparse_real(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_real(m,n,d);
};
matcl::Matrix rand_matrix_sparse::rand_sparse_float(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_float(m,n,d);
};
matcl::Matrix rand_matrix_sparse::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_compl(m,n,d);
};
matcl::Matrix rand_matrix_sparse::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_fcompl(m,n,d);
};

matcl::Matrix rand_matrix_sparse::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_int(m,n,fd,ld);
};
matcl::Matrix rand_matrix_sparse::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_real(m,n,fd,ld);
};
matcl::Matrix rand_matrix_sparse::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_float(m,n,fd,ld);
};

matcl::Matrix rand_matrix_sparse::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_compl(m,n,fd,ld);
};
matcl::Matrix rand_matrix_sparse::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_fcompl(m,n,fd,ld);
};

matcl::Matrix rand_matrix_sparse::rand_any_matrix()
{
    Real Max    = 200;
    Real Max_d  = 5./Max;
    Integer max_val = 1000;

    Integer n   = static_cast<Integer>(Max*rand());
    Integer m   = static_cast<Integer>(Max*rand());
    Real d      = Max_d * rand();
    Integer t   = static_cast<Integer>(5. * rand());

    switch(t)
    {
        case 0: return matcl::mod(isprand(m,n,d), max_val);
        case 1: return sprandn(m,n,d);
        case 2: return csprandn(m,n,d);
        case 3: return fsprandn(m,n,d);
        case 4: return fcsprandn(m,n,d);
        default: return sprandn(m,n,d);
    };    
};

//======================================================================
//                      rand_matrix_dense
//======================================================================
rand_matrix_dense::rand_matrix_dense(const rand_matrix_ptr& base_gen)
    :m_base_generator(base_gen)
{};

Integer rand_matrix_dense::rand_scalar_int()
{
    return m_base_generator->rand_scalar_int();
};
Real rand_matrix_dense::rand_scalar_real()
{
    return m_base_generator->rand_scalar_real();
};
Float rand_matrix_dense::rand_scalar_float()
{
    return m_base_generator->rand_scalar_float();
};
Complex rand_matrix_dense::rand_scalar_compl()
{
    return m_base_generator->rand_scalar_compl();
};
Float_complex rand_matrix_dense::rand_scalar_fcompl()
{
    return m_base_generator->rand_scalar_fcompl();
};

matcl::Matrix rand_matrix_dense::rand_dense_int(Integer m,Integer n)
{
    return m_base_generator->rand_dense_int(m,n);
};
matcl::Matrix rand_matrix_dense::rand_dense_real(Integer m,Integer n)
{
    return m_base_generator->rand_dense_real(m,n);
};
matcl::Matrix rand_matrix_dense::rand_dense_float(Integer m,Integer n)
{
    return m_base_generator->rand_dense_float(m,n);
};
matcl::Matrix rand_matrix_dense::rand_dense_compl(Integer m,Integer n)
{
    return m_base_generator->rand_dense_compl(m,n);
};
matcl::Matrix rand_matrix_dense::rand_dense_fcompl(Integer m,Integer n)
{
    return m_base_generator->rand_dense_fcompl(m,n);
};
matcl::Matrix rand_matrix_dense::rand_sparse_int(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_int(m,n,d);
};
matcl::Matrix rand_matrix_dense::rand_sparse_real(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_real(m,n,d);
};
matcl::Matrix rand_matrix_dense::rand_sparse_float(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_float(m,n,d);
};
matcl::Matrix rand_matrix_dense::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_compl(m,n,d);
};
matcl::Matrix rand_matrix_dense::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return m_base_generator->rand_sparse_fcompl(m,n,d);
};
matcl::Matrix rand_matrix_dense::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_int(m,n,fd,ld);
};
matcl::Matrix rand_matrix_dense::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_real(m,n,fd,ld);
};
matcl::Matrix rand_matrix_dense::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_float(m,n,fd,ld);
};
matcl::Matrix rand_matrix_dense::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_compl(m,n,fd,ld);
};
matcl::Matrix rand_matrix_dense::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return m_base_generator->rand_band_fcompl(m,n,fd,ld);
};

matcl::Matrix rand_matrix_dense::rand_any_matrix()
{
    Real Max    = 200;
    Integer n   = static_cast<Integer>(Max*rand());
    Integer m   = static_cast<Integer>(Max*rand());
    Integer t   = static_cast<Integer>(5.*rand());

    switch(t)
    {
        case 0: return matcl::mod(irand(m,n), 100);
        case 1: return randn(m,n);
        case 2: return crandn(m,n);
        case 3: return frandn(m,n);
        case 4: return fcrandn(m,n);
        default: return randn(m,n);
    };
    
};

//======================================================================
//                      rand_matrix_obj
//======================================================================
Integer rand_matrix_obj::rand_scalar_int()
{
    return irand()%max_int;
};

Real rand_matrix_obj::rand_scalar_real()
{
    return randn();
};
Float rand_matrix_obj::rand_scalar_float()
{
    return frandn();
};

Complex rand_matrix_obj::rand_scalar_compl()
{
    return Complex(randn(),randn());
};
Float_complex rand_matrix_obj::rand_scalar_fcompl()
{
    return Float_complex(frandn(),frandn());
};

matcl::Matrix rand_matrix_obj::rand_dense_int(Integer m,Integer n)
{
    return matcl::convert_object(iround(max_int*randn(m,n)), ti::predefined::get_ti_int());
};

matcl::Matrix rand_matrix_obj::rand_dense_real(Integer m,Integer n)
{
    return matcl::convert_object(randn(m,n), ti::predefined::get_ti_real());
};
matcl::Matrix rand_matrix_obj::rand_dense_float(Integer m,Integer n)
{
    return matcl::convert_object(frandn(m,n), ti::predefined::get_ti_float());
};
matcl::Matrix rand_matrix_obj::rand_dense_compl(Integer m,Integer n)
{
    return matcl::convert_object(crandn(m,n), ti::predefined::get_ti_complex());
};
matcl::Matrix rand_matrix_obj::rand_dense_fcompl(Integer m,Integer n)
{
    return matcl::convert_object(fcrandn(m,n), ti::predefined::get_ti_float_complex());
};

matcl::Matrix rand_matrix_obj::rand_sparse_int(Integer m,Integer n, Real d)
{
    Matrix tmp = sprandn(m,n,d);
    return matcl::convert_object(iround(max_int*tmp), ti::predefined::get_ti_int());
};
matcl::Matrix rand_matrix_obj::rand_sparse_real(Integer m,Integer n, Real d)
{
    return matcl::convert_object(sprandn(m,n,d), ti::predefined::get_ti_real());
};
matcl::Matrix rand_matrix_obj::rand_sparse_float(Integer m,Integer n, Real d)
{
    return matcl::convert_object(fsprandn(m,n,d), ti::predefined::get_ti_float());
};
matcl::Matrix rand_matrix_obj::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return matcl::convert_object(csprandn(m,n,d), ti::predefined::get_ti_complex());
};
matcl::Matrix rand_matrix_obj::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return matcl::convert_object(fcsprandn(m,n,d), ti::predefined::get_ti_float_complex());
};

matcl::Matrix rand_matrix_obj::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return  matcl::convert_object(iround(max_int*randn_band(m,n,fd,ld)), ti::predefined::get_ti_int());
};
matcl::Matrix rand_matrix_obj::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return matcl::convert_object(randn_band(m,n,fd,ld), ti::predefined::get_ti_real());
};
matcl::Matrix rand_matrix_obj::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return matcl::convert_object(frandn_band(m,n,fd,ld), ti::predefined::get_ti_float());
};

matcl::Matrix rand_matrix_obj::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return matcl::convert_object(crandn_band(m,n,fd,ld), ti::predefined::get_ti_complex());
};
matcl::Matrix rand_matrix_obj::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return matcl::convert_object(fcrandn_band(m,n,fd,ld), ti::predefined::get_ti_float_complex());
};

//======================================================================
//                      rand_matrix_1_sp
//======================================================================
rand_matrix_1_sp::rand_matrix_1_sp(bool with_nan)
:with_nan(with_nan)
{};

Integer rand_matrix_1_sp::rand_scalar_int()
{
    return irandn_sp(with_nan);
};

Real rand_matrix_1_sp::rand_scalar_real()
{
    return randn_sp(with_nan);
};
Float rand_matrix_1_sp::rand_scalar_float()
{
    return frandn_sp(with_nan);
};

Complex rand_matrix_1_sp::rand_scalar_compl()
{
    return Complex(randn_sp(with_nan),randn_sp(with_nan));
};
Float_complex rand_matrix_1_sp::rand_scalar_fcompl()
{
    return Float_complex(frandn_sp(with_nan),frandn_sp(with_nan));
};

matcl::Matrix rand_matrix_1_sp::rand_dense_int(Integer m,Integer n)
{
    return iround(max_int*randn_sp(m,n,false,false));
};

matcl::Matrix rand_matrix_1_sp::rand_dense_real(Integer m,Integer n)
{
    return randn_sp(m,n,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_dense_float(Integer m,Integer n)
{
    return frandn_sp(m,n,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_dense_compl(Integer m,Integer n)
{
    return crandn_sp(m,n,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_dense_fcompl(Integer m,Integer n)
{
    return fcrandn_sp(m,n,with_nan);
};

matcl::Matrix rand_matrix_1_sp::rand_sparse_int(Integer m,Integer n, Real d)
{
    return iround(max_int*sprandn_sp(m,n,d,false,false));
};
matcl::Matrix rand_matrix_1_sp::rand_sparse_real(Integer m,Integer n, Real d)
{
    return sprandn_sp(m,n,d,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_sparse_float(Integer m,Integer n, Real d)
{
    return fsprandn_sp(m,n,d,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return csprandn_sp(m,n,d,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return fcsprandn_sp(m,n,d,with_nan);
};

matcl::Matrix rand_matrix_1_sp::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return iround(max_int*randn_band_sp(m,n,fd,ld,false,false));
};
matcl::Matrix rand_matrix_1_sp::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return randn_band_sp(m,n,fd,ld,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return frandn_band_sp(m,n,fd,ld,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return crandn_band_sp(m,n,fd,ld,with_nan);
};
matcl::Matrix rand_matrix_1_sp::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return fcrandn_band_sp(m,n,fd,ld,with_nan);
};

//======================================================================
//                      rand_matrix_1_str
//======================================================================
rand_matrix_1_str::rand_matrix_1_str()
{};

Integer rand_matrix_1_str::rand_scalar_int()
{
    return irand() % max_int;
};

Real rand_matrix_1_str::rand_scalar_real()
{
    return randn();
};
Float rand_matrix_1_str::rand_scalar_float()
{
    return frandn();
};

Complex rand_matrix_1_str::rand_scalar_compl()
{
    return Complex(randn(),randn());
};
Float_complex rand_matrix_1_str::rand_scalar_fcompl()
{
    return Float_complex(frandn(),frandn());
};

matcl::Matrix rand_matrix_1_str::rand_dense_int(Integer m,Integer n)
{
    return iround(max_int*randn_str(m,n));
};

matcl::Matrix rand_matrix_1_str::rand_dense_real(Integer m,Integer n)
{
    return randn_str(m,n);
};
matcl::Matrix rand_matrix_1_str::rand_dense_float(Integer m,Integer n)
{
    return frandn_str(m,n);
};
matcl::Matrix rand_matrix_1_str::rand_dense_compl(Integer m,Integer n)
{
    return crandn_str(m,n);
};
matcl::Matrix rand_matrix_1_str::rand_dense_fcompl(Integer m,Integer n)
{
    return fcrandn_str(m,n);
};

matcl::Matrix rand_matrix_1_str::rand_sparse_int(Integer m,Integer n, Real d)
{
    return iround(max_int*sprandn_str(m,n,d));
};
matcl::Matrix rand_matrix_1_str::rand_sparse_real(Integer m,Integer n, Real d)
{
    return sprandn_str(m,n,d);
};
matcl::Matrix rand_matrix_1_str::rand_sparse_float(Integer m,Integer n, Real d)
{
    return fsprandn_str(m,n,d);
};

matcl::Matrix rand_matrix_1_str::rand_sparse_compl(Integer m,Integer n, Real d)
{
    return csprandn_str(m,n,d);
};
matcl::Matrix rand_matrix_1_str::rand_sparse_fcompl(Integer m,Integer n, Real d)
{
    return fcsprandn_str(m,n,d);
};

matcl::Matrix rand_matrix_1_str::rand_band_int(Integer m,Integer n, Integer fd, Integer ld)
{
    return iround(max_int*randn_band_str(m,n,fd,ld));
};
matcl::Matrix rand_matrix_1_str::rand_band_real(Integer m,Integer n, Integer fd, Integer ld)
{
    return randn_band_str(m,n,fd,ld);
};
matcl::Matrix rand_matrix_1_str::rand_band_float(Integer m,Integer n, Integer fd, Integer ld)
{
    return frandn_band_str(m,n,fd,ld);
};
matcl::Matrix rand_matrix_1_str::rand_band_compl(Integer m,Integer n, Integer fd, Integer ld)
{
    return crandn_band_str(m,n,fd,ld);
};
matcl::Matrix rand_matrix_1_str::rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld)
{
    return fcrandn_band_str(m,n,fd,ld);
};

};};

