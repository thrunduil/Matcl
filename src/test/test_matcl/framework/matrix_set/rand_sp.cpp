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

#include "rand_sp.h"

namespace matcl { namespace test
{

Integer irandn_sp(bool)
{
    Integer val = abs(irand())%7;
    switch (val)
    {
        case 1:
            return 1;
        case 2:
            return -1;
        case 3:
            return 2;
        case 4:
            return -2;
        case 5:
            return -2000;
        case 6:
            return 2000;
        default:
            return 0;
    };
};

Real randn_sp(bool with_nan, bool with_inf)
{
    Integer val = abs(irand())%10;
    switch (val)
    {
        case 1:
            return 1.;
        case 2:
            return -1.;
        case 3:
            return 2.;
        case 4:
            return -2.;
        case 5:
            if (with_inf)
                return constants::inf();
            else
                return 0;
        case 6:
            if (with_inf)
                return -constants::inf();
            else
                return 0;
        case 7:
        {
            if (with_nan)
            {
                return constants::nan();
            }
            else
            {
                return 0;
            };
        }
        case 8:
            if (with_inf)
                return 1.11e8;
            else
                return 1.11e2;
        case 9:
            if (with_inf)
                return -1.11e8;
            else
                return -1.11e2;
        default:
            return 0;
    };
};
Float frandn_sp(bool with_nan, bool with_inf)
{
    Integer val = abs(irand())%10;
    switch (val)
    {
        case 1:
            return 1.f;
        case 2:
            return -1.f;
        case 3:
            return 2.f;
        case 4:
            return -2.f;
        case 5:
            if (with_inf)
                return constants::f_inf();
            else
                return 0.0f;
        case 6:
            if (with_inf)
                return -constants::f_inf();
            else
                return 0.0f;
        case 7:
        {
            if (with_nan)
                return constants::f_nan();
            else
                return 0.f;
        }
        case 8:
            if (with_inf)
                return 1.11e8f;
            else
                return 1.11e2f;
        case 9:
            if (with_inf)
                return -1.11e8f;
            else
                return -1.11e2f;
        default:
            return 0.0f;
    };
};

Matrix randn_sp(Integer m, Integer n,bool with_nan, bool with_inf)
{
    Real* ptr;
    Matrix out = make_real_dense_noinit(m,n,ptr);

    for (Integer i = 0; i <m*n; i++)
        ptr[i] = randn_sp(with_nan,with_inf);

    return out;
};

Matrix frandn_sp(Integer m, Integer n,bool with_nan, bool with_inf)
{
    Float* ptr;
    Matrix out = make_float_dense_noinit(m,n,ptr);

    for (Integer i = 0; i <m*n; i++)
        ptr[i] = frandn_sp(with_nan,with_inf);

    return out;
};

Matrix crandn_sp(Integer m, Integer n,bool with_nan)
{
    Complex* ptr;
    Matrix out = make_complex_dense_noinit(m,n,ptr);
    for (Integer i = 0; i <m*n; i++)
    {
        ptr[i] = Complex(randn_sp(with_nan),randn_sp(with_nan));
    };
    return out;
};

Matrix fcrandn_sp(Integer m, Integer n,bool with_nan)
{
    Float_complex* ptr;
    Matrix out = make_float_complex_dense_noinit(m,n,ptr);
    for (Integer i = 0; i <m*n; i++)
        ptr[i] = Float_complex(frandn_sp(with_nan), frandn_sp(with_nan));

    return out;
};

Matrix sprandn_sp(Integer m, Integer n, Real d,bool with_nan, bool with_inf)
{
    Matrix out = sparse(sprandn(m,n,d));
    sparse_matrix<Real,false> rep(out);
    Real* ptr = rep.ptr_x();

    Integer nnz = out.structural_nnz();
    for (Integer i = 0; i < nnz; i++)
    {
        ptr[i] = randn_sp(with_nan,with_inf);
    };
    return out;
};
Matrix fsprandn_sp(Integer m, Integer n, Real d,bool with_nan, bool with_inf)
{
    Matrix out = sparse(fsprandn(m,n,d));
    sparse_matrix<Float,false> rep(out);
    Float* ptr = rep.ptr_x();

    Integer nnz = out.structural_nnz();
    for (Integer i = 0; i < nnz; i++)
        ptr[i] = frandn_sp(with_nan,with_inf);

    return out;
};

Matrix csprandn_sp(Integer m, Integer n, Real d,bool with_nan)
{
    Matrix out = sparse(csprandn(m,n,d));
    sparse_matrix<Complex,false> rep(out);
    Complex* ptr    = rep.ptr_x();
    Integer nnz     = out.structural_nnz();

    for (Integer i = 0; i < nnz; i++)
        ptr[i] = Complex(randn_sp(with_nan),randn_sp(with_nan));

    return out;
};
Matrix fcsprandn_sp(Integer m, Integer n, Real d,bool with_nan)
{
    Matrix out = sparse(fcsprandn(m,n,d));
    sparse_matrix<Float_complex,false> rep(out);
    Float_complex* ptr  = rep.ptr_x();
    Integer nnz         = out.structural_nnz();

    for (Integer i = 0; i < nnz; i++)
        ptr[i] = Float_complex(frandn_sp(with_nan), frandn_sp(with_nan));

    return out;
};

Matrix test::randn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan, bool with_inf)
{
    Matrix out = make_real_band_noinit(m,n,fd,ld);

    band_matrix<Real,false> rep(out);

    Real* ptr = rep.ptr();
    Integer size = rep.base_size();
    for (Integer i = 0; i < size; i++)
        ptr[i] = randn_sp(with_nan,with_inf);

    return out;
};
Matrix test::frandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan, bool with_inf)
{
    Matrix out = make_float_band_noinit(m,n,fd,ld);

    band_matrix<Float,false> rep(out);

    Float* ptr = rep.ptr();
    Integer size = rep.base_size();
    for (Integer i = 0; i < size; i++)
        ptr[i] = frandn_sp(with_nan,with_inf);

    return out;
};

Matrix test::crandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan)
{
    Matrix out = make_complex_band_noinit(m,n,fd,ld);

    band_matrix<Complex,false> rep(out);
    Complex* ptr = rep.ptr();
    Integer size = rep.base_size();
    for (Integer i = 0; i < size; i++)
    {
        ptr[i] = Complex(randn_sp(with_nan),randn_sp(with_nan));
    };
    return out;
};
Matrix test::fcrandn_band_sp(Integer m, Integer n, Integer fd, Integer ld, bool with_nan)
{
    Matrix out = make_float_complex_band_noinit(m,n,fd,ld);

    band_matrix<Float_complex,false> rep(out);
    Float_complex* ptr = rep.ptr();
    Integer size = rep.base_size();
    for (Integer i = 0; i < size; i++)
        ptr[i] = Float_complex(frandn_sp(with_nan), frandn_sp(with_nan));

    return out;
};

};};