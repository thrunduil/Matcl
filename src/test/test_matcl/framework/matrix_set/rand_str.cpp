/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "rand_str.h"
#include "rand_sp.h"

#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

template<class value_type>
struct rand_struct_impl
{
    static Matrix eval_qtril(const Matrix& A)
    {
        if (A.cols() < 2)
            return A;

        Matrix B = tril(A,1);

        dense_matrix<value_type,false> rep(B);

        Integer s       = std::min(B.rows(), B.cols()-1);
        value_type* ptr = rep.ptr() + rep.ld();

        bool is_z       = true;
        for (Integer i = 0; i < s; ++i)
        {
            if (ptr[0] == 0)
            {
                is_z = true;
            }
            else
            {
                if (is_z == true)
                {
                    is_z = false;
                }
                else
                {
                    ptr[0] = 0;
                    is_z = true;
                };
            };
            ptr += rep.ld() + 1;
        };

        return Matrix(rep);
    };

    static Matrix eval_qtriu(const Matrix& A)
    {
        if (A.rows() < 2)
            return A;

        Matrix B = triu(A,-1);

        dense_matrix<value_type,false> rep(B);

        Integer s       = std::min(B.rows()-1, B.cols());
        value_type* ptr = rep.ptr() + 1;
        bool is_z       = true;

        for (Integer i = 0; i < s; ++i)
        {
            if (ptr[0] == 0)
            {
                is_z = true;
            }
            else
            {
                if (is_z == true)
                {
                    is_z = false;
                }
                else
                {
                    ptr[0] = 0;
                    is_z = true;
                };
            };
            ptr += rep.ld() + 1;
        };

        return Matrix(rep);
    };
};

static Matrix qtril(const Matrix& A)
{    
    switch(A.get_value_code())
    {
        case value_code::v_object:
        {
            return rand_struct_impl<Integer>::eval_qtril(A);
        }
        case value_code::v_float:
        {
            return rand_struct_impl<Float>::eval_qtril(A);
        }
        case value_code::v_real:
        {
            return rand_struct_impl<Real>::eval_qtril(A);
        }
        case value_code::v_float_complex:
        {
            return rand_struct_impl<Float_complex>::eval_qtril(A);
        }
        case value_code::v_complex:
        {
            return rand_struct_impl<Complex>::eval_qtril(A);
        }
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

static Matrix qtriu(const Matrix& A)
{
    switch(A.get_value_code())
    {
        case value_code::v_object:
        {
            return rand_struct_impl<Integer>::eval_qtriu(A);
        }
        case value_code::v_float:
        {
            return rand_struct_impl<Float>::eval_qtriu(A);
        }
        case value_code::v_real:
        {
            return rand_struct_impl<Real>::eval_qtriu(A);
        }
        case value_code::v_float_complex:
        {
            return rand_struct_impl<Float_complex>::eval_qtriu(A);
        }
        case value_code::v_complex:
        {
            return rand_struct_impl<Complex>::eval_qtriu(A);
        }
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

static Matrix unitary(const Matrix& A)
{
    Matrix U;
    if (A.rows() == A.cols())
        U = svd(A,false,svd_algorithm::dc).get<1>();
    else
        U = A;

    return U;
};

Matrix rand_str(const Matrix& A)
{
    Matrix B        = A;
    Integer c       = abs(irand()) % 11;
    value_code vc_A = A.get_value_code();
    bool is_sp      = matrix_traits::is_single_precision(vc_A);

    switch(c)
    {
        case 0:
        {
            //zero
            B = is_sp ? 0.0f * A : 0 * A;
            B.set_struct(struct_flag(predefined_struct_type::diag));
            break;
        }
        case 1:
        {
            //id
            B = speye(A.rows(), A.cols(), vc_A);
            if (A.rows() == A.cols())
                B.set_struct(struct_flag(predefined_struct_type::id));

            break;
        }
        case 2:
        {
            //diag
            B = bdiags(get_diag(A),0,A.rows(),A.cols());
            B.set_struct(struct_flag(predefined_struct_type::diag));
            break;
        }
        case 3:
        {
            //tril
            B = tril(A);
            B.set_struct(struct_flag(predefined_struct_type::tril));
            break;
        }
        case 4:
        {
            //triu
            B = triu(A);
            B.set_struct(struct_flag(predefined_struct_type::triu));
            break;
        }
        case 5:
        {
            //general
            break;
        }
        case 6:
        {
            //sym
            if (A.is_square())
            {
                B = A + trans(A);
                B.set_struct(struct_flag(predefined_struct_type::sym));                
            };
            break;
        }
        case 7:
        {
            //her
            if (A.is_square())
            {
                bool is_compl   = matrix_traits::is_float_complex(B.get_value_code());
                B               = A + ctrans(A);
                if (is_compl)
                    B.set_struct(struct_flag(predefined_struct_type::her));
                else
                    B.set_struct(struct_flag(predefined_struct_type::sym));
            };
            break;
        }
        case 8:
        {
            //qtril
            B = qtril(A);
            B.set_struct(struct_flag(predefined_struct_type::qtril));
            break;
        }
        case 9:
        {
            //qtriu
            B = qtriu(A);
            B.set_struct(struct_flag(predefined_struct_type::qtriu));
            break;
        }
        case 10:
        {
            //unitary
            B = unitary(A);
            if (B.rows() == B.cols())
            {
                struct_flag sf_u;
                sf_u.set_user(unitary_flag());
                B.set_struct(sf_u);
            };
            break;
        }
    };

    check_struct(B);
    return B;
};

Matrix randn_str(Integer m, Integer n)
{
    Matrix out = randn(m,n);
    return full(rand_str(out));
};
Matrix frandn_str(Integer m, Integer n)
{
    Matrix out = frandn(m,n);
    return full(rand_str(out));
};
Matrix crandn_str(Integer m, Integer n)
{
    Matrix out = crandn(m,n);
    return full(rand_str(out));
};
Matrix fcrandn_str(Integer m, Integer n)
{
    Matrix out = fcrandn(m,n);
    return full(rand_str(out));
};
Matrix sprandn_str(Integer m, Integer n, Real d)
{
    Matrix out = sprandn(m,n,d);
    return sparse(rand_str(out));
};
Matrix fsprandn_str(Integer m, Integer n, Real d)
{
    Matrix out = fsprandn(m,n,d);
    return sparse(rand_str(out));
};
Matrix csprandn_str(Integer m, Integer n, Real d)
{
    Matrix out = csprandn(m,n,d);
    return sparse(rand_str(out));
};
Matrix fcsprandn_str(Integer m, Integer n, Real d)
{
    Matrix out = fcsprandn(m,n,d);
    return sparse(rand_str(out));
};
Matrix test::randn_band_str(Integer m, Integer n, Integer fd, Integer ld)
{
    Matrix out = randn_band(m,n,fd,ld);
    return band(rand_str(out));
};
Matrix test::frandn_band_str(Integer m, Integer n, Integer fd, Integer ld)
{
    Matrix out = frandn_band(m,n,fd,ld);
    return band(rand_str(out));
};
Matrix test::crandn_band_str(Integer m, Integer n, Integer fd, Integer ld)
{
    Matrix out = crandn_band(m,n,fd,ld);
    return band(rand_str(out));
};
Matrix test::fcrandn_band_str(Integer m, Integer n, Integer fd, Integer ld)
{
    Matrix out = fcrandn_band(m,n,fd,ld);
    return band(rand_str(out));
};

};};
