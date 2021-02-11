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

#include "matcl-matrep/algs/banded_algs.h"

namespace matcl { namespace details
{

template<class M1,class struct_type>
struct del_rows_impl
{};

template<class M1>
struct del_rows_impl<M1,struct_dense>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),1, c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        if (c_info.r_step == 1)
        {
            if (c_info.r_start == 1)
            {
                ret = Matrix(mat.make_view(c_info.r_end+1, mat.rows(), 1, mat.cols()),false);
                return;
            }
            else if (c_info.r_end == mat.rows())
            {
                ret = Matrix(mat.make_view(1,c_info.r_start-1, 1, mat.cols()),false);
                return;
            };
        };

        return algorithm::del_rows_dense(ret, mat,c_info, rvalue);
    };
};

template<class M1>
struct del_rows_impl<M1,struct_sparse>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),1, c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        return algorithm::del_rows_sparse(ret, mat,c_info, rvalue);
    };
};

template<class M1>
struct del_rows_impl<M1,struct_banded>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),1, c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        if (c_info.r_step == 1 && c_info.r_end == mat.rows())
        {
            ret = Matrix(mat.make_view(1, c_info.r_start-1, mat.cols()),false);
            return;
        };

        algorithm::del_rows_banded(ret, mat,c_info,rvalue);
        return;
    };
};

template<class M1,class struct_type>
struct del_cols_impl
{};

template<class M1>
struct del_cols_impl<M1,struct_dense>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(1,mat.cols(), c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        if (c_info.r_step == 1)
        {
            if (c_info.r_start == 1)
            {
                ret = Matrix(mat.make_view(1, mat.rows(), c_info.r_end+1, mat.cols()),false);
                return;
            }
            else if (c_info.r_end == mat.cols())
            {
                ret = Matrix(mat.make_view(1,mat.rows(), 1, c_info.r_start-1),false);
                return;
            };
        };

        return algorithm::del_cols_dense(ret,mat,c_info,rvalue);
    };
};

template<class M1>
struct del_cols_impl<M1,struct_sparse>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(1, mat.cols(), c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        if (c_info.r_step == 1)
        {
            if (c_info.r_start == 1)
            {
                ret = Matrix(mat.make_view(c_info.r_end+1, mat.cols()), false);
                return;
            }
            else if (c_info.r_end == mat.cols())
            {
                ret = Matrix(mat.make_view(1,c_info.r_start-1), false);
                return;
            };
        };

        return algorithm::del_cols_sparse(ret,mat,c_info,rvalue);
    };
};

template<class M1>
struct del_cols_impl<M1,struct_banded>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1, bool rvalue)
    {
        colon_info c_info;
        make_index(1,mat.cols(), c1, c_info);

        if (c_info.rows() == 0)
        {
            ret = Matrix(A);
            return;
        };

        if (c_info.r_step == 1 && c_info.r_end == mat.cols())
        {
             ret = Matrix(mat.make_view(1, mat.rows(), c_info.r_start-1),false);
             return;
        };

        return algorithm::del_cols_banded(ret,mat,c_info,rvalue);
    };
};

template<class M1,class struct_type>
struct del_rowscols_impl
{};

template<class M1>
struct del_rowscols_impl<M1,struct_dense>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1,
                     const colon& c2, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        if (c_info.rows() == 0)
        {
            if (c_info.cols() == 0)
            {
                ret = Matrix(A);
                return;
            }
            else
            {
                return del_cols_impl<M1,struct_dense>::eval(ret,A,mat,c2,rvalue);
            }
        }
        else if (c_info.cols() == 0)
        {
            return del_rows_impl<M1,struct_dense>::eval(ret,A,mat,c1,rvalue);
        };

        if (c_info.r_step == 1 && c_info.c_step == 1)
        {
            if (c_info.r_start == 1)
            {
                if (c_info.c_start == 1)
                {
                    ret = Matrix(mat.make_view(c_info.r_end+1, mat.rows(), c_info.c_end+1, mat.cols()),false);
                    return;
                }
                else if (c_info.c_end == mat.cols())
                {
                    ret = Matrix(mat.make_view(c_info.r_end+1, mat.rows(), 1, c_info.c_start-1),false);
                    return;
                }
            }
            else if (c_info.r_end == mat.rows())
            {
                if (c_info.c_start == 1)
                {
                    ret = Matrix(mat.make_view(1,c_info.r_start-1, c_info.c_end+1, mat.cols()),false);
                    return;
                }
                else if (c_info.c_end == mat.cols())
                {
                    ret = Matrix(mat.make_view(1,c_info.r_start-1, 1, c_info.c_start-1),false);
                    return;
                }
            };
        };

        return algorithm::del_rowscols_dense(ret,mat,c_info,rvalue);
    };
};

template<class M1>
struct del_rowscols_impl<M1,struct_sparse>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A,const M1& mat,const colon& c1, 
                     const colon& c2, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        if (c_info.rows() == 0)
        {
            if (c_info.cols() == 0)
            {
                ret = Matrix(A);
                return;
            }
            else
            {
                return del_cols_impl<M1,struct_sparse>::eval(ret,A,mat,c2,rvalue);
            }
        }
        else if (c_info.cols() == 0)
        {
            return del_rows_impl<M1,struct_sparse>::eval(ret,A,mat,c1,rvalue);
        };

        return algorithm::del_rowscols_sparse(ret,mat,c_info,rvalue);
    };
};

template<class M1>
struct del_rowscols_impl<M1,struct_banded>
{
    using value_type = typename M1::value_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1,
                     const colon& c2, bool rvalue)
    {
        colon_info c_info;
        make_index(mat.rows(),mat.cols(), c1, c2, c_info);

        if (c_info.rows() == 0)
        {
            if (c_info.cols() == 0)
            {
                ret = Matrix(A);
                return;
            }
            else
            {
                return del_cols_impl<M1,struct_banded>::eval(ret,A,mat,c2,rvalue);
            }
        }
        else if (c_info.cols() == 0)
        {
            return del_rows_impl<M1,struct_banded>::eval(ret,A,mat,c1,rvalue);
        };

        if (c_info.r_step == 1 && c_info.r_end == mat.rows())
        {
            if (c_info.c_step == 1 && c_info.c_end == mat.cols())
            {
                 ret = Matrix(mat.make_view(1, c_info.r_start-1, c_info.c_start-1),false);
                 return;
            };
        };

        return algorithm::del_rowscols_banded(ret,mat,c_info,rvalue);
    };
};

template<class M1>
struct del_functor_scal
{
    static void eval(Matrix& ret, const Matrix& A, const M1&,const colon& c, bool del_row)
    {
        colon_info c_info;
        make_index(1,1,c,c_info);

        Integer s = c_info.rows();
        if (s == 0)
        {
            ret = Matrix(A);
            return;
        };

        using FullMatrix = raw::Matrix<M1,struct_dense>;

        ti::ti_type<M1> ti_A = get_matrix_ti<M1>::eval(A);

        if (del_row)
        {
            ret = Matrix(FullMatrix(ti_A,0,1),false);
            return;
        }
        else
        {
            ret = Matrix(FullMatrix(ti_A,1,0),false);
            return;
        }
    }

    static void eval2(Matrix& ret, const Matrix& A, const M1&,const colon& c1,const colon& c2)
    {
        colon_info c_info;
        make_index(1,1,c1,c2,c_info);

        Integer sr = c_info.rows();
        Integer sc = c_info.cols();

        if (sr == 0 && sc == 0)
        {
            ret = Matrix(A);
            return;
        };

        using FullMatrix = raw::Matrix<M1,struct_dense>;

        ti::ti_type<M1> ti_A = get_matrix_ti<M1>::eval(A);

        if (sr != 0)
        {
            if (sc != 0)
                ret = Matrix(FullMatrix(ti_A,0,0),false);
            else 
                ret = Matrix(FullMatrix(ti_A,0,1),false);
            return;
        }
        else
        {
            //sr == 0, sc != 0
            ret = Matrix(FullMatrix(ti_A,1,0),false);
            return;
        };
    }
};

template<class M1>
struct del_rows_functor
{
    using value_type    = typename M1::value_type;
    using struct_type   = typename M1::struct_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1, bool rvalue)
    {
        if (c1.m_flag == colon::t_all)
        {
            using FullMatrix    = raw::Matrix<value_type,struct_dense>;
            ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);

            FullMatrix out(ti_A,0,mat.cols());
            ret = Matrix(out,false);
            return;
        };
        return del_rows_impl<M1,struct_type>::eval(ret,A,mat,c1,rvalue);
    };
};

template<class M1>
struct del_cols_functor
{
    using value_type    = typename M1::value_type;
    using struct_type   = typename M1::struct_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1, bool rvalue)
    {
        if (c1.m_flag == colon::t_all)
        {
            using FullMatrix = raw::Matrix<value_type,struct_dense>;
            ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);
            FullMatrix out(ti_A,mat.rows(),0);
            ret = Matrix(out,false);
            return;
        };
        return del_cols_impl<M1,struct_type>::eval(ret,A,mat,c1,rvalue);
    };
};

template<class M1>
struct del_rowscols_functor
{
    using value_type    = typename M1::value_type;
    using struct_type   = typename M1::struct_type;

    static void eval(Matrix& ret, const Matrix& A, const M1& mat,const colon& c1,
                     const colon& c2, bool rvalue)
    {
        if (c1.m_flag == colon::t_all && c2.m_flag == colon::t_all)
        {
            using FullMatrix = raw::Matrix<value_type,struct_dense>;
            ti::ti_type<value_type> ti_A = get_matrix_ti<value_type>::eval(A);
            FullMatrix out(ti_A,0,0);
            ret = Matrix(out,false);
            return;
        }

        return del_rowscols_impl<M1,struct_type>::eval(ret,A,mat,c1,c2,rvalue);
    };
};

};};