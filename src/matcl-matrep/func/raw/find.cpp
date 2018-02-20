/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/func/raw/find.h"
#include "matcl-matrep/base/instantiate.h"
#include <vector>
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw { namespace details 
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;
     
template<class in, class functor, class str>
struct find_impl
{};

template<class in, class functor>
struct find_impl<in,functor,struct_dense>
{
    static void eval(const in& A, functor& func)
    {
        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0)
            return;

        const typename in::value_type* ptr_A = A.ptr();
        Integer A_ld            = A.ld();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
                func.add(i+1,j+1,ptr_A[i]);

            ptr_A += A_ld;
        };

        return;
    };
};

template<class in, class functor>
struct find_impl<in,functor,struct_banded>
{
    static void eval(const in& A, functor& func)
    {
        Integer r   = A.rows();
        Integer c   = A.cols();

        if (r == 0 || c == 0)
            return;

        using VT    = typename in::value_type;

        Integer A_ld    = A.ld();
        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();

        if (fd == ld)
        {
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);
            Integer s       = A.diag_length(fd);
            Integer r0      = A.first_row_on_diag(fd);
            Integer c0      = A.first_col_on_diag(fd);

            for (Integer j = 0; j < c0; ++j)
                func.add_zero(1,r,j+1);

            for (Integer j = c0; j < c0 + s; ++j, ++r0)
            {
                if (r0 > 0)
                    func.add_zero(1, r0, j+1);

                func.add(r0 + 1,j + 1, ptr_A[0]);

                if (r0 < r-1)
                    func.add_zero(r0 + 2, r, j+1);

                ptr_A += A_ld;
            };

            for (Integer j = c0 + s; j < c; ++j)
                func.add_zero(1,r,j+1);
        }
        else
        {
            const VT* ptr_A     = A.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);

                if (first_row > 0)
                    func.add_zero(1,std::min(first_row,r),j+1);

                Integer pos = A.first_elem_pos(j);

                for (Integer i = first_row; i <= last_row; ++i, ++pos)
                    func.add(i+1,j+1,ptr_A[pos]);

                if (last_row < r-1)
                    func.add_zero(last_row + 2,r,j+1);

                ptr_A += A_ld;
            };
        };

        return;
    };
};

template<class in, class functor>
struct find_impl<in,functor,struct_sparse>
{
    static void eval(const in& A, functor& func)
    {
        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0)
            return;

        if (A.nnz() == 0)
        {
            for (Integer j = 0; j < c; ++j)
                func.add_zero(1,r,j+1);

            return;
        };

        using val_type = typename in::value_type;

        const sparse_ccs<val_type>& Ad		= A.rep();
        const Integer* Ad_c		= Ad.ptr_c();
        const Integer* Ad_r		= Ad.ptr_r();
        const val_type* Ad_x	= Ad.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            if (Ad_c[j] == Ad_c[j+1]) // empty column
            {
                func.add_zero(1, r, j + 1);
            }
            else
            {
                Integer last_checked_row = -1;
                for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
                {
                    Integer k = Ad_r[i];
                    if (k > last_checked_row + 1) // mark all elements above this non-zero as zero, if any
                    {
                        // cr + 1 because we want next element, and another + 1 for "1" based indexing 
                        // ("k" and thus "cr" are 0-based indexed)
                        func.add_zero(last_checked_row+2,k,j+1);
                    }
                    func.add(k+1,j+1,Ad_x[i]);
                    last_checked_row = k;
                }
                if (last_checked_row < r - 1) // mark all zero elements remaining in column
                {
                    // cr + 1 because we want next element, and another + 1 for "1" based indexing 
                    // ("k" and thus "cr" are 0-based indexed)
                    func.add_zero(last_checked_row + 2, r, j + 1);
                }
            }
        };

        return;
    };
};

template<class T, int type>
class find_functor_base
{};

template<class T>
class find_functor_base<T,1>
{
    protected:
        using tinfo = ti::ti_type<T>;

        Integer r;
        Integer c;
        tinfo   m_ti;

        std::vector<Integer> indices;

    public:
        find_functor_base(tinfo ti, Integer r, Integer c)
            :m_ti(ti), r(r), c(c) 
        {};

        void add_index(Integer i, Integer j, const T& )
        {
            error::check_linear_index(i,j,r);

            Integer pos = imult(j-1,r) + i;
            indices.push_back(pos);
        };

        integer_dense get_return() const
        {
            Integer rc = cast_int64(static_cast<Integer_64>(indices.size()));
            integer_dense res(ti::ti_int(),rc,1);
            Integer * ptr_res = res.ptr();
            
            for (Integer i = 0; i < rc; ++i)
            {
                ptr_res[i] = indices[i];
            };

            return res;
        };
};

template<class T>
class find_functor_base<T,2>
{
    protected:
        using tinfo = ti::ti_type<T>;

        Integer r;
        Integer c;
        tinfo   m_ti;

        std::vector<Integer> indices[2];

    public:
        find_functor_base(tinfo ti, Integer r, Integer c)
            :r(r), c(c), m_ti(ti)
        {};

        void add_index(Integer i, Integer j, const T& )
        {
            indices[0].push_back(i);
            indices[1].push_back(j);
        };

        integer_dense get_return(int type) const
        {
            Integer rc = cast_int64(static_cast<Integer_64>(indices[0].size()));
            integer_dense res(ti::ti_int(),rc,1);
            Integer * ptr_res = res.ptr();
            
            for (Integer i = 0; i < rc; ++i)
            {
                ptr_res[i] = indices[type][i];
            };

            return res;
        };
};

template<class T>
class find_functor_base<T,3>
{
    protected:
        using tinfo = ti::ti_type<T>;

        Integer r;
        Integer c;
        tinfo   m_ti;

        std::vector<Integer>            indices[2];
        std::vector<T>		            indices_x;

        using ValueMatrix = Matrix<T,struct_dense>;        

    public:
        find_functor_base(tinfo ti, Integer r, Integer c)
            :r(r), c(c), m_ti(ti)
        {};
        
        void add_index(Integer i, Integer j, const T& val)
        {
            indices[0].push_back(i);
            indices[1].push_back(j);
            indices_x.push_back(val);
        };

        integer_dense get_return(int type) const
        {
            Integer rc = cast_int64(static_cast<Integer_64>(indices[0].size()));
            integer_dense res(ti::ti_int(),rc,1);
            Integer * ptr_res = res.ptr();
            
            for (Integer i = 0; i < rc; ++i)
            {
                ptr_res[i] = indices[type][i];
            };

            return res;
        };
        
        ValueMatrix	get_return_x() const
        {
            Integer rc = cast_int64(static_cast<Integer_64>(indices[0].size()));
            ValueMatrix res(m_ti,rc,1);
            T * ptr_res = res.ptr();
            
            for (Integer i = 0; i < rc; ++i)
            {
                mrd::reset_helper(ptr_res[i],indices_x[i]);
            };

            return res;
        };
};

template<class T,int type>
class find_functor : public find_functor_base<T,type>
{
    private:
        using tinfo         = typename find_functor_base<T,type>::tinfo;
        using base_class    = find_functor_base<T,type>;

    public:
        find_functor(tinfo ti, Integer r, Integer c)
            : find_functor_base<T, type>(ti, r,c) 
        {};
        
        void add(Integer i, Integer j, const T& val)
        {
            if (mrd::is_zero(val))
            {
                return;
            };
            base_class::add_index(i,j,val);
        };
        
        void add_zero(Integer , Integer , Integer )
        {};
};

template<class T,int type>
class find_functor_t: public find_functor_base<T,type>
{
    private:
        using base_class    = find_functor_base<T,type>;
        using tinfo         = typename find_functor_base<T,type>::tinfo;

        const test_function& test;

    public:
        find_functor_t(tinfo ti, const test_function& test, Integer r, Integer c)
            :find_functor_base<T, type>(ti, r,c), test(test) 
        {};
        
        void add(Integer i, Integer j, const T& val)
        {
            if (test.eval(val) == false)
            {
                return;
            };
            base_class::add_index(i,j,val);
        };
        
        void add_zero(Integer , Integer , Integer )
        {};

    private:
        find_functor_t(const find_functor_t&);
        find_functor_t& operator=(const find_functor_t&);
};

template<class T,int type>
class find_functor_t0: public find_functor_base<T,type>
{
    private:
        using base_class    = find_functor_base<T,type>;
        using base_type     = find_functor_base<T,type>;

    private:
        const test_function& test;

    public:
        find_functor_t0(typename base_type::tinfo ti, const test_function& test, Integer r, Integer c)
            :find_functor_base<T, type>(ti,r,c), test(test)
        {};
        
        void add(Integer i, Integer j, const T& val)
        {
            if (test.eval(val) == false)
                return;

            base_class::add_index(i,j,val);
        };
        
        void add_zero(Integer r_start, Integer r_end, Integer col)
        {
            T zero = matcl::details::default_value<T>(base_type::m_ti);

            for (Integer i = r_start; i <= r_end; ++i)
                base_class::add_index(i,col,zero);
        };

    private:
        find_functor_t0(const find_functor_t0&);
        find_functor_t0& operator=(const find_functor_t0&);
};

template<class MP>
void find_helper<MP>::eval_find(matcl::Matrix& i, const MP& m)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    find_functor<val_type,1> ff(m.get_type(), m.rows(),m.cols());
    find_impl<MP,find_functor<val_type,1>,str_type>::eval(m,ff);
    
    i = matcl::Matrix(ff.get_return(), false);
};

template<class MP>
void find_helper<MP>::eval_find(matcl::Matrix& i, const MP& m, const test_function& t)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    val_type Z = matcl::details::default_value<val_type>(m.get_type());
    
    if (t.eval(Z))
    {
        find_functor_t0<val_type,1> ff(m.get_type(),t,m.rows(),m.cols());
        find_impl<MP,find_functor_t0<val_type,1>,str_type>::eval(m,ff);
        i = matcl::Matrix(ff.get_return(), false);
        return;
    }
    else
    {
        find_functor_t<val_type,1> ff(m.get_type(), t,m.rows(),m.cols());
        find_impl<MP,find_functor_t<val_type,1>,str_type>::eval(m,ff);
        i = matcl::Matrix(ff.get_return(),false);
        return;
    };
};

template<class MP>
void find_helper<MP>::eval_find_2(matcl::Matrix& i, matcl::Matrix& j, const MP& m)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    find_functor<val_type,2> ff(m.get_type(), m.rows(),m.cols());
    find_impl<MP,find_functor<val_type,2>,str_type>::eval(m,ff);

    i = matcl::Matrix(ff.get_return(0),false);
    j = matcl::Matrix(ff.get_return(1),false);
};

template<class MP>
void find_helper<MP>::eval_find_2(matcl::Matrix& i, matcl::Matrix& j, const MP& m, 
                                  const test_function& t)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    val_type Z = matcl::details::default_value<val_type>(m.get_type());

    if (t.eval(Z))
    {
        find_functor_t0<val_type,2> ff(m.get_type(),t,m.rows(),m.cols());
        find_impl<MP,find_functor_t0<val_type,2>,str_type>::eval(m,ff);

        i = matcl::Matrix(ff.get_return(0),false);
        j = matcl::Matrix(ff.get_return(1),false);
        return;
    }
    else
    {
        find_functor_t<val_type,2> ff(m.get_type(), t,m.rows(),m.cols());
        find_impl<MP,find_functor_t<val_type,2>,str_type>::eval(m,ff);

        i = matcl::Matrix(ff.get_return(0),false);
        j = matcl::Matrix(ff.get_return(1),false);
        return;
    };
};

template<class MP>
void find_helper<MP>::eval_find_3(matcl::Matrix& i, matcl::Matrix& j, matcl::Matrix& x, const MP& m)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    find_functor<val_type,3> ff(m.get_type(), m.rows(),m.cols());
    find_impl<MP,find_functor<val_type,3>,str_type>::eval(m,ff);

    i = matcl::Matrix(ff.get_return(0),false);
    j = matcl::Matrix(ff.get_return(1), false);
    x = matcl::Matrix(ff.get_return_x(),false);
};

template<class MP>
void find_helper<MP>::eval_find_3(matcl::Matrix& i, matcl::Matrix& j, matcl::Matrix& x, 
                                  const MP& m, const test_function& t)
{
    using str_type  = typename MP::struct_type;
    using val_type  = typename MP::value_type;

    val_type Z = matcl::details::default_value<val_type>(m.get_type());

    if (t.eval(Z))
    {
        find_functor_t0<val_type,3> ff(m.get_type(),t,m.rows(),m.cols());
        find_impl<MP,find_functor_t0<val_type,3>,str_type>::eval(m,ff);

        i = matcl::Matrix(ff.get_return(0),false);
        j = matcl::Matrix(ff.get_return(1),false);
        x = matcl::Matrix(ff.get_return_x(),false);
        return;
    }
    else
    {
        find_functor_t<val_type,3> ff(m.get_type(), t,m.rows(),m.cols());
        find_impl<MP,find_functor_t<val_type,3>,str_type>::eval(m,ff);

        i = matcl::Matrix(ff.get_return(0),false);
        j = matcl::Matrix(ff.get_return(1),false);
        x = matcl::Matrix(ff.get_return_x(),false);
        return;
    };
};

};};};

MACRO_INSTANTIATE_G_1(matcl::raw::details::find_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::find_helper)
