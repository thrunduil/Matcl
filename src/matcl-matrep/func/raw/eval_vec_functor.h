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

#include "matcl-internals/base/utils.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/func/converter.h"

#include <vector>

#pragma warning( push )
#pragma warning(disable:4723)    // potential divide by 0


namespace matcl { namespace raw
{

namespace md = matcl::details;

// accumulator
// value must not depend on number of zeros
template<class in_type, class out_type>
class accumulator    
{
    protected:
        Integer                 size;
        out_type                state;
        md::workspace2<out_type>state_array;    //vector of states, need to be set only if
                                                //reset_array is called

    public:
        accumulator(){};

        void set_size(Integer size) {};            //set number of elements along given dimension
        void reset() {};                            //reset state
        void reset_array(Integer s) {};             //reset state array of size s
        bool add(const in_type& value){};           //returns true is value is known
        void add(Integer p,const in_type& value){}; //add new value at pos p (zero based)
        bool add_zero(){};                          //at least one zero exists returns true is value is known
        void add_zero(Integer p){};                 //at least one zero exists at pos p, returns true is value is known

        out_type value() const {};                  //    
        out_type value(Integer p) const {};         //return value at pos p    
};

template<   class ret_type, 
            class in_type, 
            class accumulator,
            class ret_str,
            class in_str
        >
struct eval_vec_functor_impl{};

template<   class ret_type, 
            class in_type, 
            class accumulator,
            class in_str
        >
struct eval_vec_functor_vec_impl{};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_sparse
        >
struct eval_vec_functor_impl<ret_type,in_type,accumulator,struct_dense,struct_sparse>
{ 
    static void eval(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using value_type    = typename ret_type::value_type;        
        using value_type_in = typename in_type::value_type;        
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        Integer r = m.rows();
        Integer c = m.cols();

        error::check_dim(d);        

        const details::sparse_ccs<value_type_in>& Ad = m.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const Integer* Ad_r         = Ad.ptr_r();
        const value_type_in* Ad_x   = Ad.ptr_x();

        if (d == 1)
        {            
            ret_type res(ret_ti, 1,c);
            value_type* ptr_res = res.ptr();

            accum.set_size(r);

            if (c == 0)
            {
                ret = matcl::Matrix(res,false);
                return;
            };

            if (r == 0)
            {
                accum.reset();
                const value_type& val =  accum.value();

                for (Integer j = 0; j < c; ++j)
                    ptr_res[j] = val;

                ret = matcl::Matrix(res,true);
                return;
            };            

            for (Integer j = 0; j < c; ++j)
            {
                accum.reset();

                Integer nz      = Ad_c[j + 1] - Ad_c[j];
                bool ret_known  = false;

                if (r-nz > 0)
                    ret_known = accum.add_zero();

                if (ret_known)
                {
                    ptr_res[j] = accum.value();
                }
                else
                {
                    for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
                    {
                        if (accum.add(Ad_x[k]))
                            break;
                    };

                    ptr_res[j] = accum.value();
                };            
            };

            ret = matcl::Matrix(res,true);
            return;
        };
     
        ret_type res(ret_ti, r, 1);
        value_type* ptr_res = res.ptr();
        accum.set_size(c);

        if (r == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };        

        if (c == 0)
        {
            accum.reset();
            const value_type& val =  accum.value();

            for (Integer j = 0; j < r; ++j)
                ptr_res[j] = val;

            ret = matcl::Matrix(res,true);
            return;
        };
        
        accum.reset_array(r);

        matcl::pod_workspace<Integer> nz(r,0);

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                Integer p = Ad_r[k];
                ++nz[p];
                accum.add(p,Ad_x[k]);
            };
        };
        for (Integer i = 0; i < r; ++i)
        {
            if (nz[i] < c)
                accum.add_zero(i);

            ptr_res[i] = accum.value(i);
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using full_matrix = Matrix<typename in_type::value_type,struct_dense>;
        eval_vec_functor_impl<ret_type,full_matrix,accumulator,struct_dense,struct_dense>
                ::eval_cum(ret, full(m), d, accum);
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_sparse
        >
struct eval_vec_functor_vec_impl<ret_type,in_type,accumulator,struct_sparse>
{ 
    static void eval(ret_type& ret, const in_type &m, accumulator& accum)
    {
        using value_type    = ret_type;
        using value_type_in = typename in_type::value_type;        
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        Integer r = m.rows();
        Integer c = m.cols();

        const details::sparse_ccs<value_type_in>& Ad = m.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const value_type_in* Ad_x   = Ad.ptr_x();

        accum.set_size(r,c);
        accum.reset();

        Integer nz = Ad_c[c] - Ad_c[0];
        if ( (nz / r) / c == 0)
        {
            if (accum.add_zero())
            {
                ret = accum.value();
                return;
            }
        };

        for (Integer j = 0; j < c; ++j)
        {            
            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                if (accum.add(Ad_x[k]))
                    break;
            };
        };

        ret = accum.value();
        return;
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_dense
        >
struct eval_vec_functor_impl<ret_type,in_type,accumulator,struct_dense,struct_dense>
{ 
    static void eval(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using value_type    = typename ret_type::value_type;
        using value_type_in = typename in_type::value_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        ret_type res(ret_ti);        
        const value_type_in* ptr_m = m.ptr();

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();
     
        error::check_dim(d);        
     
        if (d == 1)
        {
            accum.set_size(r);
            res.reset_unique(1, c);
            value_type* ptr_res = res.ptr();

            if (r == 0)
            {
                accum.reset();
                const value_type& val =  accum.value();

                for (Integer j = 0; j < c; ++j)
                    ptr_res[j] = val;

                ret = matcl::Matrix(res,true);
                return;
            };

            for (Integer j = 0; j < c; ++j)
            {
                accum.reset();
                for (Integer i = 0; i < r; ++i)
                {
                    if (accum.add(ptr_m[i]))
                        break;
                };

                ptr_res[j]  = accum.value();
                ptr_m       += m_ld;
            };

            ret = matcl::Matrix(res,true);
            return;
        };
     
        accum.set_size(c);
        res.reset_unique(r, 1);
        value_type* ptr_res = res.ptr();
     
        if (c == 0)
        {
            accum.reset();
            const value_type& val =  accum.value();

            for (Integer j = 0; j < r; ++j)
                ptr_res[j] = val;

            ret = matcl::Matrix(res,true);
            return;
        };

        for (Integer i = 0; i < r; ++i)
        {
            accum.reset();
            ptr_m = m.ptr() + i;

            for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
            {
                if (accum.add(*ptr_m))
                    break;
            };

            ptr_res[i] = accum.value();
        }
     
        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();

        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        ret_type res(ret_ti, r,c);
        const typename in_type::value_type* ptr_m = m.ptr();
        typename ret_type::value_type* ptr_res = res.ptr();

        Integer res_ld  = res.ld();

        error::check_dim(d);

        if (d == 1)
        {
            accum.set_size(r);

            for (Integer j = 0; j < c; ++j)
            {
                accum.reset();
                for (Integer i = 0; i < r; ++i)
                {
                    accum.add(ptr_m[i]);
                    ptr_res[i] = accum.value();
                }
                ptr_m   += m_ld;
                ptr_res += res_ld;
            }

            ret = matcl::Matrix(res,true);
            return;
        }
        else
        {
            accum.set_size(c);

            for (Integer i = 0; i < r; ++i)
            {            
                accum.reset();

                ptr_m = m.ptr() + i;
                ptr_res = res.ptr() + i;

                for (Integer j = 0; j < c; ++j)
                {
                    accum.add(*ptr_m);
                    *ptr_res = accum.value();

                    ptr_m   += m_ld;
                    ptr_res += res_ld;
                }
            }

            ret = matcl::Matrix(res,true);
            return;
        }
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_dense
        >
struct eval_vec_functor_vec_impl<ret_type,in_type,accumulator,struct_dense>
{ 
    static void eval(ret_type& ret, const in_type &m, accumulator& accum)
    {
        using value_type    = ret_type;
        using value_type_in = typename in_type::value_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();
        const value_type_in* ptr_m = m.ptr(); 

        accum.set_size(r,c);
        accum.reset();

        Integer nz      = m.nnz();

        if ((nz / r) / c == 0)
        {
            if (accum.add_zero())
            {
                ret = accum.value();
                return;
            };
        };

        for (Integer j = 0; j < c; ++j)
        {            
            for (Integer i = 0; i < r; ++i)
            {
                if (accum.add(ptr_m[i]))
                    break;
            };

            ptr_m += m_ld;
        };

        ret = accum.value();
        return;
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_banded
        >
struct eval_vec_functor_impl<ret_type,in_type,accumulator,struct_dense,struct_banded>
{
    static void eval(matcl::Matrix& ret, const in_type& m, int d, accumulator& accum)
    {
        using VTR           = typename ret_type::value_type;
        using VT            = typename in_type::value_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        ret_type res(ret_ti);        

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();
        Integer fd      = m.first_diag();
        Integer ld      = m.last_diag();
     
        error::check_dim(d);              
     
        if (d == 1)
        {
            accum.set_size(r);
            res.reset_unique(1, c);
            VTR* ptr_res = res.ptr();

            if (r == 0)
            {
                accum.reset();
                const VTR& val =  accum.value();

                for (Integer j = 0; j < c; ++j)
                    ptr_res[j] = val;

                ret = matcl::Matrix(res,true);
                return;
            };

            if (fd == ld)
            {
                const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd);  
                Integer s       = m.diag_length(fd);
                Integer r0      = m.first_row_on_diag(fd);
                Integer c0      = m.first_col_on_diag(fd);

                for (Integer j = 0; j < c0; ++j)
                {
                    accum.reset();
                    accum.add_zero();
                    ptr_res[j] = accum.value();
                };

                for (Integer j = c0; j < c0 + s; ++j, ptr_m += m_ld, ++r0)
                {
                    accum.reset();

                    if (0 < r0 || r0 + 1 < r)
                    {
                        bool is_known = accum.add_zero();
                        if (is_known)
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    accum.add(ptr_m[0]);
                    ptr_res[j] = accum.value();
                };

                for (Integer j = c0 + s; j < c; ++j)
                {
                    accum.reset();
                    accum.add_zero();
                    ptr_res[j] = accum.value();
                };
            }
            else
            {
                const VT* ptr_m = m.rep_ptr();  

                for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
                {
                    accum.reset();

                    Integer row_f = m.first_row(j);
                    Integer row_l = m.last_row(j);
                    Integer row_p = m.first_elem_pos(j);

                    if (row_f > 0 || row_l < r-1)
                    {
                        bool is_known = accum.add_zero();
                        if (is_known)
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                    {
                        if (accum.add(ptr_m[row_p]))
                            break;
                    };

                    ptr_res[j] = accum.value();
                };
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        accum.set_size(c);
        res.reset_unique(r, 1);
        VTR* ptr_res = res.ptr();
     
        if (c == 0)
        {
            accum.reset();
            const VTR& val =  accum.value();

            for (Integer j = 0; j < r; ++j)
                ptr_res[j] = val;

            ret = matcl::Matrix(res,true);
            return;
        };

        if (fd == ld)
        {
            const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd);  
            Integer s       = m.diag_length(fd);
            Integer r0      = m.first_row_on_diag(fd);
            Integer c0      = m.first_col_on_diag(fd);

            for (Integer i = 0; i < r0; ++i)
            {
                accum.reset();
                accum.add_zero();
                ptr_res[i] = accum.value();
            }

            for (Integer i = r0; i < r0 + s; ++i, ptr_m += m_ld, ++c0)
            {
                accum.reset();

                if (0 < c0 || c0 + 1 < c)
                {
                    bool is_known = accum.add_zero();
                    if (is_known)
                    {
                        ptr_res[i] = accum.value();
                        continue;
                    };
                };

                accum.add(ptr_m[0]);
                ptr_res[i] = accum.value();
            }

            for (Integer i = r0 + s; i < r; ++i)
            {
                accum.reset();
                accum.add_zero();
                ptr_res[i] = accum.value();
            }
        }
        else
        {
            const VT* ptr_m = m.rep_ptr();
            Integer dpos    = m_ld - 1;

            for (Integer i = 0; i < r; ++i)
            {
                accum.reset();

                Integer col_f = m.first_col(i);
                Integer col_l = m.last_col(i);
                Integer pos   = m.first_elem_pos_row(i);

                if (col_f > 0 || col_l < c-1)
                {
                    bool is_known = accum.add_zero();
                    if (is_known)
                    {
                        ptr_res[i] = accum.value();
                        continue;
                    };
                };

                for (Integer j = col_f; j <= col_l; ++j, pos += dpos)
                {
                    if (accum.add(ptr_m[pos]))
                        break;
                };            

                ptr_res[i] = accum.value();
            }
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using full_matrix   = Matrix<typename in_type::value_type,struct_dense>;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        ret_type mp = converter<full_matrix,in_type>::eval(m);
        eval_vec_functor_impl<ret_type,full_matrix,accumulator,struct_dense,struct_dense>
                    ::eval_cum(ret,mp,d,accum);
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_banded
        >
struct eval_vec_functor_vec_impl<ret_type,in_type,accumulator,struct_banded>
{
    static void eval(ret_type& ret, const in_type& m, accumulator& accum)
    {
        using VTR           = ret_type;
        using VT            = typename in_type::value_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();
        Integer fd      = m.first_diag();
        Integer ld      = m.last_diag();        
     
        accum.set_size(r,c);
        accum.reset();

        Integer nz      = m.nnz();

        if ((nz / r) / c == 0)
        {
            if (accum.add_zero())
            {
                ret = accum.value();
                return;
            };
        };

        if (fd == ld)
        {
            const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd);
            Integer s       = m.diag_length(fd);

            for (Integer j = 0; j < s; ++j, ptr_m += m_ld)
            {
                if (accum.add(ptr_m[0]))
                    break;
            };
        }
        else
        {
            const VT* ptr_m = m.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
            {
                Integer row_f = m.first_row(j);
                Integer row_l = m.last_row(j);
                Integer row_p = m.first_elem_pos(j);

                for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                {
                    if (accum.add(ptr_m[row_p]))
                        break;
                };                
            };
        };

        ret = accum.value();
        return;
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
        >
struct eval_vec_functor
{ 
    static void eval(matcl::Matrix& ret, const in_type &x, int dim) 
    {
        using ret_str       = typename ret_type::struct_type;
        using in_str        = typename in_type::struct_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(x));
        accumulator accum(ret_ti);

        eval_vec_functor_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval(ret,x,dim,accum);
    };

    static void eval(matcl::Matrix& ret, const in_type &x, int dim, accumulator& accum) 
    {
        using ret_str   = typename ret_type::struct_type;
        using in_str    = typename in_type::struct_type;

        eval_vec_functor_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval(ret,x,dim,accum);        
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &x, int dim) 
    {
        using ret_str       = typename ret_type::struct_type;
        using in_str        = typename in_type::struct_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(x));

        accumulator accum(ret_ti);
        eval_vec_functor_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval_cum(ret,x,dim,accum);        
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &x, int dim, accumulator& accum) 
    {
        using ret_str   = typename ret_type::struct_type;
        using in_str    = typename in_type::struct_type;

        eval_vec_functor_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval_cum(ret,x,dim,accum);        
    };
};

template<   class ret_type, 
            class in_type, 
            class accumulator
        >
struct eval_vec_functor_vec
{ 
    static void eval(ret_type& ret, const in_type &x)
    {
        using in_str        = typename in_type::struct_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;

        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(x));
        accumulator accum(ret_ti);

        eval_vec_functor_vec_impl<ret_type,in_type,accumulator,in_str>::eval(ret,x,accum);
    };

    static void eval(ret_type& ret, const in_type &x, accumulator& accum)
    {
        using in_str    = typename in_type::struct_type;

        eval_vec_functor_vec_impl<ret_type,in_type,accumulator,in_str>::eval(ret,x,accum);
    };
};

}}

#pragma warning( pop )
