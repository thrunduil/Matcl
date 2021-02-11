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

#include <vector>

namespace matcl { namespace raw
{

namespace md = matcl::details;

// accumulator
// value must not depend on number of zeros
template<class in_type, class out_type>
class accumulator2	
{
    protected:
        out_type				state;
        md::workspace2<out_type>state_array;	//vector of states, need to be set only if
                                                //reset_array is called

    public:
        accumulator2(){};

        void set_size(Integer r, Integer c) {};	//set number of elements in vector and number of vectors
        void reset() {};						//reset state
        void current_vector(Integer pos) {};	//set current vector number (0-based)
        void reset_array(Integer s) {};			//reset state array of size s
        bool add(Integer v, const in_type& value){};	//returns true is value is known
                                                //v is element index in current vector (0-based)
        void add(Integer p,Integer v, const in_type& value){};	
                                                //add new value for p-th vector (zero based)
                                                //v is element index in current vector
        bool add_zero(Integer v){};				//at least one zero exists returns true is value is known
                                                //v is index of one of zero elements in current vector 
        void add_zero(Integer p,Integer v){};	//at least one zero exists at pos p, returns true is value 
                                                //is known
                                                //v is index of one of zero elements in current vector 
        out_type value() const {};				//	
        out_type value(Integer p) const {};		//return value at pos p	
};

template<	class ret_type, 
            class in_type, 
            class accumulator,
            class ret_str,
            class in_str
        >
struct eval_vec_functor2_impl{};

template<	class ret_type, 
            class in_type, 
            class accumulator,
            class in_str
        >
struct eval_vec_functor2_vec_impl{};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_sparse
        >
struct eval_vec_functor2_impl<ret_type,in_type,accumulator,struct_dense,struct_sparse>
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
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const value_type_in* Ad_x	= Ad.ptr_x();

        if (d == 1)
        {			
            ret_type res(ret_ti,1,c);
            value_type* ptr_res = res.ptr();
            accum.set_size(r,c,d);

            if (c == 0)
            {
                ret = matcl::Matrix(res,false);
                return;
            };

            if (r == 0)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    accum.reset();
                    accum.current_vector(j);

                    ptr_res[j] = accum.value();
                };

                ret = matcl::Matrix(res,true);
                return;
            };			

            for (Integer j = 0; j < c; ++j)
            {
                accum.reset();
                accum.current_vector(j);

                Integer zp = 0;
                for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
                {
                    Integer p = Ad_r[k];
                    if (zp == p)
                        ++zp;

                    if (accum.add(Ad_r[k],Ad_x[k]))
                        break;
                };

                if (zp < r)
                    accum.add_zero(zp);

                ptr_res[j] = accum.value();
            };

            ret = matcl::Matrix(res,true);
            return;
        };
     
        ret_type res(ret_ti,r, 1);
        value_type* ptr_res = res.ptr();
        accum.set_size(c,r,d);

        if (r == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        if (c == 0)
        {
            for (Integer j = 0; j < r; ++j)
            {
                accum.reset();
                accum.current_vector(j);

                ptr_res[j] = accum.value();
            };

            ret = matcl::Matrix(res,true);
            return;
        };
        
        accum.reset_array(r);

        matcl::pod_workspace<Integer> zp(r,0);

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                Integer p = Ad_r[k];
                if (zp[p] == j)
                    ++zp[p];

                accum.add(p,j,Ad_x[k]);
            };
        };

        for (Integer i = 0; i < r; ++i)
        {
            if (zp[i] < c)
                accum.add_zero(i,zp[i]);

            ptr_res[i] = accum.value(i);
        };

        ret = matcl::Matrix(res,true);
        return;	
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using full_matrix = Matrix<in_type,struct_dense>;
        return eval_vec_functor2_impl<ret_type,full_matrix,accumulator,struct_dense,struct_dense>
                        ::eval_cum(ret, full(m),d,accum);
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_sparse
        >
struct eval_vec_functor2_vec_impl<ret_type,in_type,accumulator,struct_sparse>
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
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const value_type_in* Ad_x	= Ad.ptr_x();

        accum.set_size(r,c);
        accum.reset();

        for (Integer j = 0; j < c; ++j)
        {
            accum.current_vector(j);

            Integer zp = 0;
            for (Integer k = Ad_c[j] ; k < Ad_c[j + 1] ; ++k)
            {
                Integer p = Ad_r[k];
                if (zp == p)
                    ++zp;

                if (accum.add(Ad_r[k],Ad_x[k]))
                {
                    ret = accum.value();
                    return;
                }
            };

            if (zp < r)
            {
                if (accum.add_zero(zp))
                {
                    ret = accum.value();
                    return;
                };
            }
        };

        ret = accum.value();
        return;
    }; 
};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_dense
        >
struct eval_vec_functor2_impl<ret_type,in_type,accumulator,struct_dense,struct_dense>
{ 
    static void eval(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using value_type    = typename ret_type::value_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        ret_type res(ret_ti);

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();
     
        error::check_dim(d);
        const typename in_type::value_type* ptr_m = m.ptr();
     
        if (d == 1)
        {
            accum.set_size(r,c,d);
            res.reset_unique(1, c);
            value_type* ptr_res = res.ptr();

            if (r == 0)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    accum.reset();
                    accum.current_vector(j);
                    const value_type& val =  accum.value();
                    ptr_res[j] = val;
                };

                ret = matcl::Matrix(res,true);
                return;
            };

            for (Integer j = 0; j < c; ++j)
            {
                accum.reset();
                accum.current_vector(j);

                for (Integer i = 0; i < r; ++i)
                {
                    if (accum.add(i,ptr_m[i]))
                        break;
                };

                ptr_m       += m_ld;
                ptr_res[j]  = accum.value();
            };

            ret = matcl::Matrix(res,true);
            return;
        };
     
        accum.set_size(c,r,d);
        res.reset_unique(r, 1);
        value_type* ptr_res = res.ptr();
     
        if (c == 0)
        {
            for (Integer j = 0; j < r; ++j)
            {
                accum.reset();
                accum.current_vector(j);
                const value_type& val =  accum.value();
                ptr_res[j] = val;
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        for (Integer i = 0; i < r; ++i)
        {
            accum.reset();
            accum.current_vector(i);

            ptr_m = m.ptr() + i;

            for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
            {
                if (accum.add(j,*ptr_m))
                    break;
            };

            ptr_res[i] = accum.value();
        }
     
        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using value_type = typename ret_type::value_type;

        Integer i, ii, j, r = m.rows(), c = m.cols();

        ret_type res(r,c); 		

        error::check_dim(d);
        const typename in_type::value_type* ptr_m = m.ptr();
        value_type* ptr_res = res.ptr();

        if (d == 1)
        {
            accum.set_size(r,c,d);
            accum.reset();

            for (Integer j = 0, ii = 1; j < c; ++j)
            {
                accum.current_vector(j);
                for (Integer i = 0; i < r; ++i, ++ii)
                {
                    accum.add(i,ptr_m[i]);
                    ptr_res[i] = accum.value();
                }
                ptr_m   += m_ld;
                ptr_res += res_ld;
            }

            ret = matcl::Matrix(res,true);
            return;
        }

        accum.set_size(c,r,d);
        for (Integer i = 0; i < r; ++i)
        {			
            accum.reset();
            accum.current_vector(i);

            ptr_m   = m.ptr() + i;
            ptr_res = res.ptr() + i;

            for (Integer j = 0; j < c; ++j)
            {
                accum.add(j,*ptr_m);
                *ptr_res = accum.value();

                ptr_m   += m_ld;
                ptr_res += res_ld;
            }
        }

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_dense
        >
struct eval_vec_functor2_vec_impl<ret_type,in_type,accumulator,struct_dense>
{ 
    static void eval(ret_type& ret, const in_type &m, accumulator& accum)
    {
        using value_type    = ret_type;
        using ret_ti_type   = typename ti::get_ti_type<ret_type>::type;
        ret_ti_type ret_ti  = accumulator::get_ret_ti(ti::get_ti(m));

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();

        const typename in_type::value_type* ptr_m = m.ptr();
     
        accum.set_size(r,c);
        accum.reset();

        for (Integer j = 0; j < c; ++j)
        {            
            accum.current_vector(j);

            for (Integer i = 0; i < r; ++i)
            {
                if (accum.add(i,ptr_m[i]))
                {
                    ret = accum.value();
                    return;
                }
            };

            ptr_m       += m_ld;
        };

        ret = accum.value();
        return;
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class ret_str = struct_dense,
            //class in_str = struct_banded
        >
struct eval_vec_functor2_impl<ret_type,in_type,accumulator,struct_dense,struct_banded>
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
            accum.set_size(r,c,d);
            res.reset_unique(1, c);
            
            VTR* ptr_res = res.ptr();

            if (r == 0)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    accum.reset();
                    accum.current_vector(j);
                    ptr_res[j] = accum.value();
                };

                ret = matcl::Matrix(res,true);
                return;
            };            

            if (fd == ld)
            {
                const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd); 
                Integer s       = m.diag_length(fd);
                Integer c0      = m.first_col_on_diag(fd);
                Integer r0      = m.first_row_on_diag(fd);

                for (Integer j = 0; j < c0; ++j)
                {
                    accum.reset();
                    accum.current_vector(j);

                    accum.add_zero(0);

                    ptr_res[j] = accum.value();
                };

                for (Integer j = c0; j < c0 + s; ++j, ptr_m += m_ld, ++r0)
                {
                    accum.reset();
                    accum.current_vector(j);

                    if (0 < r0)
                    {
                        if (accum.add_zero(0))
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    if (accum.add(r0,ptr_m[0]))
                    {
                        ptr_res[j] = accum.value();
                        continue;
                    };

                    if (r0 + 1 < r)
                    {
                        if (accum.add_zero(r0+1))
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    ptr_res[j] = accum.value();
                };

                for (Integer j = c0 + s; j < c; ++j)
                {
                    accum.reset();
                    accum.current_vector(j);

                    accum.add_zero(0);

                    ptr_res[j] = accum.value();
                };
            }
            else
            {
                const VT* ptr_m = m.rep_ptr(); 

                for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
                {
                    accum.reset();
                    accum.current_vector(j);

                    Integer row_f = m.first_row(j);
                    Integer row_l = m.last_row(j);
                    Integer row_p = m.first_elem_pos(j);

                    if (row_f > 0)
                    {
                        if (accum.add_zero(0))
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                    {
                        if (accum.add(i,ptr_m[row_p]))
                            break;
                    };

                    if (row_l < r-1)
                    {
                        if (accum.add_zero(row_l+1))
                        {
                            ptr_res[j] = accum.value();
                            continue;
                        };
                    };

                    ptr_res[j] = accum.value();                    
                };
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        accum.set_size(c,r,d);
        res.reset_unique(r, 1);
        VTR* ptr_res = res.ptr();
     
        if (c == 0)
        {
            for (Integer j = 0; j < r; ++j)
            {
                accum.reset();
                accum.current_vector(j);
                ptr_res[j] = accum.value();
            };

            ret = matcl::Matrix(res,true);
            return;
        };

        if (fd == ld)
        {
            const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd); 
            Integer s       = m.diag_length(fd);
            Integer c0      = m.first_col_on_diag(fd);
            Integer r0      = m.first_row_on_diag(fd);

            for (Integer i = 0; i < r0; ++i)
            {
                accum.reset();
                accum.current_vector(i);
                accum.add_zero(0);

                ptr_res[i] = accum.value();
            }

            for (Integer i = r0; i < r0 + s; ++i, ptr_m += m_ld, ++c0)
            {
                accum.reset();
                accum.current_vector(i);

                if (0 < c0 && accum.add_zero(0) )
                {
                    ptr_res[i] = accum.value();
                    continue;
                };

                if (accum.add(c0,ptr_m[0]))
                {
                    ptr_res[i] = accum.value();
                    continue;
                };

                if (c0 + 1 < c)
                    accum.add_zero(c0+1);

                ptr_res[i] = accum.value();
            }

            for (Integer i = r0 + s; i < r; ++i)
            {
                accum.reset();
                accum.current_vector(i);
                accum.add_zero(0);

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
                accum.current_vector(i);

                Integer col_f = m.first_col(i);
                Integer col_l = m.last_col(i);
                Integer pos   = m.first_elem_pos_row(i);

                bool is_known = false;

                if (col_f > 0)
                    is_known = accum.add_zero(0);

                if (is_known)
                {
                    ptr_res[i] = accum.value();
                    continue;
                };

                for (Integer j = col_f; j <= col_l; ++j, pos += dpos)
                {
                    if (accum.add(j,ptr_m[pos]))
                        break;
                };

                if (col_l < c-1)
                    accum.add_zero(col_l+1);

                ptr_res[i] = accum.value();
            }
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_cum(matcl::Matrix& ret, const in_type &m, int d, accumulator& accum)
    {
        using full_matrix = Matrix<in_type,struct_dense>;
        ret_type mp = converter<full_matrix,in_type>::eval(m);
        return eval_vec_functor2_impl<ret_type,full_matrix,accumulator,struct_dense,struct_dense>
                    ::eval_cum(ret,mp,d,accum);
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
            //class in_str = struct_banded
        >
struct eval_vec_functor2_vec_impl<ret_type,in_type,accumulator,struct_banded>
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

        if (fd == ld)
        {
            const VT* ptr_m = m.rep_ptr() + m.first_elem_diag(fd);
            Integer s       = m.diag_length(fd);
            Integer c0      = m.first_col_on_diag(fd);
            Integer r0      = m.first_row_on_diag(fd);

            if (c0 > 0)
            {
                accum.current_vector(0);
                bool is_known = accum.add_zero(0);
                
                if (is_known)
                {
                    ret = accum.value();
                    return;
                };
            };

            for (Integer j = 0; j < s; ++j, ptr_m += m_ld, ++c0, ++r0)
            {                
                accum.current_vector(c0);

                if (r0 > 0)
                {
                    bool is_known = accum.add_zero(0);
                    if (is_known)
                    {
                        ret = accum.value();
                        return;
                    };
                };

                if (accum.add(r0,ptr_m[0]))
                {
                    ret = accum.value();
                    return;
                }

                if (r0 + 1 < r-1)
                {
                    bool is_known = accum.add_zero(r0 + 1);
                    if (is_known)
                    {
                        ret = accum.value();
                        return;
                    };
                };
            };

            if (c0 < c - 1)
            {
                accum.current_vector(c0);
                bool is_known = accum.add_zero(0);
                if (is_known)
                {
                    ret = accum.value();
                    return;
                };
            };

        }
        else
        {
            const VT* ptr_m = m.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_m += m_ld)
            {
                accum.current_vector(j);

                Integer row_f = m.first_row(j);
                Integer row_l = m.last_row(j);
                Integer row_p = m.first_elem_pos(j);

                if (row_f > 0)
                {
                    bool is_known = accum.add_zero(0);
                    if (is_known)
                    {
                        ret = accum.value();
                        return;
                    };
                };

                for (Integer i = row_f; i <= row_l; ++i, ++row_p)
                {
                    if (accum.add(i,ptr_m[row_p]))
                    {
                        ret = accum.value();
                        return;
                    };
                };

                if (row_l < r-1)
                {
                    bool is_known = accum.add_zero(row_l+1);
                    if (is_known)
                    {
                        ret = accum.value();
                        return;
                    };
                };
            };
        };

        ret = accum.value();
        return;
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
        >
struct eval_vec_functor2
{ 
    static void eval(matcl::Matrix& ret, const in_type &x, int dim, accumulator& accum) 
    {
        using ret_str   = typename ret_type::struct_type;
        using in_str    = typename in_type::struct_type;

        return eval_vec_functor2_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval(ret,x,dim,accum);		
    }

    static void eval_cum(matcl::Matrix& ret, const in_type &x, int dim, accumulator& accum) 
    {
        using ret_str   = typename ret_type::struct_type;
        using in_str    = typename in_type::struct_type;

        return eval_vec_functor2_impl<ret_type,in_type,accumulator,ret_str,in_str>::eval_cum(ret,x,dim,accum);		
    };
};

template<	class ret_type, 
            class in_type, 
            class accumulator
        >
struct eval_vec_functor_vec2
{ 
    static void eval(ret_type& ret, const in_type &x, accumulator& accum)
    {
        using in_str    = typename in_type::struct_type;

        return eval_vec_functor2_vec_impl<ret_type,in_type,accumulator,in_str>::eval(ret,x,accum);		
    };
};

}}
