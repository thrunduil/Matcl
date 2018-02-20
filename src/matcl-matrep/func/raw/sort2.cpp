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

#include "matcl-matrep/func/raw/sort2.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-matrep/base/sort_iterator.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-matrep/base/colon_info.h"
#include "matcl-matrep/algs/sparse_algs.h"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/func/raw/eval_vec_functor2.h"
#include "matcl-internals/func/converter.h"

#include <vector>
#include <algorithm>
#include <boost/iterator/iterator_facade.hpp>

namespace matcl { namespace raw { namespace details
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template<class T>
struct sort_greater_nan
{
    static bool eval(const T& A, const T& B)
    {
        if (isnan_helper<T>::eval(B))
            return false;

        if (isnan_helper<T>::eval(A))
            return true;

        return (bool)gt_helper<T,T>::eval(A,B);
    };
};

template<class T>
struct sort_less_nan
{
    static bool eval(const T& A, const T& B)
    {
        if (isnan_helper<T>::eval(A))
            return false;

        if (isnan_helper<T>::eval(B))
            return true;

        return (bool)lt_helper<T,T>::eval(A,B);
    };
};

template<class T>
struct sort_leq_nan
{
    static bool eval(const T& A, const T& B)
    {
        if (isnan_helper<T>::eval(A))
        {
            if (isnan_helper<T>::eval(B))
                return true;
            else
                return false;
        };

        if (isnan_helper<T>::eval(B))
            return true;

        return (bool)leq_helper<T,T>::eval(A,B);
    };
};

template<class T>
class iterator_helper_1 : public boost::iterator_facade
                                <   iterator_helper_1<T>,
                                    T,
                                    boost::random_access_traversal_tag
                                >
{
    private:
        T*              x_ptr;
        Integer         ld;

    public:
        using parent_class      = boost::iterator_facade<iterator_helper_1<T>, T, 
                                        boost::random_access_traversal_tag>;
        using difference_type   = typename parent_class::difference_type;

        iterator_helper_1()
            :x_ptr(nullptr),ld(0)
        {};
        
        iterator_helper_1(T* x_ptr, Integer ld) 
            :x_ptr(x_ptr), ld(ld)
        {};

    private:
        friend class boost::iterator_core_access;

        void increment() 
        { 
            x_ptr += ld; 
        };

        void decrement() 
        { 
            x_ptr -= ld; 
        };
        
        bool equal(iterator_helper_1 const& other) const
        {
            return this->x_ptr == other.x_ptr;
        }
        
        T& dereference() const 
        { 
            return *x_ptr; 
        };
        
        difference_type distance_to(iterator_helper_1 j) const
        {
            return (j.x_ptr - this->x_ptr)/ld;
        }
        
        void advance(difference_type n)
        {
            this->x_ptr += imult(ld,cast_int64(static_cast<Integer_64>(n)));
        };
};

template<class value_type, bool asceding>
struct value_compare
{
    bool operator()(const value_type& A, const value_type& B)
    {
        return md::compare_nan<value_type,asceding>::eval(A,B);
    };
};

template<class ret_type, class in_type, class struct_type>
struct sort_impl{};

template<class ret_type, class in_type>
struct sort_impl<ret_type,in_type,struct_dense>
{
    using value_type_ret    = typename ret_type::value_type;
    using ret_ti_t          = ti::ti_type<value_type_ret>;

    static void eval(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& m, Integer dim, bool asceding)
    {
        using value_type = typename in_type::value_type;        

        error::check_dim(dim);

        Integer r = m.rows(), c = m.cols();

        if (r == 0 || c == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        ret_type tmp = m.copy();
        tmp.get_struct().reset();

        value_type_ret* ptr_tmp = tmp.ptr();
        Integer tmp_ld          = tmp.ld();

        if (dim == 1)
        {
            if (asceding == true)
            {
                for (Integer i = 0; i < c; ++i)
                {
                    std::stable_sort(ptr_tmp, ptr_tmp + r,value_compare<value_type,true>());
                    ptr_tmp += tmp_ld;
                };
            }
            else
            {
                for (Integer i = 0; i < c; ++i)
                {
                    std::stable_sort(ptr_tmp, ptr_tmp + r,value_compare<value_type,false>());
                    ptr_tmp += tmp_ld;
                };
            };
        }
        else
        {
            Integer lda = tmp.ld();
            if (asceding == true)
            {
                for (Integer i = 0; i < r; ++i)
                {
                    iterator_helper_1<value_type> it_begin(ptr_tmp+i,lda);
                    iterator_helper_1<value_type> it_end(ptr_tmp+i+imult(c,lda),lda);
                    std::stable_sort(it_begin, it_end,value_compare<value_type,true>());
                };
            }
            else
            {
                for (Integer i = 0; i < r; ++i)
                {
                    iterator_helper_1<value_type> it_begin(ptr_tmp+i,lda);
                    iterator_helper_1<value_type> it_end(ptr_tmp+i+imult(c,lda),lda);
                    std::stable_sort(it_begin, it_end,value_compare<value_type,false>());
                };
            };
        };

        ret = matcl::Matrix(tmp,false);
        return;
    };

    static void eval(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& m, 
                     Integer dim, bool asceding)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();
        error::check_dim(dim);

        integer_dense ind = integer_dense(ti::ti_empty());
        ind.reset_unique(r,c);

        if (r == 0 || c == 0) 
        {
            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(ret_type(ret_ti,r, c), false);
            return;
        };

        ret_type tmp = m.copy();
        tmp.get_struct().reset();

        Integer* ptr_ind = ind.ptr();
        value_type_ret* ptr_tmp = tmp.ptr();
        Integer ind_ld          = ind.ld();
        Integer tmp_ld          = tmp.ld();

        if (dim == 1)
        {
            for (Integer j = 0; j < c; ++j)
            {
                for (Integer i = 0; i < r; ++i)
                    ptr_ind[i] = i+1;

                ptr_ind += ind_ld;
            };

            ptr_ind = ind.ptr();

            for (Integer i = 0; i < c; ++i)
            {
                using iterator  = matcl::details::iterator_helper_2<value_type,Integer>;

                iterator it_begin(ptr_tmp,ptr_ind,1);
                iterator it_end(ptr_tmp + r,ptr_ind + r,1);

                if (asceding == true)
                    std::stable_sort(it_begin, it_end,
                                     matcl::details::value_ind_compare<value_type,Integer,true>());
                else
                    std::stable_sort(it_begin, it_end,
                                     matcl::details::value_ind_compare<value_type,Integer,false>());

                ptr_tmp += tmp_ld;
                ptr_ind += ind_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < c; ++j)
            {
                for (Integer i = 0; i < r; ++i)
                    ptr_ind[i] = j+1;

                ptr_ind += ind_ld;
            };
            ptr_ind = ind.ptr();

            //ind and tmp are fresh and have ld == r
            Integer ld = r;
            for (Integer i = 0; i < r; ++i)
            {
                using iterator  = matcl::details::iterator_helper_2<value_type,Integer>;

                iterator it_begin(ptr_tmp+i,ptr_ind+i,ld);
                iterator it_end(ptr_tmp+i+imult(c,ld),ptr_ind+i+imult(c,ld), ld);

                if (asceding == true)
                    std::stable_sort(it_begin, it_end,
                                     matcl::details::value_ind_compare<value_type,Integer,true>());
                else
                    std::stable_sort(it_begin, it_end,
                                     matcl::details::value_ind_compare<value_type,Integer,false>());
            };
        };

        iret    = matcl::Matrix(ind,false);
        x       = matcl::Matrix(tmp, false);
        return;
    };
};

template<class ret_type, class in_type>
struct sort_impl<ret_type,in_type,struct_sparse>
{
    using ret_ti_t = typename ti::get_ti_type<ret_type>::type;

    static void eval(matcl::Matrix& ret, ret_ti_t, const in_type& m, Integer dim, bool asceding)
    {
        Integer nz = m.nnz();

        error::check_dim(dim);

        if (nz == 0)
        {
            ret = matcl::Matrix(m,false);
            return;
        };
        
        if (dim == 2)
        {
            in_type tmp(m.get_type());
            mrd::manip_trans_helper<in_type>::eval_trans(tmp,m);

            if (tmp.is_unique() == false)
                tmp.assign_to_fresh(tmp.copy());

            if (asceding == true)
                eval_inplace<true>(tmp);
            else
                eval_inplace<false>(tmp);

            in_type ret2(tmp.get_type());
            mrd::manip_trans_helper<in_type>::eval_trans(ret2,tmp);

            ret = matcl::Matrix(ret2,false);
            return;
        };

        in_type tmp = m.copy();
        tmp.get_struct().reset();

        if (asceding == true)
            eval_inplace<true>(tmp);
        else
            eval_inplace<false>(tmp);
        
        ret = matcl::Matrix(tmp,false);
        return;
    };
    
    template<bool asceding>
    static void eval_inplace(in_type& m)
    {
        Integer r = m.rows(), c = m.cols(), nz = m.nnz();

        if (nz == 0)
            return;

        using val_type      = typename in_type::value_type;
        const val_type Z    = md::default_value<val_type>(ti::get_ti(m));

        m.get_struct().reset();

        sparse_ccs<val_type>& d  = m.rep();

        const Integer* d_c  = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        val_type* d_x       = d.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            Integer start   = d_c[j];
            Integer end     = d_c[j+1];
            Integer nz_loc  = end - start;

            if (nz_loc == 0)
                continue;

            if (nz_loc == 1)
            {
                const val_type& tmp = d_x[start];

                if(md::compare_nan<val_type,asceding>::eval(tmp,Z))
                    d_r[start] = 0;
                else if (md::compare_nan<val_type,asceding>::eval(Z,tmp))
                    d_r[start] = r-1;
                
                continue;
            };

            std::stable_sort(d_x+start,d_x+start + nz_loc,value_compare<val_type,asceding>());

            Integer k = start, pos = 0;
            for (; k < end; ++k)
            {
                if (md::compare_nan<val_type,asceding>::eval(d_x[k],Z))
                {
                    d_r[k] = pos;
                    ++pos;
                }
                else if (eeq_helper<val_type,val_type>::eval(d_x[k], Z))
                {
                    //this breaks stability of sort, solution is to remove zeroes first
                    d_r[k] = pos;
                    ++pos;
                }
                else
                {
                    break;
                };
            };

            if (k == end)
                continue;

            pos += (r-nz_loc);

            for (; k < end; ++k)
            {
                d_r[k] = pos;
                ++pos;
            };
        };

        return;
    };

    static void eval(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ,const in_type& m, Integer dim, 
                     bool asceding)
    {
        Integer r = m.rows(), c = m.cols(), nz = m.nnz();

        error::check_dim(dim);

        integer_dense ind = integer_dense(ti::ti_empty());        

        if (nz == 0)
        {
            ind.reset_unique(r,c);

            Integer ind_ld      = ind.ld();
            Integer* ptr_ind    = ind.ptr();

            if (dim == 1)
            {
                for (Integer j = 0; j < c; ++j)
                {
                    for (Integer k = 0; k < r; ++k)
                        ptr_ind[k] = k+1;

                    ptr_ind += ind_ld;
                };
            }
            else
            {
                for (Integer j = 0; j < c; ++j)
                {
                    for (Integer k = 0; k < r; ++k)
                        ptr_ind[k] = j+1;

                    ptr_ind += ind_ld;
                };
            };

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(m,false);
            return;
        };

        if (dim == 2)
        {
            in_type tmp(m.get_type());
            mrd::manip_trans_helper<in_type>::eval_trans(tmp,m);

            if (tmp.is_unique() == false)
                tmp.assign_to_fresh(tmp.copy());

            if (asceding == true)
                eval_inplace<true>(tmp,ind,true);
            else
                eval_inplace<false>(tmp,ind,true);

            in_type ret2(tmp.get_type());
            mrd::manip_trans_helper<in_type>::eval_trans(ret2,tmp);

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(ret2,false);
            return;
        };

        in_type tmp = m.copy();
        tmp.get_struct().reset();

        if (asceding == true)
            eval_inplace<true>(tmp,ind,false);
        else
            eval_inplace<false>(tmp,ind,false);

        i = matcl::Matrix(ind,false);
        x = matcl::Matrix(tmp,false);
        return;
    };

    template<bool asceding>
    static void eval_inplace(in_type& m, integer_dense& ind, bool trans)
    {
        Integer r = m.rows(), c = m.cols(), nz = m.nnz();
        Integer ind_ld  = ind.ld();

        if (nz == 0)
        {            
            if (trans == false)
            {
                ind.reset_unique(r,c);
                Integer* ptr_ind = ind.ptr();
                for (Integer j = 0; j < c; ++j)
                {
                    for (Integer i = 0; i < r; ++i)
                        ptr_ind[i] = i+1;

                    ptr_ind += ind_ld;
                };
            }
            else
            {
                ind.reset_unique(c,r);
                Integer* ptr_ind = ind.ptr();

                for (Integer j = 0; j < r; ++j)
                {
                    for (Integer i = 0; i < c; ++i)
                        ptr_ind[i] = j+1;

                    ptr_ind += ind_ld;
                };
            };
            return;
        };

        if (trans == false)
            ind.reset_unique(r,c);
        else
            ind.reset_unique(c,r);

        using val_type      = typename in_type::value_type;
        const val_type Z    = md::default_value<val_type>(ti::get_ti(m));

        sparse_ccs<val_type>& d = m.rep();

        m.get_struct().reset();

        const Integer* d_c  = d.ptr_c();
        Integer* d_r        = d.ptr_r();
        val_type* d_x       = d.ptr_x();

        matcl::pod_workspace<Integer> flags(r,-1);
        Integer* ptr_ind = ind.ptr();
        
        for (Integer j = 0; j < c; ++j)
        {
            Integer start   = d_c[j];
            Integer end     = d_c[j+1];
            Integer nz_loc  = end - start;

            if (nz_loc == 0)
            {
                if (trans == false)
                {
                    Integer pos = imult(j,r);
                    for (Integer i = 0; i < r; ++i, ++pos)
                        ptr_ind[pos] = i+1;
                }
                else
                {
                    Integer pos = j;
                    for (Integer i = 0; i < r; ++i, pos += c)
                        ptr_ind[pos] = i+1;
                };
                continue;
            }

            using iterator = matcl::details::iterator_helper_2<val_type,Integer>;

            iterator it_begin(&d_x[start],&d_r[start],1);
            iterator it_end(&d_x[start] + nz_loc,&d_r[start] + nz_loc,1);

            std::stable_sort(it_begin, it_end,matcl::details::value_ind_compare<val_type,Integer,asceding>());

            for (Integer p = start; p < end; ++p)
                flags[d_r[p]] = j;

            Integer k = start, pos = 0, pos_ind, dpos_ind;
            if (trans == false)
            {
                pos_ind     = imult(j,r);
                dpos_ind    = 1;
            }
            else
            {
                pos_ind     = j;
                dpos_ind    = c;
            };
            for (; k < end; ++k)
            {
                if (md::compare_nan<val_type,asceding>::eval(d_x[k],Z))
                {
                    ptr_ind[pos_ind]    = d_r[k]+1;
                    d_r[k]              = pos;
                    ++pos;
                    pos_ind             += dpos_ind;
                }
                else if (mrd::is_zero(d_x[k]))
                {
                    flags[d_r[k]]       = j-1;
                    d_r[k]              = pos;
                    ++pos;
                }
                else
                {
                    break;
                };
            };

            for (Integer i = 0; i < r; ++i)
            {
                if (flags[i] < j)
                {
                    ptr_ind[pos_ind]    = i+1;
                    pos_ind             += dpos_ind;
                };
            };

            pos += (r-nz_loc);

            for (; k < end; ++k)
            {
                ptr_ind[pos_ind]        = d_r[k]+1;
                d_r[k]                  = pos;
                ++pos;
                pos_ind                 += dpos_ind;
            };
        };

        return;
    };
};

template<class T>
struct row_info
{
    private:
        const T*        ptr;
        Integer         ld;
        Integer         cols;
        Integer         pos;

    public:
        row_info(const T* ptr, Integer ld, Integer cols)    
            :ptr(ptr),ld(ld),cols(cols),pos(0)
        {};

        row_info(const T* ptr, Integer ld, Integer cols, Integer pos)    
            :ptr(ptr),ld(ld),cols(cols),pos(pos)
        {};

        row_info()
            :ptr(0),ld(0),cols(0),pos(0) 
        {};

        bool operator<(const row_info& other)
        {
            for (Integer i = 0, pos_loc = 0; i < cols; ++i, pos_loc += ld)
            {
                if (sort_less_nan<T>::eval(ptr[pos_loc],other.ptr[pos_loc]))
                    return true;

                if (sort_greater_nan<T>::eval(ptr[pos_loc],other.ptr[pos_loc]))
                    return false;
            };

            return false;
        };

        const T* get_ptr() const
        {
            return ptr;
        };

        Integer get_pos() const
        {
            return pos;
        };
};

template<class T>
struct row_info_dim
{
    private:
        const T*                ptr;
        Integer                 ld;
        Integer                 pos;
        const integer_dense*    dims;   //must be continuous

    public:
        row_info_dim(const T* ptr, Integer ld, const integer_dense* dims)    
            :ptr(ptr),ld(ld),pos(0), dims(dims)
        {};

        row_info_dim(const T* ptr, Integer ld, const integer_dense* dims, Integer pos)    
            :ptr(ptr),ld(ld),pos(pos),dims(dims)
        {};

        row_info_dim()
            :ptr(0),ld(0),pos(0),dims(0)
        {};

        bool operator<(const row_info_dim& other)
        {
            Integer dims_size       = dims->size();
            const Integer* ptr_dims = dims->ptr();

            for (Integer i = 0; i < dims_size; ++i)
            {
                Integer dim = ptr_dims[i];                
                if (dim > 0)
                {
                    Integer pos_loc = imult(dim-1,ld);
                    if (sort_less_nan<T>::eval(ptr[pos_loc],other.ptr[pos_loc]))
                        return true;

                    if (sort_greater_nan<T>::eval(ptr[pos_loc],other.ptr[pos_loc]))
                        return false;
                }
                else
                {
                    Integer pos_loc = imult(-dim-1,ld);
                    if (sort_greater_nan<T>::eval(ptr[pos_loc], other.ptr[pos_loc]))
                        return true;

                    if (sort_less_nan<T>::eval(ptr[pos_loc],other.ptr[pos_loc]))
                        return false;
                };
            };

            return false;
        };

        const T* get_ptr() const
        {
            return ptr;
        };

        Integer get_pos() const
        {
            return pos;
        };
};

template<class T>
bool row_vec_comparer(row_info<T>* ptr1, row_info<T>* ptr2)
{
    return *ptr1 < * ptr2;
};

template<class T>
static bool row_vec_dim_comparer(row_info_dim<T>* ptr1, row_info_dim<T>* ptr2)
{
    return *ptr1 < * ptr2;
};

template<class T>
struct row_info_sparse
{
    private:
        const Integer * d_r;
        const T*        d_x;
        Integer         first;
        Integer         last;
        Integer         col;
        T               Z;

    public:
        row_info_sparse(ti::ti_type<T> ti, const Integer *    d_r, const T* d_x, 
                        Integer f, Integer l, Integer col )    
            :d_r(d_r),d_x(d_x), first(f), last(l), col(col)
            ,Z(md::default_value<T>(ti))
        {};

        bool operator<(const row_info_sparse& other)
        {
            Integer k1 = first, k2 = other.first;

            while (k1 < last && k2 < other.last)
            {
                if (d_r[k1] < other.d_r[k2])
                {
                    if (sort_less_nan<T>::eval(d_x[k1], Z))
                        return true;

                    if (sort_greater_nan<T>::eval(d_x[k1], Z))
                        return false;

                    ++k1;
                }
                else if (d_r[k1] > other.d_r[k2])
                {
                    if (sort_less_nan<T>::eval(Z,other.d_x[k2]))
                        return true;

                    if (sort_greater_nan<T>::eval(Z,other.d_x[k2]))
                        return false;

                    ++k2; 
                }
                else
                { 
                    if (sort_less_nan<T>::eval(d_x[k1],other.d_x[k2]))
                        return true;

                    if (sort_greater_nan<T>::eval(d_x[k1],other.d_x[k2]))
                        return false;

                    ++k1; 
                    ++k2;
                };
            };
            while (k1 < last)
            {
                if (sort_less_nan<T>::eval(d_x[k1],Z))
                    return true;

                if (sort_greater_nan<T>::eval(d_x[k1], Z))
                    return false;

                ++k1;
            };
            while (k2 < other.last)
            {                
                if (sort_less_nan<T>::eval(Z,other.d_x[k2]))
                    return true;

                if (sort_greater_nan<T>::eval(Z, other.d_x[k2]))
                    return false;

                ++k2; 
            };

            return false;
        };

        Integer get_pos() const
        {
            return col;
        };
};

template<class T>
struct row_info_sparse_dims
{
    private:
        const Integer *         d_r;
        const T*                d_x;
        Integer                 first;
        Integer                 last;
        Integer                 col;
        const integer_dense*    dims;   //must be continuous
        T                       Z;

    public:
        row_info_sparse_dims(ti::ti_type<T> ti,const Integer * d_r,const T* d_x,
                            Integer f,Integer l,Integer col,const integer_dense* dims)    
            :d_r(d_r),d_x(d_x), first(f), last(l), col(col), dims(dims)
            ,Z(md::default_value<T>(ti))
        {};

        bool operator<(const row_info_sparse_dims& other)
        {
            Integer k1 = first, k2 = other.first;

            while (k1 < last && k2 < other.last)
            {
                if (d_r[k1] < other.d_r[k2])
                {
                    if (less(d_x[k1],Z,is_increasing(d_r[k1])))
                        return true;

                    if (greater(d_x[k1],Z,is_increasing(d_r[k1])))
                        return false;

                    ++k1;
                }
                else if (d_r[k1] > other.d_r[k2])
                {
                    if (less(Z,other.d_x[k2],is_increasing(other.d_r[k2])))
                        return true;

                    if (greater(Z,other.d_x[k2],is_increasing(other.d_r[k2])))
                        return false;

                    ++k2; 
                }
                else
                { 
                    if (less(d_x[k1], other.d_x[k2],is_increasing(d_r[k1])))
                        return true;

                    if (greater(d_x[k1], other.d_x[k2],is_increasing(d_r[k1])))
                        return false;

                    ++k1; 
                    ++k2;
                };
            };
            while (k1 < last)
            {
                if (less(d_x[k1], Z,is_increasing(d_r[k1])))
                    return true;

                if (greater(d_x[k1], Z,is_increasing(d_r[k1])))
                    return false;

                ++k1;
            };
            while (k2 < other.last)
            {                
                if (less(Z, other.d_x[k2],is_increasing(other.d_r[k2])))
                    return true;

                if (greater(Z, other.d_x[k2],is_increasing(other.d_r[k2])))
                    return false;

                ++k2; 
            };

            return false;
        };

        bool is_increasing(Integer pos)
        {
            return dims->ptr()[pos] > 0;
        };

        bool less(const T& x1, const T& x2,bool is_increase)
        {
            if (is_increase)
                return sort_less_nan<T>::eval(x1,x2);
            else
                return sort_greater_nan<T>::eval(x1, x2);
        };

        bool greater(const T& x1, const T& x2,bool is_increase)
        {
            if (is_increase)
                return sort_greater_nan<T>::eval(x1,x2);
            else
                return sort_less_nan<T>::eval(x1,x2);
        };

        Integer get_pos() const
        {
            return col;
        };
};

template<class T>
static bool row_vec_sparse_comparer(row_info_sparse<T>* ptr1, row_info_sparse<T>* ptr2)
{
    return *ptr1 < * ptr2;
};

template<class T>
static bool row_vec_sparse_dims_comparer(row_info_sparse_dims<T>* ptr1, row_info_sparse_dims<T>* ptr2)
{
    return *ptr1 < * ptr2;
};

static void check_dims(const integer_dense& dims,Integer c)
{
    error::check_col_indices_sortrows(dims.size(),c);
    const Integer* ptr_dims = dims.ptr();
    Integer r   = dims.rows();
    Integer dc  = dims.cols();
    Integer ld  = dims.ld();

    for (Integer j = 0; j < dc; ++j)
    {
        for (Integer i = 0; i < r; ++i)
            error::check_col_indices_elem_sortrows(ptr_dims[i],c);

        ptr_dims += ld;
    };
};

static void check_dims_col(const integer_dense& dims,Integer c)
{
    error::check_row_indices_sortcols(dims.size(),c);

    const Integer* ptr_dims = dims.ptr();
    Integer r   = dims.rows();
    Integer dc  = dims.cols();
    Integer ld  = dims.ld();

    for (Integer j = 0; j < dc; ++j)
    {
        for (Integer i = 0; i < r; ++i)
            error::check_row_indices_elem_sortcols(ptr_dims[i],c);

        ptr_dims += ld;
    };
};

template<class ret_type, class in_type, class struct_type>
struct sortrows_impl{};

template<class ret_type, class in_type>
struct sortrows_impl<ret_type,in_type,struct_dense>
{
    using ret_2             = std::pair<ret_type,integer_dense>;
    using value_type_ret    = typename ret_type::value_type;
    using ret_ti_t          = ti::ti_type<value_type_ret>;

    static void eval(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& m)
    {
        using value_type = typename in_type::value_type;        

        Integer r = m.rows(), c = m.cols();

        if (r == 0 || c == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        std::vector<row_info<value_type>>  row_vec(r);
        std::vector<row_info<value_type>*> row_vec_ptr(r);

        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer i = 0; i < r; ++i)
        {
            row_vec[i] = row_info<value_type>(ptr_m+i,m_ld,c);
            row_vec_ptr[i] = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < r; ++i)
        {
            const value_type* ptr    = row_vec_ptr[i]->get_ptr();
            for (Integer j = 0, pos = i, pos_p = 0; j < c; ++j, pos += res_ld, pos_p += m_ld)
                mrd::reset_helper(ptr_res[pos],ptr[pos_p]);
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_col(matcl::Matrix& ret, ret_ti_t ret_ti, const in_type& m)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();

        if (r == 0 || c == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        std::vector<row_info<value_type>>  row_vec(c);
        std::vector<row_info<value_type>*> row_vec_ptr(c);
        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer i = 0; i < c; ++i)
        {
            row_vec[i] = row_info<value_type>(ptr_m,1,r);
            row_vec_ptr[i] = &row_vec[i];
            ptr_m += m_ld;
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < c; ++i)
        {
            const value_type* ptr    = row_vec_ptr[i]->get_ptr();
            for (Integer j = 0; j < r; ++j)
                mrd::reset_helper(ptr_res[j],ptr[j]);

            ptr_res += res_ld;
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval2(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& m)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();
        integer_dense ind(ti::ti_empty(),r,1);

        if (r == 0 || c == 0) 
        {
            Integer* ptr_ind = ind.ptr();
            for (Integer i = 0; i < r; ++i)
                ptr_ind[i] = i+1;

            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        std::vector<row_info<value_type>>  row_vec(r);
        std::vector<row_info<value_type>*> row_vec_ptr(r);
        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer i = 0; i < r; ++i)
        {
            row_vec[i] = row_info<value_type>(ptr_m + i,m_ld,c,i+1);
            row_vec_ptr[i] = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        Integer* ptr_ind = ind.ptr();

        for (Integer i = 0; i < r; ++i)
        {
            const value_type* ptr    = row_vec_ptr[i]->get_ptr();
            for (Integer j = 0, pos = i, pos_p = 0; j < c; ++j, pos += res_ld, pos_p += m_ld)
                mrd::reset_helper(ptr_res[pos],ptr[pos_p]);

            ptr_ind[i] = row_vec_ptr[i]->get_pos();
        };

        iret    = matcl::Matrix(ind,false);
        x       = matcl::Matrix(res,true);
        return;
    };

    static void eval_col2(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& m)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();
        integer_dense ind(ti::ti_empty(),1,c);

        Integer* ptr_ind = ind.ptr();

        if (r == 0 || c == 0) 
        {            
            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = i+1;

            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        std::vector<row_info<value_type>>  row_vec(c);
        std::vector<row_info<value_type>*> row_vec_ptr(c);

        for (Integer i = 0, pos = 0; i < c; ++i, pos += m_ld)
        {
            row_vec[i] = row_info<value_type>(ptr_m+pos,1,r,i+1);
            row_vec_ptr[i] = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < c; ++i)
        {
            const value_type* ptr    = row_vec_ptr[i]->get_ptr();
            for (Integer j = 0, pos_p = 0; j < r; ++j, ++pos_p)
                mrd::reset_helper(ptr_res[j],ptr[pos_p]);

            ptr_res     += res_ld;
            ptr_ind[i]  = row_vec_ptr[i]->get_pos();
        };

        iret    = matcl::Matrix(ind,false);
        x       = matcl::Matrix(res,true);
        return;
    };

    static void eval_dim(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& m,const integer_dense& dims)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();

        check_dims(dims,c);

        if (r == 0 || c == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        if (dims.size() == 0)
        {
            ret = matcl::Matrix(m,false);
            return;
        };

        integer_dense dims2 = dims.make_explicit();
        const value_type* ptr_m = m.ptr();        
        Integer m_ld            = m.ld();

        std::vector<row_info_dim<value_type>>  row_vec(r);
        std::vector<row_info_dim<value_type>*> row_vec_ptr(r);
        
        for (Integer i = 0; i < r; ++i)
        {
            row_vec[i]      = row_info_dim<value_type>(ptr_m+i,m_ld,&dims2);
            row_vec_ptr[i]  = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_dim_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < r; ++i)
        {
            const value_type* ptr = row_vec_ptr[i]->get_ptr();

            for (Integer j = 0, pos = i, pos_p = 0; j < c; ++j, pos += res_ld, pos_p += m_ld)
                mrd::reset_helper(ptr_res[pos],ptr[pos_p]);
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_dim_col(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& m,const integer_dense& dims)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();
        check_dims_col(dims,r);

        if (r == 0 || c == 0) 
        {
            ret = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };
        if (dims.size() == 0)
        {
            ret = matcl::Matrix(m,false);
            return;
        };

        integer_dense dims2 = dims.make_explicit();

        std::vector<row_info_dim<value_type>>  row_vec(c);
        std::vector<row_info_dim<value_type>*> row_vec_ptr(c);
        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer i = 0, pos = 0; i < c; ++i, pos += m_ld)
        {
            row_vec[i]      = row_info_dim<value_type>(ptr_m+pos,1,&dims2);
            row_vec_ptr[i]  = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_dim_comparer<value_type>);

        ret_type res(ret_ti,r,c);

        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < c; ++i)
        {
            const value_type* ptr = row_vec_ptr[i]->get_ptr();

            for (Integer j = 0; j < r; ++j)
                mrd::reset_helper(ptr_res[j],ptr[j]);

            ptr_res += res_ld;
        };

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_dim2(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& m,
                          const integer_dense& dims)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();

        check_dims(dims,c);

        integer_dense ind(ti::ti_empty(),r,1);
        Integer* ptr_ind = ind.ptr();
        const value_type* ptr_m = m.ptr();

        if (r == 0 || c == 0) 
        {
            for (Integer i = 0; i < r; ++i)
                ptr_ind[i] = i+1;

            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        if (dims.size() == 0)
        {
            for (Integer i = 0; i < r; ++i)
                ptr_ind[i] = i+1;

            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(m,false);
            return;
        };

        integer_dense dims2 = dims.make_explicit();
        Integer m_ld        = m.ld();

        std::vector<row_info_dim<value_type>>  row_vec(r);
        std::vector<row_info_dim<value_type>*> row_vec_ptr(r);
        
        for (Integer i = 0; i < r; ++i)
        {
            row_vec[i] = row_info_dim<value_type>(ptr_m+i,m_ld,&dims2,i+1);
            row_vec_ptr[i] = &row_vec[i];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_dim_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer i = 0; i < r; ++i)
        {
            const value_type* ptr = row_vec_ptr[i]->get_ptr();
            for (Integer j = 0, pos = i, pos_p = 0; j < c; ++j, pos += res_ld, pos_p += m_ld)
                mrd::reset_helper(ptr_res[pos],ptr[pos_p]);

            ptr_ind[i] = row_vec_ptr[i]->get_pos();
        };

        iret    = matcl::Matrix(ind,false);
        x       = matcl::Matrix(res,true);
        return;
    };

    static void eval_dim_col2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& m,
                              const integer_dense& dims)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();
        check_dims_col(dims,r);

        integer_dense ind(ti::ti_empty(),1,c);
        Integer* ptr_ind = ind.ptr();

        if (r == 0 || c == 0) 
        {
            for (Integer k = 0; k < c; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(ret_type(ret_ti,r, c),false);
            return;
        };

        if (dims.size() == 0)
        {
            for (Integer k = 0; k < c; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(m,false);
            return;
        };

        integer_dense dims2 = dims.make_explicit();

        std::vector<row_info_dim<value_type>>  row_vec(c);
        std::vector<row_info_dim<value_type>*> row_vec_ptr(c);

        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer k = 0, pos = 0; k < c; ++k, pos += m_ld)
        {
            row_vec[k]      = row_info_dim<value_type>(ptr_m+pos,1,&dims2,k+1);
            row_vec_ptr[k]  = &row_vec[k];
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_dim_comparer<value_type>);

        ret_type res(ret_ti,r,c);
        
        value_type_ret* ptr_res = res.ptr();
        Integer res_ld          = res.ld();

        for (Integer k = 0; k < c; ++k)
        {
            const value_type* ptr = row_vec_ptr[k]->get_ptr();

            for (Integer j = 0; j < r; ++j)
                mrd::reset_helper(ptr_res[j],ptr[j]);
            
            ptr_ind[k]  = row_vec_ptr[k]->get_pos();
            ptr_res     += res_ld;
        };

        i = matcl::Matrix(ind,false);
        x = matcl::Matrix(res,false);
        return;
    };
};

static integer_dense abs_impl(const integer_dense& mat)
{
    integer_dense ret   = mat.copy();
    Integer r           = ret.rows();
    Integer c           = ret.cols();
    Integer ld          = ret.ld();
    Integer* ptr        = ret.ptr();

    for (Integer i = 0; i < c; ++i)
    {
        for (Integer j = 0; j < r; ++j)
            ptr[j]      = abs_helper<Integer>::eval(ptr[j]);

        ptr             += ld;
    };

    return ret;
};

template<class ret_type, class in_type>
struct sortrows_impl<ret_type,in_type,struct_sparse>
{
    using ret_2     = std::pair<ret_type,integer_dense>;
    using ret_ti_t  = typename ti::get_ti_type<ret_type>::type;

    static void eval(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A)
    {
        if (A.nnz() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        in_type At(A.get_type());
        mrd::manip_trans_helper<in_type>::eval_trans(At,A);

        matcl::Matrix ret2;
        eval_impl(ret2, ret_ti, At);

        ret = matcl::trans(ret2);
        return;
    };

    static void eval_col(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A)
    {
        if (A.nnz() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        eval_impl(ret, ret_ti, A);
        return;
    };

    static void eval2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        if (A.nnz() == 0)
        {
            Integer r = A.rows();
            
            integer_dense ind(ti::ti_empty(),r,1);
            Integer* ptr_ind = ind.ptr();
            
            for (Integer k = 0; k < r; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(A,false);
            return;
        };

        in_type At(A.get_type());
        mrd::manip_trans_helper<in_type>::eval_trans(At,A);

        matcl::Matrix tmp_i, tmp_x;
        eval_impl2(tmp_i, tmp_x, ret_ti,At);

        i = matcl::trans(tmp_i);
        x = matcl::trans(tmp_x);
        return;
    };

    static void eval_col2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        if (A.nnz() == 0)
        {
            Integer c = A.cols();
            
            integer_dense ind(ti::ti_empty(),1,c);
            Integer* ptr_ind = ind.ptr();
            
            for (Integer k = 0; k < c; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(A,false);
            
            return;
        };

        eval_impl2(i,x,ret_ti,A);        
        return;
    };

    static void eval_dim(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A,const integer_dense& dims)
    {
        Integer c = A.cols(), nz = A.nnz();
        check_dims(dims,c);

        if (nz == 0 || dims.size() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        in_type A_selected = select_cols(A,dims);
        in_type At(A.get_type());
        in_type At_selected(A_selected.get_type());

        mrd::manip_trans_helper<in_type>::eval_trans(At,A);
        mrd::manip_trans_helper<in_type>::eval_trans(At_selected, A_selected);

        matcl::Matrix tmp_x;
        eval_impl(tmp_x,ret_ti,At,At_selected,dims);

        ret = matcl::trans(tmp_x);
        return;
    };

    static void eval_dim_col(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A,const integer_dense& dims)
    {
        Integer r = A.rows(), nz = A.nnz();
        check_dims_col(dims,r);

        if (nz == 0 || dims.size() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        in_type A_selected =select_rows(A,dims);
        eval_impl(ret, ret_ti,A,A_selected, dims);        
        return;
    };

    static void eval_dim2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,
                          const integer_dense& dims)
    {
        Integer c = A.cols(), nz = A.nnz();
        check_dims(dims,c);

        if (nz == 0 || dims.size() == 0)
        {
            Integer r = A.rows();
            
            integer_dense ind(ti::ti_empty(),r,1);
            Integer* ptr_ind = ind.ptr();
            
            for (Integer k = 0; k < r; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(A,false);
            
            return;
        };

        in_type A_selected = select_cols(A,dims);    
        in_type At(A.get_type());
        in_type At_selected(A.get_type());

        mrd::manip_trans_helper<in_type>::eval_trans(At,A);
        mrd::manip_trans_helper<in_type>::eval_trans(At_selected,A_selected);

        matcl::Matrix tmp_i, tmp_x;
        eval_impl2(tmp_i, tmp_x, ret_ti,At,At_selected,dims);

        i = matcl::trans(tmp_i);
        x = matcl::trans(tmp_x);

        return;
    };

    static void eval_dim_col2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,
                              const integer_dense& dims)
    {
        Integer r = A.rows(), c = A.cols();
        check_dims_col(dims,r);

        if (A.nnz() == 0 || dims.size() == 0)
        {
            integer_dense ind(ti::ti_empty(),1,c);
            Integer* ptr_ind = ind.ptr();
            
            for (Integer k = 0; k < c; ++k)
                ptr_ind[k] = k+1;

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(A,false);
            return;
        };
        
        in_type A_selected = select_rows(A,dims);
        eval_impl2(i,x,ret_ti,A,A_selected,dims);
        return;
    };

    static void eval_impl(matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        matcl::Matrix i;
        return eval_impl2(i,x,ret_ti,A);
    };

    static void eval_impl(matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,const in_type& Asel, 
                          const integer_dense& dims)
    {
        matcl::Matrix i;
        return eval_impl2(i,x,ret_ti,A,Asel,dims);
    };

    static void eval_impl2(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        using value_type = typename in_type::value_type;

        Integer r = A.rows(), c = A.cols();

        ret_type res(ret_ti,r,c,A.nnz());
        
        integer_dense ind(ti::ti_empty(),1,c);
        Integer* ptr_ind = ind.ptr();

        if (A.nnz() == 0)
        {
            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = i+1;

            iret    = matcl::Matrix(ind,false);
            x       = matcl::Matrix(res,false);

            return;
        };

        const sparse_ccs<value_type>& Ad = A.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const Integer* Ad_r         = Ad.ptr_r();
        const value_type* Ad_x      = Ad.ptr_x();

        sparse_ccs<value_type>& d        = res.rep();
        Integer* d_c                = d.ptr_c();
        Integer* d_r                = d.ptr_r();
        value_type* d_x             = d.ptr_x();

        std::vector<row_info_sparse<value_type>>  row_vec;
        std::vector<row_info_sparse<value_type>*> row_vec_ptr;

        row_vec.reserve(c);
        row_vec_ptr.reserve(c);

        for (Integer i = 0; i < c; ++i)
        {
            row_vec.push_back(row_info_sparse<value_type>(ret_ti,Ad_r,Ad_x,Ad_c[i],Ad_c[i+1],i));
            row_vec_ptr.push_back(&row_vec[i]);
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_sparse_comparer<value_type>);        

        Integer nz = 0;
        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]      = nz;
            Integer col = row_vec_ptr[j]->get_pos();

            for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
            {
                d_r[nz] = Ad_r[k];
                mrd::reset_helper(d_x[nz],Ad_x[k]);
                ++nz;
            };

            ptr_ind[j]  = col + 1;
        };

        d_c[c]          = nz;

        iret    = matcl::Matrix(ind,true);
        x       = matcl::Matrix(res,true);
        return;
    };

    static void eval_impl2(matcl::Matrix& iret, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,
                           const in_type& A_sub, const integer_dense& dims)
    {
        using value_type = typename in_type::value_type;

        Integer r = A.rows(), c = A.cols();

        const sparse_ccs<value_type>& Ad = A.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const Integer* Ad_r         = Ad.ptr_r();
        const value_type* Ad_x      = Ad.ptr_x();

        const sparse_ccs<value_type>& Ads= A_sub.rep();
        const Integer* Ads_c        = Ads.ptr_c();
        const Integer* Ads_r        = Ads.ptr_r();
        const value_type* Ads_x     = Ads.ptr_x();

        ret_type res(ret_ti,r,c,A.nnz());
        integer_dense ind(ti::ti_empty(),1,c);
        Integer* ptr_ind = ind.ptr();

        sparse_ccs<value_type>& d        = res.rep();
        Integer* d_c                = d.ptr_c();
        Integer* d_r                = d.ptr_r();
        value_type* d_x             = d.ptr_x();

        std::vector<row_info_sparse_dims<value_type>>  row_vec;
        std::vector<row_info_sparse_dims<value_type>*> row_vec_ptr;

        row_vec.reserve(c);
        row_vec_ptr.reserve(c);

        integer_dense dims2 = dims.make_explicit();

        for (Integer i = 0; i < c; ++i)
        {
            row_vec.push_back(row_info_sparse_dims<value_type>(ret_ti,Ads_r,Ads_x,Ads_c[i],Ads_c[i+1],i,&dims2));
            row_vec_ptr.push_back(&row_vec[i]);
        };

        std::stable_sort(row_vec_ptr.begin(), row_vec_ptr.end(), row_vec_sparse_dims_comparer<value_type>);        

        Integer nz = 0;
        for (Integer j = 0; j < c; ++j)
        {
            d_c[j]      = nz;
            Integer col = row_vec_ptr[j]->get_pos();

            for (Integer k = Ad_c[col]; k < Ad_c[col+1]; ++k)
            {
                d_r[nz] = Ad_r[k];
                mrd::reset_helper(d_x[nz],Ad_x[k]);
                ++nz;
            };

            ptr_ind[j]  = col + 1;
        };

        d_c[c]  = nz;

        iret    = matcl::Matrix(ind,false);
        x       = matcl::Matrix(res,false);

        return;
    };    

    static in_type select_cols(const in_type& A, const integer_dense& cols)
    {
        Integer r = A.rows();

        matcl::details::colon_info ci;

        ci.r_flag   = 1;        
        ci.r_start  = 1;
        ci.r_step   = 1;
        ci.r_end    = r;
        ci.r_size   = r;
        ci.c_flag   = 0;

        ci.set_ci(abs_impl(cols));

        return matcl::algorithm::get_submatrix(A,ci);
    };

    static in_type select_rows(const in_type& A, const integer_dense& ri)
    {
        Integer c = A.cols();

        matcl::details::colon_info ci;

        ci.r_flag   = 0;
        ci.set_ri(abs_impl(ri));
        ci.c_start  = 1;
        ci.c_step   = 1;
        ci.c_end    = c;
        ci.c_size   = c;
        ci.c_flag   = 1;

        return matcl::algorithm::get_submatrix(A,ci);
    };
};

template<class ret_type, class in_type>
struct sort_impl<ret_type,in_type,struct_banded>
{
    using ret_ti_t =  typename ti::get_ti_type<ret_type>::type;

    static void eval(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A, Integer dim, bool asceding)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);

        if (dim == 1)
        {
            if (asceding == true)
                sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_inplace<true>(tmp);
            else
                sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_inplace<false>(tmp);

            ret = matcl::Matrix(tmp,true);
            return;
        };
        return sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval(ret,ret_ti,tmp,dim,asceding);
    };

    static void eval(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A, Integer dim, 
                     bool asceding)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);

        if (dim == 1)
        {
            integer_dense ind = integer_dense(ti::ti_empty());

            if (asceding == true)
                sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_inplace<true>(tmp,ind,false);
            else
                sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_inplace<false>(tmp,ind,false);

            i = matcl::Matrix(ind,false);
            x = matcl::Matrix(tmp,false);
            return;
        };

        return sort_impl<ret_type,sparse_matrix_type,struct_sparse>::eval(i,x,ret_ti,tmp,dim,asceding);
    };
};

template<class ret_type, class in_type>
struct sortrows_impl<ret_type,in_type,struct_banded>
{
    using ret_2     = std::pair<ret_type,integer_dense>;
    using ret_ti_t  = typename ti::get_ti_type<ret_type>::type;

    static void eval(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval(ret, ret_ti,tmp);
    };

    static void eval2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval2(i,x,ret_ti,tmp);
    };

    static void eval_dim(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A,const integer_dense& dims)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_dim(ret, ret_ti,tmp,dims);
    };

    static void eval_dim2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,
                          const integer_dense& dims)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_dim2(i,x,ret_ti,tmp,dims);
    };

    static void eval_col(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>
                    ::eval_impl(ret, ret_ti,tmp);
    };

    static void eval_col2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_impl2(i,x,ret_ti,tmp);
    };

    static void eval_dim_col(matcl::Matrix& ret, ret_ti_t ret_ti,const in_type& A,const integer_dense& dims)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        Integer r = A.rows();
        check_dims_col(dims,r);

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_dim_col(ret,ret_ti,tmp,dims);
    };

    static void eval_dim_col2(matcl::Matrix& i, matcl::Matrix& x, ret_ti_t ret_ti,const in_type& A,
                              const integer_dense& dims)
    {
        using val_type              = typename ret_type::value_type;
        using sparse_matrix_type    = Matrix<val_type,struct_sparse>;

        Integer r = A.rows();
        check_dims_col(dims,r);

        sparse_matrix_type tmp = converter<sparse_matrix_type,in_type>::eval(A,ret_ti);
        return sortrows_impl<ret_type,sparse_matrix_type,struct_sparse>::eval_dim_col2(i,x,ret_ti,tmp,dims);
    };
};

template<class MP>    
void sort_helper<MP>::eval_sort(matcl::Matrix& ret, const MP& m,Integer dim, bool asceding)
{
    using str_type = typename MP::struct_type;
    return sort_impl<ret_type_sort,MP,str_type>::eval(ret, ti::get_ti(m),m,dim,asceding);
};

template<class MP>
void sort_helper<MP>::eval_sort_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m,Integer dim, 
                                  bool asceding)
{
    using str_type = typename MP::struct_type;
    return sort_impl<ret_type_sort,MP,str_type>::eval(i,x,ti::get_ti(m),m,dim,asceding);
};

template<class MP>
void sort_helper<MP>::eval_sortrows(matcl::Matrix& ret, const MP& m)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval(ret, ti::get_ti(m),m);
};

template<class MP>
void sort_helper<MP>::eval_sortrows_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval2(i,x,ti::get_ti(m),m);
};

template<class MP>
void sort_helper<MP>::eval_sortrows(matcl::Matrix& ret, const MP& m, const integer_dense& dims)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_dim(ret, ti::get_ti(m),m,dims);
};

template<class MP>
void sort_helper<MP>::eval_sortrows_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m, 
                                      const integer_dense& dims)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_dim2(i,x,ti::get_ti(m),m,dims);
};

template<class MP>
void sort_helper<MP>::eval_sortcols(matcl::Matrix& ret, const MP& m)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_col(ret, ti::get_ti(m),m);
};

template<class MP>
void sort_helper<MP>::eval_sortcols_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_col2(i,x,ti::get_ti(m),m);
};

template<class MP>
void sort_helper<MP>::eval_sortcols(matcl::Matrix& i, const MP& m, const integer_dense& dims)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_dim_col(i, ti::get_ti(m),m,dims);
};

template<class MP>
void sort_helper<MP>::eval_sortcols_2(matcl::Matrix& i, matcl::Matrix& x, const MP& m, 
                                      const integer_dense& dims)
{
    using str_type = typename MP::struct_type;
    return sortrows_impl<ret_type_sort,MP,str_type>::eval_dim_col2(i,x,ti::get_ti(m),m,dims);
};

template<class M, class struct_type>
struct is_sortedcols_helper {};

template<class M, class struct_type>
struct is_sortedrows_helper {};

template<class in_type>
struct is_sortedcols_helper<in_type,struct_dense> 
{
    static bool eval(const in_type& m)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();

        if (r == 0 || c == 0) 
            return true;

        Integer pos_c           = m.ld();
        Integer pos_cl          = 0;
        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer j = 1; j < c; ++j)
        {
            Integer pos     = pos_c;
            Integer pos_l   = pos_cl;

            for (Integer i = 0; i < r; ++i)
            {
                if (sort_greater_nan<value_type>::eval(ptr_m[pos],ptr_m[pos_l]))
                    break;
                else if (sort_less_nan<value_type>::eval(ptr_m[pos],ptr_m[pos_l]))
                    return false;

                ++pos;
                ++pos_l;
            };

            pos_c   += m_ld;
            pos_cl  += m_ld;

        };
        return true;
    };
};

template<class in_type>
struct is_sortedrows_helper<in_type,struct_dense> 
{
    static bool eval(const in_type& m)
    {
        using value_type = typename in_type::value_type;

        Integer r = m.rows(), c = m.cols();

        if (r == 0 || c == 0) 
            return true;

        const value_type* ptr_m = m.ptr();
        Integer m_ld            = m.ld();

        for (Integer i = 1; i < r; ++i)
        {
            Integer pos     = i;
            Integer pos_l   = i - 1;

            for (Integer j = 0; j < c; ++j)
            {
                if (sort_greater_nan<value_type>::eval(ptr_m[pos],ptr_m[pos_l]))
                    break;
                else if (sort_less_nan<value_type>::eval(ptr_m[pos],ptr_m[pos_l]))
                    return false;

                pos     += m_ld;
                pos_l   += m_ld;
            };
        };
        return true;
    };
};

template<class in_type>
struct is_sortedcols_helper<in_type,struct_sparse> 
{
    static bool eval(const in_type& A)
    {
        using value_type = typename in_type::value_type;

        const value_type Z = md::default_value<value_type>(ti::get_ti(A));

        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0 || A.nnz() == 0)
            return true;

        const sparse_ccs<value_type>& Ad = A.rep();
        const Integer* Ad_c         = Ad.ptr_c();
        const Integer* Ad_r         = Ad.ptr_r();
        const value_type* Ad_x      = Ad.ptr_x();

        for (Integer i = 1; i < c; ++i)
        {
            Integer k1 = Ad_c[i], k2 = Ad_c[i-1];
            Integer l1 = Ad_c[i+1], l2 = k1;

            while (k1 < l1 && k2 < l2)
            {
                if (Ad_r[k1] < Ad_r[k2])
                {
                    if (sort_less_nan<value_type>::eval(Ad_x[k1],Z))
                        return false;

                    if (sort_greater_nan<value_type>::eval(Ad_x[k1],Z))
                        goto label_exis_row;

                    ++k1;
                }
                else if (Ad_r[k1] > Ad_r[k2])
                {
                    if (sort_less_nan<value_type>::eval(Z, Ad_x[k2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(Z,Ad_x[k2]))
                        goto label_exis_row;

                    ++k2; 
                }
                else
                { 
                    if (sort_less_nan<value_type>::eval(Ad_x[k1], Ad_x[k2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(Ad_x[k1], Ad_x[k2]))
                        goto label_exis_row;

                    ++k1; 
                    ++k2;
                };
            };
            while (k1 < l1)
            {
                if (sort_less_nan<value_type>::eval(Ad_x[k1], Z))
                    return false;

                if (sort_greater_nan<value_type>::eval(Ad_x[k1], Z))
                    goto label_exis_row;

                ++k1;
            };

            while (k2 < l2)
            {                
                if (sort_less_nan<value_type>::eval(Z, Ad_x[k2]))
                    return false;

                if (sort_greater_nan<value_type>::eval(Z, Ad_x[k2]))
                    goto label_exis_row;

                ++k2; 
            };

            label_exis_row:
                (void)0;
        };

        return true;
    };
};

template<class in_type>
struct is_sortedrows_helper<in_type,struct_sparse> 
{
    static bool eval(const in_type& A)
    {
        in_type At(A.get_type());
        mrd::manip_trans_helper<in_type>::eval_trans(At, A);
        return is_sortedcols_helper<in_type,struct_sparse>::eval(At);
    };
};

template<class in_type>
struct is_sortedcols_helper<in_type,struct_banded> 
{
    static bool eval(const in_type& A)
    {
        using value_type    = typename in_type::value_type;
        const value_type Z  = md::default_value<value_type>(ti::get_ti(A));

        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0) 
            return true;

        Integer pos_c   = A.ld();
        Integer A_ld    = A.ld();
        Integer pos_cl  = 0;
        const value_type* ptr_A = A.rep_ptr();

        for (Integer i = 1; i < c; ++i)
        {
            Integer k1 = A.first_row(i),    k2 = A.first_row(i-1);
            Integer l1 = A.last_row(i)+1,   l2 = A.last_row(i-1)+1;

            Integer pos_1 = A.first_elem_pos(i) + pos_c;
            Integer pos_2 = A.first_elem_pos(i-1) + pos_cl;

            while (k1 < l1 && k2 < l2)
            {
                if (k1 < k2)
                {
                    if (sort_less_nan<value_type>::eval(ptr_A[pos_1], Z))
                        return false;

                    if (sort_greater_nan<value_type>::eval(ptr_A[pos_1], Z))
                        goto label_exit_col;

                    ++k1;
                    ++pos_1;
                }
                else if (k1 > k2)
                {
                    if (sort_less_nan<value_type>::eval(Z, ptr_A[pos_2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(Z, ptr_A[pos_2]))
                        goto label_exit_col;

                    ++k2; 
                    ++pos_2;
                }
                else
                { 
                    if (sort_less_nan<value_type>::eval(ptr_A[pos_1], ptr_A[pos_2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(ptr_A[pos_1], ptr_A[pos_2]))
                        goto label_exit_col;

                    ++k1; 
                    ++k2;
                    ++pos_1;
                    ++pos_2;
                };
            };

            while (k1 < l1)
            {
                if (sort_less_nan<value_type>::eval(ptr_A[pos_1],Z))
                    return false;

                if (sort_greater_nan<value_type>::eval(ptr_A[pos_1], Z))
                    goto label_exit_col;

                ++k1;
                ++pos_1;
            };

            while (k2 < l2)
            {                
                if (sort_less_nan<value_type>::eval(Z,ptr_A[pos_2]))
                    return false;

                if (sort_greater_nan<value_type>::eval(Z,ptr_A[pos_2]))
                    goto label_exit_col;

                ++k2; 
                ++pos_2;
            };

            label_exit_col:

            pos_c   += A_ld;
            pos_cl  += A_ld;
        };

        return true;
    };
};

template<class in_type>
struct is_sortedrows_helper<in_type,struct_banded> 
{
    static bool eval(const in_type& A)
    {
        using value_type    = typename in_type::value_type;
        const value_type Z  = md::default_value<value_type>(ti::get_ti(A));

        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0) 
            return true;

        Integer dpos = A.ld() - 1;

        const value_type* ptr_A = A.rep_ptr();

        for (Integer i = 1; i < r; ++i)
        {
            Integer k1 = A.first_col(i),    k2 = A.first_col(i-1);
            Integer l1 = A.last_col(i)+1,   l2 = A.last_col(i-1)+1;

            Integer pos_1 = A.first_elem_pos_row(i);
            Integer pos_2 = A.first_elem_pos_row(i-1);

            while (k1 < l1 && k2 < l2)
            {
                if (k1 < k2)
                {
                    if (sort_less_nan<value_type>::eval(ptr_A[pos_1],Z))
                        return false;

                    if (sort_greater_nan<value_type>::eval(ptr_A[pos_1],Z))
                        goto label_exit_row;

                    ++k1;
                    pos_1+= dpos;
                }
                else if (k1 > k2)
                {
                    if (sort_less_nan<value_type>::eval(Z,ptr_A[pos_2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(Z,ptr_A[pos_2]))
                        goto label_exit_row;

                    ++k2; 
                    pos_2+= dpos;
                }
                else
                { 
                    if (sort_less_nan<value_type>::eval(ptr_A[pos_1],ptr_A[pos_2]))
                        return false;

                    if (sort_greater_nan<value_type>::eval(ptr_A[pos_1],ptr_A[pos_2]))
                        goto label_exit_row;

                    ++k1; 
                    ++k2;
                    pos_1+= dpos;
                    pos_2+= dpos;
                };
            };

            while (k1 < l1)
            {
                if (sort_less_nan<value_type>::eval(ptr_A[pos_1],Z))
                    return false;

                if (sort_greater_nan<value_type>::eval(ptr_A[pos_1],Z))
                    goto label_exit_row;

                ++k1;
                pos_1+= dpos;
            };

            while (k2 < l2)
            {                
                if (sort_less_nan<value_type>::eval(Z,ptr_A[pos_2]))
                    return false;

                if (sort_greater_nan<value_type>::eval(Z,ptr_A[pos_2]))
                    goto label_exit_row;

                ++k2; 
                pos_2+= dpos;
            };

            label_exit_row:
                (void)0;
        };

        return true;
    };
};

template<class in_type, bool asceding>
struct issorted_accumulator
{
    private:
        using in_ti             = ti::ti_type<in_type>;
        using ret_ti            = ti::ti_type<Integer>;

    private:
        enum state_t { false_s, true_s,  uninitialized };

        state_t                 state;
        in_type                 current_value;
        const in_type           Z;
        Integer                 pos;
        std::vector<state_t>    st_array;
        std::vector<Integer>    st_pos;
        md::workspace2<in_type> st_current_value;
        ti::ti_type<in_type>    m_ti;

        Integer                 size;        

        issorted_accumulator&   operator=(const issorted_accumulator&) = delete;

    public:
        issorted_accumulator(ti::ti_type<in_type> ti)    
            :current_value(md::default_value<in_type>(ti))
            ,Z(md::default_value<in_type>(ti))
            ,m_ti(ti),st_current_value(ti)
        {};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer r, Integer , Integer )    
        { 
            size = r;
        };

        void reset()                        
        { 
            state = uninitialized;
            pos   = 0;
        };

        void reset_array(Integer s)            
        { 
            st_array.clear(); 
            st_pos.clear(); 

            st_array.resize(s,uninitialized);    
            st_pos.resize(s,0);    
            st_current_value.resize(s);
        };

        void current_vector(Integer )
        {};

        bool add(Integer v, const in_type& val)
        { 
            if (state == uninitialized)
            {                
                pos = v + 1;
                mrd::reset_helper(current_value,val);

                if ( v > 0 && md::compare_nan<in_type,asceding>::eval(val,Z))
                    state = false_s;
                else
                    state = true_s;

                return false;
            };

            if (state == true_s)
            {
                if ( v == pos )
                {
                    if ( md::compare_nan<in_type,asceding>::eval(val,current_value))
                        state = false_s;
                }
                else
                {
                    if ( md::compare_nan<in_type,asceding>::eval(Z,current_value))
                        state = false_s;

                    if ( md::compare_nan<in_type,asceding>::eval(val,Z))
                        state = false_s;
                };

                pos = v + 1;
                mrd::reset_helper(current_value,val);
            };

            return !(state==true_s);
        };

        void add(Integer p,Integer v, const in_type& val)        
        { 
            if (st_array[p] == uninitialized)
            {                
                st_pos[p] = v + 1;
                mrd::reset_helper(st_current_value[p],val);

                if ( v > 0 && md::compare_nan<in_type,asceding>::eval(val,Z))
                    st_array[p] = false_s;
                else
                    st_array[p] = true_s;

                return;
            };

            if (st_array[p] == true_s)
            {
                if ( v == st_pos[p] )
                {
                    if ( md::compare_nan<in_type,asceding>::eval(val,st_current_value[p]))
                        st_array[p] = false_s;
                }
                else
                {
                    if ( md::compare_nan<in_type,asceding>::eval(Z,st_current_value[p]))
                        st_array[p] = false_s;

                    if ( md::compare_nan<in_type,asceding>::eval(val,Z))
                        st_array[p] = false_s;
                };

                st_pos[p] = v + 1;
                mrd::reset_helper(st_current_value[p],val);
            };
        };

        bool add_zero(Integer )
        { 
            return false;
        };

        void add_zero(Integer ,Integer )
        {};

        bool value()                    
        { 
            if (state == uninitialized)
                return true;

            if (pos < size)
            {
                if ( md::compare_nan<in_type,asceding>::eval(Z,current_value))
                    return false;
            };

            return state == true_s;
        };

        bool value(Integer p)            
        { 
            if (st_array[p] == uninitialized)
                return true;

            if (st_pos[p] < size)
            {
                if ( md::compare_nan<in_type,asceding>::eval(Z,st_current_value[p]))
                    return false;
            };

            return st_array[p]==true_s;
        };
};

template<class MP>
struct is_sorted_helper
{
    static void eval(matcl::Matrix& ret, const MP& A, Integer dim, bool asceding)
    {
        using in_type = typename MP::value_type;

        if (asceding == true)
        {
            using accum = issorted_accumulator<in_type,true>;

            accum ac(ti::get_ti(A));
            eval_vec_functor2<integer_dense,MP, accum>::eval(ret,A,dim,ac);
            return;
        }
        else
        {
            using accum = issorted_accumulator<in_type,false>;

            accum ac(ti::get_ti(A));
            eval_vec_functor2<integer_dense,MP, accum>::eval(ret,A,dim,ac);
            return;
        }
    };
};

template<class MP>
void sort_helper<MP>::eval_issorted(matcl::Matrix& ret, const MP& m, Integer dim, bool asceding)
{
    return is_sorted_helper<MP>::eval(ret, m,dim,asceding);
};

template<class MP>
bool sort_helper<MP>::eval_issorted_rows(const MP& m)
{
    return is_sortedrows_helper<MP,typename MP::struct_type>::eval(m);
};        

template<class MP>
bool sort_helper<MP>::eval_issorted_cols(const MP& m)
{
    return is_sortedcols_helper<MP,typename MP::struct_type>::eval(m);
};        

}}}

MACRO_INSTANTIATE_G_1(matcl::raw::details::sort_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::sort_helper)
