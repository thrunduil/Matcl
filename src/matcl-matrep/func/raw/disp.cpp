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

#include "matcl-matrep/func/raw/disp.h"
#include <iomanip>
#include "matcl-matrep/func/raw/io.h"
#include "matcl-internals/func/converter.h"
#include "matcl-core/details/integer.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/IO/disp_stream_impl.h"
#include "matcl-core/details/IO/printer.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/objects/omatrix.h"
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-scalar/IO/scalar_io.h"

#pragma warning(push) 
#pragma warning(disable:4702) //unreachable code

namespace matcl { namespace raw
{
 
namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

namespace details
{
    template<class matrix_type, class V>
    struct matrix_provider_dense : public matrix_provider_base<V>
    {
        const matrix_type&  m_mat;

        matrix_provider_dense(const matrix_type& mat)
            :m_mat(mat)
        {};
        virtual ~matrix_provider_dense(){};

        virtual Integer rows() const override
        {
            return m_mat.rows();
        };

        virtual Integer cols() const override
        {
            return m_mat.cols();
        };

        virtual V get_value(const disp_stream* user, Integer width, 
                            align_type at, Integer r, Integer c, bool& is_zero) const override
        {
            (void)user;
            (void)width;
            (void)at;
            V val   = m_mat.ptr()[r + c * m_mat.ld()];
            is_zero = mrd::is_zero(val);
            return val;
        };

        virtual bool is_symher() const override
        {
            struct_flag sf = m_mat.get_struct();

            return sf.has_sym_flag() || sf.has_her_flag();
        };
    };

    template<class matrix_type, class V>
    struct matrix_provider_band : public matrix_provider_base<V>
    {
        const matrix_type&  m_mat;
        V                   m_zero;

        matrix_provider_band(const matrix_type& mat)
            :m_mat(mat), m_zero(md::default_value<V>(ti::get_ti(mat)))
        {};

        virtual ~matrix_provider_band(){};

        virtual Integer rows() const override
        {
            return m_mat.rows();
        };

        virtual Integer cols() const override
        {
            return m_mat.cols();
        };

        virtual V get_value(const disp_stream* user, Integer width, align_type at, 
                           Integer r, Integer c, bool& is_zero) const override
        {
            (void)user;
            (void)width;
            (void)at;

            Integer fr  = m_mat.first_row(c);
            Integer lr  = m_mat.last_row(c);
            Integer row = r;

            if (row < fr || row > lr)
            {
                is_zero = true;
                return m_zero;
            }

            V ret   = m_mat(r+1,c+1);
            is_zero = mrd::is_zero(ret);
            return ret;
        };

        virtual bool is_symher() const override
        {
            struct_flag sf = m_mat.get_struct();

            return sf.has_sym_flag() || sf.has_her_flag();
        };

    };

    template<class matrix_type, class V>
    struct matrix_provider_sparse : public matrix_provider_base<V>
    {
        const matrix_type&  m_mat;
        V                   m_zero;

        matrix_provider_sparse(const matrix_type& mat)
            :m_mat(mat), m_zero(md::default_value<V>(ti::get_ti(mat)))
        {};

        virtual ~matrix_provider_sparse(){};

        virtual Integer rows() const override
        {
            return m_mat.rows();
        };

        virtual Integer cols() const override
        {
            return m_mat.cols();
        };

        virtual V get_value(const disp_stream* user, Integer width, align_type at, 
                            Integer r, Integer c, bool& is_zero) const override
        {
            (void)user;
            (void)width;
            (void)at;

            V ret   = m_mat(r+1,c+1);
            is_zero = mrd::is_zero(ret);
            return ret;
        };

        virtual bool is_symher() const override
        {
            struct_flag sf = m_mat.get_struct();

            return sf.has_sym_flag() || sf.has_her_flag();
        };
    };

    template<class V>
    struct disp_matrix<V,struct_banded>
    {
        using matrix_type =  raw::Matrix<V,struct_banded>;

        static void eval(md::disp_stream_impl& os, const disp_stream* user, const matrix_type& mat)
        {
            os.do_init_display(user,matrix_traits::value_code<V>::value);
            os.do_start_display(user);	
            
            os.do_displaying_banded_matrix2(user,mat.rows(),mat.cols(), mat.first_diag(), mat.last_diag(), 
                        matrix_traits::value_code<V>::value, mat.get_struct().to_string());

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            matrix_provider_band<matrix_type,V> mdp(mat);
            return disp_matrix<V,struct_dense>::eval_matrix_body(os,user,mdp);
        };
    };

    template<class V>
    struct disp_matrix<V,struct_sparse>
    {
        using matrix_type = raw::Matrix<V,struct_sparse>;

        struct measure_params
        {
            Integer                 vcvw;
            std::vector<Integer>    min_value_width;
            Integer                 min_width;
            Integer                 width_row_header;
            Integer                 width_col_header;
            Integer                 width_header_addin;
        };

        static void measure_elem(md::disp_stream_impl& os, const disp_stream* user, 
                                 const V& x, measure_params& m)
        {
            Integer n_subv                  = os.do_get_number_subvalues(user);

            for(Integer v = 0; v < n_subv; ++v)
            {
                Integer vw			        = os.do_get_value_min_width(user, x, v);
                vw					        = std::min(vw,os.do_get_max_value_width(user,v) + m.vcvw);
                m.min_value_width[v]	    = std::max(m.min_value_width[v],vw);
        
                if(os.do_show_values_labels(user))
                {
                    m.min_value_width[v]	= std::max(m.min_value_width[v],
                                                       os.do_get_subvalue_label_width(user,v+1));
                    m.min_value_width[v]	= std::min(m.min_value_width[v],
                                                       os.do_get_max_value_width(user,v) + m.vcvw);
                }
            }

            m.min_width = m.min_value_width[0];

            for(int j = 1; j < n_subv; ++j)
            {
                m.min_width += m.min_value_width[j];
            }
        };

        static void measure_sparse(md::disp_stream_impl& os, const disp_stream* user, 
                                   const matrix_type& m, Integer mr, Integer mc, Integer nnz, 
                                   measure_params& mp)
        {
            Integer c = m.cols();

            const details::sparse_ccs<V>& tmp = m.rep();

            const Integer* tmp_c    = tmp.ptr_c();
            const Integer* tmp_r    = tmp.ptr_r();
            const V* tmp_x          = tmp.ptr_x();

            Integer pos         = 0;
            Integer tc          = 0;
            
            Integer wr_min;
            Integer wr_max;

            os.do_get_row_header_width(user, wr_min, wr_max);

            Integer width_c     = 1;
            Integer width_r     = wr_min;

            std::set<Integer>   added_rows;

            for (Integer i = 0; i < c; ++i)
            {
                Integer w_col_h     = os.do_get_column_label_width_sparse(user, i);
                Integer pos_old     = pos;                

                Integer tr          = 0;

                if (tc > mc)
                {
                    for (Integer k = tmp_c[i]; k < tmp_c[i + 1]; ++k)
                    {
                        Integer row = tmp_r[k];

                        if (added_rows.find(row) == added_rows.end())                        
                        {
                            //measure row header
                            added_rows.insert(row);

                            Integer w_col_r = os.do_get_row_header_width_sparse(user, row);
                            width_r         = std::max(width_r, w_col_r);
                        };

                        measure_elem(os,user,tmp_x[k],mp);
                        ++pos;
                        break;
                    }

                    if (pos > pos_old)
                    {
                        //some elements in col are printed
                        width_c     = std::max(width_c, w_col_h);
                    };

                    break;
                }
                else
                {
                    for (Integer k = tmp_c[i]; k < tmp_c[i + 1]; ++k)
                    {
                        Integer row     = tmp_r[k];

                        if (added_rows.find(row) == added_rows.end())                        
                        {
                            //measure row header
                            added_rows.insert(row);

                            Integer w_col_r = os.do_get_row_header_width_sparse(user, row);
                            width_r         = std::max(width_r, w_col_r);
                        };

                        measure_elem(os,user,tmp_x[k],mp);

                        ++pos;

                        if (tr > mr)
                        {                            
                            break;
                        };

                        ++tr;

                        if (pos >= nnz)
                            break;
                    }
                    if (tr > 0)
                    {
                        ++tc;
                    };

                    if (pos > pos_old)
                    {
                        //some elements in col are printed
                        width_c     = std::max(width_c, w_col_h);
                    };
                };

                if (pos >= nnz)
                    break;
            };

            width_r     = std::min(width_r, wr_max);

            if (user->show_column_header_row() == false)
            {
                mp.width_row_header     = 0;
                mp.width_col_header     = 0;
                mp.width_header_addin   = 0;
            }
            else
            {
                Integer width_addin     = 4;

                mp.width_row_header     = std::max(mp.width_row_header, width_r);
                mp.width_col_header     = std::max(mp.width_col_header, width_c);
                mp.width_header_addin   = width_addin;
            };
        };

        static void eval(md::disp_stream_impl& os, const disp_stream* user, const matrix_type& m)
        {
            Integer r = m.rows(), c = m.cols(), nnz = m.nnz();
            
            matcl::value_code vt =  matrix_traits::value_code<V>::value;
            os.do_init_display(user,vt);

            os.do_start_display(user);	            
            os.do_displaying_sparse_matrix(user,r,c,nnz,vt,m.get_struct().to_string());            

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            disp_mode dm = os.do_get_display_mode(user);

            switch(dm)
            {
                case disp_mode::all_dense:
                case disp_mode::matrix_dense:
                {
                    matrix_provider_sparse<matrix_type,V> mdp(m);
                    return disp_matrix<V,struct_dense>::eval_matrix_body(os,user,mdp);
                }
            };

            if (nnz == 0)
            {
                bool short_print    = os.do_short_print_empty_matrix(user);

                if (short_print == true)
                {
                    os.do_display_empty_matrix(user,r,c);
                    os.do_end_display(user);
                    return;
                };
            };

            const details::sparse_ccs<V>& tmp = m.rep();

            const Integer* tmp_c    = tmp.ptr_c();
            const Integer* tmp_r    = tmp.ptr_r();
            const V* tmp_x          = tmp.ptr_x();

            Integer w_sep           = os.do_get_sparse_separator_width(user);

            Integer mc              = os.do_max_matrix_cols_sparse(user, c);
            Integer mr              = os.do_max_matrix_rows_sparse(user, r);
            Integer max_nnz         = os.do_max_matrix_nnz_sparse(user, nnz);

            measure_params mp;
            mp.min_width            = 1;
            mp.width_row_header     = 0;
            mp.width_col_header     = 0;
            mp.width_header_addin   = 0;

            Integer n_subv          = os.do_get_number_subvalues(user);

            mp.min_value_width      = std::vector<Integer>(n_subv, 0);

            measure_sparse(os,user,m,mr,mc,max_nnz, mp);

            Integer max_matrix_width= os.do_get_max_matrix_width(user);

            os.do_adjust_sparse_row_header(user,mp.width_row_header, mp.width_col_header, mp.width_header_addin,
                                           max_matrix_width/2);
            Integer rhw                 = mp.width_row_header + mp.width_col_header + mp.width_header_addin;

            Integer pos0                = os.do_get_line_prefix_width(user) + rhw + w_sep;
            Integer max_allowed_width   = max_matrix_width - pos0;

            os.do_adjust_sparse_column_min_width(user,mp.min_value_width, max_allowed_width);
            Integer chw     = os.do_get_values_labels_width(user,mp.min_value_width);

            os.do_start_display_matrix_block_sparse(user, rhw + w_sep + chw);

            if(os.do_show_values_labels(user))
            {
                os.do_disp_first_row_separator(user,rhw + w_sep + chw);

                os.do_disp_labels_row_id(user,rhw);
                os.do_disp_sparse_separator(user,w_sep);

                os.do_disp_values_labels(user,mp.min_value_width,0);
                os.do_disp_end_line(user,0);
                os.do_disp_first_row_separator(user,rhw + w_sep + chw);
            }

            Integer pos             = 0;
            Integer tc              = 0;
            bool need_separaror     = false;
            bool too_many_nnz       = false;

            for (Integer i = 0; i < c; ++i)
            {
                Integer col = i+1;

                Integer tr                  = 0;
                bool add_col_separator      = need_separaror;

                if (tc > mc || too_many_nnz == true)
                {
                    for (Integer k = tmp_c[i]; k < tmp_c[i + 1]; ++k)
                    {
                        ++pos;

                        os.do_disp_new_line(user,pos);
                        os.do_sparse_row(user,-1,-1,mp.width_row_header,mp.width_col_header);
                        os.do_disp_sparse_separator(user,w_sep);
                        os.do_disp_continuation_value(user,mp.min_width);
                        os.do_disp_end_line(user,pos);
                        break;
                    }
                    break;
                }
                else
                {
                    for (Integer k = tmp_c[i]; k < tmp_c[i + 1]; ++k)
                    {
                        if (add_col_separator)
                        {
                            os.do_disp_begin_next_column_sparse(user,rhw + w_sep + chw);
                            add_col_separator = false;
                        };
                        
                        Integer row = tmp_r[k] + 1;
                        ++pos;

                        if (tr > mr || too_many_nnz == true)
                        {
                            os.do_disp_new_line(user,pos);
                            os.do_sparse_row(user,-1,col,mp.width_row_header,mp.width_col_header);
                            os.do_disp_sparse_separator(user,w_sep);
                            os.do_disp_continuation_value(user,mp.min_width);
                            os.do_disp_end_line(user,pos);
                            break;
                        }
                        else
                        {			            
                            os.do_disp_new_line(user,pos);
                            os.do_sparse_row(user,row,col,mp.width_row_header,mp.width_col_header);
                            os.do_disp_sparse_separator(user,w_sep);

                            disp_elem_helper<V,struct_sparse>
                                    ::eval(os, user, tmp_x[k], mp.min_value_width,0);

                            os.do_disp_end_line(user,pos);
                        };

                        ++tr;
                        
                        if (pos >= max_nnz)
                            too_many_nnz = true;
                    }
                    if (tr > 0)
                    {
                        ++tc;
                        need_separaror = true;
                    };
                    if (pos >= max_nnz)
                        too_many_nnz = true;
                };
            };
            if (pos < nnz)
            {
            };

            os.do_end_display_matrix_block_sparse(user, rhw + w_sep + chw);

            os.do_end_display(user);
        };
    };

    template<class V>
    struct disp_matrix_dense
    {
        static void eval(md::disp_stream_impl& os, const disp_stream* user, 
                            const raw::Matrix<V,struct_dense>& m)
        {     
            using matrix_type = raw::Matrix<V,struct_dense>;            

            os.do_init_display(user, matrix_traits::value_code<V>::value);
            os.do_start_display(user);	
            os.do_displaying_dense_matrix(user,m.rows(),m.cols(),matrix_traits::value_code<V>::value, 
                                m.get_struct().to_string());

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            matrix_provider_dense<matrix_type,V> mdp(m);
            return disp_matrix<V,struct_dense>::eval_matrix_body(os,user,mdp);
        }
    };

};

template<class value_type>
void raw::disp(const disp_stream_ptr& os, const Matrix<value_type,struct_dense>& m)
{
    details::disp_matrix_dense<value_type>::eval(*os->impl(),os.get(),m);
};

template<class value_type>
void raw::disp(const disp_stream_ptr& os, const Matrix<value_type,struct_sparse>& m)
{
    details::disp_matrix<value_type,struct_sparse>::eval(*os->impl(),os.get(),m);
};

template<class value_type>
void raw::disp(const disp_stream_ptr& os, const Matrix<value_type,struct_banded>& m)
{
    details::disp_matrix<value_type,struct_banded>::eval(*os->impl(),os.get(),m);
};

template void raw::disp(const disp_stream_ptr& os, const Matrix<Integer,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Real,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Complex,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float_complex,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Object,struct_dense>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Integer,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Real,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Complex,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float_complex,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Object,struct_sparse>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Integer,struct_banded>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Real,struct_banded>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float,struct_banded>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Complex,struct_banded>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Float_complex,struct_banded>& m);
template void raw::disp(const disp_stream_ptr& os, const Matrix<Object,struct_banded>& m);

};};

#pragma warning(pop) 