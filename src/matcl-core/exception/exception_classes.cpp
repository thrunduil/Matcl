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

#include "matcl-core/error/exception_classes.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/general/thread.h"

#include <iostream>

#pragma warning( push )
#pragma warning(disable:4702)	// unreachable code
#include <boost/lexical_cast.hpp>
#pragma warning( pop )

namespace matcl { namespace error
{

void error::memory_corrupted()
{
    std::cerr << "memory is corrupted \n";
    std::terminate();
};

const char* error_general::what(exception_message& em) const
{
    return em.error_general(msg);
};

const char* overflow_int_mult::what(exception_message& em) const
{
    return em.overflow_int_mult(m_rows,m_cols);
};

const char* overflow_sizet_mult::what(exception_message& em) const
{
    return em.overflow_sizet_mult(m_rows,m_cols);
};

const char* overflow_int_cast::what(exception_message& em) const
{
    return em.overflow_int_cast(m_val);
};

const char* invalid_diag::what(exception_message& em) const
{
    return em.invalid_diag(m_d,m_r,m_c);
};

const char* invalid_return_type::what(exception_message& em) const
{
    return em.invalid_return_type();
};

const char* invalid_index::what(exception_message& em) const
{
    return em.invalid_index(m_i,m_j,m_r,m_c);
};

const char* invalid_colon_too_many_mat::what(exception_message& em) const
{
    return em.invalid_colon_too_many_mat();
};


const char* invalid_reshape::what(exception_message& em) const
{
    return em.invalid_reshape(m_r1,m_c1,m_r2, m_c2);
};

const char* bspdiags_nonconf::what(exception_message& em) const
{
    return em.bspdiags_nonconf();
};

const char* alloc_ext::what(exception_message& em) const
{
    return em.alloc_ext();
};

const char* scalar_required::what(exception_message& em) const
{
    return em.scalar_required(m_r,m_c);
};
const char* vector_required::what(exception_message& em) const
{
    return em.vector_required(m_r,m_c);
};

const char* invalid_type_get::what(exception_message& em) const
{
    return em.invalid_type_get(m_ret,m_in);
};
const char* unable_to_convert::what(exception_message& em) const
{
    return em.unable_to_convert(m_ret,m_in);
};
const char* unable_to_read_matrix::what(exception_message& em) const
{
    return em.unable_to_read_matrix();
};
const char* unable_to_read_mm_matrix::what(exception_message& em) const
{
    return em.unable_to_read_mm_matrix(m_reason);
};
const char* invalid_vectors_spmat::what(exception_message& em) const
{
    return em.invalid_vectors_spmat(m_sr,m_sc,m_sx);
};
const char* unable_to_convert_invalid_code::what(exception_message& em) const
{
    return em.unable_to_convert_invalid_code();
};

const char* object_value_type_not_allowed::what(exception_message& em) const
{
    return em.object_value_type_not_allowed(m_func_name);
};
const char* integer_value_type_not_allowed::what(exception_message& em) const
{
    return em.integer_value_type_not_allowed(m_func_name);
};

const char* invalid_struct::what(exception_message& em) const
{
    return em.invalid_struct(struct_flag_string);
};

const char* invalid_perm_vectors_composition::what(exception_message& em) const
{
    return em.invalid_permvec_composition(m_length_1,m_length_2);
};

const char* alloc::what(exception_message& em) const
{
    return em.alloc(m_size);
};

const char* invalid_size::what(exception_message& em) const
{
    return em.invalid_size(m_rows,m_cols);
};
const char* invalid_size2::what(exception_message& em) const
{
    return em.invalid_size2(m_rows,m_cols,m_exp_rows,m_exp_cols);
};

const char* invalid_resize::what(exception_message& em) const
{
    return em.invalid_resize(m_rows,m_cols);
};

const char* invalid_single_index::what(exception_message& em) const
{
    return em.invalid_single_index(m_pos,m_size);
};

const char* invalid_double_index::what(exception_message& em) const
{
    return em.invalid_double_index(m_i,m_j,m_r,m_c);
};
const char* invalid_index_band::what(exception_message& em) const
{
    return em.invalid_index_band(m_i,m_j,m_r,m_c,m_l,m_u);
};

const char* invalid_size_band::what(exception_message& em) const
{
    return em.invalid_size_band(m_r,m_c,m_l,m_u);
};

const char* invalid_size_sp::what(exception_message& em) const
{
    return em.invalid_size_sp(m_r,m_c);
};
const char* invalid_horzcat::what(exception_message& em) const
{
    return em.invalid_horzcat(m_r1,m_c1,m_r2,m_c2);
};

const char* invalid_vertcat::what(exception_message& em) const
{
    return em.invalid_vertcat(m_r1,m_c1,m_r2,m_c2);
};
const char* invalid_assign_1::what(exception_message& em) const
{
    return em.invalid_assign_1(m_s1,m_r2,m_c2);
};
const char* invalid_assign_2::what(exception_message& em) const
{
    return em.invalid_assign_2(m_r1,m_r2,m_r2,m_c2);
};
const char* invalid_row::what(exception_message& em) const
{
    return em.invalid_row(m_i,m_r,m_c);
};

const char* invalid_col::what(exception_message& em) const
{
    return em.invalid_col(m_j,m_r,m_c);
};
const char* invalid_dim::what(exception_message& em) const
{
    return em.invalid_dim(m_i,m_d);
};

const char* randperm_arg_neg::what(exception_message& em) const
{
    return em.randperm_arg_neg(m_n);
};

const char* bspdiag_1starg_not_vec::what(exception_message& em) const
{
    return em.bspdiag_1starg_not_vec(m_r,m_c);
};
const char* bspdiags_2ndarg_not_vec::what(exception_message& em) const
{
    return em.bspdiags_2ndarg_not_vec(m_r,m_c);
};

const char* diag_arg_not_vec::what(exception_message& em) const
{
    return em.diag_arg_not_vec(m_r,m_c);
};
const char* linear_index_too_large::what(exception_message& em) const
{
    return em.linear_index_too_large(m_r,m_c,m_rows);
};
const char* invalid_row_indices_sortcols::what(exception_message& em) const
{
    return em.invalid_row_indices_sortcols(m_s,m_c);
};
const char* invalid_cols_indices_sortrows::what(exception_message& em) const
{
    return em.invalid_cols_indices_sortrows(m_s,m_c);
};
const char* invalid_row_indices_elem_sortcols::what(exception_message& em) const
{
    return em.invalid_row_indices_elem_sortcols(m_elem,m_c);
};
const char* invalid_col_indices_elem_sortrows::what(exception_message& em) const
{
    return em.invalid_col_indices_elem_sortrows(m_elem,m_c);
};

const char* invalid_permvec::what(exception_message& em) const
{
    return em.invalid_permvec();
};
const char* invalid_permvec_length::what(exception_message& em) const
{
    return em.invalid_permvec_length(m_length, m_exp_length);
};
const char* invalid_eeop::what(exception_message& em) const
{
    return em.invalid_eeop(m_r1,m_c1,m_r2, m_c2);
};
const char* invalid_mul::what(exception_message& em) const
{
    return em.invalid_mul(m_r1,m_c1,m_r2, m_c2, m_t1, m_t2);
};

const char* function_not_defined_for_complex::what(exception_message& em) const
{
    return em.function_not_defined_for_complex(m_name);
};

const char* invalid_gemm_C::what(exception_message& em) const
{
    return em.invalid_gemm_C(m_vc_C, m_sc_C, m_vc_exp);
};

const char* unable_to_create_view::what(exception_message& em) const
{
    return em.unable_to_create_view(m_reason);
};

const char* dense_matrix_required::what(exception_message& em) const
{
    return em.dense_matrix_required(m_struct);
};

const char* invalid_extract_from_unique_matrix::what(exception_message& em) const
{
    return em.invalid_extract_from_unique_matrix();
};

const char* square_matrix_required::what(exception_message& em) const
{
    return em.square_matrix_required(m_rows, m_cols);
};

const char* option_validator_error::what(exception_message& em) const
{
    return em.option_validator_error(m_reason);
};

invalid_option_type::invalid_option_type(const std::string& opt_name, const std::string& req_type,
                const std::string& opt_type)
    :m_opt_name(opt_name), m_req_type(req_type), m_opt_type(opt_type)
{};

const char* invalid_option_type::what(exception_message& em) const
{
    return em.invalid_option_type(m_opt_name, m_req_type, m_opt_type);
};

const char* optional_value_not_set::what(exception_message& em) const
{
    return em.optional_value_not_set();
};

option_unregistered::option_unregistered(const std::string& opt_name)
    :m_opt_name(opt_name)
{};

const char* option_unregistered::what(exception_message& em) const
{
    return em.option_unregistered(m_opt_name);
};

const char* uninitialized_disp_stream::what(exception_message& em) const
{
    return em.uninitialized_disp_stream();
};

const char* uninitialized_output_stream::what(exception_message& em) const
{
    return em.uninitialized_output_stream();
};

const char* formatted_disp_invalid_column::what(exception_message& em) const
{
    return em.formatted_disp_invalid_column(m_col, m_num_cols);
};

const char* formatted_disp_invalid_row_size::what(exception_message& em) const
{
    return em.formatted_disp_invalid_row_size(m_size, m_req_size);
};

const char* value_not_in_cache::what(exception_message& em) const
{
    return em.value_not_in_cache(m_name, m_prec, m_cache_prec);
};

const char* limest_error_valid_initial_point_not_found::what(exception_message& em) const
{
    return em.limest_error_valid_initial_point_not_found(m_point);
};

const char* limest_error_nonfinite_starting_point::what(exception_message& em) const
{
    return em.limest_error_nonfinite_starting_point();
};

const char* seq_error_omega_is_zero::what(exception_message& em) const
{
    return em.seq_error_omega_is_zero();
};

const char* seq_error_two_equal_points::what(exception_message& em) const
{
    return em.seq_error_two_equal_points(m_point);
};

};};