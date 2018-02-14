/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-core/error/safe_string_message.h"
#include "matcl-core/memory/alloc.h"

#include <string>
#pragma warning(push)
#pragma warning(disable:4251)   // class needs to have dll-interface

namespace matcl { namespace error
{

// shared pointer to exception message class
using exception_message_ptr = std::shared_ptr<exception_message>;

// set global exception message
MATCL_CORE_EXPORT 
void                    set_global_messanger(exception_message_ptr msg);

// ge current exception message
MATCL_CORE_EXPORT 
exception_message_ptr   get_global_messanger();

// enum used by unable_create_view class
enum class unable_create_view_reason
{
    invalid_colon_type, step_not_one
};

// class responsible for creating string message for exception classes
class MATCL_CORE_EXPORT exception_message : public matcl_new_delete
{
    public:
        virtual ~exception_message(){};

    public:
        virtual void        warning(const std::string& msg) = 0;
        virtual void        warning_precision_lost_real_to_int(Real val) = 0;
        virtual void        warning_precision_lost_real_to_float(Real val) = 0;
        virtual void        warning_precision_lost_float_compl_to_float
                                (const Float_complex& val) = 0;
        virtual void        warning_precision_lost_float_to_int(Float val) = 0;        
        virtual void        warning_precision_lost_compl_to_real(Complex val) = 0;
        virtual void        warning_precision_lost_int_to_float(Integer val) = 0;
        virtual void        possibly_inaccurate_result(const std::string& func_name) = 0;

        virtual const char* error_general(const std::string& msg) = 0;
        virtual const char* overflow_int_mult(Integer r, Integer c) = 0;
        virtual const char* overflow_sizet_mult(size_t r, size_t c) = 0;
        virtual const char* overflow_int_cast(Real val) = 0;
        virtual const char* invalid_diag(Integer d,Integer r,Integer c) = 0;
        virtual const char* invalid_return_type() = 0;
        virtual const char* invalid_index(Integer i, Integer j, Integer r, 
                                Integer c) = 0;
        virtual const char* invalid_colon_too_many_mat() = 0;
        virtual const char* invalid_reshape(Integer r1, Integer c1, Integer r2, 
                                Integer c2) = 0;
        virtual const char* bspdiags_nonconf() = 0;
        virtual const char* alloc_ext() = 0;
        virtual const char* scalar_required(Integer r, Integer c) = 0;
        virtual const char* vector_required(Integer r, Integer c) = 0;
        virtual const char* invalid_type_get(matcl::mat_code ret, matcl::mat_code in)= 0;
        virtual const char* unable_to_convert(matcl::mat_code ret, matcl::mat_code in)= 0;
        virtual const char* unable_to_read_matrix()= 0;
        virtual const char* unable_to_read_mm_matrix(const std::string& reason)= 0;
        virtual const char* invalid_vectors_spmat(Integer sr, Integer sc, Integer sx)= 0;
        virtual const char* unable_to_convert_invalid_code()= 0;
        virtual const char* invalid_struct(const std::string& struct_flag_string) = 0;        
        virtual const char* alloc(size_t size) = 0;
        virtual const char* invalid_size(Integer r, Integer c) = 0;
        virtual const char* invalid_size2(Integer r, Integer c, Integer exp_r, 
                                Integer exp_c) = 0;
        virtual const char* invalid_resize(Integer r, Integer c) = 0;
        virtual const char* invalid_single_index(Integer pos,Integer size) = 0;
        virtual const char* invalid_double_index(Integer i,Integer j,Integer r,Integer c) = 0;
        virtual const char* invalid_index_band(Integer i,Integer j, Integer r, Integer c, 
                                Integer l, Integer u) = 0;
        virtual const char* invalid_size_band(Integer r, Integer c, Integer l, Integer u) = 0;
        virtual const char* invalid_size_sp(Integer r, Integer c) = 0;
        virtual const char* invalid_horzcat(Integer r1, Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* invalid_vertcat(Integer r1, Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* invalid_assign_1(Integer s1, Integer r2, Integer c2) = 0;
        virtual const char* invalid_assign_2(Integer r1, Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* invalid_row(Integer i, Integer r, Integer c) = 0;
        virtual const char* invalid_col(Integer j, Integer r, Integer c) = 0;
        virtual const char* invalid_dim(Integer i, Integer d)	= 0;
        virtual const char* randperm_arg_neg(Integer n) = 0;
        virtual const char* bspdiag_1starg_not_vec(Integer r, Integer c) = 0;
        virtual const char* bspdiags_2ndarg_not_vec(Integer r, Integer c) = 0;
        virtual const char* diag_arg_not_vec(Integer r, Integer c) = 0;
        virtual const char* linear_index_too_large(Integer r, Integer c, Integer rows) = 0;
        virtual const char* invalid_row_indices_sortcols(Integer s, Integer c)= 0;
        virtual const char* invalid_cols_indices_sortrows(Integer s, Integer c)= 0;
        virtual const char* invalid_row_indices_elem_sortcols(Integer elem, Integer c)= 0;
        virtual const char* invalid_col_indices_elem_sortrows(Integer elem, Integer c)= 0;
        virtual const char* invalid_permvec_composition(Integer l1, Integer l2) = 0;
        virtual const char* invalid_permvec() = 0;
        virtual const char* invalid_permvec_length(Integer len, Integer exp_len) = 0;
        virtual const char* invalid_eeop(Integer r1, Integer c1, Integer r2, Integer c2) = 0;
        virtual const char* invalid_mul(Integer r1, Integer c1, Integer r2, Integer c2, 
                                trans_type t1, trans_type t2) = 0;
        virtual const char* invalid_gemm_C(value_code vc_C, struct_code sc_C, 
                                value_code vc_exp) = 0;
        virtual const char* unable_to_create_view(unable_create_view_reason r) = 0;
        virtual const char* dense_matrix_required(struct_code s) = 0;
        virtual const char* invalid_extract_from_unique_matrix() = 0;
        virtual const char* square_matrix_required(Integer rows, Integer cols) = 0;

        virtual const char* option_validator_error(const std::string& reason) = 0;
        virtual const char* invalid_option_type(const std::string& opt_name, const std::string& req_type,
                                const std::string& opt_type) = 0;
        virtual const char* optional_value_not_set() = 0;
        virtual const char* option_unregistered(const std::string& opt_name) = 0;
        virtual const char* uninitialized_disp_stream() = 0;
        virtual const char* uninitialized_output_stream() = 0;

        virtual const char* formatted_disp_invalid_column(Integer col, Integer num_cols) = 0;
        virtual const char* formatted_disp_invalid_row_size(Integer size, Integer req_size) = 0;
        virtual const char* value_not_in_cache(const std::string& name, Integer prec, 
                                Integer cache_prec) = 0;

        virtual const char* object_value_type_not_allowed(const std::string& func) = 0;
        virtual const char* integer_value_type_not_allowed(const std::string& func) = 0;
        virtual const char* function_not_defined_for_complex(const std::string& name) = 0;

        virtual const char* limest_error_valid_initial_point_not_found(double x) = 0;
        virtual const char* limest_error_nonfinite_starting_point() = 0;
        virtual const char* seq_error_omega_is_zero() = 0;
        virtual const char* seq_error_two_equal_points(double x) = 0;
};

// default strings associated with exception classes
class MATCL_CORE_EXPORT default_exception_message : public exception_message
{
    private:
        safe_string_message current_message;

    public:
        static std::string  val_to_string(matcl::value_code vt);
        static std::string  struct_to_string(matcl::struct_code st);

    public:
        virtual void        warning(const std::string& msg) override;
        virtual void        warning_precision_lost_float_to_int(Float val) override;  
        virtual void        warning_precision_lost_real_to_int(Real val) override;
        virtual void        warning_precision_lost_real_to_float(Real val) override;
        virtual void        warning_precision_lost_compl_to_real(Complex val) override;
        virtual void        warning_precision_lost_float_compl_to_float(const 
                                Float_complex& val) override;
        virtual void        warning_precision_lost_int_to_float(Integer val) override;
        virtual void        possibly_inaccurate_result(const std::string& func_name) override;

        virtual const char* error_general(const std::string& msg) override;
        virtual const char* overflow_int_mult(Integer r, Integer c) override;
        virtual const char* overflow_sizet_mult(size_t r, size_t c) override;
        virtual const char* overflow_int_cast(Real val) override;
        virtual const char* invalid_diag(Integer d,Integer r,Integer c) override;
        virtual const char* invalid_return_type() override;
        virtual const char* invalid_index(Integer i, Integer j, Integer r, Integer c) override;
        virtual const char* invalid_colon_too_many_mat() override;
        virtual const char* invalid_reshape(Integer r1, Integer c1, Integer r2, 
                                Integer c2) override;
        virtual const char* bspdiags_nonconf() override;
        virtual const char* alloc_ext() override;
        virtual const char* scalar_required(Integer r, Integer c) override;
        virtual const char* vector_required(Integer r, Integer c) override;
        virtual const char* invalid_type_get(matcl::mat_code ret, matcl::mat_code in) override;
        virtual const char* unable_to_convert(matcl::mat_code ret, matcl::mat_code in) override;
        virtual const char* unable_to_read_matrix() override;
        virtual const char* unable_to_read_mm_matrix(const std::string& reason) override;
        virtual const char* invalid_vectors_spmat(Integer sr, Integer sc, Integer sx) override;
        virtual const char* unable_to_convert_invalid_code() override;
        virtual const char* invalid_struct(const std::string& struct_flag_string) override;
        virtual const char* invalid_permvec_composition(Integer l1, Integer l2) override;
        virtual const char* alloc(size_t size) override;
        virtual const char* invalid_size(Integer r, Integer c) override;
        virtual const char* invalid_size2(Integer r, Integer c, Integer exp_r, 
                                Integer exp_c) override;
        virtual const char* invalid_resize(Integer r, Integer c) override;
        virtual const char* invalid_single_index(Integer pos,Integer size) override;
        virtual const char* invalid_double_index(Integer i,Integer j,Integer r,Integer c) override;
        virtual const char* invalid_index_band(Integer i,Integer j, Integer r, Integer c, 
                                Integer l, Integer u) override;
        virtual const char* invalid_size_band(Integer r, Integer c, Integer l, Integer u) override;
        virtual const char* invalid_size_sp(Integer r, Integer c) override;
        virtual const char* invalid_horzcat(Integer r1, Integer c1, Integer r2, Integer c2) override;
        virtual const char* invalid_vertcat(Integer r1, Integer c1, Integer r2, Integer c2) override;
        virtual const char* invalid_assign_1(Integer s1, Integer r2, Integer c2) override;
        virtual const char* invalid_assign_2(Integer r1, Integer c1, Integer r2, Integer c2) override;
        virtual const char* invalid_row(Integer i, Integer r, Integer c) override;
        virtual const char* invalid_col(Integer j, Integer r, Integer c) override;
        virtual const char* invalid_dim(Integer i, Integer d) override;
        virtual const char* randperm_arg_neg(Integer n) override;
        virtual const char* bspdiag_1starg_not_vec(Integer r, Integer c) override;
        virtual const char* bspdiags_2ndarg_not_vec(Integer r, Integer c) override;
        virtual const char* diag_arg_not_vec(Integer r, Integer c) override;
        virtual const char* linear_index_too_large(Integer r, Integer c, Integer rows) override;
        virtual const char* invalid_row_indices_sortcols(Integer s, Integer c) override;
        virtual const char* invalid_cols_indices_sortrows(Integer s, Integer c) override;
        virtual const char* invalid_row_indices_elem_sortcols(Integer elem, Integer c) override;
        virtual const char* invalid_col_indices_elem_sortrows(Integer elem, Integer c) override;
        virtual const char* invalid_permvec() override;
        virtual const char* invalid_permvec_length(Integer len, Integer exp_len) override;
        virtual const char* invalid_eeop(Integer r1, Integer c1, Integer r2, Integer c2) override;
        virtual const char* invalid_mul(Integer r1, Integer c1, Integer r2, Integer c2,
                                trans_type t1, trans_type t2) override;
        virtual const char* invalid_gemm_C(value_code vc_C, struct_code sc_C, value_code vc_exp) override;
        virtual const char* unable_to_create_view(unable_create_view_reason r) override;
        virtual const char* dense_matrix_required(struct_code s) override;
        virtual const char* invalid_extract_from_unique_matrix() override;
        virtual const char* square_matrix_required(Integer rows, Integer cols) override;

        virtual const char* option_validator_error(const std::string& reason) override;
        virtual const char* invalid_option_type(const std::string& opt_name, const std::string& req_type,
                                const std::string& opt_type) override;
        virtual const char* optional_value_not_set() override;
        virtual const char* option_unregistered(const std::string& opt_name) override;
        virtual const char* uninitialized_disp_stream() override;
        virtual const char* uninitialized_output_stream() override;

        virtual const char* formatted_disp_invalid_column(Integer col, Integer num_cols) override;
        virtual const char* formatted_disp_invalid_row_size(Integer size, Integer req_size) override;
        virtual const char* value_not_in_cache(const std::string& name, Integer prec, 
                                Integer cache_prec) override;

        virtual const char* object_value_type_not_allowed(const std::string& func) override;
        virtual const char* integer_value_type_not_allowed(const std::string& func) override;
        virtual const char* function_not_defined_for_complex(const std::string& name) override;

        virtual const char* limest_error_valid_initial_point_not_found(double x) override;
        virtual const char* limest_error_nonfinite_starting_point() override;
        virtual const char* seq_error_omega_is_zero() override;
        virtual const char* seq_error_two_equal_points(double x) override;

        static std::string  matrix_type_string(matcl::mat_code mt);        
};

};};

#pragma warning(pop)