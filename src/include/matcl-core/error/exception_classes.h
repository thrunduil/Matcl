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
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/error/exception.h"
#include "matcl-core/error/exception_message.h"
#include "matcl-core/details/exception_details.h"

#include <exception>
#include <string>

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 
#pragma warning(disable:4275)	//non dll-interface class used as base for dll-interface class

namespace matcl { namespace error
{

// classes representing exceptions thrown in matcl_basic
// these classes belongs to private interface and should not be used directly

// this function is called when memory corruption is detected
// error message is printed and std::terminate() is called
void MATCL_CORE_EXPORT memory_corrupted();

//=====================================================================
//
// classes representing exceptions throws in matcl
// these classes should not be used directly, since they can change
// in next versions of matcl and new exception classes may appear;
// therefore exception classes will not be documented
//
// matcl may throw matcl_exception (derived from std::exception),
// std::exception or assert_exception (not derived from std::exception)
//=====================================================================

// this class is throws if no additional information is available
class MATCL_CORE_EXPORT error_general : public matcl_exception
{
    public:
        std::string msg;

    public:
        error_general(const std::string& msg)	: msg(msg) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT overflow_int_cast : public matcl_exception
{
    public:		
        Real m_val;

    public:
        overflow_int_cast(Real val)	: m_val(val){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT overflow_int_mult : public matcl_exception
{
    public:
        Integer m_rows;
        Integer m_cols;

    public:
        overflow_int_mult(Integer r, Integer c)	: m_rows(r), m_cols(c) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT overflow_sizet_mult : public matcl_exception
{
    public:
        size_t  m_rows;
        size_t  m_cols;

    public:
        overflow_sizet_mult(size_t r, size_t c)	: m_rows(r), m_cols(c) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_perm_vectors_composition : public matcl_exception
{
    private:
        Integer m_length_1;
        Integer m_length_2;

    public:
        invalid_perm_vectors_composition(Integer p1_length, Integer p2_lenght)
            :m_length_1(p1_length), m_length_2(p2_lenght)
        {};
        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_struct: public matcl_exception
{
    private:
        std::string struct_flag_string;

    public:
        invalid_struct(const std::string& struct_flag_string_) 
            : struct_flag_string(struct_flag_string_)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_index : public matcl_exception
{
    public:
        Integer m_i;
        Integer m_j;
        Integer m_r;
        Integer m_c;

    public:
        invalid_index(Integer i, Integer j, Integer r, Integer c)
            : m_i(i),m_j(j),m_r(r),m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_colon_too_many_mat : public matcl_exception
{
    public:
        invalid_colon_too_many_mat(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_diag : public matcl_exception
{
    public:
        Integer m_d;
        Integer m_r;
        Integer m_c;

    public:
        invalid_diag(Integer d, Integer r, Integer c)
            : m_d(d), m_r(r), m_c(c)	
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT scalar_required : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        scalar_required(Integer r, Integer c)
            : m_r(r), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT vector_required : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        vector_required(Integer r, Integer c)
            : m_r(r), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_type_get : public matcl_exception
{
    public:
        matcl::mat_code m_ret;
        matcl::mat_code m_in;

    public:
        invalid_type_get(matcl::mat_code ret, matcl::mat_code in)
            : m_ret(ret), m_in(in)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_return_type : public matcl_exception
{
    public:
        invalid_return_type()	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT integer_value_type_not_allowed : public matcl_exception
{
    private:
        std::string m_func_name;

    public:
        integer_value_type_not_allowed(const std::string& func_name)
            :m_func_name(func_name)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT object_value_type_not_allowed : public matcl_exception
{
    private:
        std::string m_func_name;

    public:
        object_value_type_not_allowed(const std::string& func_name)
            :m_func_name(func_name)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT unable_to_read_matrix : public matcl_exception
{
    public:
        unable_to_read_matrix(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT unable_to_read_mm_matrix : public matcl_exception
{
    private:
        std::string m_reason;

    public:
        unable_to_read_mm_matrix(const std::string& reason)
            :m_reason(reason)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT unable_to_convert_invalid_code : public matcl_exception
{
    public:
        unable_to_convert_invalid_code() {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT alloc_ext : public matcl_exception
{
    public:
        alloc_ext() {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT bspdiags_nonconf : public matcl_exception
{
    public:
        bspdiags_nonconf() {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_reshape : public matcl_exception
{
    public:
        Integer m_r1;
        Integer m_c1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_reshape(Integer r1, Integer c1, Integer r2, Integer c2)		
            : m_r1(r1), m_c1(c1), m_r2(r2), m_c2(c2)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT unable_to_convert : public matcl_exception
{
    public:
        matcl::mat_code m_ret;
        matcl::mat_code m_in;

    public:
        unable_to_convert(matcl::mat_code ret, matcl::mat_code in)
            : m_ret(ret), m_in(in)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_vectors_spmat : public matcl_exception
{
    public:
        Integer m_sr;
        Integer m_sc;
        Integer m_sx;

    public:
        invalid_vectors_spmat(Integer sr, Integer sc, Integer sx)
            : m_sr(sr), m_sc(sc), m_sx(sx)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_double_index : public matcl_exception
{
    public:
        Integer m_i;
        Integer m_j;
        Integer m_r;
        Integer m_c;

    public:
        invalid_double_index(Integer i, Integer j, Integer r, Integer c)	
            : m_i(i), m_j(j), m_r(r), m_c(c)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_resize : public matcl_exception
{
    public:
        Integer m_rows;
        Integer m_cols;

    public:
        invalid_resize(Integer r, Integer c)	: m_rows(r), m_cols(c) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_index_band : public matcl_exception
{
    public:
        Integer m_i;
        Integer m_j;
        Integer m_r;
        Integer m_c;
        Integer m_l;
        Integer m_u;

    public:
        invalid_index_band(Integer i,Integer j, Integer r, Integer c,
                           Integer l, Integer u)
            : m_i(i), m_j(j), m_r(r), m_c(c), m_l(l), m_u(u)	
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_size_band : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;
        Integer m_l;
        Integer m_u;

    public:
        invalid_size_band(Integer r, Integer c, Integer l, Integer u)
            : m_r(r),m_c(c),m_l(l),m_u(u)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_size_sp : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        invalid_size_sp(Integer r, Integer c)
            : m_r(r),m_c(c){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_horzcat : public matcl_exception
{
    public:
        Integer m_r1;
        Integer m_c1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_horzcat(Integer r1, Integer c1, Integer r2, Integer c2)
            : m_r1(r1),m_c1(c1),m_r2(r2),m_c2(c2)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_vertcat : public matcl_exception
{
    public:
        Integer m_r1;
        Integer m_c1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_vertcat(Integer r1, Integer c1, Integer r2, Integer c2)
            : m_r1(r1),m_c1(c1),m_r2(r2),m_c2(c2)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_assign_2 : public matcl_exception
{
    public:
        Integer m_r1;
        Integer m_c1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_assign_2(Integer r1, Integer c1, Integer r2, Integer c2)
            : m_r1(r1),m_c1(c1),m_r2(r2),m_c2(c2)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_assign_1 : public matcl_exception
{
    public:
        Integer m_s1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_assign_1(Integer s1, Integer r2, Integer c2)
            : m_s1(s1),m_r2(r2),m_c2(c2)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_dim : public matcl_exception
{
    public:
        Integer m_i;
        Integer m_d;

    public:
        invalid_dim(Integer i, Integer d)		
            : m_i(i), m_d(d)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT randperm_arg_neg : public matcl_exception
{
    public:
        Integer m_n;

    public:
        randperm_arg_neg(Integer n)
            : m_n(n)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT bspdiag_1starg_not_vec : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        bspdiag_1starg_not_vec(Integer r, Integer c)
            : m_r(r), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT bspdiags_2ndarg_not_vec : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        bspdiags_2ndarg_not_vec(Integer r, Integer c)
            : m_r(r), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT diag_arg_not_vec : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;

    public:
        diag_arg_not_vec(Integer r, Integer c)
            : m_r(r), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT linear_index_too_large : public matcl_exception
{
    public:
        Integer m_r;
        Integer m_c;
        Integer m_rows;

    public:
        linear_index_too_large(Integer r, Integer c, Integer rows)
            : m_r(r), m_c(c),m_rows(rows)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_row_indices_sortcols : public matcl_exception
{
    public:
        Integer m_s;
        Integer m_c;

    public:
        invalid_row_indices_sortcols(Integer s, Integer c)
            : m_s(s), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_row_indices_elem_sortcols : public matcl_exception
{
    public:
        Integer m_elem;
        Integer m_c;

    public:
        invalid_row_indices_elem_sortcols(Integer elem, Integer c)
            : m_elem(elem), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_cols_indices_sortrows : public matcl_exception
{
    public:
        Integer m_s;
        Integer m_c;

    public:
        invalid_cols_indices_sortrows(Integer s, Integer c)
            : m_s(s), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_col_indices_elem_sortrows : public matcl_exception
{
    public:
        Integer m_elem;
        Integer m_c;

    public:
        invalid_col_indices_elem_sortrows(Integer elem, Integer c)
            : m_elem(elem), m_c(c)	{};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_permvec : public matcl_exception
{
    public:
        invalid_permvec() {};
        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_permvec_length : public matcl_exception
{
    private:
        Integer m_length;
        Integer m_exp_length;

    public:
        invalid_permvec_length(Integer length, Integer exp_length)
            :m_length(length), m_exp_length(exp_length)
        {};
        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_single_index : public matcl_exception
{
    public:
        Integer m_pos;
        Integer m_size;

    public:
        invalid_single_index(Integer i, Integer size)
            : m_pos(i), m_size(size) 
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT alloc : public matcl_exception
{
    public:
        size_t m_size;

    public:
        alloc(size_t size)  : m_size(size) {};

        //unknown memory error or unknown size
        alloc()             : m_size(0) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_size : public matcl_exception
{
    public:
        Integer m_rows;
        Integer m_cols;

    public:
        invalid_size(Integer r, Integer c)	: m_rows(r), m_cols(c) {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_size2 : public matcl_exception
{
    public:
        Integer m_rows;
        Integer m_cols;
        Integer m_exp_rows;
        Integer m_exp_cols;

    public:
        invalid_size2(Integer r, Integer c, Integer exp_r, Integer exp_c)
            : m_rows(r), m_cols(c), m_exp_rows(exp_r), m_exp_cols(exp_c)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_row : public matcl_exception
{
    public:
        Integer m_i;
        Integer m_r;
        Integer m_c;

    public:
        invalid_row(Integer i, Integer r, Integer c)
            : m_i(i), m_r(r), m_c(c)	
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_col : public matcl_exception
{
    public:
        Integer m_j;
        Integer m_r;
        Integer m_c;

    public:
        invalid_col(Integer j, Integer r, Integer c)
            : m_j(j), m_r(r), m_c(c)	
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_eeop : public matcl_exception
{
    public:
        Integer m_r1;
        Integer m_c1;
        Integer m_r2;
        Integer m_c2;

    public:
        invalid_eeop(Integer r1, Integer c1, Integer r2, Integer c2)		
            : m_r1(r1), m_c1(c1), m_r2(r2), m_c2(c2)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_mul : public matcl_exception
{
    public:
        Integer     m_r1;
        Integer     m_c1;
        Integer     m_r2;
        Integer     m_c2;
        trans_type  m_t1;
        trans_type  m_t2;

    public:
        invalid_mul(Integer r1, Integer c1, Integer r2, Integer c2, 
                    trans_type t1, trans_type t2)		
            : m_r1(r1), m_c1(c1), m_r2(r2), m_c2(c2), m_t1(t1), m_t2(t2)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT function_not_defined_for_complex : public matcl_exception
{
    private:
        std::string m_name;

    public:
        function_not_defined_for_complex(const std::string& s)
            :m_name(s)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_gemm_C : public matcl_exception
{
    private:
        value_code  m_vc_C;
        struct_code m_sc_C;
        value_code  m_vc_exp;

    public:
        invalid_gemm_C(value_code vc_C, struct_code sc_C, value_code vc_exp)
            :m_vc_C(vc_C), m_sc_C(sc_C), m_vc_exp(vc_exp)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT unable_to_create_view : public matcl_exception
{
    private:
        unable_create_view_reason   m_reason;

    public:
        unable_to_create_view(unable_create_view_reason r) :m_reason(r){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT dense_matrix_required : public matcl_exception
{
    private:
        struct_code     m_struct;

    public:
        dense_matrix_required(struct_code s) 
            :m_struct(s)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_extract_from_unique_matrix : public matcl_exception
{
    public:
        invalid_extract_from_unique_matrix(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT square_matrix_required: public matcl_exception
{
    private:
        Integer m_rows;
        Integer m_cols;

    public:
        square_matrix_required(Integer r, Integer c)
            :m_rows(r), m_cols(c)
        {};
        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT option_validator_error : public matcl_exception
{
    private:
        std::string m_reason;

    public:
        option_validator_error(const std::string& reason)
            :m_reason(reason)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT invalid_option_type  : public matcl_exception
{
    private:
        std::string m_opt_name;
        std::string m_req_type;
        std::string m_opt_type;

    public:
        invalid_option_type(const std::string& opt_name, const std::string& req_type,
                        const std::string& opt_type);

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT optional_value_not_set : public matcl_exception
{
    public:
        optional_value_not_set(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT option_unregistered : public matcl_exception
{
    private:
        std::string m_opt_name;

    public:
        option_unregistered(const std::string& opt_name);

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT uninitialized_disp_stream : public matcl_exception
{
    public:
        uninitialized_disp_stream(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT uninitialized_output_stream : public matcl_exception
{
    public:
        uninitialized_output_stream(){};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT formatted_disp_invalid_column : public matcl_exception
{
    private:
        Integer     m_col;
        Integer     m_num_cols;

    public:
        formatted_disp_invalid_column(Integer col, Integer num_cols)
            :m_col(col), m_num_cols(num_cols)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT formatted_disp_invalid_row_size : public matcl_exception
{
    private:
        Integer     m_size;
        Integer     m_req_size;

    public:
        formatted_disp_invalid_row_size(Integer size, Integer req_size)
            :m_size(size), m_req_size(req_size)
        {};

        virtual const char* what(exception_message& em) const override;
};

class MATCL_CORE_EXPORT value_not_in_cache : public matcl_exception
{
    private:
        std::string m_name;
        Integer     m_prec;
        Integer     m_cache_prec;

    public:
        // set cache_prec = -1, if value not in cache
        value_not_in_cache(const std::string& name, Integer prec, 
                            Integer cache_prec)
            :m_name(name), m_prec(prec), m_cache_prec(cache_prec)
        {};

        virtual const char* what(exception_message& em) const override;
};

};};

#pragma warning(pop)