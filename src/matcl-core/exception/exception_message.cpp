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

#include "matcl-core/error/exception_message.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-core/error/exception.h"
#include "matcl-core/IO/output_stream.h"
#include "matcl-core/memory/alloc.h"

#include <sstream>
#include <iostream>

namespace matcl { namespace error
{

static exception_message_ptr get_global_message()
{
    exception_message_ptr em(new default_exception_message());
    return em;
};
static exception_message_ptr messanger = get_global_message();

std::string default_exception_message::val_to_string(matcl::value_code vt)
{
    switch (vt)
    {
        case value_code::v_integer:         return "integer";
        case value_code::v_float:           return "float";
        case value_code::v_real:            return "real";
        case value_code::v_float_complex:   return "float complex";
        case value_code::v_complex:         return "complex";
        case value_code::v_object:          return "object";

        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

std::string default_exception_message::struct_to_string(matcl::struct_code st)
{
    switch (st)
    {
        case struct_code::struct_dense:
            return "dense matrix";
        case struct_code::struct_banded:
            return "banded matrix";
        case struct_code::struct_sparse:
            return "sparse matrix";
        case struct_code::struct_scalar:
            return "scalar";
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

const char* default_exception_message:: error_general(const std::string& msg)
{
    current_message = msg;
    return current_message.c_str();
};

std::string default_exception_message::matrix_type_string(matcl::mat_code mt)
{
    matcl::value_code val_type = matrix_traits::get_value_type(mt);
    matcl::struct_code struct_type = matrix_traits::get_struct_type(mt);

    std::string val_name = val_to_string(val_type);
    std::string struct_name = struct_to_string(struct_type);

    return val_name + " " + struct_name; 
};

const char* default_exception_message::alloc(size_t size)
{
    std::ostringstream buf;

    buf << "failed to allocate ";
    if (size == 0)
    {
        buf << "memory";
    }
    else
    {
        if (size >= 1048576) buf << size/1048576 << "MB";
        else
        {
            if (size >= 1024) buf << size/1024 << "kB";
            else buf << size << "B";
        }
    };

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_size(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix size: " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_size2(Integer r, Integer c, 
                                            Integer exp_r, Integer exp_c)
{
    std::ostringstream buf;

    buf << "invalid matrix size: " << r << "x" << c 
        << "; expecting matrix of size: " << exp_r << "x" << exp_c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_resize(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix size: " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::overflow_int_mult(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "unable to multply two integers: " << r << " and " << c 
        << ", resulting value is too large";

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::overflow_sizet_mult(size_t r, size_t c)
{
    std::ostringstream buf;

    buf << "unable to multply two integers: " << r << " and " << c 
        << ", resulting value is too large";

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::overflow_int_cast(Real val)
{
    std::ostringstream buf;

    buf << "unable to cast real to integer: value " << val << "is too large";

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_single_index(Integer i,Integer size)
{
    std::ostringstream buf;

    buf << "invalid single index; index: " << i << ", object size: " << size;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_double_index(Integer i,
                                        Integer j,Integer r,Integer c)
{
    std::ostringstream buf;

    buf << "invalid double index; index: (" << i << "," << j
        << "); matrix size: " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_diag(Integer d,Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix diagonal; diagonal: "<< d << ", matrix size: " 
        << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_index_band(Integer i, Integer j,
                                  Integer r, Integer c, Integer l, Integer u)
{
    std::ostringstream buf;

    buf << "invalid band matrix index; index: (" << i << "," << j
        << "); matrix size: " << r << "x" << c << " with " << l 
        << " subdiagonals and " << u << " superdiagonals";

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_size_band(Integer r, 
                                   Integer c, Integer l, Integer u)
{
    std::ostringstream buf;

    buf << "invalid band matrix size: " << r << "x" << c << " with " 
        << l << " subdiagonals and " << u << " superdiagonals";

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_size_sp(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid sparse matrix size: " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_return_type()
{
    current_message = "invalid return type";
    return current_message.c_str();
};
const char* default_exception_message::invalid_horzcat(Integer r1, 
                               Integer c1, Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "cannot horizontally concatenate " << r1 << "x" << c1
        << " and " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();

};
const char* default_exception_message::invalid_vertcat(Integer r1,
                               Integer c1, Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "cannot vertically concatenate " << r1 << "x" << c1
        << " and " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_index(Integer i, 
                               Integer j, Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid double index; index: (" << i << "," << j
        << "); matrix size: " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_colon_too_many_mat()
{
    current_message = "too many matrices in colon expression";
    return current_message.c_str();
};

const char* default_exception_message::invalid_assign_1(Integer s1,
                                            Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "nonconformant sizes in assigment; LHS size: " << s1
        << ", RHS is " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_assign_2(Integer r1,
                                 Integer c1, Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "nonconformant sizes in assigment; LHS is " << r1 << "x" << c1 
        << ", RHS is " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();
};

const char*
default_exception_message::invalid_row(Integer i, Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid row; row: "<< i << ", object size: " << r 
        << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_col(Integer j, Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "invalid column; column: "<< j << ", object size: " << r 
        << "x" << c;


    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_reshape(Integer r1,
                                  Integer c1, Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "cannot reshape " << r1 << "x" << c1 << " to " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_dim(Integer i, Integer d)
{
    std::ostringstream buf;

    buf << "invalid dimension; dimension: " << i 
        << ", object dimension: " << d;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::bspdiags_nonconf()
{
    current_message = "nonconformant arguments supplied to `diags'"
                      ", `bdiags' or `spdiags'";
    return current_message.c_str();
};

const char* default_exception_message::randperm_arg_neg(Integer n)
{
    std::ostringstream buf;

    buf << "argument supplied to `randperm' function is negative; it is " << n;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::bspdiag_1starg_not_vec(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "1st argument supplied to `bdiag' or `spdiag' function is "
        << "not vector; its size is " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::bspdiags_2ndarg_not_vec(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "2nd argument supplied to `diags', `bdiags' or `spdiags' function is "
        << "not vector; its size is " << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::diag_arg_not_vec(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "argument supplied to `diag' function is not vector; its size is "
        << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};		
const char* default_exception_message::linear_index_too_large(Integer r,
                                                Integer c, Integer rows)
{
    std::ostringstream buf;

    buf << "unable to create linear index, index is too large; row: " 
            << r << ", column : " << c << ", number of rows:" << rows;

    current_message = buf.str();
    return current_message.c_str();
};	
const char* default_exception_message::alloc_ext()
{
    current_message = "memory allocation error in external routine";
    return current_message.c_str();
};	

const char* default_exception_message::scalar_required(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "scalar or 1x1 matrix is required, matrix size is: " 
            << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::vector_required(Integer r, Integer c)
{
    std::ostringstream buf;

    buf << "vector is required, matrix size is: " 
            << r << "x" << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::unable_to_read_matrix()
{
    current_message = "unable to read matrix";
    return current_message.c_str();
};
const char*
default_exception_message::unable_to_read_mm_matrix(const std::string& reason)
{
    current_message = "unable to read Matrix Marker data; resonon: " + reason;
    return current_message.c_str();
};


const char* default_exception_message::unable_to_convert_invalid_code()
{
    current_message = "unable to convert matrix, invalid new matrix type code";
    return current_message.c_str();	
};

const char* 
default_exception_message::invalid_type_get(matcl::mat_code ret, matcl::mat_code in)
{
    std::ostringstream buf;

    std::string ret_type = matrix_type_string(ret);
    std::string in_type = matrix_type_string(in);

    buf << "unable to return " << ret_type << ", matrix type is " << in_type;

    current_message = buf.str();
    return current_message.c_str();
};
const char*
default_exception_message::unable_to_convert(matcl::mat_code ret, matcl::mat_code in)
{
    std::ostringstream buf;

    std::string ret_type = matrix_type_string(ret);
    std::string in_type = matrix_type_string(in);

    buf << "unable to convert " << in_type << " to " << ret_type;

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_vectors_spmat(Integer sr, Integer sc, Integer sx)
{
    std::ostringstream buf;

    buf << "invalid triplet representation, vectors size not match" << '\n' 
        << "row indices size: " << sr << ", columns indiced size: " 
        << sc << ", value size: " << sx;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_row_indices_sortcols(Integer s, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix of row indices supplied to function sortcols"
           ", matrix size exceeded number of rows\n"
        << "matrix size: " << s << ", number of rows: " << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_cols_indices_sortrows(Integer s, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix of row indices supplied to function sortcols"
        << ", matrix size exceeded number of rows\n"
        << "matrix size: " << s << ", number of rows: " << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_row_indices_elem_sortcols(Integer elem, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix of row indices supplied to function sortcols"
        << ", invalid row index\n"
        << "row index: " << elem << ", number of rows: " << c;

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_col_indices_elem_sortrows(Integer elem, Integer c)
{
    std::ostringstream buf;

    buf << "invalid matrix of column indices supplied to function sortrows"
        << ", invalid column index\n"
        << "column index: " << elem << ", number of columns: " << c;

    current_message = buf.str();
    return current_message.c_str();
};

const char*
default_exception_message::invalid_struct(const std::string& struct_flag_string)
{
    std::ostringstream buf;

    buf << "invalid struct flag: ";
    buf << struct_flag_string;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_permvec()
{
    current_message = "invalid permutation vector";
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_permvec_length(Integer len, Integer exp_len)
{
    std::ostringstream os;
    os << "invalid permutation vector length: " << len 
       << "; expected length: " << exp_len;

    current_message = os.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_permvec_composition(Integer l1, Integer l2)
{
    std::ostringstream os;
    os << "invalid permutation vector composition: first vector length: " 
       << l1 << ", second vector length: " << l2;

    current_message = os.str();
    return current_message.c_str();
};
const char* default_exception_message::invalid_eeop(Integer r1, Integer c1,
                                                   Integer r2, Integer c2)
{
    std::ostringstream buf;

    buf << "nonconfomant sizes in element by element operation; "
        << "first operand: " << r1 << "x" << c1
        << ", second operand: " << r2 << "x" << c2;

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::invalid_mul(Integer r1, Integer c1, Integer r2,
                                      Integer c2, trans_type t1, trans_type t2)
{
    std::ostringstream buf;

    buf << "nonconfomant sizes in multiplication; "
        << "first operand: " << r1 << "x" << c1
        << ", second operand: " << r2 << "x" << c2;

    if (t1 != trans_type::no_trans && t2 == trans_type::no_trans)
        buf << " , with transposition of the first argument";
    else if (t1 != trans_type::no_trans && t2 != trans_type::no_trans)
        buf << " , with transposition of the first and second argument";
    else if (t1 == trans_type::no_trans && t2 != trans_type::no_trans)
        buf << " , with transposition of the second argument";

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_gemm_C(value_code vc_C, 
                                   struct_code sc_C, value_code vc_exp)
{
    std::ostringstream buf;
    buf << "invalid argument passed to gemm function";

    if (vc_C == vc_exp)
    {
        //invalid struct
        buf << ", matrix C must be a dense matrix, C is: " 
            << struct_to_string(sc_C);
    }
    else
    {
        //invalid value type
        buf << ", matrix C has invalid value_code, expecting: " 
            << val_to_string(vc_exp)
            << ", C stores elements of type: " << val_to_string(vc_C);
    };

    current_message = buf.str();
    return current_message.c_str();    
};

const char* 
default_exception_message::unable_to_create_view(unable_create_view_reason r) 
{
    std::ostringstream buf;
    buf << "unable to create rectangle view, reason: ";
    switch(r)
    {
        case unable_create_view_reason::invalid_colon_type:
        {
            buf << "invalid subsref type";
            break;
        }
        case unable_create_view_reason::step_not_one:
        {
            buf << "colon step is not 1";
            break;
        };
        default:
        {
            buf << "unknown";
            break;
        };
    };

    current_message = buf.str();
    return current_message.c_str();    
};

const char* default_exception_message::dense_matrix_required(struct_code s)
{
    std::ostringstream buf;
    buf << "dense matrix is required; matrix type is: " << struct_to_string(s);

    current_message = buf.str();
    return current_message.c_str();  
};

const char* default_exception_message::invalid_extract_from_unique_matrix()
{
    current_message = "attempt to get invalid data from unique_matrix";
    return current_message.c_str();  
};

const char*
default_exception_message::object_value_type_not_allowed(const std::string& func)
{
    std::ostringstream buf;
    buf <<  "function " << func << " is not defined for object values";

    current_message = buf.str();
    return current_message.c_str();
};
const char* 
default_exception_message::integer_value_type_not_allowed(const std::string& name)
{
    std::ostringstream buf;
    buf <<  "function " << name << " is not defined for integer values";

    current_message = buf.str();
    return current_message.c_str();
};

const char* 
default_exception_message::function_not_defined_for_complex(const std::string& name)
{
    std::ostringstream buf;
    buf <<  "function " << name << " is not defined for complex values";

    current_message = buf.str();
    return current_message.c_str();
};

void default_exception_message::warning_precision_lost_float_to_int(Float val)
{
    (void)val;
    warning("conversion of float scalar to integer scalar: precision is lost");
};
void default_exception_message::warning_precision_lost_real_to_int(Real)
{
    warning("conversion of real scalar to integer scalar: precision is lost");
};
void default_exception_message::warning_precision_lost_compl_to_real(Complex val)
{
    (void)val;
    warning("conversion of complex scalar to real scalar: imaginary part ignored");
};
void default_exception_message::warning_precision_lost_real_to_float(Real val)
{
    (void)val;
    warning("conversion of real scalar to float scalar: precision is lost");
};
void default_exception_message::warning_precision_lost_float_compl_to_float
                                                (const Float_complex& val)
{
    (void)val;
    warning("conversion of float complex scalar to float scalar: imaginary part ignored");
};
void default_exception_message::warning_precision_lost_int_to_float(Integer val)
{
    (void)val;
    warning("conversion of integer scalar to float scalar: precision is lost");
};

void default_exception_message::possibly_inaccurate_result(const std::string& func_name)
{
    std::ostringstream os;
    os << "possibly inaccurate result from " << func_name;
    warning(os.str());
};

const char* default_exception_message::square_matrix_required(Integer rows, Integer cols)
{
    std::ostringstream msg;
    msg << "square matrix is required, matrix size is: " << rows << "x" << cols;
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::option_validator_error(const std::string& reason)
{
    std::ostringstream msg;
    msg << "invalid option value: " << reason;
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::invalid_option_type(const std::string& opt_name, const std::string& req_type,
                const std::string& opt_type)
{
    std::ostringstream msg;
    msg << "invalid option type " << req_type
        << ", option " << opt_name << " has type " << opt_type;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::optional_value_not_set()
{
    std::ostringstream msg;
    msg << "optional value not set";
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::option_unregistered(const std::string& opt_name)
{
    std::ostringstream msg;
    msg << "unknown option: " << opt_name;
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::uninitialized_disp_stream()
{
    std::ostringstream msg;
    msg << "uninitialized disp stream used";
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::uninitialized_output_stream()
{
    std::ostringstream msg;
    msg << "uninitialized output stream used";
    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::value_not_in_cache(const std::string& name, Integer prec, 
                Integer cache_prec)
{
    std::ostringstream msg;
    if (cache_prec < 0)
    {
        msg << "value with id " << name << " is not stored in cache";
    }
    else
    {
        msg << "value with id " << name << " and precision " << prec 
            << " is not stored in cache; there is a value with precision " << cache_prec;
    };

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::formatted_disp_invalid_column(Integer col, Integer num_cols)
{
    std::ostringstream msg;

    msg << "access to invalid column " << col << "; formatted disp has " << num_cols << " columns";

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::formatted_disp_invalid_row_size(Integer size, Integer req_size)
{
    std::ostringstream msg;

    msg << "invalid call to formatted_disp::disp_row function"
        << "; expecting " << req_size << " values"
        << "; number of supplied values is " << size;

    current_message = msg.str();
    return current_message.c_str();
};

const char* default_exception_message::limest_error_valid_initial_point_not_found(double x)
{
    std::ostringstream msg;

    msg << "limest error: could not find initial point with well defined value of supplied function"
        << "; point is: " << x;

    current_message = msg.str();
    return current_message.c_str();
}

const char* default_exception_message::limest_error_nonfinite_starting_point()
{
    std::ostringstream msg;

    msg << "limest error: unable to evaluate limit; point is not finite";

    current_message = msg.str();
    return current_message.c_str();
}

const char* default_exception_message::seq_error_omega_is_zero()
{
    std::ostringstream msg;

    msg << "sequence transformation error: omega is zero";

    current_message = msg.str();
    return current_message.c_str();
}

const char* default_exception_message::seq_error_two_equal_points(double x)
{
    std::ostringstream msg;

    msg << "sequence transformation error: two interpolation points are equal or too close to each other"
        << "; invalid point is: " << x;

    current_message = msg.str();
    return current_message.c_str();
}

void default_exception_message::warning(const std::string& msg)
{
    std::stringstream of;

    bool show_warning   = check_warnings_enabled();
    if (show_warning == true)
    {
        of << "warning: " << msg << "\n";
        global_output_stream()->disp(of);
    };
};
void set_global_messanger(exception_message_ptr msg)
{
    if (msg)
    {
        messanger = msg;
    };
};

exception_message_ptr get_global_messanger()
{
    return messanger;
};

};};
