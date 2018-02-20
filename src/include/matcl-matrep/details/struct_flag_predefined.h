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

#pragma once

#include "matcl-matrep/general/config.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-core/matrix/enums.h"

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl { namespace details
{

class MATCL_MATREP_EXPORT user_struct_flags
{
    public:
        static void         check_struct(size_t code, const matcl::Matrix& mat);

        static bool         test_flags(size_t code);
        static void         to_string(size_t code, std::ostringstream& str, bool& add);
        static void         save_as_string(size_t code, std::ostringstream& str, bool& add);
        static bool         load_from_string(struct_flag& flags, const std::string& struct_name);

        static void         get_abs(size_t code, struct_flag sf_mat, bool is_square, struct_flag& flags);
        static void         get_conj(size_t code, struct_flag sf_mat, struct_flag& flags);
        static void         get_real(size_t code, struct_flag sf_mat, struct_flag& flags);
        static void         get_trans(size_t code, struct_flag sf_mat, struct_flag& flags);
        static void         get_ctrans(size_t code, struct_flag sf_mat, struct_flag& flags);
        static void         get_uminus(size_t code_user, struct_flag sf_mat, struct_flag& sfu);
        static void         get_warnings(size_t code, struct_flag sf_mat, struct_flag& flags);
        static void         precision_increased(size_t code, struct_flag sf_mat, struct_flag& flags);

        static void         get_seft_mult(size_t code, struct_flag sf_mat, trans_type t1, 
                                trans_type t2, value_code vc, bool is_square, struct_flag& flags);
        static void         get_kron(size_t f1, size_t f2, struct_flag sf_left, struct_flag sf_right,
                                struct_flag& flags);
        static void         get_mult_scal(size_t code, struct_flag sf_mat, value_struct_class vc, 
                                struct_flag& flags);
        static void         get_mult(size_t f1, size_t f2, struct_flag sf_left, struct_flag sf_right, 
                                bool is_square, struct_flag& flags);
        static void         get_plus(size_t f1, size_t f2, struct_flag sf_left, struct_flag sf_right,
                                struct_flag& flags);
        static void         get_minus(size_t f1, size_t f2, struct_flag sf_left, struct_flag sf_right,
                                struct_flag& flags);
};

class MATCL_MATREP_EXPORT predefined_struct
{
    public:
        using value_class   = value_struct_class;

    public:        
        static struct_flag  get_set_diag(const struct_flag& sf, Integer d, value_struct_class vt, 
                                    bool is_real_mat, bool is_square);
        static struct_flag  get_rectangle_view(const struct_flag& sf, bool is_square);
        static struct_flag  get_resize(const struct_flag& sf, bool is_sym, bool is_real_mat);
        
        static struct_flag  get_tril(const struct_flag& sf, Integer d, bool is_real_mat);
        static struct_flag  get_triu(const struct_flag& sf, Integer d, bool is_real_mat);
        static struct_flag  get_qtril(const struct_flag& sf);
        static struct_flag  get_qtriu(const struct_flag& sf);

        static struct_flag  get_trans(const struct_flag& sf, trans_type t);
        static struct_flag  get_trans(const struct_flag& sf);
        static struct_flag  get_ctrans(const struct_flag& sf);
        
        static struct_flag  get_abs(const struct_flag& sf, bool is_square);
        static struct_flag  get_conj(const struct_flag& sf);
        static struct_flag  get_real(const struct_flag& sf);

        // set additional structures after calling to eval_struct
        static struct_flag  uminus_cont(const struct_flag& sf_ret, const struct_flag& sf_mat);
        static struct_flag  inv_cont(const struct_flag& sf_ret, const struct_flag& sf_mat);
        
        static struct_flag  set_warnings(const struct_flag& s1, bool prec_lost, bool compl_to_real);
        static struct_flag  precision_increased(struct_flag f1);

        static value_class  get_value_type(Integer val, bool test_zero);
        static value_class  get_value_type(Real val, bool test_zero);
        static value_class  get_value_type(Float val, bool test_zero);
        static value_class  get_value_type(Complex val, bool test_zero);
        static value_class  get_value_type(Float_complex val, bool test_zero);
        static value_class  get_value_type(const Object& val, bool test_zero);        

        static struct_flag  eval_struct(struct_flag f1, bool is_zero_id);        
};

};};

#pragma warning(pop)
