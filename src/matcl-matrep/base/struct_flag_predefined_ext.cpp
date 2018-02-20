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

#include "matcl-matrep/matrix/struct_flag_ext.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/manip.h"

#include <algorithm>

namespace matcl 
{

namespace md = matcl::details;

struct_flag predefined_struct_ext::seft_mult_struct(struct_flag f, trans_type t1, trans_type t2, 
                                                    value_code vc, bool is_square_X, bool is_square_ret)
{
    bool is_re      = matrix_traits::is_real(vc);
    struct_flag sf  = mult_struct(f, f, t1, t2, is_re, is_re, is_square_X, is_square_X, is_square_ret);
    bool sym        = false;
    bool her        = false;

    if (t1 == trans_type::no_trans)
    {
        if (t2 == trans_type::trans)
        {
            sym     = true;
        }
        else if (t2 == trans_type::conj_trans)
        {
            if (is_re == true)
                sym = true;
            else
                her = true;
        }
        else
        {
            return sf;
        };
    }
    else if (t2 == trans_type::no_trans)
    {
        if (t1 == trans_type::trans)
        {
            sym     = true;
        }
        else if (t1 == trans_type::conj_trans)
        {
            if (is_re == true)
                sym = true;
            else
                her = true;
        }
        else
        {
            return sf;
        };
    }
    else
    {
        //somethis is wrong, one of t1, t2 should be no_trans
        return sf;
    };

    if (sym)
        sf.add(predefined_struct_type::sym);
    if (her)
        sf.add(predefined_struct_type::her);

    struct_flag su;
    md::user_struct_flags::get_seft_mult(f.to_int_user(), f, t1, t2, vc, is_square_ret, su);

    //no reset of previous flags; this was done by mult_struct
    sf.add_user(su);

    return sf;
};

struct_flag predefined_struct_ext::kron_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                               bool is_square_1, bool is_square_2, bool is_square)
{
    return op_struct(f1, f2, op_kron, is_re_1, is_re_2, is_square_1, is_square_2, is_square);
};

struct_flag predefined_struct_ext::mult_struct(struct_flag f1, struct_flag f2, trans_type t1, trans_type t2,
                                               bool is_re_1, bool is_re_2, bool is_square_1, bool is_square_2,
                                               bool is_square_ret)
{
    if (f1.is_general() || f2.is_general())
        return struct_flag();

    struct_flag f1_t, f2_t;

    switch(t1)
    {
        case trans_type::no_trans:
            f1_t = f1; break;
        case trans_type::trans:
            f1_t = predefined_struct::get_trans(f1); break;
        case trans_type::conj_trans:
            f1_t = predefined_struct::get_ctrans(f1); break;
        default:
            assert(0);
            throw;
    };
    switch(t2)
    {
        case trans_type::no_trans:
            f2_t = f2; break;
        case trans_type::trans:
            f2_t = predefined_struct::get_trans(f2); break;
        case trans_type::conj_trans:
            f2_t = predefined_struct::get_ctrans(f2); break;
        default:
            assert(0);
            throw;
    };

    return op_struct(f1_t,f2_t,op_mult, is_re_1, is_re_2, is_square_1, is_square_2, is_square_ret);
};

struct_flag predefined_struct_ext::plus_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                               bool is_square)
{
    if (f1.is_general() || f2.is_general())
        return struct_flag();

    return op_struct(f1,f2,op_plus, is_re_1, is_re_2, is_square, is_square, is_square);
};

struct_flag predefined_struct_ext::dmult_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                                bool is_square)
{
    return op_struct(f1,f2,op_dmult, is_re_1, is_re_2, is_square, is_square, is_square);
};

struct_flag predefined_struct_ext::minus_struct(struct_flag f1, struct_flag f2, bool is_re_1, bool is_re_2,
                                                bool is_square)
{
    if (f1.is_general() || f2.is_general())
        return struct_flag();

    return op_struct(f1,f2,op_minus, is_re_1, is_re_2, is_square, is_square, is_square);
};

struct_flag predefined_struct_ext::get_scal_rowcol(const struct_flag& sf)
{
    struct_flag out;
    out.set_ldiags(sf.get_ldiags());
    out.set_udiags(sf.get_udiags());

    //sym her flags are invalid
    //user flags are invalid
    return out;
};

struct_flag predefined_struct_ext::get_scal_mult(const struct_flag& sf, value_struct_class vt)
{
    bool is_zero    = (vt == value_struct_class::vc_zero);
    bool is_one     = (vt == value_struct_class::vc_one);

    if (is_zero == true)
    {
        struct_flag out;
        out.set_ldiags(struct_flag::diag_type::zero);
        out.set_udiags(struct_flag::diag_type::zero);
        return out;
    }

    if (sf.is_general())
        return sf;

    struct_flag out = sf;

    if (is_one == true)
        return out;

    if (out.is_id())
    {
        out = predefined_struct_type::diag;
    };

    if (out.has_her_flag() == true)
    {
        if (vt == value_struct_class::vc_general)
            out.set_her(false);
    };

    out.reset_user();

    struct_flag su;
    md::user_struct_flags::get_mult_scal(sf.to_int_user(), sf, vt, su);

    out.reset_user();
    out.add_user(su);

    return out;
};

struct_flag predefined_struct_ext::op_struct(struct_flag f1, struct_flag f2, op_type op,
                            bool is_re_1, bool is_re_2, bool is_square_1, bool is_square_2, 
                            bool is_square_ret)
{
    struct_flag out;

    if (op == op_mult)
    {
        if (f1.is_id()) return f2;
        if (f2.is_id()) return f1;
    };       
    if (op == op_kron)
    {
        if (f1.is_id() && f2.is_id())
            return predefined_struct_type::id;
    };

    struct_flag lu  = op_diags(f1, f2, op, is_square_2);
    struct_flag sh;
    
    if (is_square_ret == true)
        sh = op_symher(f1, f2, op, is_square_1, is_square_2, is_re_1, is_re_2);

    struct_flag su;
    switch(op)
    {
        case op_mult:
            md::user_struct_flags::get_mult(f1.to_int_user(), f2.to_int_user(), f1, f2, is_square_ret, su);
            break;
        case op_plus:
            md::user_struct_flags::get_plus(f1.to_int_user(), f2.to_int_user(), f1, f2, su);
            break;
        case op_minus:
            md::user_struct_flags::get_minus(f1.to_int_user(), f2.to_int_user(), f1, f2, su);
            break;
        case op_dmult:
            break;
        case op_kron:
            md::user_struct_flags::get_kron(f1.to_int_user(), f2.to_int_user(), f1, f2, su);
            break;
    };

    out.add(lu);
    out.add(sh);

    out.reset_user();
    out.add_user(su);

    return out;
};

struct_flag predefined_struct_ext::op_diags(struct_flag f1, struct_flag f2, op_type op, 
                                            bool is_square_2)
{
    struct_flag::diag_type ld = struct_flag::general;
    struct_flag::diag_type ud = struct_flag::general;;

    switch (op)
    {
        case op_plus:
        case op_minus:
        {
            ld  = std::min(f1.get_ldiags(), f2.get_ldiags());
            ud  = std::min(f1.get_udiags(), f2.get_udiags());

            if (f1.get_ldiags() == struct_flag::qtriang && f2.get_ldiags() == struct_flag::qtriang)
                ld  = struct_flag::one;

            if (f1.get_udiags() == struct_flag::qtriang && f2.get_udiags() == struct_flag::qtriang)
                ud  = struct_flag::one;

            break;
        }
        case op_dmult:
        {
            ld  = std::max(f1.get_ldiags(), f2.get_ldiags());
            ud  = std::max(f1.get_udiags(), f2.get_udiags());
            break;
        }
        case op_mult:
        {
            if (f1.get_ldiags() == struct_flag::zero || f2.get_ldiags() == struct_flag::zero)
                ld  = std::min(f1.get_ldiags(), f2.get_ldiags());
            else
                ld  = struct_flag::general;

            if (f1.get_udiags() == struct_flag::zero || f2.get_udiags() == struct_flag::zero)
                ud  = std::min(f1.get_udiags(), f2.get_udiags());
            else
                ud  = struct_flag::general;
            break;
        }
        case op_kron:
        {
            if (is_square_2 == true )
            {
                if (f1.get_ldiags() == struct_flag::zero)
                    ld  = f2.get_ldiags();
                if (f1.get_udiags() == struct_flag::zero)
                    ud  = f2.get_udiags();
            };
            break;
        }
    };

    struct_flag out;
    out.set_ldiags(ld);
    out.set_udiags(ud);

    return out;
};

struct_flag predefined_struct_ext::op_symher(struct_flag f1, struct_flag f2, op_type op, 
                                            bool is_square_1, bool is_square_2, 
                                            bool is_re_1, bool is_re_2)
{
    if (op == op_mult)
        return struct_flag();

    //dont set sym/her flag if neither f1 nor f2 has one; in this case
    //resulting struct can have sym/her flag only if f1 and f2 are diagonal
    //we dont set sym/her flags for diagonal matrices
    bool has_symher = f1.has_symher_flag() || f2.has_symher_flag();

    if (has_symher == false)
        return struct_flag();

    bool ret_re = is_re_1 && is_re_2;

    bool is_sym = f1.is_symmetric(is_square_1, is_re_1) && f2.is_symmetric(is_square_2, is_re_2);
    bool is_her = f1.is_hermitian(is_square_1, is_re_1) && f2.is_hermitian(is_square_2, is_re_2) 
                    && (ret_re == false);    

    struct_flag out;
    out.set_sym(is_sym);
    out.set_her(is_her);

    return out;
};

};
