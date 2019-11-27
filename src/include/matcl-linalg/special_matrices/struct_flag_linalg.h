/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl
{

/// structure flag for positive definite matrices (not necessary symmetric/hermitian)
class MATCL_LINALG_EXPORT posdef_flag : public user_flag_config, public user_flag
{
    public:
        posdef_flag() : user_flag(user_flag::get<posdef_flag>()) {};

    private:
        virtual bool        test(const matcl::Matrix& mat) const override;
        
        virtual std::string tag() const override;

        virtual user_flag   abs(struct_flag sf, bool is_square) const override;
        virtual user_flag   real(struct_flag sf) const override;
        virtual user_flag   conj(struct_flag sf) const override;

        virtual user_flag   trans(struct_flag sf) const override;
        virtual user_flag   ctrans(struct_flag sf) const override;

        virtual user_flag   seft_mult(struct_flag sf, trans_type t1, trans_type t2, 
                                value_code vc, bool is_square) const override;
        virtual user_flag   scal(struct_flag sf, value_struct_class vc) const override;

        virtual user_flag   plus_both(struct_flag sf_left, struct_flag sf_right) const override;
        virtual user_flag   plus_1(struct_flag sf_left, struct_flag sf_right) const override;
        virtual user_flag   plus_2(struct_flag sf_left, struct_flag sf_right) const override;

        virtual user_flag   mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;
        virtual user_flag   mult_1(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;
        virtual user_flag   mult_2(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;

        virtual user_flag   kron_both(struct_flag sf_left, struct_flag sf_right) const;
        virtual user_flag   kron_1(struct_flag sf_left, struct_flag sf_right) const;
        virtual user_flag   kron_2(struct_flag sf_left, struct_flag sf_right) const;

        virtual user_flag   precision_increased(struct_flag sf_mat) const override;

        virtual bool        register_visit_scal() const override;
        virtual bool        visit_scal(struct_flag mat_flags, value_struct_class vc) const override;

        virtual bool        register_visit_plus() const;
        virtual bool        visit_plus(struct_flag sf_X, struct_flag sf_Y) const;
};

/// structure flag for positive semi-definite matrices (not necessary symmetric/hermitian)
class MATCL_LINALG_EXPORT semi_posdef_flag : public user_flag_config, public user_flag
{
    public:
        semi_posdef_flag() : user_flag(user_flag::get<semi_posdef_flag>()) {};

    private:
        virtual bool        test(const matcl::Matrix& mat) const override;
        virtual std::string tag() const override;

        virtual user_flag   abs(struct_flag sf, bool is_square) const override;
        virtual user_flag   real(struct_flag sf) const override;
        virtual user_flag   conj(struct_flag sf) const override;

        virtual user_flag   trans(struct_flag sf) const override;
        virtual user_flag   ctrans(struct_flag sf) const override;

        virtual user_flag   seft_mult(struct_flag sf, trans_type t1, trans_type t2, 
                                value_code vc, bool is_square) const override;
        virtual user_flag   scal(struct_flag sf, value_struct_class vc) const override;

        virtual user_flag   plus_both(struct_flag sf_left, struct_flag sf_right) const override;
        virtual user_flag   plus_1(struct_flag sf_left, struct_flag sf_right) const override;
        virtual user_flag   plus_2(struct_flag sf_left, struct_flag sf_right) const override;

        virtual user_flag   kron_both(struct_flag sf_left, struct_flag sf_right) const;
        virtual user_flag   kron_1(struct_flag sf_left, struct_flag sf_right) const;
        virtual user_flag   kron_2(struct_flag sf_left, struct_flag sf_right) const;

        virtual user_flag   mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;
        virtual user_flag   mult_1(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;
        virtual user_flag   mult_2(struct_flag sf_left, struct_flag sf_right, bool is_square) const override;

        virtual user_flag   precision_increased(struct_flag sf_mat) const override;

        virtual bool        register_visit_abs() const override;
        virtual bool        visit_abs(struct_flag mat_flags, bool is_square) const override;

        virtual bool        register_visit_self_mult() const override;
        virtual bool        visit_self_mult(struct_flag mat_flags, trans_type t1, trans_type t2, 
                                value_code vc, bool is_square) const override;
};

bool MATCL_LINALG_EXPORT is_posdef(const struct_flag& sf);
bool MATCL_LINALG_EXPORT is_semi_posdef(const struct_flag& sf);

};