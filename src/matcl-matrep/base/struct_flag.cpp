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

#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-core/IO/archive.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/func/raw/raw_manip.h"

#include <algorithm>
#include <boost/tokenizer.hpp>

namespace matcl 
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;

user_flag::user_flag(const size_t* code)
    :m_code(code)
{};

user_flag::user_flag()
    :m_code(nullptr)
{};

size_t user_flag::get_code() const
{
    return m_code? *m_code : 0;
};

bool user_flag::is_empty() const
{
    return m_code == nullptr;
};

user_flag user_flag_config::conj(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::abs(struct_flag, bool) const
{
    return user_flag();
};

user_flag user_flag_config::real(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::uminus(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::trans(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::ctrans(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::precision_lost(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::precision_increased(struct_flag) const
{
    return user_flag();
};

user_flag user_flag_config::scal(struct_flag, value_struct_class) const
{
    return user_flag();
};

user_flag user_flag_config::seft_mult(struct_flag, trans_type, trans_type, value_code, bool) const
{
    return user_flag();
}

user_flag user_flag_config::kron_both(struct_flag, struct_flag) const
{
    return user_flag();
}

user_flag user_flag_config::kron_1(struct_flag, struct_flag) const
{
    return user_flag();
}

user_flag user_flag_config::kron_2(struct_flag, struct_flag) const
{
    return user_flag();
}

user_flag user_flag_config::mult_both(struct_flag sf_left, struct_flag sf_right, bool) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::mult_1(struct_flag sf_left, struct_flag sf_right, bool) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::mult_2(struct_flag sf_left, struct_flag sf_right, bool) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::plus_both(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::plus_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::plus_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::minus_both(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::minus_1(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

user_flag user_flag_config::minus_2(struct_flag sf_left, struct_flag sf_right) const
{
    (void)sf_left;
    (void)sf_right;
    return user_flag();
};

bool user_flag_config::register_visit_conj() const
{
    return false;
}

bool user_flag_config::visit_conj(struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_abs() const
{
    return false;
}

bool user_flag_config::visit_abs(struct_flag, bool) const
{
    return false;
}

bool user_flag_config::register_visit_real() const
{
    return false;
}

bool user_flag_config::visit_real(struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_uminus() const
{
    return false;
}

bool user_flag_config::visit_uminus(struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_trans() const
{
    return false;
}

bool user_flag_config::visit_trans(struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_ctrans() const
{
    return false;
}

bool user_flag_config::visit_ctrans(struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_scal() const
{
    return false;
}

bool user_flag_config::visit_scal(struct_flag, value_struct_class) const
{
    return false;
}

bool user_flag_config::register_visit_self_mult() const
{
    return false;
}

bool user_flag_config::visit_self_mult(struct_flag, trans_type, trans_type, value_code, bool) const
{
    return false;
}

bool user_flag_config::register_visit_kron() const
{
    return false;
}

bool user_flag_config::visit_kron(struct_flag, struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_mmul() const
{
    return false;
}

bool user_flag_config::visit_mmul(struct_flag, struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_plus() const
{
    return false;
}

bool user_flag_config::visit_plus(struct_flag, struct_flag) const
{
    return false;
}

bool user_flag_config::register_visit_minus() const
{
    return false;
}

bool user_flag_config::visit_minus(struct_flag, struct_flag) const
{
    return false;
}

struct_flag::impl_type::impl_type()
{
    m_ldiags    = struct_flag::diag_type::general;
    m_udiags    = struct_flag::diag_type::general;
    m_sym       = false;
    m_her       = false;
    m_id        = false;
    m_user      = 0;
};

struct_flag::impl_type::impl_type(id_tag)
{
    m_ldiags    = struct_flag::diag_type::zero;
    m_udiags    = struct_flag::diag_type::zero;
    m_sym       = true;
    m_her       = true;
    m_id        = true;
    m_user      = 0;
}

struct_flag::impl_type::impl_type(diag_type ld, diag_type ud, bool sym, bool her)
{
    m_ldiags    = ld;
    m_udiags    = ud;
    m_sym       = sym;
    m_her       = her;
    m_id        = false;
    m_user      = 0;
}

struct_flag::impl_type::impl_type(size_t code)
{
    impl_type other = *reinterpret_cast<impl_type*>(&code);

    if (md::user_struct_flags::test_flags(other.m_user) == true)
    {
        *this = other;
        return;
    };

    //code is not valid, set general structure
    *this = impl_type();
    return;
};

struct_flag::impl_type struct_flag::impl_type::add(impl_type other) const
{
    impl_type ret;

    ret.m_ldiags    = std::max(this->m_ldiags, other.m_ldiags);
    ret.m_udiags    = std::max(this->m_udiags, other.m_udiags);
    ret.m_sym       = this->m_sym || other.m_sym;
    ret.m_her       = this->m_her || other.m_her;
    ret.m_id        = this->m_id || other.m_id;
    ret.m_user      = this->m_user | other.m_user;

    return ret;
}

struct_flag::impl_type struct_flag::impl_type::add_user(impl_type other) const
{
    impl_type ret   = *this;
    ret.m_user      = ret.m_user | other.m_user;

    return ret;
}

struct_flag::struct_flag()
{};

struct_flag::struct_flag(predefined_struct_type type) 
{
    switch(type)
    {        
        case predefined_struct_type::diag:
            m_flag = impl_type(diag_type::zero, diag_type::zero, false, false); return;
        case predefined_struct_type::tril:
            m_flag = impl_type(diag_type::general, diag_type::zero, false, false); return;
        case predefined_struct_type::triu:
            m_flag = impl_type(diag_type::zero, diag_type::general, false, false); return;
        case predefined_struct_type::qtril:  
            m_flag = impl_type(diag_type::general, diag_type::qtriang, false, false); return;
        case predefined_struct_type::qtriu:
            m_flag = impl_type(diag_type::qtriang, diag_type::general, false, false); return;
        case predefined_struct_type::hessl:
            m_flag = impl_type(diag_type::general, diag_type::one, false, false); return;
        case predefined_struct_type::hessu:
            m_flag = impl_type(diag_type::one, diag_type::general, false, false); return;

        case predefined_struct_type::sym:
            m_flag = impl_type(diag_type::general, diag_type::general, true, false); return;
        case predefined_struct_type::her:
            m_flag = impl_type(diag_type::general, diag_type::general, false, true); return;

        case predefined_struct_type::id:
            m_flag = impl_type(impl_type::id_tag()); return;

        default:
            //unknown case, ignore
            return;
    };
};

struct_flag::struct_flag(user_flag uf) 
    :m_flag()
{
    set_user(uf);
};

bool struct_flag::is_symmetric(bool is_square, bool is_real) const
{
    if (is_real)
        return has_sym_flag() || has_her_flag() || is_id() || is_square && is_diag();
    else
        return has_sym_flag();
};

bool struct_flag::is_hermitian(bool is_square, bool is_real) const
{
    if (is_real)
        return has_sym_flag() || has_her_flag() || is_id() || is_square && is_diag();
    else
        return has_her_flag();
};

void struct_flag::set_user(user_flag uf)
{
    if (uf.is_empty() == true)
        return;

    size_t pos      = uf.get_code() - 1;
    size_t flags    = m_flag.m_user;
    flags           |= (size_t(1) << pos);
    m_flag.m_user   = flags;
};

bool struct_flag::get_user(user_flag uf) const
{
    if (uf.is_empty() == true)
        return false;

    size_t pos          = uf.get_code() - 1;
    size_t test_res     = m_flag.m_user & (size_t(1) << pos);

    return test_res != 0;
};

void struct_flag::reset()
{
    m_flag = impl_type();
};

void struct_flag::reset_value()
{
    m_flag.m_id         = 0;
    m_flag.m_user       = 0;
};

void struct_flag::add(const struct_flag& t) const
{
    m_flag  = m_flag.add(t.m_flag);
}

void struct_flag::add_user(const struct_flag& t) const
{
    m_flag  = m_flag.add_user(t.m_flag);
}

bool struct_flag::operator==(const struct_flag& other) const
{
    return to_int() == other.to_int();
}

bool struct_flag::operator!=(const struct_flag& other) const
{
    return !operator==(other);
};

size_t struct_flag::to_int() const
{
    return *reinterpret_cast<const size_t*>(&this->m_flag);
}

size_t struct_flag::to_int_user() const
{
    return this->m_flag.m_user;
}

struct_flag struct_flag::from_int(size_t code)
{
    struct_flag ret;
    ret.m_flag = impl_type(code);
    return ret;
};

struct_flag struct_flag::from_int_user(size_t code)
{
    if (md::user_struct_flags::test_flags(code) == false)
        return struct_flag();

    struct_flag ret;
    ret.m_flag.m_user = code;
    return ret;
};

std::string struct_flag::to_string() const
{
    std::ostringstream msg;
    bool add = false;

    if (this->is_id())
    {
        msg << "identity";
        add = true;
    }
    else
    {
        if (this->is_diag())
        {
            msg << "diagonal";
            add = true;
        }
        else
        {
            if (this->is_tril())
            {
                msg << "lower triangular";
                add = true;
            }
            else if (this->is_qtril())
            {
                msg << "quasi lower triangular";
                add = true;
            }
            else if (this->is_hessl())
            {
                msg << "lower hessenberg";
                add = true;
            }

            if (this->is_triu())
            {
                if (add)
                    msg << ", ";

                msg << "upper triangular";
                add = true;
            }
            else if (this->is_qtriu())
            {
                if (add)
                    msg << ", ";

                add = true;
                msg << "quasi upper triangular";
            }
            else if (this->is_hessu())
            {
                if (add)
                    msg << ", ";

                add = true;
                msg << "upper hessenberg";
            }
        }

        if (this->has_sym_flag())
        {
            if (add)
                msg << ", ";

            add = true;
            msg << "symmetric";
        }
        if (this->has_her_flag())
        {
            if (add)
                msg << ", ";

            add = true;
            msg << "hermitian";
        }
    };

    md::user_struct_flags::to_string(m_flag.m_user, msg, add);

    if (add == false)
    {
        msg << "general";
    };

    return msg.str();
}

void struct_flag::save(oarchive_impl & ar, const unsigned int ver) const
{
    (void) ver;

    struct_flag sf = *this;
    sf.reset_user();

    bool has_user   = false;

    if (this->to_int_user() != 0)
    {
        has_user            = true;
        sf.m_flag.m_user    = 1;
    };

    size_t code = sf.to_int();
    ar << code;

    if (has_user == true)
    {
        bool add = false;
        std::ostringstream msg;

        md::user_struct_flags::save_as_string(this->m_flag.m_user, msg, add);

        ar << msg.str();
    };

    return;
};

void struct_flag::load(iarchive_impl & ar, const unsigned int ver)
{
    (void) ver;
    size_t code;
    ar >> code;

    struct_flag sf;
    sf = struct_flag::from_int(code);

    if (sf.to_int_user() != 0)
    {
        sf.reset_user();

        std::string flags;
        ar >> flags;

        struct_flag sf_u;
        sf_u.load_from_string(flags);

        sf.add_user(sf_u);
    };

    *this = sf;
}

std::string struct_flag::save_as_string() const
{
    std::ostringstream msg;
    bool add = false;

    if (this->is_id())
    {
        msg << "identity";
        add = true;
    }
    else
    {
        if (this->is_diag())
        {
            msg << "diagonal";
            add = true;
        }
        else
        {
            if (this->is_tril())
            {
                msg << "ltriang";
                add = true;
            }
            else if (this->is_qtril())
            {
                msg << "qltriang";
                add = true;
            }
            else if (this->is_hessl())
            {
                msg << "lhess";
                add = true;
            }

            if (this->is_triu())
            {
                if (add)
                    msg << "_";

                msg << "utriang";
                add = true;
            }
            else if (this->is_qtriu())
            {
                if (add)
                    msg << "_";

                add = true;
                msg << "qutriang";
            }
            else if (this->is_hessu())
            {
                if (add)
                    msg << "_";

                add = true;
                msg << "uhess";
            }
        }

        if (this->has_sym_flag())
        {
            if (add)
                msg << "_";

            add = true;
            msg << "symmetric";
        }
        if (this->has_her_flag())
        {
            if (add)
                msg << "_";

            add = true;
            msg << "hermitian";
        }
    };

    md::user_struct_flags::save_as_string(m_flag.m_user, msg, add);

    if (add == false)
    {
        msg << "general";
    };

    return msg.str();
}

static inline 
char tolower_impl(char c)
{
    return (char)::tolower(c);
};

void struct_flag::load_from_string(const std::string& name)
{
    using tokenizer     = boost::tokenizer<boost:: char_separator<char> >;
    using tok_iterator  = tokenizer::iterator;

    boost::char_separator<char> sep("_");
                
    tokenizer tok(name,sep);

    tok_iterator beg = tok.begin(); 
    tok_iterator end = tok.end();

    struct_flag out;

    for(; beg != end; ++beg)
    {        
        std::string cur = *beg;

        std::transform(cur.begin(), cur.end(), cur.begin(), &tolower_impl);

        if (cur == "identity")          out.add(predefined_struct_type::id);
        else if (cur == "diagonal")     out.add(predefined_struct_type::diag);
        else if (cur == "ltriang")      out.add(predefined_struct_type::tril);
        else if (cur == "qltriang")     out.add(predefined_struct_type::qtril);
        else if (cur == "lhess")        out.add(predefined_struct_type::hessl);
        else if (cur == "utriang")      out.add(predefined_struct_type::triu);
        else if (cur == "qutriang")     out.add(predefined_struct_type::qtriu);
        else if (cur == "uhess")        out.add(predefined_struct_type::hessu);
        else if (cur == "symmetric")    out.add(predefined_struct_type::sym);
        else if (cur == "hermitian")    out.add(predefined_struct_type::her);
        else if (cur == "general")      continue;
        else                            md::user_struct_flags::load_from_string(out, cur);
    };

    *this = out;
    return;
};

template<class Val, class Struct>
struct eval_struct_impl
{};

template<class Val>
struct eval_struct_impl<Val, struct_dense>
{
    using Mat = mr::Matrix<Val,struct_dense>;

    static void eval(const Mat& mat)
    {
        Integer ld = raw::get_ld(mat, 1, false);
        Integer ud = raw::get_ud(mat, 1, false);

        using diag_type = struct_flag::diag_type;

        diag_type str_ld    = mat.get_struct().get_ldiags();
        diag_type str_ud    = mat.get_struct().get_udiags();

        if (str_ld == struct_flag::zero && ld > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::qtriang && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::one && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::zero && ud > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::qtriang && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::one && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        Integer this_ld  = mat.ld();

        if (mat.get_struct().is_qtril() && ud != 0)
        {
            Integer s       = std::min(mat.rows(), mat.cols()-1);
            const Val* ptr  = mat.ptr() + this_ld;
            bool z          = false;

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_zero(val) == false)
                {
                    if (z)  throw error::invalid_struct(mat.get_struct().to_string());
                    else    z = true;
                }
                else
                {
                    z       = false;
                };

                ptr         += this_ld + 1;
            };
        };

        if (mat.get_struct().is_qtriu() && ld != 0)
        {
            Integer s       = std::min(mat.rows()-1, mat.cols());
            const Val* ptr  = mat.ptr() + 1;
            bool z          = false;

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_zero(val) == false)
                {
                    if (z)  throw error::invalid_struct(mat.get_struct().to_string());
                    else    z = true;
                }
                else
                {
                    z       = false;
                };

                ptr         += this_ld + 1;
            };
        }

        if (mat.get_struct().is_id())
        {
            if (ld > 0 || ud > 0)
                throw error::invalid_struct(mat.get_struct().to_string());

            if (mat.rows() != mat.cols())
                throw error::invalid_struct(mat.get_struct().to_string());

            Integer s       = std::min(mat.rows(), mat.cols());
            const Val* ptr  = mat.ptr();

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_one(val) == false)
                    throw error::invalid_struct(mat.get_struct().to_string());

                ptr         += this_ld + 1;
            };
        }

        if (mat.get_struct().has_sym_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_sym(mat,tol, false);

            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }

        if (mat.get_struct().has_her_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_her(mat,tol,false);
            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }
    };
};

template<class Val>
struct eval_struct_impl<Val, struct_sparse>
{
    using Mat = mr::Matrix<Val,struct_sparse>;
    static void eval(const Mat& mat)
    {
        Integer ld = raw::get_ld(mat, 1, false);
        Integer ud = raw::get_ud(mat, 1, false);

        using diag_type = struct_flag::diag_type;

        diag_type str_ld    = mat.get_struct().get_ldiags();
        diag_type str_ud    = mat.get_struct().get_udiags();

        if (str_ld == struct_flag::zero && ld > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::qtriang && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::one && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::zero && ud > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::qtriang && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::one && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (mat.get_struct().has_sym_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_sym(mat,tol,false);
            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }

        if (mat.get_struct().has_her_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_her(mat,tol,false);

            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }

        if (mat.get_struct().is_qtril() && ud != 0 && mat.nnz() != 0)
        {
            bool is_nz              = false;
            const Integer* ptr_c    = mat.rep().ptr_c();
            const Integer* ptr_r    = mat.rep().ptr_r();
            const Val* ptr_x        = mat.rep().ptr_x();
            Integer c               = mat.cols();

            for (Integer i = 1; i < c; ++i)
            {            
                if (ptr_c[i+1] == ptr_c[i])
                {
                    is_nz           = false;
                    continue;
                };

                bool is1            = false;
                for (Integer k = ptr_c[i]; k < ptr_c[i+1]; ++k)
                {
                    if (mrd::is_zero(ptr_x[k]) == true)
                        continue;

                    if (ptr_r[k] == i-1)
                    {
                        if (is_nz == true)
                        {
                            throw error::invalid_struct(mat.get_struct().to_string());
                        }
                        else
                        {
                            is_nz   = true;
                            is1     = true;
                        };
                    }
                    else
                    {
                        break;
                    };
                };
                if (is1 == false)
                {
                    is_nz           = false;
                };
            };
        };

        if (mat.get_struct().is_qtriu() && ld != 0 && mat.nnz() != 0)
        {
            bool is_nz = false;
            const Integer* ptr_c    = mat.rep().ptr_c();
            const Integer* ptr_r    = mat.rep().ptr_r();
            const Val* ptr_x        = mat.rep().ptr_x();
            Integer c               = mat.cols();

            for (Integer i = 1; i < c; ++i)
            {            
                if (ptr_c[i+1] == ptr_c[i])
                {
                    is_nz           = false;
                    continue;
                };

                bool is1            = false;
                for (Integer k = ptr_c[i+1]-1; k >= ptr_c[i]; --k)
                {
                    if (mrd::is_zero(ptr_x[k]) == true)
                        continue;

                    if (ptr_r[k] == i+1)
                    {
                        if (is_nz == true)
                            throw error::invalid_struct(mat.get_struct().to_string());
                        else
                        {
                            is_nz   = true;
                            is1     = true;
                        };
                    }
                    else
                    {
                        break;
                    };
                };

                if (is1 == false)
                    is_nz           = false;
            };
        }

        if (mat.get_struct().is_id())
        {
            if (ld > 0 || ud > 0)
                throw error::invalid_struct(mat.get_struct().to_string());

            if (mat.rows() != mat.cols())
                throw error::invalid_struct(mat.get_struct().to_string());

            const Integer* ptr_c    = mat.rep().ptr_c();
            const Integer* ptr_r    = mat.rep().ptr_r();
            const Val* ptr_x        = mat.rep().ptr_x();
            Integer c               = mat.cols();

            for (Integer i = 0; i < c; ++i)
            {
                bool is1            = false;
                for (Integer k = ptr_c[i]; k < ptr_c[i+1]; ++k)
                {
                    if (mrd::is_zero(ptr_x[k]) == true)
                        continue;

                    if (ptr_r[k] == i)
                        is1         = true;
                    else
                        throw error::invalid_struct(mat.get_struct().to_string());

                    if (mrd::is_one(ptr_x[k]) == false)
                        throw error::invalid_struct(mat.get_struct().to_string());
                };

                if (is1 == false)
                    throw error::invalid_struct(mat.get_struct().to_string());
            };
        }
    };
};

template<class Val>
struct eval_struct_impl<Val, struct_banded>
{
    using Mat = mr::Matrix<Val,struct_banded>;

    static void eval(const Mat& mat)
    {
        Integer ld = raw::get_ld(mat, 1, false);
        Integer ud = raw::get_ud(mat, 1, false);

        using diag_type = struct_flag::diag_type;

        diag_type str_ld    = mat.get_struct().get_ldiags();
        diag_type str_ud    = mat.get_struct().get_udiags();

        if (str_ld == struct_flag::zero && ld > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::qtriang && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ld == struct_flag::one && ld > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::zero && ud > 0)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::qtriang && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        if (str_ud == struct_flag::one && ud > 1)
            throw error::invalid_struct(mat.get_struct().to_string());

        Integer this_ld  = mat.ld();

        if (mat.get_struct().is_qtril() && ud != 0)
        {
            Integer s       = mat.diag_length(1);
            const Val* ptr  = mat.rep_ptr() + mat.first_elem_diag(1);
            bool z          = false;

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_zero(val) == false)
                {
                    if (z)  
                        throw error::invalid_struct(mat.get_struct().to_string());
                    else
                        z = true;
                }
                else
                {
                    z       = false;
                };

                ptr         += this_ld;
            };
        };

        if (mat.get_struct().is_qtriu() && ld != 0)
        {
            Integer s       = mat.diag_length(-1);
            const Val* ptr  = mat.rep_ptr() + mat.first_elem_diag(-1);
            bool z          = false;

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_zero(val) == false)
                {
                    if (z)  
                        throw error::invalid_struct(mat.get_struct().to_string());
                    else
                        z = true;
                }
                else
                {
                    z       = false;
                };

                ptr         += this_ld;
            };
        }

        if (mat.get_struct().is_id())
        {
            if (ld > 0 || ud > 0)
                throw error::invalid_struct(mat.get_struct().to_string());

            if (mat.rows() != mat.cols())
                throw error::invalid_struct(mat.get_struct().to_string());

            Integer s       = mat.diag_length(0);
            const Val* ptr  = mat.rep_ptr() + mat.first_elem_diag(0);

            for (Integer i = 0; i < s; ++i)
            {
                const Val& val = *ptr;

                if (mrd::is_one(val) == false)
                    throw error::invalid_struct(mat.get_struct().to_string());

                ptr         += this_ld;
            };
        };

        if (mat.get_struct().has_sym_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_sym(mat, tol, false);
            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }

        if (mat.get_struct().has_her_flag())
        {
            Real tol    = 1.0e2;
            bool is_sym = raw::is_her(mat, tol, false);

            if (is_sym == false)
                throw error::invalid_struct(mat.get_struct().to_string());
        }
    };
};

struct eval_check_struct : public details::extract_type_switch<void, eval_check_struct, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat)
    {
        using value_type    = typename T::value_type;
        using struct_type   = typename T::struct_type;
        return eval_struct_impl<value_type, struct_type>::eval(mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat)
    {
        (void)handle;
        (void)mat;
        return;
    };
};

void matcl::check_struct(const matcl::Matrix& A)
{
    eval_check_struct::make<const Matrix&>(A);
    md::user_struct_flags::check_struct(A.get_struct().to_int_user(), A);    
};

};
