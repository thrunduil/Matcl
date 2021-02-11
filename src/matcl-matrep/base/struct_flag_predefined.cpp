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

#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-internals/base/utils.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-core/IO/archive.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-matrep/details/bit_manip.h"

#include <algorithm>
#include <map>

namespace matcl { namespace details
{

namespace mrd = matcl::raw::details;

//------------------------------------------------------------------------
//                  registered_user_flags
//------------------------------------------------------------------------        
class registered_user_flags
{
    private:
        enum reg_ops
        {
            op_conj = 0, op_abs, op_real, op_uminus, op_trans, op_ctrans,
            op_scal, op_self_mult, op_kron, op_mult, op_plus, op_minus,

            op_last
        };

    private:
        using config        = std::shared_ptr<user_flag_config>;
        using code_name     = std::pair<size_t,std::string>;
        using string_map    = std::map<std::string, code_name>;
        using vec_confing   = std::vector<config>;
        using vec_string    = std::vector<std::string>;
        using vec_sizet     = std::vector<const size_t*>;
        using vec_visit     = std::vector<vec_sizet>;

    private:
        string_map          m_string_map;
        vec_confing         m_configs;
        vec_string          m_flag_names;
        vec_visit           m_visit;

    public:        
        registered_user_flags();

        size_t              number_registered_flags() const;
        user_flag           flag_from_string(const std::string& str) const;
        std::string         to_string_one(size_t struct_code);
        const config&       get_user_flag_config(size_t struct_code);
        const size_t*       make_register(const config& conf, const char* cl_name);

        static registered_user_flags*
                            get();        

    private:
        void                get_registered_vis(Integer& N, const size_t* const*& reg_flags, reg_ops op) const;
        void                register_visitors(const size_t* code, const config& conf);
        void                register_visitor(const size_t* code, reg_ops op);

    public:
        void                visit_conj_impl(struct_flag mat_flags, struct_flag& flags);
        void                visit_abs_impl(struct_flag mat_flags, bool is_square, struct_flag& flags);
        void                visit_real_impl(struct_flag mat_flags, struct_flag& flags);
        void                visit_uminus_impl(struct_flag mat_flags, struct_flag& flags);
        void                visit_trans_impl(struct_flag mat_flags, struct_flag& flags);
        void                visit_ctrans_impl(struct_flag mat_flags, struct_flag& flags);
        void                visit_scal_impl(struct_flag mat_flags, value_struct_class vc,
                                    struct_flag& flags);
        void                visit_self_mult_impl(struct_flag mat_flags, trans_type t1, trans_type t2,
                                    value_code vc, bool is_square, struct_flag& flags);
        void                visit_kron_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags);
        void                visit_mmul_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags);
        void                visit_plus_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags);
        void                visit_minus_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags);

    public:
        static void         check_struct_one(size_t struct_code, const matcl::Matrix& mat);                 

        static user_flag    get_conj_one(size_t struct_code, struct_flag sf);
        static user_flag    get_abs_one(size_t struct_code, struct_flag sf, bool is_square);
        static user_flag    get_uminus_one(size_t struct_code, struct_flag sf);
        static user_flag    get_real_one(size_t struct_code, struct_flag sf);
        static user_flag    get_trans_one(size_t struct_code, struct_flag sf);
        static user_flag    get_ctrans_one(size_t struct_code, struct_flag sf);
        static user_flag    get_warnings_one(size_t struct_code, struct_flag sf);
        static user_flag    precision_increased_one(size_t struct_code, struct_flag sf);

        static user_flag    get_mult_scal(size_t struct_code, struct_flag sf, value_struct_class vc);
        static user_flag    get_seft_mult(size_t struct_code, struct_flag sf, trans_type t1, trans_type t2, 
                                value_code vc, bool is_square);
        static user_flag    get_kron_one_both(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        static user_flag    get_kron_one_1(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        static user_flag    get_kron_one_2(size_t struct_code, struct_flag sf_left, struct_flag sf_right);

        static user_flag    get_mult_one_both(size_t struct_code, struct_flag sf_left, struct_flag sf_right, 
                                bool is_square);
        // first has given structure, second does not
        static user_flag    get_mult_one_1(size_t struct_code, struct_flag sf_left, struct_flag sf_right, 
                                bool is_square);
        // second has given structure, first does not
        static user_flag    get_mult_one_2(size_t struct_code, struct_flag sf_left, struct_flag sf_right, 
                                bool is_square);

        static user_flag    get_plus_one_both(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        // first has given structure, second does not
        static user_flag    get_plus_one_1(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        // second has given structure, first does not
        static user_flag    get_plus_one_2(size_t struct_code, struct_flag sf_left, struct_flag sf_right);

        static user_flag    get_minus_one_both(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        // first has given structure, second does not
        static user_flag    get_minus_one_1(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
        // second has given structure, first does not
        static user_flag    get_minus_one_2(size_t struct_code, struct_flag sf_left, struct_flag sf_right);
};

registered_user_flags g_instance;

registered_user_flags* registered_user_flags::get()
{
    return &g_instance;
};

registered_user_flags::registered_user_flags()
    :m_visit(op_last)
{};

void registered_user_flags::get_registered_vis(Integer& N, const size_t* const*& reg_flags, reg_ops op) const
{
    const vec_sizet& v = m_visit[(size_t)op];
    N           = (Integer)v.size();
    reg_flags   = v.data();
};

const size_t* struct_flag_register_impl::get_code(const std::shared_ptr<user_flag_config>& conf,
                                           const type_info& cl_name)
{
    return registered_user_flags::get()->make_register(conf, cl_name.name());
};

static inline 
char tolower_impl(char c)
{
    return (char)::tolower(c);
};

const size_t* registered_user_flags::make_register(const std::shared_ptr<user_flag_config>& conf,
                                            const char* cl_name)
{
    //make_register is called during dynamic phase of static object initialization
    std::string name = conf->tag();
    std::transform(name.begin(), name.end(), name.begin(), &tolower_impl);
    auto pos = m_string_map.find(name);

    std::string class_name = cl_name? std::string(cl_name) : std::string();

    if (pos != m_string_map.end())
    {
        if (pos->second.second == class_name)
            return &pos->second.first;

        std::ostringstream msg;
        msg << "user-defined structure " << name << " is already defined";
        throw std::runtime_error(msg.str());
    }

    size_t new_code = m_string_map.size();
    auto pos2 = m_string_map.insert(string_map::value_type(name,code_name(new_code+1,class_name)));
    m_flag_names.push_back(name);
    m_configs.push_back(conf);

    const size_t* code = &pos2.first->second.first;
    register_visitors(code, conf);

    return code;
};

void registered_user_flags::register_visitors(const size_t* code, const config& conf)
{
    if (conf->register_visit_conj())        register_visitor(code, op_conj);
    if (conf->register_visit_abs())         register_visitor(code, op_abs);
    if (conf->register_visit_real())        register_visitor(code, op_real);
    if (conf->register_visit_uminus())      register_visitor(code, op_uminus);
    if (conf->register_visit_ctrans())      register_visitor(code, op_ctrans);
    if (conf->register_visit_trans())       register_visitor(code, op_trans);
    if (conf->register_visit_scal())        register_visitor(code, op_scal);
    if (conf->register_visit_self_mult())   register_visitor(code, op_self_mult);
    if (conf->register_visit_kron())        register_visitor(code, op_kron);
    if (conf->register_visit_mmul())        register_visitor(code, op_mult);
    if (conf->register_visit_plus())        register_visitor(code, op_plus);
    if (conf->register_visit_minus())       register_visitor(code, op_minus);
};

void registered_user_flags::register_visitor(const size_t* code, reg_ops op)
{
    m_visit[op].push_back(code);
};

size_t registered_user_flags::number_registered_flags() const
{
    return m_configs.size();
};

user_flag registered_user_flags::flag_from_string(const std::string& str) const
{
    auto pos = m_string_map.find(str);
    if (pos == m_string_map.end())
        return user_flag();

    const size_t* code = &pos->second.first;
    return user_flag(code);
};

std::string registered_user_flags::to_string_one(size_t struct_code)
{
    if (struct_code == 0 || struct_code > number_registered_flags())
        return "";
    return m_flag_names[struct_code - 1];
}

const registered_user_flags::config& registered_user_flags::get_user_flag_config(size_t struct_code)
{
    return m_configs[struct_code - 1];
};

void registered_user_flags::visit_conj_impl(struct_flag mat_flags, struct_flag& flags)
{
    Integer N;
    const size_t* const* reg_flags;
    get_registered_vis(N, reg_flags, op_conj);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_conj(mat_flags);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_abs_impl(struct_flag mat_flags, bool is_square, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_abs);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_abs(mat_flags, is_square);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_real_impl(struct_flag mat_flags, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_real);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_real(mat_flags);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_uminus_impl(struct_flag mat_flags, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_uminus);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_uminus(mat_flags);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_trans_impl(struct_flag mat_flags, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_trans);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_trans(mat_flags);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_ctrans_impl(struct_flag mat_flags, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_ctrans);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_ctrans(mat_flags);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_scal_impl(struct_flag mat_flags, value_struct_class vc, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_scal);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_scal(mat_flags, vc);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
}

void registered_user_flags::visit_self_mult_impl(struct_flag mat_flags, trans_type t1, trans_type t2, 
                                                 value_code vc, bool is_square, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_self_mult);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_self_mult(mat_flags, t1, t2, vc, is_square);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
};

void registered_user_flags::visit_kron_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_kron);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_kron(sf_X, sf_Y);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
};

void registered_user_flags::visit_mmul_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_mult);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_mmul(sf_X, sf_Y);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
};

void registered_user_flags::visit_plus_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_plus);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_plus(sf_X, sf_Y);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
};

void registered_user_flags::visit_minus_impl(struct_flag sf_X, struct_flag sf_Y, struct_flag& flags)
{
    Integer N;
    const size_t* const * reg_flags;
    get_registered_vis(N, reg_flags, op_minus);

    for (Integer i = 0; i < N; ++i)
    {
        const size_t* code  = reg_flags[i];
        bool set            = get_user_flag_config(*code)->visit_minus(sf_X, sf_Y);

        if (set == true)
        {
            user_flag uf(code);
            flags.set_user(uf);
        };
    };
};

void registered_user_flags::check_struct_one(size_t struct_code, const matcl::Matrix& mat)
{
    bool valid = get()->get_user_flag_config(struct_code)->test(mat);

    if (valid == false)
        throw error::invalid_struct(mat.get_struct().to_string());
};

user_flag registered_user_flags::get_conj_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->conj(sf);
};

user_flag registered_user_flags::get_abs_one(size_t struct_code, struct_flag sf, bool is_square)
{
    return get()->get_user_flag_config(struct_code)->abs(sf, is_square);
};

user_flag registered_user_flags::get_uminus_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->uminus(sf);
};

user_flag registered_user_flags::get_real_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->real(sf);
};

user_flag registered_user_flags::get_trans_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->trans(sf);
};

user_flag registered_user_flags::get_ctrans_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->ctrans(sf);
};

user_flag registered_user_flags::get_warnings_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->precision_lost(sf);
};

user_flag registered_user_flags::precision_increased_one(size_t struct_code, struct_flag sf)
{
    return get()->get_user_flag_config(struct_code)->precision_increased(sf);
};

user_flag registered_user_flags::get_mult_scal(size_t struct_code, struct_flag sf, value_struct_class vc)
{
    return get()->get_user_flag_config(struct_code)->scal(sf,vc);
};

user_flag registered_user_flags::get_seft_mult(size_t struct_code, struct_flag sf, trans_type t1, 
                                               trans_type t2, value_code vc, bool is_square)
{
    return get()->get_user_flag_config(struct_code)->seft_mult(sf, t1, t2, vc, is_square);
};

user_flag registered_user_flags::get_kron_one_both(size_t struct_code, struct_flag sf_left, struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->kron_both(sf_left, sf_right);
};

user_flag registered_user_flags::get_kron_one_1(size_t struct_code, struct_flag sf_left, struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->kron_1(sf_left, sf_right);
};

user_flag registered_user_flags::get_kron_one_2(size_t struct_code, struct_flag sf_left, struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->kron_2(sf_left, sf_right);
};

user_flag registered_user_flags::get_mult_one_both(size_t struct_code, struct_flag sf_left, 
                                                   struct_flag sf_right, bool is_square)
{
    return get()->get_user_flag_config(struct_code)->mult_both(sf_left, sf_right, is_square);
};

user_flag registered_user_flags::get_mult_one_1(size_t struct_code, struct_flag sf_left, 
                                                struct_flag sf_right, bool is_square)
{
    return get()->get_user_flag_config(struct_code)->mult_1(sf_left, sf_right, is_square);
};

user_flag registered_user_flags::get_mult_one_2(size_t struct_code, struct_flag sf_left, 
                                               struct_flag sf_right, bool is_square)
{
    return get()->get_user_flag_config(struct_code)->mult_2(sf_left, sf_right, is_square);
};

user_flag registered_user_flags::get_plus_one_both(size_t struct_code, struct_flag sf_left, 
                                                   struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->plus_both(sf_left, sf_right);
};

user_flag registered_user_flags::get_plus_one_1(size_t struct_code, struct_flag sf_left, 
                                                struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->plus_1(sf_left, sf_right);
};

user_flag registered_user_flags::get_plus_one_2(size_t struct_code, struct_flag sf_left, 
                                               struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->plus_2(sf_left, sf_right);
};

user_flag registered_user_flags::get_minus_one_both(size_t struct_code, struct_flag sf_left, 
                                                   struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->minus_both(sf_left, sf_right);
};

user_flag registered_user_flags::get_minus_one_1(size_t struct_code, struct_flag sf_left, 
                                                struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->minus_1(sf_left, sf_right);
};

user_flag registered_user_flags::get_minus_one_2(size_t struct_code, struct_flag sf_left, 
                                               struct_flag sf_right)
{
    return get()->get_user_flag_config(struct_code)->minus_2(sf_left, sf_right);
};

//------------------------------------------------------------------------
//                  user_struct_flags
//------------------------------------------------------------------------
void user_struct_flags::check_struct(size_t code, const matcl::Matrix& mat)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
            registered_user_flags::check_struct_one(flag, mat);

        ++flag;
        code = code/2;
    };
};

bool user_struct_flags::test_flags(size_t code)
{
    if (code == 0)
        return true;

    size_t last_flag = bit_manip<size_t>::most_significant_bit_pos(code) + 1;

    if (last_flag > registered_user_flags::get()->number_registered_flags())
        return false;
    else
        return true;
};

void user_struct_flags::to_string(size_t code, std::ostringstream& str, bool& add)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            if (add)
                str << ", ";

            str << registered_user_flags::get()->to_string_one(flag);
            add = true;
        };

        ++flag;
        code = code/2;
    };
}

void user_struct_flags::save_as_string(size_t code, std::ostringstream& str, bool& add)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {

            if (add)
                str << "_";

            str << registered_user_flags::get()->to_string_one(flag);
            add = true;
        };

        ++flag;
        code = code/2;
    };
};

bool user_struct_flags::load_from_string(struct_flag& flags, const std::string& struct_name)
{
    user_flag uf = registered_user_flags::get()->flag_from_string(struct_name);
    flags.set_user(uf);

    return (uf.is_empty() == false);
};

void user_struct_flags::get_abs(size_t code, struct_flag mat_sf, bool is_square, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_abs_one(flag,mat_sf, is_square);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_abs_impl(mat_sf, is_square, flags);
};

void user_struct_flags::get_uminus(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_uminus_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_uminus_impl(sf_mat, flags);
};

void user_struct_flags::get_conj(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_conj_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_conj_impl(sf_mat, flags);
};

void user_struct_flags::get_real(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_real_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_real_impl(sf_mat, flags);
};

void user_struct_flags::get_trans(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_trans_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_trans_impl(sf_mat, flags);
};

void user_struct_flags::get_ctrans(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_ctrans_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };

    registered_user_flags::get()->visit_ctrans_impl(sf_mat, flags);
};

void user_struct_flags::get_warnings(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::get_warnings_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };
};

void user_struct_flags::precision_increased(size_t code, struct_flag sf_mat, struct_flag& flags)
{
    size_t flag = 1;

    while(code != 0)
    {
        if (code % 2 != 0)
        {
            user_flag uf = registered_user_flags::precision_increased_one(flag,sf_mat);
            flags.set_user(uf);
        };

        ++flag;
        code = code/2;
    };
};

void user_struct_flags::get_seft_mult(size_t f1, struct_flag sf_mat, trans_type t1, trans_type t2,
                                      value_code vc, bool is_square, struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0)
    {
        bool has_1  = (f1 % 2 != 0);

        if (has_1)
        {
            user_flag uf = registered_user_flags::get_seft_mult(flag,sf_mat, t1, t2, vc, is_square);
            flags.set_user(uf);
        };

        ++flag;
        f1 = f1/2;
    };

    registered_user_flags::get()->visit_self_mult_impl(sf_mat, t1, t2, vc, is_square, flags);
};

void user_struct_flags::get_mult_scal(size_t f1, struct_flag mat_sf, value_struct_class vc, struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0)
    {
        bool has_1  = (f1 % 2 != 0);

        if (has_1)
        {
            user_flag uf = registered_user_flags::get_mult_scal(flag, mat_sf,vc);
            flags.set_user(uf);
        };

        ++flag;
        f1 = f1/2;
    };

    registered_user_flags::get()->visit_scal_impl(mat_sf, vc, flags);
};

void user_struct_flags::get_kron(size_t f1, size_t f2, struct_flag sf_l, struct_flag sf_r,
                 struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0 || f2 != 0)
    {
        bool has_1  = (f1 % 2 != 0);
        bool has_2  = (f2 % 2 != 0);

        if (has_1 || has_2)
        {
            if (has_1 && has_2)
            {
                user_flag uf = registered_user_flags::get_kron_one_both(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_1)
            {
                user_flag uf = registered_user_flags::get_kron_one_1(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_2)
            {
                user_flag uf = registered_user_flags::get_kron_one_2(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
        };

        ++flag;
        f1 = f1/2;
        f2 = f2/2;
    };

    registered_user_flags::get()->visit_kron_impl(sf_l, sf_r, flags);
};

void user_struct_flags::get_mult(size_t f1, size_t f2, struct_flag sf_l, struct_flag sf_r, bool is_square,
                                 struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0 || f2 != 0)
    {
        bool has_1  = (f1 % 2 != 0);
        bool has_2  = (f2 % 2 != 0);

        if (has_1 || has_2)
        {
            if (has_1 && has_2)
            {
                user_flag uf = registered_user_flags::get_mult_one_both(flag, sf_l, sf_r, is_square);
                flags.set_user(uf);
            }
            else if (has_1)
            {
                user_flag uf = registered_user_flags::get_mult_one_1(flag, sf_l, sf_r, is_square);
                flags.set_user(uf);
            }
            else if (has_2)
            {
                user_flag uf = registered_user_flags::get_mult_one_2(flag, sf_l, sf_r, is_square);
                flags.set_user(uf);
            }
        };

        ++flag;
        f1 = f1/2;
        f2 = f2/2;
    };

    registered_user_flags::get()->visit_mmul_impl(sf_l, sf_r, flags);
};

void user_struct_flags::get_plus(size_t f1, size_t f2, struct_flag sf_l, struct_flag sf_r,
                                 struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0 || f2 != 0)
    {
        bool has_1  = (f1 % 2 != 0);
        bool has_2  = (f2 % 2 != 0);

        if (has_1 || has_2)
        {
            if (has_1 && has_2)
            {
                user_flag uf = registered_user_flags::get_plus_one_both(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_1)
            {
                user_flag uf = registered_user_flags::get_plus_one_1(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_2)
            {
                user_flag uf = registered_user_flags::get_plus_one_2(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
        };

        ++flag;
        f1 = f1/2;
        f2 = f2/2;
    };

    registered_user_flags::get()->visit_plus_impl(sf_l, sf_r, flags);
};

void user_struct_flags::get_minus(size_t f1, size_t f2, struct_flag sf_l, struct_flag sf_r,
                                 struct_flag& flags)
{
    size_t flag = 1;

    while(f1 != 0 || f2 != 0)
    {
        bool has_1  = (f1 % 2 != 0);
        bool has_2  = (f2 % 2 != 0);

        if (has_1 || has_2)
        {
            if (has_1 && has_2)
            {
                user_flag uf = registered_user_flags::get_minus_one_both(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_1)
            {
                user_flag uf = registered_user_flags::get_minus_one_1(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
            else if (has_2)
            {
                user_flag uf = registered_user_flags::get_minus_one_2(flag, sf_l, sf_r);
                flags.set_user(uf);
            }
        };

        ++flag;
        f1 = f1/2;
        f2 = f2/2;
    };

    registered_user_flags::get()->visit_minus_impl(sf_l, sf_r, flags);
};

//------------------------------------------------------------------------
//                  predefined_struct
//------------------------------------------------------------------------
value_struct_class predefined_struct::get_value_type(Integer val, bool test_zero)
{
    //if test_zero == false, then we are assigning to Object;
    //resulting value is unknown
    if (test_zero == false)             return value_struct_class::vc_general;

    if (mrd::is_zero(val))              return value_struct_class::vc_zero;
    else if (mrd::is_one(val))          return value_struct_class::vc_one;
    else if (val < 0)                   return value_struct_class::vc_neg_real;
    else if (val > 0)                   return value_struct_class::vc_pos_real;
    else                                return value_struct_class::vc_general;
};

value_struct_class predefined_struct::get_value_type(Real val, bool test_zero)
{
    //if test_zero == false, then we are assigning to Object;
    //resulting value is unknown
    if (test_zero == false)             return value_struct_class::vc_general;

    if (mrd::is_zero(val))              return value_struct_class::vc_zero;
    else if (mrd::is_one(val))          return value_struct_class::vc_one;
    else if (!mrd::isfinite_helper<Real>::eval(val))
                                        return value_struct_class::vc_real;
    else if (val < 0)                   return value_struct_class::vc_neg_real;
    else if (val > 0)                   return value_struct_class::vc_pos_real;
    else                                return value_struct_class::vc_general;
};

value_struct_class predefined_struct::get_value_type(Float val, bool test_zero)
{
    //if test_zero == false, then we are assigning to Object;
    //resulting value is unknown
    if (test_zero == false)             return value_struct_class::vc_general;

    if (mrd::is_zero(val))              return value_struct_class::vc_zero;
    else if (mrd::is_one(val))          return value_struct_class::vc_one;
    else if (!mrd::isfinite_helper<Float>::eval(val))
                                        return value_struct_class::vc_real;
    else if (val < 0)                   return value_struct_class::vc_neg_real;
    else if (val > 0)                   return value_struct_class::vc_pos_real;
    else                                return value_struct_class::vc_general;
};

value_struct_class predefined_struct::get_value_type(Complex val, bool test_zero)
{
    //if test_zero == false, then we are assigning to Object;
    //resulting value is unknown
    if (test_zero == false)             return value_struct_class::vc_general;
    if (mrd::is_zero(imag(val)))        return get_value_type(real(val), test_zero);

    return value_struct_class::vc_general;
};

value_struct_class predefined_struct::get_value_type(Float_complex val, bool test_zero)
{
    //if test_zero == false, then we are assigning to Object;
    //resulting value is unknown
    if (test_zero == false)             return value_struct_class::vc_general;
    if (mrd::is_zero(imag(val)))        return get_value_type(real(val), test_zero);

    return value_struct_class::vc_general;
};

value_struct_class predefined_struct::get_value_type(const Object& val, bool test_zero)
{
    if (test_zero && mrd::is_zero(val)) return value_struct_class::vc_zero;
    else                                return value_struct_class::vc_general;
};

struct_flag predefined_struct::get_set_diag(const struct_flag& sf, Integer d, value_struct_class vt,
                                            bool is_real_mat, bool is_square)
{
    bool is_zero    = (vt == value_struct_class::vc_zero);
    bool is_one     = (vt == value_struct_class::vc_one);
    bool is_gen     = (vt == value_struct_class::vc_general);

    struct_flag out = sf;

    out.reset_user();

    if (out.is_id())
    {
        out         = predefined_struct_type::diag;

        if (is_gen == true)
            out.set_sym(true);
        else if (is_real_mat == true)
            out.set_sym(true);
        else
            out.set_her(true);
    };

    if (d == 0)
    {
        if (out.is_diag() == true && is_one == true && is_square == true)
            out.add(predefined_struct_type::id);

        if (out.has_her_flag() == true)
        {
            if (vt == value_struct_class::vc_general)
                out.set_her(false);
        };
    }
    else if (d > 0)
    {
        if (is_zero == false)
            out.set_udiags(struct_flag::diag_type::general);

        out.set_sym(false);
        out.set_her(false);
    }
    else
    {
        if (is_zero == false)
            out.set_ldiags(struct_flag::diag_type::general);

        out.set_sym(false);
        out.set_her(false);
    };    

    return out;
};

struct_flag predefined_struct::get_resize(const struct_flag& sf, bool is_sym, bool is_real_mat)
{
    struct_flag out = sf;
    
    out.reset_user();

    if (out.is_id())
    {
        out = predefined_struct_type::diag;

        if (is_sym == true)
        {
            if (is_real_mat)
                out.set_sym(true);
            else
                out.set_her(true);
        };

        return out;
    };

    if (is_sym == false)
    {
        out.set_sym(false);
        out.set_her(false);
    };

    return out;
};

struct_flag predefined_struct::get_rectangle_view(const struct_flag& sf, bool is_square)
{
    struct_flag out = sf;

    out.reset_user();

    if (out.is_id())
    {
        if (is_square)
            return out;
        else
        {
            out = predefined_struct_type::diag;
            return out;
        };
    };

    if (is_square == false)
    {
        out.set_sym(false);
        out.set_her(false);
    };

    return out;
};

struct_flag predefined_struct::get_qtril(const struct_flag& sf)
{
    struct_flag out = sf;    

    if (out.is_qtril() == false)
    {
        out.set_udiags(std::max(out.get_udiags(), struct_flag::diag_type::qtriang));

        out.set_sym(false);
        out.set_her(false);
        out.reset_user();
    };

    return out;
}

struct_flag predefined_struct::get_qtriu(const struct_flag& sf)
{
    struct_flag out = sf;    

    if (out.is_qtriu() == false)
    {
        out.set_ldiags(std::max(out.get_ldiags(), struct_flag::diag_type::qtriang));

        out.set_sym(false);
        out.set_her(false);
        out.reset_user();
    };

    return out;
}

struct_flag predefined_struct::get_tril(const struct_flag& sf, Integer d, bool is_real_mat)
{
    struct_flag out = sf;    

    if (out.is_tril() == false)
    {
        if (d > 1)
        {
        }
        else if (d == 1)
        {
            out.set_udiags(std::max(struct_flag::diag_type::one, out.get_udiags()));
        }
        else if (d <= 0)
        {
            out.set_udiags(struct_flag::diag_type::zero);
        }

        out.set_sym(false);
        out.set_her(false);
        out.reset_user();
    }
    else if (sf.is_id())
    {
        if (d < 0)
        {
            out = predefined_struct_type::diag;

            if (is_real_mat)
                out.set_sym(true);
            else
                out.set_her(true);

            out.reset_user();
        };
    }
    else
    {
        if (d < 0)
        {
            out.set_sym(false);
            out.set_her(false);
            out.reset_user();
        };
    };

    return out;
};

struct_flag predefined_struct::get_triu(const struct_flag& sf, Integer d, bool is_real_mat)
{
    struct_flag out = sf;    

    if (out.is_triu() == false)
    {
        if (d < -1)
        {
        }
        else if (d == -1)
        {
            out.set_ldiags(std::max(struct_flag::diag_type::one, out.get_ldiags()));
        }
        else if (d >= 0)
        {
            out.set_ldiags(struct_flag::diag_type::zero);
        }

        out.set_sym(false);
        out.set_her(false);
        out.reset_user();
    }
    else if (sf.is_id())
    {
        if (d > 0)
        {
            out = predefined_struct_type::diag;

            if (is_real_mat)
                out.set_sym(true);
            else
                out.set_her(true);

            out.reset_user();
        };
    }
    else
    {
        if (d > 0)
        {
            out.set_sym(false);
            out.set_her(false);
            out.reset_user();
        };
    };

    return out;
};

struct_flag predefined_struct::get_abs(const struct_flag& sf, bool is_square)
{
    struct_flag out     = sf;    

    size_t code_user    = sf.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_abs(code_user, sf, is_square, sfu);

    out.reset_user();
    out.add_user(sfu);

    if (out.has_her_flag() == true)
        out.set_sym(true);

    return out;
}

struct_flag predefined_struct::get_conj(const struct_flag& sf)
{
    struct_flag out     = sf;    

    size_t code_user    = out.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_conj(code_user, sf, sfu);

    out.reset_user();
    out.add_user(sfu);

    return out;
};

struct_flag predefined_struct::get_real(const struct_flag& sf)
{
    struct_flag out = sf;    

    size_t code_user    = out.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_real(code_user, sf, sfu);

    out.reset_user();
    out.add_user(sfu);

    if (out.has_her_flag() == true)
        out.set_sym(true);

    return out;
}

struct_flag predefined_struct::get_trans(const struct_flag& sf, trans_type t)
{
    switch(t)
    {
        case trans_type::no_trans:
            return sf;
        case trans_type::trans:
            return get_trans(sf);
        case trans_type::conj_trans:
            return get_ctrans(sf);
        default:
            //impossible case
            return sf;
    }
}

struct_flag predefined_struct::get_trans(const struct_flag& sf)
{
    struct_flag out = sf;

    struct_flag::diag_type ld   = out.get_ldiags();
    struct_flag::diag_type ud   = out.get_udiags();

    out.set_ldiags(ud);
    out.set_udiags(ld);

    //trans(sym) -> sym
    //trans(her) -> her

    size_t code_user    = out.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_trans(code_user,sf, sfu);

    out.reset_user();
    out.add_user(sfu);

    return out;
};

struct_flag predefined_struct::get_ctrans(const struct_flag& sf)
{
    struct_flag out = sf;

    struct_flag::diag_type ld   = out.get_ldiags();
    struct_flag::diag_type ud   = out.get_udiags();

    out.set_ldiags(ud);
    out.set_udiags(ld);

    //ctrans(sym) -> sym
    //ctrans(her) -> her

    size_t code_user    = out.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_ctrans(code_user,sf, sfu);

    out.reset_user();
    out.add_user(sfu);

    return out;
};

struct_flag predefined_struct::precision_increased(struct_flag f1)
{
    size_t code_user    = f1.to_int_user();
    struct_flag sfu;
    user_struct_flags::precision_increased(code_user,f1,sfu);

    f1.reset_user();
    f1.add_user(sfu);

    return f1;
};

struct_flag predefined_struct::set_warnings(const struct_flag& s1, bool prec_lost, bool compl_to_real)
{
    struct_flag out = compl_to_real == false? s1 : get_real(s1);

    if (prec_lost)
    {
        size_t code_user    = out.to_int_user();
        struct_flag sfu;
        user_struct_flags::get_warnings(code_user,s1,sfu);

        out.reset_user();
        out.add_user(sfu);
    };

    return out;
}

struct_flag predefined_struct::eval_struct(struct_flag f1, bool is_zero_id)
{
    struct_flag out;

    if (is_zero_id)
    {    
        if (out.is_id())
            out = predefined_struct_type::diag;
    }
    else
    {
        out.set_ldiags(struct_flag::diag_type::general);
        out.set_udiags(struct_flag::diag_type::general);
    };

    out.set_her(false);
    out.reset_user();

    //for a symmetric matrix X f(X) is also symmetric
    if (f1.has_sym_flag())
        out.set_sym(true);

    return out;
};

struct_flag predefined_struct::uminus_cont(const struct_flag& sf_ret, const struct_flag& sf_mat)
{
    struct_flag ret = sf_ret;
    ret.set_her(sf_mat.has_her_flag());

    size_t code_user    = sf_mat.to_int_user();
    struct_flag sfu;
    user_struct_flags::get_uminus(code_user,sf_mat,sfu);

    ret.reset_user();
    ret.add_user(sfu);

    return ret;
};

struct_flag predefined_struct::inv_cont(const struct_flag& sf_ret, const struct_flag& sf_mat)
{
    struct_flag ret = sf_ret;
    ret.set_her(sf_mat.has_her_flag());
    return ret;
};

};};
