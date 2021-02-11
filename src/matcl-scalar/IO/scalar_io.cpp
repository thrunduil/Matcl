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

#include "matcl-scalar/IO/scalar_io.h"
#include "matcl-core/details/IO/disp_stream_options.h"
#include "matcl-core/details/IO/printer.h"
#include "matcl-core/details/IO/io_impl.h"
#include "matcl-core/details/IO/disp_impl.h"

#include "matcl-core/details/IO/disp_stream_impl.h"
#include "matcl-dynamic/details/object.inl"
#include "matcl-dynamic/object_type.h"
#include "matcl-scalar/object.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/details/object_interface.h"

#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-core/IO/disp_data_provider.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

namespace mr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

namespace matcl { namespace raw { namespace details
{
    struct struct_scalar{};

    template<class V>
    struct matrix_provider_scalar : public matrix_provider_base<V>
    {
        const V&            m_mat;
        V                   m_zero;

        matrix_provider_scalar(const V& mat)
            :m_mat(mat), m_zero(0)
        {};

        virtual ~matrix_provider_scalar(){};

        virtual Integer rows() const override
        {
            return 1;
        };
        virtual Integer cols() const override
        {
            return 1;
        };

        virtual V get_value(const disp_stream* user, Integer width, align_type at, 
                           Integer r, Integer c, bool& is_zero) const override
        {
            (void)user;
            (void)width;
            (void)at;
            (void)r;
            (void)c;

            is_zero = matcl::is_zero(m_mat);
            return m_mat;
        };

        virtual bool is_symher() const override
        {
            return false;   
        };
    };

    template<>
    struct disp_matrix<std::string,struct_scalar>
    {
        static void eval(md::disp_stream_impl& os, const disp_stream* user, const std::string& str)
        {
            os.do_init_display(user, value_code::v_object);
            os.do_start_display(user);	
                
            os.do_displaying_string(user);

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            if (str.empty() == false)
            {
                std::istringstream iss(str);
                std::string line;

                while (std::getline(iss, line))
                {
                    os.do_disp_new_line(user,1);
                    os.get_printer().disp_elem(-os.do_get_max_matrix_width(user), line, align_type::left, 0);
                    os.do_disp_end_line(user,1);
                };
            };

            os.do_end_display(user);
        };
    };

    template<class V>
    struct disp_matrix<V,struct_scalar>
    {
        static void eval(md::disp_stream_impl& os, const disp_stream* user, const V& i)
        {
            matcl::value_code vt =  matrix_traits::value_code<V>::value;

            os.do_init_display(user,vt);
            os.do_start_display(user);	

            os.do_displaying_scalar(user,vt, dynamic::object_type<V>::get_static_type().to_string());

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            disp_mode dm = os.do_get_display_mode(user);

            switch(dm)
            {
                case disp_mode::standard:
                default:
                {
                    os.do_disp_new_line(user,1);
                    os.get_printer().disp_elem(-os.do_get_max_matrix_width(user),i,align_type::left,0);
                    os.do_disp_end_line(user,1);

                    os.do_end_display(user);
                    return;
                }
                case disp_mode::all_dense:
                case disp_mode::scalar_dense:
                {
                    matrix_provider_scalar<V> mdp(i);
                    return disp_matrix<V,struct_dense>::eval_matrix_body(os,user,mdp);
                }
            };
        };
    };

    template<>
    struct disp_matrix<Object,struct_scalar>
    {
        static void eval(md::disp_stream_impl& os, const disp_stream* user, const Object& i)
        {
            matcl::value_code vt    =  value_code::v_object;

            os.do_init_display(user, vt);
            os.do_start_display(user);	

            os.do_displaying_scalar(user,vt, i.get_type().to_string());

            if (os.do_display_data(user) == false)
            {
                os.do_end_display(user);
                return;
            };

            disp_mode dm = os.do_get_display_mode(user);

            switch(dm)
            {
                case disp_mode::standard:
                default:
                {
                    os.do_disp_new_line(user,1);
                    os.get_printer().disp_elem(-os.do_get_max_matrix_width(user), i,
                                               align_type::left, 0);
                    os.do_disp_end_line(user,1);

                    os.do_end_display(user);
                    return;
                }
                case disp_mode::all_dense:
                case disp_mode::scalar_dense:
                {
                    matrix_provider_scalar<Object> mdp(i);
                    return disp_matrix<Object,struct_dense>::eval_matrix_body(os,user,mdp);
                }
            };
        };
    };
}}};

namespace matcl { namespace raw
{

void disp(const disp_stream_ptr& os,Integer i)
{	
    mrd::disp_matrix<Integer, mrd::struct_scalar>::eval(*os->impl(),os.get(), i);
}

void disp(const disp_stream_ptr& os, Real val)
{
    mrd::disp_matrix<Real,mrd::struct_scalar>::eval(*os->impl(),os.get(), val);
};

void disp(const disp_stream_ptr& os, Float val)
{
    mrd::disp_matrix<Float, mrd::struct_scalar>::eval(*os->impl(),os.get(), val);
};

void disp(const disp_stream_ptr& os, const Complex& val)
{
    mrd::disp_matrix<Complex, mrd::struct_scalar>::eval(*os->impl(),os.get(), val);
};

void disp(const disp_stream_ptr& os, const Float_complex& val)
{
    mrd::disp_matrix<Float_complex, mrd::struct_scalar>::eval(*os->impl(),os.get(), val);
};

void disp(const disp_stream_ptr& os, const Object& c)
{
    md::matrix_object_interface oi(&c);

    if (oi.displayed_object_matrix(os) == true)
        return;

    mrd::disp_matrix<Object, mrd::struct_scalar>::eval(*os->impl(),os.get(), c);
};

void disp(const disp_stream_ptr& os, const char* str)
{
    mrd::disp_matrix<std::string, mrd::struct_scalar>::eval(*os->impl(),os.get(), str);
};

void disp(const disp_stream_ptr& os, const std::string& str)
{
    mrd::disp_matrix<std::string, mrd::struct_scalar>::eval(*os->impl(),os.get(), str);
};

}};

namespace matcl
{

void details::disp_impl(Integer val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

void details::disp_impl(const Real& val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

void details::disp_impl(const Complex& val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

void details::disp_impl(const Object& val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

void matcl::disp_header(const char* val, const disp_stream_ptr& os, const options& opts)
{
    options opts2 = opts;
    opts2.set(matcl::opt::disp::header_only(true));
    return disp(val, os, opts2);
};

void matcl::disp_header(const std::string& val, const disp_stream_ptr& os, const options& opts)
{
    options opts2 = opts;
    opts2.set(matcl::opt::disp::header_only(true));
    return disp(val, os, opts2);
};

void matcl::disp(const char* val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os, val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

void matcl::disp(const std::string& val, const disp_stream_ptr& os, const options& opts)
{
    if (opts.size() == 0)
        return raw::disp(os,val);
    else
    {
        return raw::disp(std::make_shared<disp_stream_options>(os, opts), val);
    }
};

};
