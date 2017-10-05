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
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/options/optional.h"
#include "matcl-core/memory/alloc.h"

#include <ostream>
#include <stdexcept>
#include <map>
#include <memory>

namespace matcl
{
    class option;
    class option_visitor;

    template<class T, class Derived>
    class option_base;
};

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients

namespace matcl { namespace details
{

class MATCL_CORE_EXPORT option_impl : public std::enable_shared_from_this<option_impl>
{
    public:
        virtual ~option_impl(){};

        virtual std::string name() const = 0;
        virtual std::string description() const = 0;

        virtual void        accept(option_visitor& vis, const option& opt) const = 0;
        virtual bool        has_value() const = 0;
        virtual std::string get_option_type() const = 0;
        void                help(const disp_stream_ptr& ds, const option& opt,
                                const options& disp_options) const;

    public:
        static std::string  pretty_type_name(const std::string& type_name);
};

template<class T>
class option_impl_type : public option_impl
{
    private:
        using impl_type     = optional<T>;

    private:
        impl_type           m_value;
        std::string         m_name;
        std::string         m_description;

    public:
        option_impl_type(const std::string& name, const std::string& descr,
                         const impl_type& val);

        virtual ~option_impl_type() override;

        option_impl_type(const option_impl_type&)             = delete;
        option_impl_type& operator=(const option_impl_type&)  = delete;

    public:
        virtual std::string name() const override                       { return m_name; }
        virtual std::string description() const override                { return m_description;}
        virtual bool        has_value() const override                  { return (bool)m_value; }
        impl_type           get() const                                 { return m_value; }

        virtual void        accept(option_visitor& vis, const option& opt) const override;
        virtual std::string get_option_type() const override;
};

class MATCL_CORE_EXPORT options_impl : public matcl_new_delete
{
    private:
        using option_map    = std::map<std::string, option>;
        using ptr_type      = std::shared_ptr<options_impl>;
        using notifier_ptr  = std::shared_ptr<options_notifier>;

    public:
        using mutex_type    = matcl::default_spinlock_mutex;

    private:
        option_map          m_options;
        notifier_ptr        m_notifier;
        mutable mutex_type  m_mutex_local;        

    public:
        options_impl();
        ~options_impl();

        void                set(const option& option_value);
        void                remove(const option& option_value);
        void                set(const options_impl& other);
        void                remove(const options_impl& other);
        void                disp(const disp_stream_ptr& ds, const options& disp_options) const;
        static void         help(const disp_stream_ptr& ds, const options& disp_options);
        ptr_type            copy() const;

        void                clear();
        Integer             size() const;

        template<class T>
        optional<T>         get(const option& option_value) const;

        template <class T>
        static T            get_predefined(const option& option_value);

        bool                has_option(const option& option_value) const;

        void                install_notifier(const notifier_ptr& notif);

    private:
        options_impl(const options_impl&) = delete;
        options_impl& operator=(const options_impl&) = delete;

        static void         error_unregistered_option(const std::string& opt_name);
        static void         register_option(const option& opt);
        static option_map&  get_options_predefined();
        static mutex_type&  get_mutex_global();

        template<class T, class Derived>
        friend class matcl::option_base;
};

}};

#pragma warning( pop )