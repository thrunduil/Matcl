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
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/details/options/matcl_options_details.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/options/option_validators.h"
#include "matcl-core/memory/alloc.h"

#include <map>
#include <type_traits>
#include <vector>

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients

namespace matcl 
{

//--------------------------------------------------------------------------
//                              OPTION
//--------------------------------------------------------------------------
// visitor class for inspecting options' properties
class MATCL_CORE_EXPORT option_visitor
{
    public:
        virtual ~option_visitor(){};

        // visit option, which holds one of recognized type, opt argument
        // represents the option currently visited
        virtual void    visit(const optional<Integer>& val, const option& opt) = 0;
        virtual void    visit(const optional<Real>& val, const option& opt) = 0;
        virtual void    visit(const optional<Complex>& val, const option& opt) = 0;
        virtual void    visit(const optional<std::string>& val, const option& opt) = 0;
        virtual void    visit(const optional<bool>& val, const option& opt) = 0;

        virtual void    visit(const optional<std::vector<Integer>>& val, const option& opt) = 0;
        virtual void    visit(const optional<std::vector<Real>>& val, const option& opt) = 0;
        virtual void    visit(const optional<std::vector<Complex>>& val, const option& opt) = 0;
        virtual void    visit(const optional<std::vector<std::string>>& val, const option& opt) = 0;
        virtual void    visit(const optional<std::vector<bool>>& val, const option& opt) = 0;
};

//base class for concrete options; this class should never be used directly.
class MATCL_CORE_EXPORT option : public matcl_new_delete
{
    private:
        using impl_type     = std::shared_ptr<details::option_impl>;

    private:
        impl_type           m_impl;

    protected:
        option() = delete;
        option(const impl_type&);

        static std::string  make_name(const std::string& type_name);

    public:
        ~option();

        // name of the option.
        std::string         name() const;
        
        // description of the options
        std::string         description() const;

        // get option value. Function returns optional value, which 
        // contains option value of the option is present.
        template <class T>
        optional<T>         get() const;

        // return true if value of this option is set.
        bool                has_value() const;

        // visitor is entering given option
        void                accept(option_visitor& vi) const;

        // return name of value type stored in given option
        std::string         get_option_type() const;
        
        // display information about this option
        void                help(const disp_stream_ptr& ds,const options& disp_options) const;

    private:
        void                error_invalid_option_type(const std::string& get_type) const;
};

// internal use only
template<class T>
class option_type : public option
{
    private:
        using option_impl_type  = details::option_impl_type<T>;
        using impl_type         = std::shared_ptr<option_impl_type>;
        using base_type         = option;

    protected:
        option_type(const std::string& name, const std::string& descr,
                    const optional<T>& val);

        ~option_type();
};

// base class for user defined options. User defined options must implement a
// function static void config(); which sets value of all public
// static members of this class. Template arguments T: type of value stored in
// given class, Derived: type of option class.
// allowed types T: bool, Integer, Real, Complex, std::string, and
// std::vector of these types
template<class T, class Derived>
class option_base : public option_type<T>
{    
    private:
        using base_type         = option_type<T>;
        using mutex_type        = matcl::default_spinlock_mutex;

    public:
        // optional type
        using opt_type          = optional<T>;

        // type of stored value
        using value_type        = T;
        
        // type of validator function, function may throw exception
        // or correct supplied value
        using validator_type    = std::function<opt_type (opt_type)>;        

    private:
        static bool             m_initialized;
        static mutex_type       m_mutex;
        static std::string      m_name;

    public:        
        // these static members must be initialized in Derived::config() function
        // Default value of given option
        static T                m_default_value;
        
        // short description of given option
        static std::string      m_description;
        
        // a function, that checks if option parameter is valid. If not set, then
        // no checks are performed
        static validator_type   m_validator;        

    private:
        static Derived          m_default_opt;

    protected:
        explicit option_base(const opt_type& val = opt_type());

        ~option_base();

    private:
        static void             initialize();        
        static std::string      get_name();
        static std::string      get_description();
        static validator_type   get_validator();  
};

// option storing optional value of type T identified by name
// allowed types T: bool, Integer, Real, Complex, std::string, and
// std::vector of these types
template<class T>
class string_option : public option_type<T>
{    
    private:
        using base_type         = option_type<T>;

    public:
        // create option of name opt with undefined value
        string_option(const std::string& opt);

        // create option of name opt with given value
        string_option(const std::string& opt, const optional<T>& value);
        string_option(const std::string& opt, const T& value);

        // create option of name opt with given value and description
        // description is used only for options printing purpose
        string_option(const std::string& opt, const optional<T>& value, 
                const std::string& descr);
        string_option(const std::string& opt, const T& value, 
                const std::string& descr);
};

//--------------------------------------------------------------------------
//                              OPTIONS
//--------------------------------------------------------------------------

// helper class to inspect tested options
class MATCL_CORE_EXPORT options_notifier
{
    public:
        virtual ~options_notifier(){};

        // option with given name is tested
        virtual void    report(const std::string& name) = 0;
};

// class to handle a set of arbitrary options with values. There are three 
// levels of options: local user defined options, global user defined default
// options, and global predefined options, that cannot change;
// options implements copy-on-write type behavior, i.e. modification of non
// unique object will create a copy
class MATCL_CORE_EXPORT options
{
    public:
        using notifier_ptr  = std::shared_ptr<options_notifier>;

    private:
        using impl_type     = std::shared_ptr<details::options_impl>;

    private:
        impl_type           m_options;

    public: 
        // create empty set of options.
        options();
        ~options();

        // create options set from initializer list
        options(std::initializer_list<option> options_vec);        
        explicit options(const std::vector<option>& options_vec);

        // convert on option to option set
        options(const option& one_option);

        // standard copy and move constructors
        options(const options& other);
        options(options&& other);

        // standard assign and move assign operators
        options&            operator=(const options& other) &;
        options&            operator=(options&& other) &;

        // get option value. Function returns optional value, which 
        // contains concrete value of the option is present, value stored
        // in passed argument if present, or undefined value otherwise.
        template <class T>
        optional<T>         get(const option& option_type) const;

        // return true if given option is set
        bool                has_option(const option& option_type) const;

        // get option value. If option is not present in given set, value stored
        // in passed argument is taken, then user default value is taken. 
        // If user default is not set, then predefined value is returned.
        template <class T>
        T                   get_option(const option& option_type) const;
        
        // get default option value if is set, then check if value is stored in 
        // passed argument, finally predefined value is taken.
        template <class T>
        static T            get_default(const option& option_type);
        
        // get predefine option value
        template <class T>
        static T            get_predefined(const option& option_type);

        // set value of existing option or add a new one. 
        options&            set(const option& other);

        // set options stored in other object, replacing options if they exist
        // in this set.
        options&            set(const options& other);

        // remove option opt
        options&            remove(const option& opt);

        // remove all options existing in other options set
        options&            remove(const options& other);

        // return default option set
        static options&     default_options();

        // remove all options from options set
        void                clear();
        
        // number of options in options set
        Integer             size() const;

        // make this object unique
        void                make_unique();

        // return true is this object is unique
        bool                is_unique() const;

        // display options stored in given object
        void                disp(const disp_stream_ptr& ds = default_disp_stream(),
                                 const options& disp_options = options()) const;
        // display information about all available options
        static void         help(const disp_stream_ptr& ds = default_disp_stream(),
                                 const options& disp_options = options());

        // install options_notifier to check options that was tested
        void                install_notifier(const notifier_ptr& notif);

    private:
        template<class Iterator>
        void                init(Iterator beg, Iterator end);
};

// display options stored in given object
MATCL_CORE_EXPORT void  disp(const options& opts, const disp_stream_ptr& ds = default_disp_stream(),
                            const options& print_options = options());

// display information about all available options
MATCL_CORE_EXPORT void  help(const options& opts, const disp_stream_ptr& ds = default_disp_stream(),
                            const options& print_options = options());

// display information about this option
MATCL_CORE_EXPORT void  help(const option& opts, const disp_stream_ptr& ds = default_disp_stream(),
                            const options& print_options = options());
};

#pragma warning( pop )

#include "matcl-core/details/options/matcl_options_impl.inl"