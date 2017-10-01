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
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/special_types.h"
#include "matcl-dynamic/object_type.h"

#include "matcl-dynamic/dynamic_function/dynamic_function.h"
#include "matcl-dynamic/details/data_constructor.inl"
#include "matcl-dynamic/details/func_evaler_conv.h"

#include <string>

namespace matcl { namespace dynamic { namespace details
{

struct MATCL_DYN_EXPORT unifier_tag { static function_name eval(); };
struct MATCL_DYN_EXPORT assign_tag  { static function_name eval(); };

struct enum_tag;

struct delayed_function_register
{
    virtual function            execute_registration() = 0;
    virtual function_name       get_fun_name() const = 0;
};

using make_return_fptr  = Type (*)(int n_template, const Type templates[],
                                           int n_args, const Type args[]);

struct delayed_function_template_register
{
    virtual function            execute_registration(std::vector<Type>& types, 
                                                     make_return_fptr& ret) = 0;
    virtual function_name       get_fun_name() const = 0;
};

struct MATCL_DYN_EXPORT fun_templ_register_helper
{
    fun_templ_register_helper(delayed_function_template_register* fr_int);
    void initialize();
};

struct MATCL_DYN_EXPORT fun_register_helper
{
    fun_register_helper(delayed_function_register* fr_int);
    void initialize();
};

template<class Fun>
class function_register : public delayed_function_register
{
    private:
        Fun                 m_fun;
        function_name       m_fun_name;        

    public:
        function_register(function_name fun_name, Fun f)
            : m_fun_name(fun_name), m_fun(f)
        {};

        static function         process_function(Fun fun);

        virtual function        execute_registration()      { return process_function(m_fun); };
        virtual function_name   get_fun_name() const        { return m_fun_name;    }
};

template<class Fun, class ... Templates>
class function_template_register : public delayed_function_template_register
{
    private:
        using type_vec      = std::vector<Type>;
        using ret_func      = make_return_fptr;

    private:
        Fun                 m_fun;
        function_name       m_fun_name;
        ret_func            m_ret;

    public:
        function_template_register(function_name fun_name, Fun f, ret_func ret)
            : m_fun_name(fun_name), m_fun(f), m_ret(ret)
        {};

        static function         process_function(Fun fun, type_vec& types);

        virtual function        execute_registration(type_vec& types, ret_func& ret);
        virtual function_name   get_fun_name() const        { return m_fun_name;    }
};

template<class Fun>
delayed_function_register* reg_fun(function_name fn, Fun f)
{
    return new function_register<Fun>(fn, f);
};

template<class Fun, class ... Templates>
delayed_function_template_register* reg_fun_templ(function_name fn, Fun f, 
                                                  make_return_fptr ret)
{
    return new function_template_register<Fun, Templates...>(fn, f, ret);
};

//add reference if required
template<class T> struct make_obj_reference_type
{
    using type = T;
};
template<class T> struct make_obj_reference_type<const T&>
{
    using T0    = typename std::decay<T>::type;
    static_assert(std::is_same<T0, T>::value, "unsupported type");
    using type = T;
};
template<class T> struct make_obj_reference_type<T&>
{
    using type = object_reference<T>;
};
template<class T> struct make_obj_reference_type<T&&>
{
    static_assert(matcl::details::dependent_false<T>::value, 
                  "rvalue reference not allowed");
};
template<class T> struct make_obj_reference_type<const T&&>
{
    static_assert(matcl::details::dependent_false<T>::value, 
                  "rvalue reference not allowed");
};

//helper for registering functions
template<class function_class, class Name, class Function_Type, Function_Type F>
struct make_function_type
{
    private:
        using this_type     = make_function_type;
        using function_type = Function_Type;

    private:
        static this_type    g_instance;
        fun_register_helper m_register;

        make_function_type();

    public:
        static this_type&   get()           { return g_instance; };
        void                initialize()    { return m_register.initialize(); };
};

template<class Derived, class Name, class Function_Type, Function_Type F>
make_function_type<Derived, Name, Function_Type,F> 
make_function_type<Derived, Name, Function_Type,F>::g_instance;

template<class Derived, class Name, class Function_type, Function_type F>
make_function_type<Derived, Name, Function_type,F>::make_function_type()
    :m_register(reg_fun<Function_type>(Name::eval(), F))
{};

//helper for registering template functions
template<class function_class, class Name, class Function_Type, Function_Type F,
        make_return_fptr Ret, class ... Templates>
struct make_function_template_type
{
    private:
        using this_type     = make_function_template_type;
        using function_type = Function_Type;
        using registerer    = fun_templ_register_helper;

    private:
        static this_type    g_instance;
        registerer          m_register;

        make_function_template_type();

    public:
        static this_type&   get()           { return g_instance; };
        void                initialize()    { return m_register.initialize(); };
};

template<class Derived, class Name, class Function_Type, Function_Type F, 
    make_return_fptr Ret, class ... Templ>
make_function_template_type<Derived, Name, Function_Type,F, Ret, Templ...> 
make_function_template_type<Derived, Name, Function_Type,F, Ret, Templ ...>::g_instance;

template<class Derived, class Name, class Function_type, Function_type F, 
    make_return_fptr Ret, class ... Templ>
make_function_template_type<Derived, Name, Function_type,F, Ret, Templ ...>
        ::make_function_template_type()
    :m_register(reg_fun_templ<Function_type, Templ...>(Name::eval(), F, Ret))
{};

template<int nargs, class args_list>
struct init_args
{
    static void init(Type* args_ti)
    {
        using obj_type_0                = typename get_elem<args_list, 0>::type;
        using obj_type_base             = typename std::decay<obj_type_0>::type;

        static const bool is_obj        = is_object<obj_type_base>::value;
        static_assert(is_obj == true, "argument must have object type");

        using obj_type                  = typename make_obj_reference_type<obj_type_0>::type;
        using next_type                 = typename tail<args_list>::type;

        args_ti[0]                      = obj_type::get_static_type();
        init_args<nargs - 1, next_type>::init(args_ti + 1);
    };
};

template<class args_list>
struct init_args<0, args_list>
{
    static void init(Type*){};
};

template<class T>
struct init_ret
{
    using obj_type  = typename std::decay<T>::type;

    static void init(Type& ret_ti)
    {    
        ret_ti = obj_type::get_static_type();
    };
};
template<>
struct init_ret<void>
{
    static void init(Type& ret_ti)
    {    
        ret_ti = predefined::type_unit();
    };
};

template<class T>
struct init_template_arg
{
    using obj_type  = typename std::decay<T>::type;

    static void init(Type& ret_ti)
    {    
        ret_ti = obj_type::get_static_type();
    };
};
template<>
struct init_template_arg<void>
{
    static void init(Type& ret_ti)
    {    
        ret_ti = predefined::type_unit();
    };
};

template<class ... Templates>
struct init_template_types;

template<class T, class ... Templates>
struct init_template_types<T, Templates...>
{
    static void eval(std::vector<Type>& types)
    {
        Type t;
        init_template_arg<T>::init(t);
        types.push_back(t);
        init_template_types<Templates...>::eval(types);
    }
};

template<>
struct init_template_types<>
{
    static void eval(std::vector<Type>&){};
};

template<class data_constructor, typename Fun, class base>
class matcl_dynamic_function : public dynamic_function<Fun,base>
{
    private:
        using base_type         = typename dynamic_function<Fun,base>;
        using function_traits   = typename matcl::dynamic::details::function_traits<Fun>;

    public:
        matcl_dynamic_function(Fun f) : base_type(f)
        {
            static const int n_inputs   = function_traits::n_inputs;
            static const bool is_obj    = is_object<return_type>::value;
            static const bool is_void   = std::is_same<void, return_type>::value;

            static_assert(is_obj == true || is_void == true, "return must have object or void type");

            base_type::m_args_size      = n_inputs;
            base_type::m_arg_ti         = new Type[n_inputs];            

            using input_type            = typename function_traits::input_type;
            using return_type           = typename base_type::return_type;
            init_args<n_inputs, input_type>::init(base_type::m_arg_ti);

            init_ret<return_type>::init(base_type::m_ret_ti);
        };

        virtual dynamic::function make_converter(int n_deduced, const Type deduced[], Type ded_ret,
                                    const std::vector<dynamic::function>& conv_vec) const override
        {
            static const int N = function_traits::n_inputs;
            return new details::fun_evaler_conv(this, n_deduced, deduced, ded_ret, conv_vec);
        };

        virtual bool make_eval(const dynamic::object** args, object& ret) const override
        {
            using ret_t = typename function_traits::return_type;

            data_constructor ds;
            ds.m_args   = args;
            ds.m_ret    = &ret;
            this->eval(&ds);
            return true;
        };
};

template<class Fun>
function function_register<Fun>::process_function(Fun fun)
{
    using fun_traits        = matcl::dynamic::details::function_traits<Fun>;
    using data_cons_type    = data_constructor<fun_traits>;

    return function(new matcl_dynamic_function<data_cons_type, Fun, 
                            matcl::dynamic::details::evaler>(fun));
};

template<class Fun, class ... Templates>
function function_template_register<Fun, Templates...>::process_function(Fun fun, type_vec& types)
{
    using fun_traits        = matcl::dynamic::details::function_traits<Fun>;
    using data_cons_type    = data_constructor<fun_traits>;

    types.clear();
    init_template_types<Templates...>::eval(types);
    
    return function(new matcl_dynamic_function<data_cons_type, Fun, 
                            matcl::dynamic::details::evaler>(fun));
};

template<class Fun, class ... Templates>
function function_template_register<Fun, Templates...>
        ::execute_registration(type_vec& types, ret_func& ret)
{ 
    ret = m_ret; 
    return process_function(m_fun, types); 
};

template<class T>
struct make_object_type
{
    using T0    = typename std::decay<T>::type;
    using type  = const object_type<T0>&;
};
template<class T>
struct make_object_type<T&>
{
    using T0    = typename std::decay<T>::type;
    using type  = object_type<T0>&;
};
template<class T>
struct make_object_type<const T&>
{
    using T0    = typename std::decay<T>::type;
    using type  = const object_type<T0>&;
};
template<class T>
struct make_object_type<T&&>
{
    using T0    = typename std::decay<T>::type;
    using type  = const object_type<T0>&;
};
template<class T>
struct make_object_type<const T&&>
{
    using T0    = typename std::decay<T>::type;
    using type  = const object_type<T0>&;
};

template<class List>
struct make_object_type_list{};

template<class ... Args>
struct make_object_type_list<make_list<Args...>>
{
    using type = make_list<typename make_object_type<Args>::type ...>;
};

template<class T>
struct get_from_object 
{ 
    using T0    = typename std::decay<T>::type;
    using Ret   = typename T0::value_type;
    static const Ret& eval(const T& val) { return val.get(); }; 
};

template<class T>
struct get_from_object<T&>
{
    using T0    = typename std::decay<T>::type;
    using Ret   = typename T0::value_type;
    static Ret& eval(T& val) { return val.get_unique(); }; 
};

template<class T>
struct get_from_object<T&&>
{
    using T0    = typename std::decay<T>::type;
    using Ret   = typename T0::value_type;
    static const Ret& eval(const T& val) { return val.get(); }; 
};

template<class T>
struct get_from_object<const T&&>
{
    using T0    = typename std::decay<T>::type;
    using Ret   = typename T0::value_type;
    static const Ret& eval(const T& val) { return val.get(); }; 
};

template<class T>
struct get_from_object<const T&>
{
    using T0    = typename std::decay<T>::type;
    using Ret   = typename T0::value_type;
    static const Ret& eval(const T& val) { return val.get(); }; 
};

template<class T, bool Is_enum = std::is_enum<T>::value>
struct make_object_return_type
{
    using T0    = typename std::decay<T>::type;
    using type  = object_type<T0>;
};
template<class T>
struct make_object_return_type<T, true>
{
    using type  = enum_tag;
};
template<>
struct make_object_return_type<void, false>
{
    using type  = void;
};

template<class Func>
struct make_function_args
{
    using func_traits   = function_traits<Func>;
    using return_type   = typename func_traits::return_type;
    using input_type    = typename func_traits::input_type;

    using return_object = typename make_object_return_type<return_type>::type;
    using input_object  = typename make_object_type_list<input_type>::type;
};

template<class Ret, class In, class Func, Func f>
struct register_function_eval_arg{};

template<class Ret, class ... In, class Func, Func f>
struct register_function_eval_arg<Ret, make_list<In...>,Func,f>
{
	static Ret eval(In ... val)
	{
		return Ret(f(get_from_object<In>::eval(val) ...));
	};
};

template<class ... In, class Func, Func f>
struct register_function_eval_arg<void, make_list<In...>,Func,f>
{
    using Unit = object_type<unit_type>;

	static Unit eval(In ... val)
	{
		f(get_from_object<In>::eval(val) ...);
        return Unit();
	};
};

template<class ... In, class Func, Func f>
struct register_function_eval_arg<enum_tag, make_list<In...>,Func,f>
{
    using Ret = object_type<Integer>;

	static Ret eval(In ... val)
	{
		return Ret((Integer)f(get_from_object<In>::eval(val) ...));
	};
};

template<   class Func, Func f, 
            class Ret = typename make_function_args<Func>::return_object,
            class In = typename make_function_args<Func>::input_object
        >
struct register_function_eval : register_function_eval_arg<Ret,In,Func,f>
{};

};};};

