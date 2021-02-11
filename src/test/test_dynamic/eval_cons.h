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

#include "rand_scalars.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-scalar/objects/typed_object_functions.h"
#include "matcl-scalar/lib_functions/manip.h"
#include "matcl-dynamic/result_of.h"

#pragma warning(push)
#pragma warning(disable:4800) // forcing value to bool 'true' or 'false' (performance warning)
#pragma warning(disable:4127) // conditional expression is constant

namespace matcl { namespace result_of
{
template<>
struct result_of_uminus<mp_int>
{
    using T             = mp_int;
    using type          = decltype(-(std::declval<T>()));
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;
};
template<>
struct result_of_uminus<mp_float>
{
    using T             = mp_float;
    using type          = decltype(-(std::declval<T>()));
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;
};
template<>
struct result_of_uminus<mp_complex>
{
    using T             = mp_complex;
    using type          = decltype(-(std::declval<T>()));
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;
};
template<>
struct result_of_uminus<mp_rational>
{
    using T             = mp_rational;
    using type          = decltype(-(std::declval<T>()));
    using type_object   = dynamic::object_type<typename std::decay<type>::type>;
};
}}

namespace matcl
{
    inline fp_type fpclassify(const Float_complex&)
    {
        return fp_type();
    };
    inline fp_type fpclassify(const Complex&)
    {
        return fp_type();
    };

    inline fp_type fpclassify(const mp_complex&)
    {
        return fp_type();
    };

    inline fp_type fpclassify(const OFloat_complex&)
    {
        return fp_type();
    };
    inline fp_type fpclassify(const OComplex&)
    {
        return fp_type();
    };

    inline fp_type fpclassify(const MP_complex&)
    {
        return fp_type();
    };

    inline bool is_finite(bool)
    {
        return true;
    };

    inline bool is_nan(bool)
    {
        return false;
    };

    inline Real eps(bool)
    {
        return 0.0;
    };

    inline Real abs(bool v)
    {
        return Integer(v);
    };
};

namespace matcl { namespace test
{

namespace mdd = matcl::dynamic::details;

template<class T1, class T2, bool Can_cons>
struct construct
{
    static T1 eval(const T2& val)
    {
        return T1(val);
    };
};
template<class T1, class T2>
struct construct<T1, T2, false>
{
    static T1 eval(const T2& val)
    {
        (void)val;
        return T1();
    };
};

template<class T1, class T2, bool Can_cons>
struct construct_convert_func
{
    using Ret = object_type<T1>;
    static Ret eval(const T2& val)
    {
        return dynamic::convert<T1>(val);
    };
};
template<class T1, class T2>
struct construct_convert_func<T1, T2, false>
{
    using Ret = object_type<T1>;

    static Ret eval(const T2& val)
    {
        (void)val;
        return Ret();
    };
};

template<class T1, class T2, bool Can_cons>
struct construct_convert_member
{
    using Ret = object_type<T1>;

    static Ret eval(const T2& val)
    {
        return val.convert<T1>();
    };
};
template<class T1, class T2>
struct construct_convert_member<T1, T2, false>
{
    using Ret = object_type<T1>;

    static Ret eval(const T2&)
    {
        return Ret();
    };
};

struct eval_cons_val : eval_scalars<eval_cons_val>
{
    Integer code;

    eval_cons_val(Integer c) : code(c){};

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        using type_1    = mdy::object_type<T1>;
        using type_2    = mdy::object_type<T2>;

        static const bool can_construct  =  matcl::dynamic::details::is_convertible<type_2, type_1, true>::value
                                || matcl::dynamic::details::is_convertible<type_2, type_1, false>::value;
        bool check_construct    = matcl::details::is_mp_scalar<type_1>::value;

        double res = 0;

        type_1 v1(s1);
        type_2 v2(s2);

        type_1 val_11;
        type_1 val_12;
        type_1 val_13;
        type_1 val_14;

        try
        {
            val_11 = type_1(mdy::object(v2), dynamic::from_object());            

            if (can_construct == false && check_construct == true)
            {
                res += 1;
                return res;
            };

            val_12  = construct<type_1, type_2, can_construct>::eval(v2);
            val_13  = construct_convert_func<T1, type_2, can_construct>::eval(v2);
            val_14  = construct_convert_member<T1, type_2, can_construct>::eval(v2);
        }
        catch(error::matcl_dynamic_exception& ex)
        {
            if (can_construct)
            {
                if (type_1::get_static_type() == MP_int::get_static_type()
                    && (type_2::get_static_type() == OFloat::get_static_type()
                        ||type_2::get_static_type() == OReal::get_static_type()))
                {
                    //conversion Float->mp_int is allowed in C++ and cannot be disabled
                    //this conversion is not allowed in dynamic
                }
                else
                {
                    disp(ex.what());
                    res += 1;                    
                };                
            }     
            return res;
        }
        catch(...)
        {
            res += 1;
            return res;
        }

        mdy::object rhs(v2);
        mdy::object lhs(v1);

        mdy::object val_21(type_1::get_static_type(), rhs);
        mdy::object val_22 = convert(type_1::get_static_type(), rhs);

        if (val_21.get_type() != type_1::get_static_type())
            res += 1;
        if (val_22.get_type() != type_1::get_static_type())
            res += 1;

        type_1 val_31(val_21, dynamic::from_object());
        type_1 val_32(val_22, dynamic::from_object());

        if (val_11.get() != val_12.get() && matcl::is_nan(val_11.get()) == false
            && check_construct == true)
        {
            res += 1;
        }
        if (val_11.get() != val_13.get() && matcl::is_nan(val_11.get()) == false
            && check_construct == true)
        {
            res += 1;
        }
        if (val_11.get() != val_14.get() && matcl::is_nan(val_11.get()) == false
            && check_construct == true)
        {
            res += 1;
        }

        if (val_31.get() != val_11.get() && matcl::is_nan(val_11.get()) == false)
            res += 1;
        if (val_32.get() != val_11.get() && matcl::is_nan(val_11.get()) == false)
            res += 1;

        if (res != 0)
            disp(code);

        return res;
    };
};

template<class T1, class T2, bool Has_assign = mdd::is_assignable<T2,T1>::value>
struct eval_assign_val_impl
{
    static double eval(Integer code, const T1& s1, const T2& s2)
    {
        //if (code == 22)
        //    disp("break");

        using type_1    = mdy::object_type<T1>;
        using type_2    = mdy::object_type<T2>;

        type_1 v1(s1);
        type_2 v2(s2);

        mdy::Type t1 = type_1::get_static_type();
        mdy::Type t2 = type_2::get_static_type();

        mdy::object rhs(v2);
        mdy::object lhs(v1);

        v1 = v2;

        mdy::object val_2(type_1::get_static_type());

        double res = 0;

        try
        {
            val_2 = rhs;
        }
        catch(std::exception& ex)
        {
            disp(ex.what());
            res += 1;
        };        

        if (val_2.get_type() != type_1::get_static_type())
            res += 1;

        type_1 val_3(val_2, dynamic::from_object());

        if (val_3.get() != v1.get() && matcl::is_nan(v1.get()) == false)
            res += 1;

        if (res != 0)
        {
            //out_stream << s1 << " " << s2 << "\n";
            //out_stream << val_3 << " " << v1 << "\n";
            disp(code);
        }

        return res;
    };
};

template<class T1, class T2>
struct eval_assign_val_impl<T1, T2, false>
{
    static double eval(Integer code, const T1& s1, const T2& s2)
    {
        (void)code;

        using type_1    = mdy::object_type<T1>;
        using type_2    = mdy::object_type<T2>;

        type_1 v1(s1);
        type_2 v2(s2);

        mdy::object rhs(v2);
        mdy::object lhs(v1);

        mdy::object val_2(type_1::get_static_type());

        double res = 0;

        try
        {
            val_2 = rhs;
            res += 1;
        }
        catch(...){};
        
        return res;
    };
};

template<class T1, class T2>
struct eval_cast_val_impl
{
    static double eval(Integer code, const T1& s1, const T2& s2)
    {
        if (code == -1)
            disp("break");

        using type_1    = mdy::object_type<T1>;
        using type_2    = mdy::object_type<T2>;

        type_1 v1(s1);
        type_2 v2(s2);

        mdy::Type t1 = type_1::get_static_type();
        mdy::Type t2 = type_2::get_static_type();
        
        mdy::object lhs;
        mdy::object rhs(v2);

        T1 ret_val  = convert_scalar<T1>(s2);

        type_2 v3   = v2;
        type_1 res2 = v3.cast<T1>();
        type_1 res3 = matcl::dynamic::cast<T1>(v3);

        double res = 0;

        try
        {
            lhs.reset(cast(t1, rhs));
        }
        catch(std::exception& ex)
        {
            disp(ex.what());
            res += 1;
        };        

        if (lhs.get_type() != type_1::get_static_type())
            res += 1;

        type_1 val_3(lhs, dynamic::from_object());

        if (val_3.get() != ret_val && matcl::is_nan(ret_val) == false)
            res += 1;
        if (res2.get() != ret_val && matcl::is_nan(ret_val) == false)
            res += 1;
        if (res3.get() != ret_val && matcl::is_nan(ret_val) == false)
            res += 1;

        if (res != 0)
        {
            out_stream << s1 << " " << s2 << "\n";
            out_stream << val_3 << " " << v1 << "\n";
            out_stream << code << "\n";
        }

        return res;
    };
};

struct eval_assign_val : eval_scalars<eval_assign_val>
{
    Integer             code;

    eval_assign_val(Integer c)  : code(c) {};

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        return eval_assign_val_impl<T1, T2>::eval(code, s1, s2);
    };
};

struct eval_cast_val : eval_scalars<eval_cast_val>
{
    Integer             code;

    eval_cast_val(Integer c)  : code(c) {};

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        return eval_cast_val_impl<T1, T2>::eval(code, s1, s2);
    };
};

struct eval_object_func_val : eval_scalars_1<eval_object_func_val>
{
    Integer code;

    eval_object_func_val(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        //if (code == 4)
        //    out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1 v1(s1);

        bool res1_1         = (bool)s1;
        bool res2_1         = !s1;

        bool res1_21        = (bool)v1;
        bool res1_22        = cast_bool(v1);
        bool res2_2         = !v1;        

        mdy::object lhs(v1);

        bool res1_31        = (bool)lhs;
        bool res1_32        = cast_bool(lhs);
        bool res2_3         = !lhs;

        type_1 res4_2       = v1.clone();
        mdy::object res4_3   = lhs.clone();

        if (res1_1 != res1_21)
            res += 1;
        if (res1_1 != res1_22)
            res += 1;
        if (res1_1 != res1_31)
            res += 1;
        if (res1_1 != res1_32)
            res += 1;
        if (res2_1 != res2_2)
            res += 1;
        if (res2_1 != res2_3)
            res += 1;

        if (res4_2.get() != type_1(res4_3, dynamic::from_object()).get() && matcl::is_nan(res4_2.get()) == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_uminus_val : eval_scalars_1<eval_uminus_val>
{
    Integer code;

    eval_uminus_val(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        //if (code == 12)
        //    out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        using T0    = typename result_of::result_of_uminus<T1>::type_object;
        using Type = decltype(-v1);

        T1 res1_1           = -s1;
        type_1 res1_2       = -v1;
        mdy::object res1_3   = -lhs;

        if (res1_1 != res1_2.get() && matcl::is_nan(s1) == false)
            res += 1;
        if (res1_3.get_type() != type_1::get_static_type())
            res += 1;
        if (res1_1 != type_1(res1_3, dynamic::from_object()).get() && matcl::is_nan(s1) == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_reim_val : eval_scalars_1<eval_reim_val>
{
    Integer code;

    eval_reim_val(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        using T1_re     = typename matcl::details::real_type<T1>::type;
        using T1_fre    = decltype(arg(s1));
        using type_1    = mdy::object_type<T1>;
        using type_1_re = mdy::object_type<T1_re>;
        using type_1_fre= mdy::object_type<T1_fre>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        T1_re res1_1        = real(s1);
        T1_re res2_1        = imag(s1);
        T1    res3_1        = conj(s1);
        T1_fre res4_1       = arg(s1);
        T1_fre res5_1       = angle(s1);
        T1_re res6_1        = abs(s1);
        T1_re res7_1        = abs2(s1);

        type_1_re res1_2    = mdy::real(v1);
        type_1_re res2_2    = mdy::imag(v1);
        type_1    res3_2    = mdy::conj(v1);
        type_1_fre res4_2   = mdy::arg(v1);
        type_1_fre res5_2   = mdy::angle(v1);
        type_1_re res6_2    = mdy::abs(v1);
        type_1_re res7_2    = mdy::abs2(v1);

        mdy::object res1_3   = real(lhs);
        mdy::object res2_3   = imag(lhs);
        mdy::object res3_3   = conj(lhs);
        mdy::object res4_3   = arg(lhs);
        mdy::object res5_3   = angle(lhs);
        mdy::object res6_3   = abs(lhs);
        mdy::object res7_3   = abs2(lhs);

        bool is_nan_1       = matcl::is_nan(s1);

        if (res1_1 != res1_2.get() && is_nan_1 == false)
            res += 1;
        if (res2_1 != res2_2.get() && is_nan_1 == false)
            res += 1;
        if (res3_1 != res3_2.get() && is_nan_1 == false)
            res += 1;
        if (res4_1 != res4_2.get() && is_nan_1 == false)
            res += 1;
        if (res5_1 != res5_2.get() && is_nan_1 == false)
            res += 1;
        if (res6_1 != res6_2.get() && is_nan_1 == false)
            res += 1;
        if (res7_1 != res7_2.get() && is_nan_1 == false)
            res += 1;

        if (res1_3.get_type() != type_1_re::get_static_type())
            res += 1;
        if (res2_3.get_type() != type_1_re::get_static_type())
            res += 1;
        if (res3_3.get_type() != type_1::get_static_type())
            res += 1;
        if (res4_3.get_type() != type_1_fre::get_static_type())
            res += 1;
        if (res5_3.get_type() != type_1_fre::get_static_type())
            res += 1;
        if (res6_3.get_type() != type_1_re::get_static_type())
            res += 1;
        if (res7_3.get_type() != type_1_re::get_static_type())
            res += 1;

        if (res1_1 != type_1_re(res1_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res2_1 != type_1_re(res2_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res3_1 != type_1(res3_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res4_1 != type_1_fre(res4_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res5_1 != type_1_fre(res5_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res6_1 != type_1_re(res6_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res7_1 != type_1_re(res7_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_is : eval_scalars_1<eval_is>
{
    Integer code;

    eval_is(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        //if (code == 10)
        //    out_stream << "break" << "\n";

        using T1_re     = typename matcl::details::real_type<T1>::type;
        using type_1    = mdy::object_type<T1>;
        using type_1_re = mdy::object_type<T1_re>;

        bool is_compl   = md::is_complex<T1>::value;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        bool res1_1     = is_nan(s1);
        bool res2_1     = is_inf(s1);
        bool res3_1     = is_finite(s1);
        bool res4_1     = is_regular(s1);
        bool res5_1     = is_int(s1);
        bool res6_1     = is_real(s1);
        bool res7_1     = is_zero(s1);
        bool res8_1     = is_one(s1);
        bool res9_1     = is_normal(s1);        

        bool res1_2     = (bool)is_nan(v1);
        bool res2_2     = (bool)is_inf(v1);
        bool res3_2     = (bool)is_finite(v1);
        bool res4_2     = (bool)is_regular(v1);
        bool res5_2     = (bool)is_int(v1);
        bool res6_2     = (bool)is_real(v1);
        bool res7_2     = is_zero(v1);
        bool res8_2     = is_one(v1);
        bool res9_2     = (bool)is_normal(v1);        

        bool res1_3     = (bool)is_nan(lhs);
        bool res2_3     = (bool)is_inf(lhs);
        bool res3_3     = (bool)is_finite(lhs);
        bool res4_3     = (bool)is_regular(lhs);
        bool res5_3     = (bool)is_int(lhs);
        bool res6_3     = (bool)is_real(lhs);
        bool res7_3     = is_zero(lhs);
        bool res8_3     = is_one(lhs);
        bool res9_3     = (bool)is_normal(lhs);        

        if (res1_1 != res1_2 || res1_1 != res1_3)
            res += 1;
        if (res2_1 != res2_2 || res2_1 != res2_3)
            res += 1;
        if (res3_1 != res3_2 || res3_1 != res3_3)
            res += 1;
        if (res4_1 != res4_2 || res4_1 != res4_3)
            res += 1;
        if (res5_1 != res5_2 || res5_1 != res5_3)
            res += 1;
        if (res6_1 != res6_2 || res6_1 != res6_3)
            res += 1;
        if (res7_1 != res7_2 || res7_1 != res7_3)
            res += 1;
        if (res8_1 != res8_2 || res8_1 != res8_3)
            res += 1;
        if (res9_1 != res9_2 || res9_1 != res9_3)
            res += 1;

        if (is_compl == false)
        {
            fp_type res10_1 = fpclassify(s1);
            fp_type res10_2 = fpclassify(v1);
            fp_type res10_3 = fpclassify(lhs);

            if (res10_1 != res10_2 || res10_1 != res10_3)
                res += 1;
        };

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_next : eval_scalars_1<eval_next>
{
    Integer code;

    eval_next(Integer c) : code(c){};

    double eval_scal_func(const mp_complex&)    { return 0.0 ;};
    double eval_scal_func(const Complex&)       { return 0.0 ;};
    double eval_scal_func(const Float_complex&) { return 0.0 ;};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        //using T1_re     = typename matcl::details::real_type<T1>::type;
        //using T1_fre    = decltype(arg(s1));        
        //using type_1_re = mdy::object_type<T1_re>;
        //using type_1_fre= mdy::object_type<T1_fre>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        auto res1_1         = nextabove(s1);
        auto res2_1         = nextbelow(s1);

        auto res1_2         = nextabove(v1);
        auto res2_2         = nextbelow(v1);

        mdy::object res1_3   = nextabove(lhs);
        mdy::object res2_3   = nextbelow(lhs);

        bool is_nan_1       = matcl::is_nan(s1);

        if (res1_1 != res1_2.get() && is_nan_1 == false)
            res += 1;
        if (res2_1 != res2_2.get() && is_nan_1 == false)
            res += 1;

        using type_res_1    = decltype(res1_2);
        using type_res_2    = decltype(res2_2);

        if (res1_3.get_type() != type_res_1::get_static_type())
            res += 1;
        if (res2_3.get_type() != type_res_2::get_static_type())
            res += 1;

        if (res1_1 != type_res_1(res1_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;
        if (res2_1 != type_res_2(res2_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct test_eval_eps : eval_scalars_1<test_eval_eps>
{
    Integer code;

    test_eval_eps(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        auto res1_1         = eps(s1);

        auto res1_2         = eps(v1);

        mdy::object res1_3   = eps(lhs);

        bool is_nan_1       = matcl::is_nan(res1_1);

        if (res1_1 != res1_2.get() && is_nan_1 == false)
            res += 1;

        using type_res_1    = decltype(res1_2);

        if (res1_3.get_type() != type_res_1::get_static_type())
            res += 1;

        if (res1_1 != type_res_1(res1_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_signbit : eval_scalars_1<eval_signbit>
{
    Integer code;

    eval_signbit(Integer c) : code(c){};

    double eval_scal_func(const mp_complex&)    { return 0.0 ;};
    double eval_scal_func(const Complex&)       { return 0.0 ;};
    double eval_scal_func(const Float_complex&) { return 0.0 ;};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        bool res1_1         = signbit(s1);
        bool res1_2         = (bool)signbit(v1);
        bool res1_3         = (bool)signbit(lhs);

        bool is_nan_1       = matcl::is_nan(s1);
        (void)is_nan_1;

        if (res1_1 != res1_2)
            res += 1;

        if (res1_1 != res1_3)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

struct eval_isign : eval_scalars_1<eval_isign>
{
    Integer code;

    eval_isign(Integer c) : code(c){};

    double eval_scal_func(const mp_complex&)    { return 0.0 ;};
    double eval_scal_func(const Complex&)       { return 0.0 ;};
    double eval_scal_func(const Float_complex&) { return 0.0 ;};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        auto res1_1         = isign(s1);
        auto res1_2         = isign(v1);
        auto res1_3         = isign(lhs);

        bool is_nan_1       = matcl::is_nan(s1);
        (void)is_nan_1;

        if (res1_1 != res1_2)
            res += 1;

        if (res1_1 != res1_3)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

template<class Func, bool With_mp>
struct eval_scalar_func_templ_obj : eval_scalars_1<eval_scalar_func_templ_obj<Func,With_mp>>
{
    Integer code;

    eval_scalar_func_templ_obj(Integer c) : code(c){};

    template<class T1>
    static auto eval(const T1& a1) -> decltype(Func::eval(a1))
    {
        return Func::eval(a1);
    };

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        if (code == -1)
            out_stream << "break" << "\n";

        if (md::is_complex<T1>::value && Func::is_complex_allowed == false)
            return 0.0;

        using type_1    = mdy::object_type<T1>;

        double res = 0;

        type_1      v1(s1);
        mdy::object  lhs(v1);

        auto res1_1         = eval(s1);

        auto res1_2         = eval(v1);

        mdy::object res1_3   = eval(lhs);

        bool is_nan_1       = matcl::is_nan(res1_1);

        if (res1_1 != res1_2.get() && is_nan_1 == false)
            res += 1;

        using type_res_1    = decltype(res1_2);

        if (res1_3.get_type() != type_res_1::get_static_type())
            res += 1;

        if (res1_1 != type_res_1(res1_3, dynamic::from_object()).get() && is_nan_1 == false)
            res += 1;

        if (res != 0)
            out_stream << code << "\n";

        return res;
    };
};

template<class Func, class T1, bool Is_mp>
struct eval_scal_func_impl
{
    static double eval(const T1& s, Integer code)
    {
        return eval_scalar_func_templ_obj<Func,true>(code).eval_scal_func(s);
    };
};

template<class Func, class T1>
struct eval_scal_func_impl<Func, T1, true>
{
    static double eval(const T1&, Integer)
    {
        return 0.0;
    };
};

template<class Func>
struct eval_scalar_func_templ_obj<Func,false> 
    : eval_scalars_1<eval_scalar_func_templ_obj<Func,false>>
{
    Integer code;

    eval_scalar_func_templ_obj(Integer c) : code(c){};

    template<class T1>
    double eval_scal_func(const T1& s1)
    {
        static const bool is_mp = matcl::details::is_mp_scalar<T1>::value;
        return eval_scal_func_impl<Func, T1, is_mp>::eval(s1, code);
    };
};

}};

#pragma warning(pop)