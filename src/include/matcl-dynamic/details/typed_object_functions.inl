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

#pragma once

#include "matcl-dynamic/typed_object_functions.h"

#pragma warning(push)
#pragma warning(disable:4800) //forcing value to bool 'true' or 'false' (performance warning)

namespace matcl { namespace unqualified_call
{

//evaluators;
//evaluations must be performed in the matcl namespace

template<class Ret, class T>
struct eval_real
{
    static Ret eval(const T& a)
    {
        return Ret(real(a));
    }
};

template<class Ret, class T>
struct eval_imag
{
    static Ret eval(const T& a)
    {
        return Ret(imag(a));
    }
};

template<class Ret, class T, class S>
struct eval_eeq
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a == b);
    }
};

template<class Ret, class T, class S>
struct eval_neq
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a != b);
    }
};

template<class Ret, class T, class S>
struct eval_geq
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a >= b);
    }
};

template<class Ret, class T, class S>
struct eval_leq
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a <= b);
    }
};

template<class Ret, class T, class S>
struct eval_gt
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a > b);
    }
};

template<class Ret, class T, class S>
struct eval_lt
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a < b);
    }
};

template<class T>
struct eval_cast_bool
{
    static bool eval(const T& a)
    {
        return (bool)a;
    }
};

template<class T>
struct eval_not
{
    static bool eval(const T& a)
    {
        return !a;
    }
};

template<class Ret, class T, class S>
struct eval_plus
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a + b);
    }
};

template<class Ret, class T, class S>
struct eval_minus
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a - b);
    }
};

template<class Ret, class T, class S>
struct eval_op_mul
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a * b);
    }
};

template<class Ret, class T, class S>
struct eval_div
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(a / b);
    }
};
template<class Ret>
struct eval_div<Ret, Integer, Integer>
{
    static Ret eval(Integer a, Integer b)
    {
        return Ret(Real(a) / Real(b));
    }
};
template<class Ret>
struct eval_div<Ret, Integer, Float>
{
    static Ret eval(const Integer& a, const Float& b)
    {
        return Ret(Real(a) / Real(b));
    }
};
template<class Ret>
struct eval_div<Ret, Float, Integer>
{
    static Ret eval(const Float& a, const Integer& b)
    {
        return Ret(Real(a) / Real(b));
    }
};


template<class Ret, class T>
struct eval_uminus
{
    static Ret eval(const T& a)
    {
        return Ret( -a );
    }
};

template<class Ret, class T, class S>
struct eval_idiv
{
    static Ret eval(const T& a, const S& b)
    {
        return Ret(idiv(a,b));
    }
};

template<class Ret, class To, class From, bool Use_cast>
struct eval_cast
{
    static Ret eval(const From& a)
    {
        return Ret(matcl::convert_scalar<To, From>(a));
    }
};

template<class Ret, class To, class From>
struct eval_cast<Ret, To, From, false>
{
    static Ret eval(const From& a)
    {
        return Ret(To(a));
    }
};

}};

namespace matcl { namespace dynamic
{

template<class To, class From, class Enable>
object_type<To> dynamic::cast(const object_type<From>& a)
{
    using info  = result_of::result_of_cast<To,From>;
    using ret   = typename info::type_object;
    static const bool use_cast = info::use_cast;

    return matcl::unqualified_call::eval_cast<ret, To,From,use_cast>::eval(a.get());
};

template<class To, class From>
object_type<To> dynamic::convert(const object_type<From>& a, 
                                  typename details::enable_convert<From, To, void*>::type)
{
    return object_type<To>(To(a.get()));
};

template<class T>
typename result_of::result_of_real<T>::type_object
dynamic::real(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_real<T>::type_object;
    return matcl::unqualified_call::eval_real<return_type, T>::eval(a.get());
};

template<class T>
typename result_of::result_of_imag<T>::type_object
dynamic::imag(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_imag<T>::type_object;
    return matcl::unqualified_call::eval_imag<return_type, T>::eval(a.get());
};

template<class T>
bool dynamic::is_zero(const object_type<T>& a)
{
    return a.is_zero();
};

template<class T>
bool dynamic::is_one(const object_type<T>& a)
{
    return a.is_one();
};

template<class T, class Enable>
bool dynamic::cast_bool(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_cast_bool<T>::eval(a.get());
}

template<class T, class Enable>
bool dynamic::operator!(const object_type<T>& a)
{
    return matcl::unqualified_call::eval_not<T>::eval(a.get());
};

template<class T>
typename result_of::result_of_uminus<T>::type_object
dynamic::operator-(const object_type<T>& a)
{
    using return_type = typename result_of::result_of_uminus<T>::type_object;
    return matcl::unqualified_call::eval_uminus<return_type, T>::eval(a.get());
};

template<class T, class S, class Ret>
Ret dynamic::operator==(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator==(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator==(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_eeq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator!=(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator!=(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator!=(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_neq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator>=(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator>=(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator>=(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_geq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator<=(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator<=(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator<=(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_leq<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator>(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator>(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator>(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_gt<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator<(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator<(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator<(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_lt<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator+(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator+(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator+(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_plus<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator-(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator-(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator-(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_minus<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator*(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator*(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator*(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_op_mul<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::operator/(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::operator/(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::operator/(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_div<Ret,T,S>::eval(a,b.get());
}

template<class T, class S, class Ret>
Ret dynamic::idiv(const object_type<T>& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_idiv<Ret,T,S>::eval(a.get(),b.get());
}
template<class T, class S, class Ret>
Ret dynamic::idiv(const object_type<T>& a, const S& b)
{
    return matcl::unqualified_call::eval_idiv<Ret,T,S>::eval(a.get(),b);
}
template<class T, class S, class Ret>
Ret dynamic::idiv(const T& a, const object_type<S>& b)
{
    return matcl::unqualified_call::eval_idiv<Ret,T,S>::eval(a,b.get());
}

};};

#pragma warning(pop)