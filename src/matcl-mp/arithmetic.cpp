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

#include "matcl-mp/details/arithmetic.h"
#include "matcl-mp/func_unary.h"
#include "matcl-mp/func_binary.h"
#include "matcl-mp/constants.h"
#include "utils/impl_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "utils/utils.h"
#include "utils/extend_precision.h"
#include "matcl-mp/details/prec_utils.h"
#include "matcl-mp/error_flags.h"
#include "impl_functions.h"

namespace matcl { namespace mp { namespace details
{

namespace bm = boost::multiprecision;
namespace mmd = matcl::mp::details;

//---------------------------------------------------------------
//                  plus_impl
//---------------------------------------------------------------
mp_int details::plus_impl(const mp_int& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_int details::plus_impl(const mp_int& a, Integer b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_int details::plus_impl(Integer a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_float details::plus_impl(const mp_int& a, Real b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_float details::plus_impl(Real a, const mp_int& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_int& a, const Complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_complex details::plus_impl(const Complex& a, const mp_int& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_int& a, const Float_complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_complex details::plus_impl(const Float_complex& a, const mp_int& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}

mp_rational details::plus_impl(const mp_rational& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_rational details::plus_impl(const mp_rational& a, Integer b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_rational details::plus_impl(Integer a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_float details::plus_impl(const mp_rational& a, Real b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_float details::plus_impl(Real a, const mp_rational& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_rational& a, const Complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_complex details::plus_impl(const Complex& a, const mp_rational& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_rational& a, const Float_complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_complex details::plus_impl(const Float_complex& a, const mp_rational& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_rational details::plus_impl(const mp_rational& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}
mp_rational details::plus_impl(const mp_int& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)+mmd::impl_value(b)));
}

mp_float details::plus_impl(const mp_float& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());
	mp_float ret(0, p);
	mpfr_add(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::plus_impl(const mp_float& a, Integer b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_si(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::plus_impl(Integer a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_si(mmd::impl_value(ret), mmd::impl_value(b), a, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::plus_impl(const mp_float& a, Real b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_d(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::plus_impl(Real a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_d(mmd::impl_value(ret), mmd::impl_value(b), a, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_complex details::plus_impl(const mp_float& a, const Complex& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_d(mmd::impl_value(ret), mmd::impl_value(a), real(b), MPFR_RNDN);
    ret.update_debug();
	return mp_complex(ret, mp_float(imag(b), p));
};
mp_complex details::plus_impl(const Complex& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_d(mmd::impl_value(ret), mmd::impl_value(b), real(a), MPFR_RNDN);
    ret.update_debug();
	return mp_complex(ret, mp_float(imag(a), p));
};
mp_complex details::plus_impl(const mp_float& a, const Float_complex& b, precision p)
{
    return plus_impl(a,Complex(b), p);
}
mp_complex details::plus_impl(const Float_complex& a, const mp_float& b, precision p)
{
    return plus_impl(Complex(a), b, p);
}
mp_float details::plus_impl(const mp_float& a, const mp_int& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_z(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::plus_impl(const mp_int& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_z(mmd::impl_value(ret), mmd::impl_value(b), 
               mmd::impl_value(a).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::plus_impl(const mp_float& a, const mp_rational& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_q(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::plus_impl(const mp_rational& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_add_q(mmd::impl_value(ret), mmd::impl_value(b), 
               mmd::impl_value(a).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};

mp_complex details::plus_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(real(a),real(b),p),plus_impl(imag(a),imag(b),p));
}
mp_complex details::plus_impl(const mp_complex& a, Integer b, precision p)
{
    return mp_complex(plus_impl(real(a),Real(b),p),imag(a), p);
}
mp_complex details::plus_impl(Integer a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(Real(a),real(b),p), imag(b), p);
}
mp_complex details::plus_impl(const mp_complex& a, Real b, precision p)
{
    return mp_complex(plus_impl(real(a),b,p),imag(a), p);
}
mp_complex details::plus_impl(Real a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(Real(a),real(b),p),imag(b), p);
}
mp_complex details::plus_impl(const mp_complex& a, const Complex& b, precision p)
{
    return mp_complex(plus_impl(real(a), real(b), p), plus_impl(imag(a), imag(b), p));
}
mp_complex details::plus_impl(const Complex& a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(real(a),real(b),p), plus_impl(imag(a),imag(b),p));
}
mp_complex details::plus_impl(const mp_complex& a, const Float_complex& b, precision p)
{
    return mp_complex(plus_impl(real(a),real(b),p), plus_impl(imag(a),imag(b),p));
}
mp_complex details::plus_impl(const Float_complex& a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(real(a),real(b),p), plus_impl(imag(a),imag(b),p));
}
mp_complex details::plus_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}
mp_complex details::plus_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return mp_complex(plus_impl(real(a),b,p), imag(a), p);
}
mp_complex details::plus_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return mp_complex(plus_impl(a,real(b),p), imag(b), p);
}
mp_complex details::plus_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return plus_impl(a, mp_float(b), p);
}
mp_complex details::plus_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return plus_impl(mp_float(a), b, p);
}

//---------------------------------------------------------------
//                  minus_impl
//---------------------------------------------------------------
mp_int details::minus_impl(const mp_int& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_int details::minus_impl(const mp_int& a, Integer b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_int details::minus_impl(Integer a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_float details::minus_impl(const mp_int& a, Real b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_float details::minus_impl(Real a, const mp_int& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_int& a, const Complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_complex details::minus_impl(const Complex& a, const mp_int& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_int& a, const Float_complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_complex details::minus_impl(const Float_complex& a, const mp_int& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}

mp_rational details::minus_impl(const mp_rational& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_rational details::minus_impl(const mp_rational& a, Integer b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_rational details::minus_impl(Integer a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_float details::minus_impl(const mp_rational& a, Real b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_float details::minus_impl(Real a, const mp_rational& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_rational& a, const Complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_complex details::minus_impl(const Complex& a, const mp_rational& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_rational& a, const Float_complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_complex details::minus_impl(const Float_complex& a, const mp_rational& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_rational details::minus_impl(const mp_rational& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}
mp_rational details::minus_impl(const mp_int& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)-mmd::impl_value(b)));
}

mp_float details::minus_impl(const mp_float& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());
	mp_float ret(0, p);
	mpfr_sub(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::minus_impl(const mp_float& a, Integer b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_sub_si(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::minus_impl(Integer a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_si_sub(mmd::impl_value(ret), a, mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::minus_impl(const mp_float& a, Real b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_sub_d(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::minus_impl(Real a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_d_sub(mmd::impl_value(ret), a, mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_complex details::minus_impl(const mp_float& a, const Complex& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_sub_d(mmd::impl_value(ret), mmd::impl_value(a), real(b), MPFR_RNDN);
    ret.update_debug();
	return mp_complex(ret, -imag(b));
};
mp_complex details::minus_impl(const Complex& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_d_sub(mmd::impl_value(ret), real(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return mp_complex(ret, imag(a), p);
};
mp_complex details::minus_impl(const mp_float& a, const Float_complex& b, precision p)
{
    return minus_impl(a,Complex(b), p);
}
mp_complex details::minus_impl(const Float_complex& a, const mp_float& b, precision p)
{
    return minus_impl(Complex(a), b, p);
}
mp_float details::minus_impl(const mp_float& a, const mp_int& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_sub_z(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::minus_impl(const mp_int& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_z_sub(mmd::impl_value(ret), mmd::impl_value(a).backend().data(), 
               mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::minus_impl(const mp_float& a, const mp_rational& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_sub_q(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::minus_impl(const mp_rational& a, const mp_float& b, precision p)
{
    mp_float ret = minus_impl(b, a, p);
	return -ret;
};

mp_complex details::minus_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(real(a),real(b),p),minus_impl(imag(a),imag(b),p));
}
mp_complex details::minus_impl(const mp_complex& a, Integer b, precision p)
{
    return mp_complex(minus_impl(real(a),Real(b),p),imag(a), p);
}
mp_complex details::minus_impl(Integer a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(Real(a),real(b),p), -imag(b), p);
}
mp_complex details::minus_impl(const mp_complex& a, Real b, precision p)
{
    return mp_complex(minus_impl(real(a),b,p), imag(a), p);
}
mp_complex details::minus_impl(Real a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(Real(a),real(b),p), -imag(b), p);
}
mp_complex details::minus_impl(const mp_complex& a, const Complex& b, precision p)
{
    return mp_complex(minus_impl(real(a),real(b),p), minus_impl(imag(a),imag(b),p));
}
mp_complex details::minus_impl(const Complex& a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(real(a),real(b),p), minus_impl(imag(a),imag(b),p));
}
mp_complex details::minus_impl(const mp_complex& a, const Float_complex& b, precision p)
{
    return mp_complex(minus_impl(real(a),real(b),p), minus_impl(imag(a),imag(b),p));
}
mp_complex details::minus_impl(const Float_complex& a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(real(a),real(b),p), minus_impl(imag(a),imag(b),p));
}
mp_complex details::minus_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}
mp_complex details::minus_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return mp_complex(minus_impl(real(a),b,p), imag(a), p);
}
mp_complex details::minus_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return mp_complex(minus_impl(a,real(b),p), -imag(b), p);
}
mp_complex details::minus_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return minus_impl(a, mp_float(b), p);
}
mp_complex details::minus_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return minus_impl(mp_float(a), b, p);
}

//---------------------------------------------------------------
//                  mul_impl
//---------------------------------------------------------------
mp_int details::mul_impl(const mp_int& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_int details::mul_impl(const mp_int& a, Integer b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_int details::mul_impl(Integer a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_float details::mul_impl(const mp_int& a, Real b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_float details::mul_impl(Real a, const mp_int& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_int& a, const Complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_complex details::mul_impl(const Complex& a, const mp_int& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_int& a, const Float_complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_complex details::mul_impl(const Float_complex& a, const mp_int& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}

mp_rational details::mul_impl(const mp_rational& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_rational details::mul_impl(const mp_rational& a, Integer b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_rational details::mul_impl(Integer a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_float details::mul_impl(const mp_rational& a, Real b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_float details::mul_impl(Real a, const mp_rational& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_rational& a, const Complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_complex details::mul_impl(const Complex& a, const mp_rational& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_rational& a, const Float_complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_complex details::mul_impl(const Float_complex& a, const mp_rational& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_rational details::mul_impl(const mp_rational& a, const mp_int& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}
mp_rational details::mul_impl(const mp_int& a, const mp_rational& b, precision)
{
    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)*mmd::impl_value(b)));
}

mp_float details::mul_impl(const mp_float& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());
	mp_float ret(0, p);
	mpfr_mul(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::mul_impl(const mp_float& a, Integer b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_si(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::mul_impl(Integer a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_si(mmd::impl_value(ret), mmd::impl_value(b), a, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::mul_impl(const mp_float& a, Real b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_d(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::mul_impl(Real a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_d(mmd::impl_value(ret), mmd::impl_value(b), a, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_complex details::mul_impl(const mp_float& a, const Complex& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret_re(0, p);
    mp_float ret_im(0, p);
	mpfr_mul_d(mmd::impl_value(ret_re), mmd::impl_value(a), real(b), MPFR_RNDN);
    mpfr_mul_d(mmd::impl_value(ret_im), mmd::impl_value(a), imag(b), MPFR_RNDN);
    ret_re.update_debug();
    ret_im.update_debug();
	return mp_complex(ret_re, ret_im);
};
mp_complex details::mul_impl(const Complex& a, const mp_float& b, precision p)
{
	return mul_impl(b,a,p);
};
mp_complex details::mul_impl(const mp_float& a, const Float_complex& b, precision p)
{
    return mul_impl(a,Complex(b),p);
}
mp_complex details::mul_impl(const Float_complex& a, const mp_float& b, precision p)
{
    return mul_impl(Complex(a), b,p);
}
mp_float details::mul_impl(const mp_float& a, const mp_int& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_z(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::mul_impl(const mp_int& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_z(mmd::impl_value(ret), mmd::impl_value(b), 
               mmd::impl_value(a).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::mul_impl(const mp_float& a, const mp_rational& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_mul_q(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::mul_impl(const mp_rational& a, const mp_float& b, precision p)
{
    return mul_impl(b, a, p);
};

template<class Compl_2>
struct recover_nan_mult
{
    static mp_complex eval(const mp_complex& a, const Compl_2& b,
                           const mp_float& r_re, const mp_float& r_im)
    {
        //code taken from C99 standard; ISO/IEC 9899:TC2; WG14/N1124
        //notice that correct inf/nan handling is not required neither
        //by C not by C++ standard

        bool recalc     = false;
        precision sp    = precision(10);
        precision p     = r_re.get_precision();

        mp_float a_re   = real(a);        
        mp_float b_re   = real(b);        
        mp_float a_im   = imag(a);
        mp_float b_im   = imag(b);

        if (is_inf(a_re) || is_inf(b_re) ) 
        { 
            // z is infinite
            a_re        = copysign(mp_float(is_inf(a_re) ? 1.0 : 0.0, sp), a_re);
            a_im        = copysign(mp_float(is_inf(a_im) ? 1.0 : 0.0, sp), a_im);

            if (is_nan(b_re)) 
                b_re    = copysign(mp_float(0.0,sp), b_re);
            if (is_nan(b_im)) 
                b_im    = copysign(mp_float(0.0,sp), b_im);

            recalc = true;
        }

        if (is_inf(b_re) || is_inf(b_im) ) 
        { 
            // w is infinite
            b_re        = copysign(mp_float(is_inf(b_re) ? 1.0 : 0.0, sp), b_re);
            b_im        = copysign(mp_float(is_inf(b_im) ? 1.0 : 0.0, sp), b_im);

            if (is_nan(a_re)) 
                a_re    = copysign(mp_float(0.0,sp), a_re);
            if (is_nan(a_im)) 
                a_im    = copysign(mp_float(0.0,sp), a_im);

            recalc = true;
        }

        if (recalc == false)
        {
            if (is_inf(mul(a_re, b_re, sp)) == true || is_inf(mul(a_im, b_im, sp)) ==  true 
                || is_inf(mul(a_re, b_im, sp)) == true || is_inf(mul(a_im, b_re, sp)) == true)
            {
                // recover infinities from overflow by changing NaNs to 0
                if (is_nan(a_re)) 
                    a_re    = copysign(mp_float(0.0, sp), a_re);
                if (is_nan(a_im)) 
                    a_im    = copysign(mp_float(0.0, sp), a_im);
                if (is_nan(b_re))
                    b_re    = copysign(mp_float(0.0, sp), b_re);
                if (is_nan(b_im))
                    b_im    = copysign(mp_float(0.0, sp), b_im);

                recalc = true;
            }
        }

        if (recalc) 
        {
            mp_float r_re2  = dot2_a(a_re, b_re, -a_im, b_im, sp);
            mp_float r_im2  = dot2_a(a_re, b_im,  a_im, b_re, sp);
            r_re2           = mul(constants::mp_inf(sp), r_re2, p);
            r_im2           = mul(constants::mp_inf(sp), r_im2, p);

            return mp_complex(r_re2, r_im2);
        }
        else
        {
            return mp_complex(r_re, r_im, p);
        }
    };
};

mp_complex details::mul_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    p               = mmd::result_prec(p, get_precision(a), get_precision(b));

    mp_float r_re   = dot2_a(real(a), real(b), -imag(a), imag(b), p);
    mp_float r_im   = dot2_a(real(a), imag(b),  imag(a), real(b), p);

    if (is_nan(r_re) == true && is_nan(r_im) == true)
        return recover_nan_mult<mp_complex>::eval(a, b, r_re, r_im);
    else
        return mp_complex(r_re, r_im);
}
mp_complex details::mul_impl(const mp_complex& a, Integer b, precision p)
{
    return mul_impl(a, Real(b), p);
}
mp_complex details::mul_impl(Integer a, const mp_complex& b, precision p)
{
    return mul_impl(Real(a), b, p);
}
mp_complex details::mul_impl(const mp_complex& a, Real b, precision p)
{
    return mp_complex(mul_impl(real(a),b,p), mul_impl(imag(a),b,p));
}
mp_complex details::mul_impl(Real a, const mp_complex& b, precision p)
{
    return mul_impl(b,a,p);
}
mp_complex details::mul_impl(const mp_complex& a, const Complex& b, precision p)
{
    p               = mmd::result_prec(p, get_precision(a), get_precision(b));
    mp_float r_re   = dot2_a(real(a), real(b), -imag(a), imag(b), p);
    mp_float r_im   = dot2_a(real(a), imag(b),  imag(a), real(b), p);

    if (is_nan(r_re) == true && is_nan(r_im) == true)
        return recover_nan_mult<Complex>::eval(a, b, r_re, r_im);
    else
        return mp_complex(r_re, r_im);
}
mp_complex details::mul_impl(const Complex& a, const mp_complex& b, precision p)
{
    return mul_impl(b,a,p);
}
mp_complex details::mul_impl(const mp_complex& a, const Float_complex& b, precision p)
{
    return mul_impl(a, Complex(b), p);
}
mp_complex details::mul_impl(const Float_complex& a, const mp_complex& b, precision p)
{
    return mul_impl(Complex(a), b, p);
}
mp_complex details::mul_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}
mp_complex details::mul_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return mp_complex(mul_impl(real(a),b,p), mul_impl(imag(a),b,p));
}
mp_complex details::mul_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return mul_impl(b, a, p);
}
mp_complex details::mul_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return mul_impl(a, mp_float(b), p);
}
mp_complex details::mul_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return mul_impl(mp_float(a), b, p);
}

//---------------------------------------------------------------
//                  div_impl
//---------------------------------------------------------------
mp_rational details::div_impl(const mp_int& a, const mp_int& b, precision)
{    
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    return mp_rational(a, b);
}
mp_rational details::div_impl(const mp_int& a, Integer b, precision)
{
    if (b == 0)
        return mp_rational(0);
    
    return mp_rational(a, b);
}
mp_rational details::div_impl(Integer a, const mp_int& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    return mp_rational(a, b);
}
mp_float details::div_impl(const mp_int& a, Real b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_float details::div_impl(Real a, const mp_int& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_int& a, const Complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_int& a, const Float_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Float_complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a,b,p);
}

mp_rational details::div_impl(const mp_rational& a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_rational details::div_impl(const mp_rational& a, Integer b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_rational details::div_impl(Integer a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_float details::div_impl(const mp_rational& a, Real b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_float details::div_impl(Real a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_rational& a, const Complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_rational& a, const Float_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Float_complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_rational details::div_impl(const mp_rational& a, const mp_int& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_rational details::div_impl(const mp_int& a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}

mp_float details::div_impl(const mp_float& a, const mp_float& b, precision p)
{
	return idiv_impl(a,b,p);
};
mp_float details::div_impl(const mp_float& a, Integer b, precision p)
{
	return idiv_impl(a,b,p);
}
mp_float details::div_impl(Integer a, const mp_float& b, precision p)
{
	return idiv_impl(a,b,p);
}
mp_float details::div_impl(const mp_float& a, Real b, precision p)
{
	return idiv_impl(a,b,p);
};
mp_float details::div_impl(Real a, const mp_float& b, precision p)
{
	return idiv_impl(a,b,p);
};
mp_complex details::div_impl(const mp_float& a, const Complex& b, precision p)
{
    return idiv_impl(a,b,p);
};
mp_complex details::div_impl(const Complex& a, const mp_float& b, precision p)
{
    return idiv_impl(a,b,p);
};
mp_complex details::div_impl(const mp_float& a, const Float_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Float_complex& a, const mp_float& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_float details::div_impl(const mp_float& a, const mp_int& b, precision p)
{
	return idiv_impl(a,b,p);
};
mp_float details::div_impl(const mp_int& a, const mp_float& b, precision p)
{
    return idiv_impl(a,b,p);
};
mp_float details::div_impl(const mp_float& a, const mp_rational& b, precision p)
{
	return idiv_impl(a,b,p);
};
mp_float details::div_impl(const mp_rational& a, const mp_float& b, precision p)
{
	return idiv_impl(a,b,p);
};

mp_complex details::div_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, Integer b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(Integer a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, Real b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(Real a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, const Complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Complex& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, const Float_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const Float_complex& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a,b,p);
}
mp_complex details::div_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return idiv_impl(a,b,p);
}

//---------------------------------------------------------------
//                  idiv_impl
//---------------------------------------------------------------
namespace impl
{

struct less_abs
{
    static bool eval(const mp_float& a, const mp_float& b)
    {
        return compare_abs(a, b) == cmp_type::less;
    };
};

inline mp_complex make_nan(precision prec)
{
    mp_float nan = constants::mp_nan(prec);
    return mp_complex(nan, nan);
};

template<class T1, class T2>
struct div_complex_complex_helper
{
    static mp_complex recover_nan_div(const mp_float& a_re, const mp_float& a_im, 
                            const mp_float& b_re, const mp_float& b_im, 
                            const mp_float& r_re, const mp_float& r_im)
    {
        precision p     = r_re.get_precision();
        precision sp    = precision(10);

        if (is_zero(b_re) && is_zero(b_im) && (is_nan(a_re) == false || is_nan(a_im) == false)) 
        {
            mp_float inf    = constants::mp_inf(p);
            mp_float ret_re = mul(copysign(inf, b_re), a_re, p);
            mp_float ret_im = mul(copysign(inf, b_re), a_im, p);

            return mp_complex(ret_re, ret_im, p);
        }
        else if ((is_inf(a_re) == true || is_inf(a_im) == true) 
                 && is_finite(b_re) && is_finite(b_im)) 
        {
            mp_float one    = mp_float(1.0, sp);
            mp_float zero   = mp_float(0.0, sp);
            mp_float inf    = constants::mp_inf(sp);

            mp_float a_re2  = copysign(is_inf(a_re) ? one : zero, a_re);
            mp_float a_im2  = copysign(is_inf(a_im) ? one : zero, a_im);            

            mp_float ret_re = inf * dot2_a(a_re2, b_re, a_im2, b_im, sp);
            mp_float ret_im = inf * dot2_a(a_im2, b_re, -a_re2, b_im, sp);

            return mp_complex(ret_re, ret_im, p);
        }
        else if ((is_inf(b_re) == true || is_inf(b_im) == true) 
                 && is_finite(a_re) == true && is_finite(a_im) == true 
                 && is_nan(b_re) == false && is_nan(b_im) == false) 
        {
            mp_float one    = mp_float(1.0, sp);
            mp_float zero   = mp_float(0.0, sp);
            mp_float inf    = constants::mp_inf(sp);

            mp_float b_re2  = copysign(is_inf(b_re) ? one : zero, b_re);
            mp_float b_im2  = copysign(is_inf(b_im) ? one : zero, b_im);
            mp_float ret_re = zero * dot2_a(a_re, b_re2, a_im, b_im2, sp);
            mp_float ret_im = zero * dot2_a(a_im, b_re2, -a_re, b_im2, sp);

            return mp_complex(ret_re, ret_im, p);
        }

        return mp_complex(r_re, r_im, p);
    };

    static mp_complex eval(const T1& a, const T2& b, precision p0)
    {
        precision prec      = mmd::result_prec(p0, get_precision(a), get_precision(b));
        precision ip        = extend_prec_div_cc(prec);

        const mp_float& a_re    = real(a);
        const mp_float& a_im    = imag(a);
        mp_float b_re           = real(b);
        mp_float b_im           = imag(b);

        int sing_a_re           = signbit(a_re) ? -1 : 1;
        int sing_a_im           = signbit(a_im) ? -1 : 1;
        int sing_b_re           = signbit(b_re) ? -1 : 1;
        int sing_b_im           = signbit(b_im) ? -1 : 1;
        int sing_test           = sing_a_re * sing_a_im * sing_b_re * sing_b_im;

        mp_float r_re, r_im;

        //scale in order to avoid overflows
        mp_float scale          = logb(std::max(abs(b_re), abs(b_im)));
        Integer iscale          = 0;

        if (is_finite(scale) == true) 
        {
            iscale              = scale.cast_int();
            b_re                = ldexp(b_re, -iscale); 
            b_im                = ldexp(b_im, -iscale);
        }

        mp_float d              = abs2(b_re, ip) + abs2(b_im, ip);  //1.0 ulp
        mp_float di             = div(1, d, ip);                    //1.5 ulp

        if (sing_test == 1)
        {
            //in a_re * b_re + a_im * b_im there are no cancelations
            r_re                = mul(a_re, b_re, ip) + mul(a_im, b_im, ip);    //1 ulp
            r_im                = dot2_a(a_re, -b_im, a_im, b_re, ip);          //1 ulp            
        }
        else
        {
            //in -a_re * b_im + a_im * b_re there are no cancelations
            r_re                = dot2_a(a_re, b_re, a_im, b_im, ip);           //1 ulp
            r_im                = mul(a_im, b_re, ip) - mul(a_re, b_im, ip);    //1 ulp
        };
        
        r_re                    = mul(r_re, di, prec);  //3 ulp
        r_im                    = mul(r_im, di, prec);  //3 ulp

        if (iscale != 0)
        {
            r_re                = ldexp(r_re, -iscale);
            r_im                = ldexp(r_im, -iscale);
        };

        if (is_nan(r_re) == true && is_nan(r_im) == true)
            return recover_nan_div(a_re, a_im, b_re, b_im, r_re, r_im);
        else
            return mp_complex(r_re, r_im);
    };
};

template<class T1, class T2>
struct div_real_complex_helper
{
    static mp_complex recover_nan(const mp_float& a_re, const mp_float& b_re, const mp_float& b_im,
                                  const mp_float& r_re, const mp_float& r_im)
    {
        precision p     = r_re.get_precision();
        precision sp    = precision(10);

        if ((is_zero(b_re) == true && is_zero(b_im) == true) && is_nan(a_re) == false)
        {
            mp_float inf    = constants::mp_inf(sp);
            mp_float ret_re = mul(copysign(inf, b_re), a_re, p);
            mp_float ret_im = mp_float(p);

            return mp_complex(ret_re, ret_im, p);
        }
        else if (is_inf(a_re) == true && is_finite(b_re) == true && is_finite(b_im) == true) 
        {
            mp_float a_re2  = copysign(mp_float(is_inf(a_re) ? 1.0 : 0.0, sp), a_re);
            mp_float inf    = constants::mp_inf(sp);

            mp_float ret_re = inf * mul(a_re2, b_re, sp);
            mp_float ret_im = inf * mul(-a_re2, b_im, sp);

            return mp_complex(ret_re, ret_im, p);
        }
        else if ((is_inf(b_re) == true || is_inf(b_im) == true) 
                 && is_finite(a_re) == true
                 && is_nan(b_re) == false && is_nan(b_im) == false) 
        {
            mp_float b_re2  = copysign(mp_float(is_inf(b_re) ? 1.0 : 0.0, sp), b_re);
            mp_float b_im2  = copysign(mp_float(is_inf(b_im) ? 1.0 : 0.0, sp), b_im);
            
            mp_float ret_re = mp_float(0.0, sp) * mul(a_re, b_re2, sp);
            mp_float ret_im = mp_float(0.0, sp) * mul(-a_re, b_im2, sp);

            return mp_complex(ret_re, ret_im, p);
        }

        return mp_complex(r_re, r_im, p);
    };

    static mp_complex eval(const T1& a, const T2& b, precision p)
    {
        precision a_prec    = get_precision(a);
        precision b_prec    = get_precision(b);
        precision prec      = mmd::result_prec(p, a_prec, b_prec);
        precision int_prec  = extend_prec_div_rc(prec);

	    const mp_float& b_real = real(b);
        const mp_float& b_imag = imag(b);

        mp_float ret_re, ret_im;

	    if (less_abs::eval(b_imag, b_real) == true)
	    {
		    // |imag(b)| < |real(b)|

		    mp_float Wr = div(b_imag, b_real, int_prec);            // 0.5 ulp
		    mp_float Wd = mul(Wr, b_imag, int_prec);                // 1.0 ulp
            Wd          = plus(b_real, std::move(Wd), int_prec);    // 1.5 ulp

		    // compute representable result
            mp_float aW = div(a, Wd, int_prec);                     // 2.0 ulp
			ret_re      = mp_float(aW, prec);
			ret_im      = -mul(aW, Wr, prec);                       // 3.5 ulp
	    }
        else
        {
		    // |real(b)| <= |imag(b)|
		    mp_float Wr = div(b_real, b_imag, int_prec);
		    mp_float Wd = mul(Wr, b_real, int_prec);
            Wd          = plus(b_imag, std::move(Wd), int_prec);

            // compute representable result
            mp_float aW = div(a, Wd, int_prec);
		    ret_re      = mul(aW, Wr, prec);
		    ret_im      = uminus(aW, prec);
        };

        if (is_nan(ret_re) == true && is_nan(ret_im) == true)
            return recover_nan(a, b_real, b_imag, ret_re, ret_im);
        else
            return mp_complex(ret_re, ret_im, prec);
    };

    static mp_complex eval_inv(const T2& b, precision p)
    {
        return eval(mp_float(1), b, p);
    };
};

};

mp_int details::idiv_impl(const mp_int& a, const mp_int& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_int(0);
    }

    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_int details::idiv_impl(const mp_int& a, Integer b, precision)
{
    if (b == 0)
        return mp_int(0);
    
    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_int details::idiv_impl(Integer a, const mp_int& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_int(0);
    }

    using impl_type = mmd::impl_int;
    return mp_int(mp_int::impl_tag(), impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_float details::idiv_impl(const mp_int& a, Real b, precision p)
{
    return div_impl(mp_float(a), b, p);
}
mp_float details::idiv_impl(Real a, const mp_int& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_int& a, const Complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_complex details::idiv_impl(const Complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_int& a, const Float_complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_complex details::idiv_impl(const Float_complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}

mp_rational details::idiv_impl(const mp_rational& a, const mp_rational& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_rational details::idiv_impl(const mp_rational& a, Integer b, precision)
{
    if (b == 0)
        return mp_rational(0);

    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_rational details::idiv_impl(Integer a, const mp_rational& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_float details::idiv_impl(const mp_rational& a, Real b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_float details::idiv_impl(Real a, const mp_rational& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_rational& a, const Complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_complex details::idiv_impl(const Complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_rational& a, const Float_complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_complex details::idiv_impl(const Float_complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_rational details::idiv_impl(const mp_rational& a, const mp_int& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}
mp_rational details::idiv_impl(const mp_int& a, const mp_rational& b, precision)
{
    if (is_zero(b) == true)
    {
        error_flags::set_integer_overflow();
        return mp_rational(0);
    }

    using impl_type = mmd::impl_rat;
    return mp_rational(mp_rational::impl_tag(), 
                impl_type(mmd::impl_value(a)/mmd::impl_value(b)));
}

mp_float details::idiv_impl(const mp_float& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());
	mp_float ret(0, p);
	mpfr_div(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::idiv_impl(const mp_float& a, Integer b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_div_si(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::idiv_impl(Integer a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_si_div(mmd::impl_value(ret), a, mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
}
mp_float details::idiv_impl(const mp_float& a, Real b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_div_d(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::idiv_impl(Real a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_d_div(mmd::impl_value(ret), a, mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_complex details::idiv_impl(const mp_float& a, const Complex& b, precision p)
{
    return impl::div_real_complex_helper<mp_float, Complex>::eval(a,b,p);
};
mp_complex details::idiv_impl(const Complex& a, const mp_float& b, precision p)
{
    return mp_complex(idiv_impl(real(a),b,p), idiv_impl(imag(a),b,p));
};
mp_complex details::idiv_impl(const mp_float& a, const Float_complex& b, precision p)
{
    return idiv_impl(a,Complex(b),p);
}
mp_complex details::idiv_impl(const Float_complex& a, const mp_float& b, precision p)
{
    return idiv_impl(Complex(a), b,p);
}
mp_float details::idiv_impl(const mp_float& a, const mp_int& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_div_z(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::idiv_impl(const mp_int& a, const mp_float& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
};
mp_float details::idiv_impl(const mp_float& a, const mp_rational& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));
	mp_float ret(0, p);
	mpfr_div_q(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();
	return ret;
};
mp_float details::idiv_impl(const mp_rational& a, const mp_float& b, precision p)
{
	return idiv_impl(mp_float(a), b, p);
};

mp_complex details::idiv_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    return impl::div_complex_complex_helper<mp_complex, mp_complex>::eval(a,b,p);
}
mp_complex details::idiv_impl(const mp_complex& a, Integer b, precision p)
{
    return idiv_impl(a, Real(b), p);
}
mp_complex details::idiv_impl(Integer a, const mp_complex& b, precision p)
{
    return idiv_impl(Real(a), b, p);
}
mp_complex details::idiv_impl(const mp_complex& a, Real b, precision p)
{
    return mp_complex(idiv_impl(real(a),b,p), idiv_impl(imag(a),b,p));
}
mp_complex details::idiv_impl(Real a, const mp_complex& b, precision p)
{
    return impl::div_real_complex_helper<Real, mp_complex>::eval(a,b, p);
}
mp_complex details::idiv_impl(const mp_complex& a, const Complex& b, precision p)
{
    return impl::div_complex_complex_helper<mp_complex, Complex>::eval(a,b, p);
}
mp_complex details::idiv_impl(const Complex& a, const mp_complex& b, precision p)
{
    return impl::div_complex_complex_helper<Complex, mp_complex>::eval(a,b, p);
}
mp_complex details::idiv_impl(const mp_complex& a, const Float_complex& b, precision p)
{
    return idiv_impl(a, Complex(b), p);
}
mp_complex details::idiv_impl(const Float_complex& a, const mp_complex& b, precision p)
{
    return idiv_impl(Complex(a), b, p);
}
mp_complex details::idiv_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}
mp_complex details::idiv_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return mp_complex(idiv_impl(real(a),b,p), idiv_impl(imag(a),b,p));
}
mp_complex details::idiv_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return impl::div_real_complex_helper<mp_float, mp_complex>::eval(a,b,p);
}
mp_complex details::idiv_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return idiv_impl(a, mp_float(b), p);
}
mp_complex details::idiv_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return idiv_impl(mp_float(a), b, p);
}

mp_complex details::inv_impl(const mp_complex& a, precision p)
{
    return impl::div_real_complex_helper<mp_float, mp_complex>::eval_inv(a,p);
};

//---------------------------------------------------------------
//                  pow_impl
//---------------------------------------------------------------
mp_float details::pow_impl(const mp_int& a, const mp_int& b, precision p)
{
    //integer powers overflow too easily, therefore not allowed
    return pow_impl(mp_float(a), mp_float(b), p);
}

mp_float details::pow_impl(const mp_int& a, Integer b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
};
mp_float details::pow_impl(Integer a, const mp_int& b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}
mp_float details::pow_impl(const mp_rational& a, Integer b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}
mp_float details::pow_impl(const mp_rational& a, const mp_int& b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}

mp_float details::pow_impl(const mp_rational& a, const mp_rational& b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}
mp_float details::pow_impl(Integer a, const mp_rational& b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}
mp_float details::pow_impl(const mp_int& a, const mp_rational& b, precision p)
{
    return pow_impl(mp_float(a), mp_float(b), p);
}

mp_float details::pow_impl(const mp_float& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());

    mp_float ret(p);
    mpfr_pow(mmd::impl_value(ret), mmd::impl_value(a), mmd::impl_value(b), MPFR_RNDN);
    ret.update_debug();

    return ret;
}
mp_float details::pow_impl(const mp_float& a, Integer b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));

    mp_float ret(p);
    mpfr_pow_si(mmd::impl_value(ret), mmd::impl_value(a), b, MPFR_RNDN);
    ret.update_debug();

    return ret;
}
mp_float details::pow_impl(Integer a, const mp_float& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}
mp_float details::pow_impl(const mp_float& a, const mp_int& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));

    mp_float ret(p);
    mpfr_pow_z(mmd::impl_value(ret), mmd::impl_value(a), 
               mmd::impl_value(b).backend().data(), MPFR_RNDN);
    ret.update_debug();

    return ret;
}
mp_float details::pow_impl(const mp_int& a, const mp_float& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}
mp_float details::pow_impl(const mp_float& a, const mp_rational& b, precision p)
{
    return pow_impl(a, mp_float(b), p);
}
mp_float details::pow_impl(const mp_rational& a, const mp_float& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}

mp_complex details::pow_impl(Integer a, const mp_complex& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}
mp_complex details::pow_impl(const mp_int& a, const mp_complex& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}
mp_complex details::pow_impl(const mp_complex& a, const mp_rational& b, precision p)
{
    return pow_impl(a, mp_float(b), p);
}
mp_complex details::pow_impl(const mp_rational& a, const mp_complex& b, precision p)
{
    return pow_impl(mp_float(a), b, p);
}
mp_complex details::pow_impl(const mp_complex& a, Integer b, precision p)
{
    return pow_impl(a, mp_float(b), p);
}
mp_complex details::pow_impl(const mp_complex& a, const mp_int& b, precision p)
{
    return pow_impl(a, mp_float(b), p);
}
mp_complex details::pow_impl(const mp_float& a, const mp_complex& b, precision p)
{
    return pow_real_compl(a, b, p);
}
mp_complex details::pow_impl(const mp_complex& a, const mp_float& b, precision p)
{
    return pow_compl_real(a, b, p);
}
mp_complex details::pow_impl(const mp_complex& a, const mp_complex& b, precision p)
{
    return pow_compl_compl(a, b, p);
}

}}}
