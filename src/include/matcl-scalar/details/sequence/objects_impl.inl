/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017-2018
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

#include "matcl-scalar/lib_functions/sequence/objects_impl.h"

namespace matcl
{

//-----------------------------------------------------------------------
//                   extrapolator_wynn_epsilon
//-----------------------------------------------------------------------

template<class Float>
extrapolator_wynn_epsilon<Float>::extrapolator_wynn_epsilon(int precision)
    : m_impl(precision)
{};

template<class Float>
void extrapolator_wynn_epsilon<Float>::eval(const Float& x, const Float& val, Float& lim, 
                                            Float& err)
{
    (void)x;
    return m_impl.eval(val, lim, err);
}

template<class Float>
void extrapolator_wynn_epsilon<Float>::eval(const Float& x, const Float& val, 
                            const Float& val_abs_err, Float& lim, Float& err)
{
    (void)x;
    return m_impl.eval(val, val_abs_err, lim, err);
}

template<class Float>
const Float& extrapolator_wynn_epsilon<Float>::last_result() const
{
    return m_impl.last_result();
}

template<class Float>
const Float& extrapolator_wynn_epsilon<Float>::numerical_error() const
{
    return m_impl.numerical_error();
}

template<class Float>
const Float& extrapolator_wynn_epsilon<Float>::limit_error() const
{
    return m_impl.limit_error();
}

template<class Float>
Float extrapolator_wynn_epsilon<Float>::total_error() const
{
    return m_impl.total_error();
}

template<class Float>
void extrapolator_wynn_epsilon<Float>::clear(int precision)
{
    return m_impl.clear(precision);
};

//-----------------------------------------------------------------------
//                   extrapolator_aitken_delta
//-----------------------------------------------------------------------

template<class Float>
extrapolator_aitken_delta<Float>::extrapolator_aitken_delta(int precision)
    : m_impl(precision)
{};

template<class Float>
void extrapolator_aitken_delta<Float>::eval(const Float& x, const Float& val, Float& lim, 
                                            Float& err)
{
    (void)x;
    return m_impl.eval(val, lim, err);
}

template<class Float>
void extrapolator_aitken_delta<Float>::eval(const Float& x, const Float& val, 
                            const Float& val_abs_err, Float& lim, Float& err)
{
    (void)x;
    return m_impl.eval(val, val_abs_err, lim, err);
}

template<class Float>
const Float& extrapolator_aitken_delta<Float>::last_result() const
{
    return m_impl.last_result();
}

template<class Float>
const Float& extrapolator_aitken_delta<Float>::numerical_error() const
{
    return m_impl.numerical_error();
}

template<class Float>
const Float& extrapolator_aitken_delta<Float>::limit_error() const
{
    return m_impl.limit_error();
}

template<class Float>
Float extrapolator_aitken_delta<Float>::total_error() const
{
    return m_impl.total_error();
}

template<class Float>
void extrapolator_aitken_delta<Float>::clear(int precision)
{
    return m_impl.clear(precision);
};

//-----------------------------------------------------------------------
//                   extrapolator_richardson
//-----------------------------------------------------------------------

template<class Float>
extrapolator_richardson<Float>::extrapolator_richardson(int precision)
    : m_impl(precision)
{};

template<class Float>
void extrapolator_richardson<Float>::eval(const Float& x, const Float& val, Float& lim, 
                                            Float& err)
{
    return m_impl.eval(x, val, lim, err);
}

template<class Float>
void extrapolator_richardson<Float>::eval(const Float& x, const Float& val, 
                            const Float& val_abs_err, Float& lim, Float& err)
{
    return m_impl.eval(x, val, val_abs_err, lim, err);
}

template<class Float>
const Float& extrapolator_richardson<Float>::last_result() const
{
    return m_impl.last_result();
}

template<class Float>
const Float& extrapolator_richardson<Float>::numerical_error() const
{
    return m_impl.numerical_error();
}

template<class Float>
const Float& extrapolator_richardson<Float>::limit_error() const
{
    return m_impl.limit_error();
}

template<class Float>
Float extrapolator_richardson<Float>::total_error() const
{
    return m_impl.total_error();
}

template<class Float>
void extrapolator_richardson<Float>::clear(int precision)
{
    return m_impl.clear(precision);
};

};

namespace matcl
{

template<class Float>
extrapolator_ptr<Float> matcl::make_extrapolator_wynn_epsilon(int precision)
{
    return extrapolator_ptr<Float>(new extrapolator_wynn_epsilon<Float>(precision));
}

template<class Float>
extrapolator_ptr<Float> matcl::make_extrapolator_aitken_delta(int precision)
{
    return extrapolator_ptr<Float>(new extrapolator_aitken_delta<Float>(precision));
}

template<class Float>
extrapolator_ptr<Float> matcl::make_extrapolator_richardson(int precision)
{
    return extrapolator_ptr<Float>(new extrapolator_richardson<Float>(precision));
}

}