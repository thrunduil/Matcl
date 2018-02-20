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

#include "matcl-blas/level1/level1_basic.h"

namespace matcl { namespace level1 { namespace details
{

template<class T1>
struct func_set_val_zero
{    
    template<Integer N>
    force_inline
    void eval(T1* Y1, Integer rows) const
    {
        return set_val<T1,N>::eval(Y1, rows, T1(0.0));
    };
};

template<class T1>
struct func_set_val
{    
    const T1& a;

    func_set_val(const T1& a_)  :a(a_){};

    template<Integer N>
    force_inline
    void eval(T1* Y1, Integer rows) const
    {
        return set_val<T1,N>::eval(Y1, rows, a);
    };
};

template<class TY, class TX, class TA, class TB>
struct func_axpby
{ 
    const TA& a;
    const TB& b;

    func_axpby(const TA& a_, const TB& b_) :a(a_), b(b_){};

    template<Integer N>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        //Y = a*X + bY
        return axpby<TY, TX, TA, TB, N>::eval(Y1, X1, rows, a, b);
    };
};
template<class T1, class T2>
struct func_copy
{ 
    template<Integer N>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return copy<T1, T2, N>::eval(Y1, X1, rows);
    };
};
template<class T1, class T2>
struct func_mx
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return mx<T1, T2, Rows>::eval(Y1, X1, rows);
    };
};
template<class TY, class TX, class TA>
struct func_ax
{ 
    const TA& a;

    func_ax(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return ax<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class T1, class T2>
struct func_ypx
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return ypx<T1, T2, Rows>::eval(Y1, X1, rows);
    };
};

template<class T1, class T2>
struct func_ymx
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return ymx<T1, T2, Rows>::eval(Y1, X1, rows);
    };
};
template<class TY, class TX, class TA>
struct func_axpy
{ 
    const TA& a;

    func_axpy(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return axpy<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};
template<class T1>
struct func_my
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, Integer rows) const
    {
        return my<T1, Rows>::eval(Y1, rows);
    };
};

template<class T1, class T2>
struct func_xmy
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return xmy<T1, T2, Rows>::eval(Y1, X1, rows);
    };
};

template<class T1, class T2>
struct func_ypxm
{ 
    template<Integer Rows>
    force_inline
    void eval(T1* Y1, const T2* X1, Integer rows) const
    {
        return ypxm<T1, T2, Rows>::eval(Y1, X1, rows);
    };
};

template<class TY, class TX, class TA>
struct func_axmy
{ 
    const TA& a;

    func_axmy(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return axmy<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class TY, class TA>
struct func_ay
{ 
    const TA& a;

    func_ay(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, Integer rows) const
    {
        return ay<TY, TA, Rows>::eval(Y1, rows, a);
    };
};

template<class TY, class TX, class TA>
struct func_xpby
{ 
    const TA& a;

    func_xpby(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return xpby<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class TY, class TX, class TA>
struct func_mxpby
{ 
    const TA& a;

    func_mxpby(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return mxpby<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class TY, class TX, class TA>
struct func_xpya
{ 
    const TA& a;

    func_xpya(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return xpya<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class TY, class TX, class TA>
struct func_xmya
{ 
    const TA& a;

    func_xmya(const TA& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return xmya<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class T1>
struct func_apy
{ 
    const T1& a;

    func_apy(const T1& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(T1* Y1, Integer rows) const
    {
        return apy<T1, Rows>::eval(Y1, rows, a);
    };
};

template<class T1>
struct func_amy
{ 
    const T1& a;

    func_amy(const T1& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(T1* Y1, Integer rows) const
    {
        return amy<T1, Rows>::eval(Y1, rows, a);
    };
};

template<class TY, class TB>
struct func_apby
{ 
    const TY& a;
    const TB& b;

    func_apby(const TY& a_, const TB& b_) : a(a_), b(b_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, Integer rows) const
    {
        return apby<TY, TB, Rows>::eval(Y1, rows, a, b);
    };
};

template<class TY, class TX>
struct func_apx_abs
{ 
    const TY& a;

    func_apx_abs(const TY& a_) : a(a_){};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return apx_abs<TY, TX, Rows>::eval(Y1, X1, rows, a);
    };
};

template<class TY, class TX>
struct func_x_abs
{ 
    func_x_abs() {};

    template<Integer Rows>
    force_inline
    void eval(TY* Y1, const TX* X1, Integer rows) const
    {
        return x_abs<TY, TX, Rows>::eval(Y1, X1, rows);
    };
};

}}}