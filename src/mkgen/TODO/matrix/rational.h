#pragma once

#include "mkgen/mkgen_fwd.h"

namespace matcl { namespace mkgen
{

template<Integer A, Integer B>
struct gdc
{
    static const Integer value          = gdc<B, A % B>::value;
};

template<Integer A>
struct gdc<A,0>
{
    static const Integer value          = A;
};

template<Integer A>
struct gdc<A,1>
{
    static const Integer value          = 1;
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct rational_plus
{
    static const Integer nom            = N1*D2 + N2*D1;
    static const Integer den            = D1*D2;
    static const Integer d              = gdc<nom,den>::value;

    static const Integer nominator      = nom/d;
    static const Integer denominator    = den/d;
};

template<Integer N1, Integer D1, Integer N2, Integer D2>
struct rational_mult
{
    static const Integer nom            = N1*N2;
    static const Integer den            = D1*D2;
    static const Integer d              = gdc<nom,den>::value;

    static const Integer nominator      = nom/d;
    static const Integer denominator    = den/d;
};

}}