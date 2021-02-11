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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

real_dense              range(Real s, Real i, Real e);
inline real_dense       range(Real s, Real e)                       { return range(s, 1., e); }

integer_dense           irange(Integer s, Integer i, Integer e);
inline integer_dense    irange(Integer s, Integer e)                { return irange(s, 1, e); }

float_dense             frange(Float s, Float i, Float e);
inline float_dense      frange(Float s, Float e)                    { return frange(s, 1.0f, e); }

real_dense              linspace(Real s, Real e, Integer n);
float_dense             flinspace(Float s, Float e, Integer n);
real_dense              logspace(Real s, Real e, Integer n);
float_dense             flogspace(Float s, Float e, Integer n);

inline real_dense zeros(Integer r, Integer c)             
{ 
    real_dense out(ti::ti_empty(), 0., r, c); 
    out.get_struct().set(predefined_struct_type::diag); 
    return out; 
}

inline real_dense zeros(Integer n)
{
    real_dense out(ti::ti_empty(),0., n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline object_dense zeros(ti::ti_object ti,Integer r, Integer c)
{
    object_dense out(ti, md::default_value<Object>(ti), r, c);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline object_dense zeros(ti::ti_object ti,Integer n)
{
    object_dense out(ti,md::default_value<Object>(ti), n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline integer_dense izeros(Integer r, Integer c)               
{ 
    integer_dense out(ti::ti_empty(),0, r, c);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline integer_dense izeros(Integer n)
{
    integer_dense out(ti::ti_empty(),0, n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline float_dense fzeros(Integer r, Integer c)               
{ 
    float_dense out(ti::ti_empty(), 0.0f, r, c);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline float_dense fzeros(Integer n)
{
    float_dense out(ti::ti_empty(),0.0f, n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline complex_dense czeros(Integer r, Integer c)
{
    complex_dense out(ti::ti_empty(),0., r, c);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline complex_dense czeros(Integer n)
{
    complex_dense out(ti::ti_empty(), 0., n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline float_complex_dense fczeros(Integer r, Integer c)
{
    float_complex_dense out(ti::ti_empty(),0., r, c);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline float_complex_dense fczeros(Integer n)
{
    float_complex_dense out(ti::ti_empty(), 0., n, n);
    out.get_struct().set(predefined_struct_type::diag);
    return out;
}

inline real_dense       ones(Integer r, Integer c)      { return real_dense(ti::ti_empty(),1., r, c); }
inline real_dense       ones(Integer n)                 { return real_dense(ti::ti_empty(),1., n, n); }
inline integer_dense    iones(Integer r, Integer c)     { return integer_dense(ti::ti_empty(),1, r, c); }
inline integer_dense    iones(Integer n)                { return integer_dense(ti::ti_empty(),1, n, n); }
inline complex_dense    cones(Integer r, Integer c)     { return complex_dense(ti::ti_empty(), 1., r, c); }
inline complex_dense    cones(Integer n)                { return complex_dense(ti::ti_empty(), 1., n, n); }
inline float_dense      fones(Integer r, Integer c)     { return float_dense(ti::ti_empty(), 1.f, r, c); }
inline float_dense      fones(Integer n)                { return float_dense(ti::ti_empty(), 1.f, n, n); }
inline float_complex_dense
                        fcones(Integer r, Integer c)    { return float_complex_dense(ti::ti_empty(), 1.f, r, c); }
inline float_complex_dense
                        fcones(Integer n)               { return float_complex_dense(ti::ti_empty(), 1.f, n, n); }

inline object_dense ones(ti::ti_object ti,Integer r, Integer c)
{ 
    return object_dense(ti, md::one_value<Object>(ti), r, c); 
}
inline object_dense ones(ti::ti_object ti,Integer n)
{ 
    return object_dense(ti, md::one_value<Object>(ti), n, n); 
}

real_dense              eye(Integer m, Integer n);
object_dense            eye(ti::ti_object ti, Integer m, Integer n);
inline real_dense       eye(Integer n)                              { return raw::eye(n, n);}
inline object_dense     eye(ti::ti_object ti,Integer n)             { return raw::eye(ti, n, n); }
integer_dense           ieye(Integer m, Integer n);
inline integer_dense    ieye(Integer n)                             { return ieye(n, n); }
complex_dense           ceye(Integer m, Integer n); 
inline complex_dense    ceye(Integer n)                             { return ceye(n, n); }
float_dense             feye(Integer m, Integer n); 
inline float_dense      feye(Integer n)                             { return feye(n, n); }
float_complex_dense     fceye(Integer m, Integer n); 
inline float_complex_dense
                        fceye(Integer n)                            { return fceye(n, n); }

real_sparse             speye(Integer m, Integer n);
real_sparse             speye(Integer n);
object_sparse           speye(ti::ti_object ti,Integer m, Integer n);
object_sparse           speye(ti::ti_object ti,Integer n);
integer_sparse          ispeye(Integer m, Integer n);
integer_sparse          ispeye(Integer n);
complex_sparse          cspeye(Integer m, Integer n);
complex_sparse          cspeye(Integer n);
float_sparse            fspeye(Integer m, Integer n);
float_sparse            fspeye(Integer n);
float_complex_sparse    fcspeye(Integer m, Integer n);
float_complex_sparse    fcspeye(Integer n);

real_band               beye(Integer m, Integer n, Integer fd, Integer ld);
real_band               beye(Integer n, Integer fd, Integer ld);
object_band             beye(ti::ti_object ti,Integer m, Integer n, Integer fd, Integer ld);
object_band             beye(ti::ti_object ti,Integer n, Integer fd, Integer ld);
integer_band            ibeye(Integer m, Integer n, Integer fd, Integer ld);
integer_band            ibeye(Integer n, Integer fd, Integer ld);
complex_band            cbeye(Integer m, Integer n, Integer fd, Integer ld);
complex_band            cbeye(Integer n, Integer fd, Integer ld);
float_band              fbeye(Integer m, Integer n, Integer fd, Integer ld);
float_band              fbeye(Integer n, Integer fd, Integer ld);
float_complex_band      fcbeye(Integer m, Integer n, Integer fd, Integer ld);
float_complex_band      fcbeye(Integer n, Integer fd, Integer ld);

integer_dense           diag(const integer_dense &v, Integer i = 0);
real_dense              diag(const real_dense &v, Integer i = 0);
float_dense             diag(const float_dense &v, Integer i = 0);
complex_dense           diag(const complex_dense &v, Integer i = 0);
float_complex_dense     diag(const float_complex_dense &v, Integer i = 0);
object_dense            diag(const object_dense &v, Integer i = 0);

integer_band            bdiag(const integer_dense &v, Integer d = 0);
integer_band            bdiag(const integer_sparse &v, Integer d = 0);
real_band               bdiag(const real_dense &v, Integer d = 0);
real_band               bdiag(const real_sparse &v, Integer d = 0);
float_band              bdiag(const float_dense &v, Integer d = 0);
float_band              bdiag(const float_sparse &v, Integer d = 0);
complex_band            bdiag(const complex_dense &v, Integer d = 0);
complex_band            bdiag(const complex_sparse &v, Integer d = 0);
float_complex_band      bdiag(const float_complex_dense &v, Integer d = 0);
float_complex_band      bdiag(const float_complex_sparse &v, Integer d = 0);
object_band             bdiag(const object_dense &v, Integer d = 0);
object_band             bdiag(const object_sparse &v, Integer d = 0);

integer_sparse          spdiag(const integer_dense &v, Integer d = 0);
real_sparse             spdiag(const real_dense &v, Integer d = 0);
float_sparse            spdiag(const float_dense &v, Integer d = 0);
complex_sparse          spdiag(const complex_dense &v, Integer d = 0);
float_complex_sparse    spdiag(const float_complex_dense &v, Integer d = 0);
object_sparse           spdiag(const object_dense &v, Integer d = 0);

integer_dense           diags(const integer_dense &B, const integer_dense &d, Integer m, Integer n);
real_dense              diags(const real_dense &B, const integer_dense &d, Integer m, Integer n);
float_dense             diags(const float_dense &B, const integer_dense &d, Integer m, Integer n);
complex_dense           diags(const complex_dense &B, const integer_dense &d, Integer m, Integer n);
float_complex_dense     diags(const float_complex_dense &B, const integer_dense &d, Integer m, Integer n);
object_dense            diags(const object_dense &B, const integer_dense &d, Integer m, Integer n);

integer_sparse          spdiags(const integer_dense &B, const integer_dense &d, Integer m, Integer n);
real_sparse             spdiags(const real_dense &B, const integer_dense &d, Integer m, Integer n);
float_sparse            spdiags(const float_dense &B, const integer_dense &d, Integer m, Integer n);
complex_sparse          spdiags(const complex_dense &B, const integer_dense &d, Integer m, Integer n);
float_complex_sparse    spdiags(const float_complex_dense &B, const integer_dense &d, Integer m, Integer n);
object_sparse           spdiags(const object_dense &B, const integer_dense &d, Integer m, Integer n);

integer_band            bdiags(const integer_dense &B, const integer_dense &d, Integer m, Integer n);
real_band               bdiags(const real_dense &B, const integer_dense &d, Integer m, Integer n);
float_band              bdiags(const float_dense &B, const integer_dense &d, Integer m, Integer n);
complex_band            bdiags(const complex_dense &B, const integer_dense &d, Integer m, Integer n);
float_complex_band      bdiags(const float_complex_dense &B, const integer_dense &d, Integer m, Integer n);
object_band             bdiags(const object_dense &B, const integer_dense &d, Integer m, Integer n);

real_dense              rand(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline real_dense       rand(const matcl::rand_state& rand_ptr, Integer n) 
                                                                    { return rand(rand_ptr, n, n); }
integer_dense           irand(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline integer_dense    irand(const matcl::rand_state& rand_ptr, Integer n)
                                                                    { return irand(rand_ptr, n, n); }
complex_dense           crand(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline complex_dense    cranc(const matcl::rand_state& rand_ptr, Integer n)
                                                                    { return crand(rand_ptr, n, n); }
float_dense             frand(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline float_dense      franc(const matcl::rand_state& rand_ptr, Integer n)
                                                                    { return frand(rand_ptr, n, n); }
float_complex_dense     fcrand(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline float_complex_dense
                        fcranc(const matcl::rand_state& rand_ptr, Integer n)
                                                                    { return fcrand(rand_ptr, n, n); }

real_dense              randn(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline real_dense       randn(const matcl::rand_state& rand_ptr, Integer n) 
                                                                    { return randn(rand_ptr, n, n); }
complex_dense           crandn(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline complex_dense    crandn(const matcl::rand_state& rand_ptr, Integer n)  
                                                                    { return crandn(rand_ptr, n, n); }
float_dense             frandn(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline float_dense      frandn(const matcl::rand_state& rand_ptr, Integer n)  
                                                                    { return frandn(rand_ptr, n, n); }
float_complex_dense     fcrandn(const matcl::rand_state& rand_ptr, Integer m, Integer n);
inline float_complex_dense
                        fcrandn(const matcl::rand_state& rand_ptr, Integer n)  
                                                                    { return fcrandn(rand_ptr, n, n); }

real_sparse             sprand(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
integer_sparse          isprand(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
complex_sparse          csprand(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
float_sparse            fsprand(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
float_complex_sparse    fcsprand(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
real_sparse             sprandn(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
complex_sparse          csprandn(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
float_sparse            fsprandn(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);
float_complex_sparse    fcsprandn(const matcl::rand_state& rand_ptr, Integer m, Integer n, Real d);

integer_dense           randperm(const matcl::rand_state& rand_ptr, Integer n);

};};