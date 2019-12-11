#pragma once

#include "mkgen/mkgen_fwd.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              is_simd_type
//----------------------------------------------------------------------------------

// is_simd_type<V>::value is true if V is a matcl::simd type

template<class V>
struct is_simd_type
{
    static const bool value = false;
};

template<class V, int Bits, class Tag>
struct is_simd_type<matcl::simd::simd<V, Bits, Tag>>
{
    static const bool value = true;
};

//----------------------------------------------------------------------------------
//                              simd helpers
//----------------------------------------------------------------------------------

// simd_vector_size<V>::value returns number of scalars stored in simd type V
// or return 1 if V is a scalar type

template<class V, bool Is_simd = is_simd_type<V>::value>
struct simd_vector_size;

template<class V>
struct simd_vector_size<V, true>
{
    static const int value  = V::vector_size;
};

template<class V>
struct simd_vector_size<V, false>
{
    static const int value  = 1;
};

//----------------------------------------------------------------------------------
//                              store reverse
//----------------------------------------------------------------------------------

template<class Simd, class Aligned, 
            bool Is_simd = is_simd_type<Simd>::value>
struct store_reverse
{    
    static_assert(Is_simd == true, "simd type required");

    using Scalar    = typename Simd::value_type;

    static void eval(Scalar* arr, const Simd& v)
    {
        static const Integer off        = Simd::vector_size - 1;

        Simd vr = matcl::simd::reverse(v);
        v.store(arr - off, Aligned());
    };
};

//----------------------------------------------------------------------------------
//                              load reverse
//----------------------------------------------------------------------------------

template<class Simd, class Aligned, 
            bool Is_simd = is_simd_type<Simd>::value>
struct load_reverse
{    
    static_assert(Is_simd == true, "simd type required");

    using Scalar    = typename Simd::value_type;

    static Simd eval(const Scalar* arr)
    {
        static const Integer off    = Simd::vector_size - 1;

        Simd v;
        v.load(arr - off, Aligned());
        
        Simd vr     = matcl::simd::reverse(v);
        return vr;
    };
};

}}}