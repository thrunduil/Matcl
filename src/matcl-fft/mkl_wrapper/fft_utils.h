#pragma once

#include "fft_wrapper.h"
#include "matcl-matrep/matcl_matrep.h"

namespace matcl {namespace fft { namespace details { namespace mkl_fft
{

inline fft_int* m2fft(matcl::Integer* ptr)
{
    static_assert(sizeof(fft_int) == sizeof(matcl::Integer), "incompatible types");
    return reinterpret_cast<fft_int*>(ptr);
}
inline fft_real* m2fft(matcl::Real* ptr)
{
    static_assert(sizeof(fft_real) == sizeof(matcl::Real), "incompatible types");
    return reinterpret_cast<fft_real*>(ptr);
}
inline fft_cplx* m2fft(matcl::Complex* ptr)
{
    static_assert(sizeof(fft_cplx) == sizeof(matcl::Complex), "incompatible types");
    return reinterpret_cast<fft_cplx*>(ptr);
}

inline const fft_int* m2fft(const matcl::Integer* ptr)
{
    static_assert(sizeof(fft_int) == sizeof(matcl::Integer), "incompatible types");
    return reinterpret_cast<const fft_int*>(ptr);
}
inline const fft_real* m2fft(const matcl::Real* ptr)
{
    static_assert(sizeof(fft_real) == sizeof(matcl::Real), "incompatible types");
    return reinterpret_cast<const fft_real*>(ptr);
}
inline const fft_cplx* m2fft(const matcl::Complex* ptr)
{
    static_assert(sizeof(fft_cplx) == sizeof(matcl::Complex), "incompatible types");
    return reinterpret_cast<const fft_cplx*>(ptr);
}

}}}}