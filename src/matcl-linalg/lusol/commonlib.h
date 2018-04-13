#pragma once

#include <string.h>
#include <complex>

namespace lusol
{

using REAL      = double;
using FLOAT     = float;
using COMPLEX   = std::complex<double>;
using FCOMPLEX  = std::complex<float>;
using Integer   = int;
using INT       = int;

static const int LINEARSEARCH   = 5;

template<class T>
inline void lusol_free(T*& ptr)
{
    if(ptr != nullptr) 
    {
        free(ptr);
        ptr = nullptr;
    };
};

template<class T>
inline void memcopy(T* nptr, const T* optr, size_t nr)
{
    memcpy(nptr, optr, nr * sizeof(T));
}

template<class T>
inline void mem_move(T* nptr, const T* optr, size_t nr)
{
    memmove(nptr, optr, nr * sizeof(T));
}
template<class T>
inline void memclear(T* ptr, size_t nr)
{
    memset(ptr, '\0', nr * sizeof(T));
};

template<class T> inline T minimum(T x, T y) { return x<y? x : y; };
template<class T> inline T maximum(T x, T y) { return x>y? x : y; };

typedef int findCompare_func(const void *current, const void *candidate);

};