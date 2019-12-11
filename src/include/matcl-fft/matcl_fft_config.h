#pragma once

#include "matcl_config.h"

#ifdef MATCL_FFT_EXPORTS
    #define MATCL_FFT_EXPORT _declspec(dllexport)
#else
    #define MATCL_FFT_EXPORT _declspec(dllimport)
#endif

#define restricted __restrict
