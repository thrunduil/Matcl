#pragma once 

#include "mkgen/mkgen.h"
#include "matcl-fft/matcl_fft_config.h"
#include "matcl-fft/matcl_dct.h"
#include "static_dct.h"

namespace matcl { namespace fft
{

namespace mk = mkgen;

#define DCT_SIZE        8

struct basic_codegen
{
    static const bool simd_enable                   = true;
    static const bool simd_allow_negative_step      = true;
    static const bool simd_allow_unaligned          = true;

    static const bool simd_half_allow               = true;
    static const bool simd_half_allow_negative_step = true;
    static const bool simd_half_allow_unaligned     = true;

    using simd_type                                 = simd::maximum_tag;
    static const int simd_bits                      = simd::maximum_bits;
};

struct dct_config
{
    static const Integer seed                       = -1;

    static const Integer min_rec_length_dct1        = 2;
    static const Integer min_rec_length_dct2        = 4;
    static const Integer min_rec_length_dct3        = 4;
    static const Integer min_rec_length_dct4        = 4;

    static const Integer max_genmat_dct1            = 0;
    static const Integer max_genmat_dct2            = 4;    
    static const Integer max_genmat_dct3            = 4;
    static const Integer max_genmat_dct4            = 4;
};

struct dct_config_1
{
    static const Integer seed                       = -1;

    static const Integer min_rec_length_dct1        = 2;
    static const Integer min_rec_length_dct2        = 2;
    static const Integer min_rec_length_dct3        = 2;
    static const Integer min_rec_length_dct4        = 2;

    static const Integer max_genmat_dct1            = 0;
    static const Integer max_genmat_dct2            = 0;    
    static const Integer max_genmat_dct3            = 0;
    static const Integer max_genmat_dct4            = 0;
};

template<Integer M>
using dct1_evaler_type   = dct1_impl<dct_config_1, M,1,basic_codegen,mk::two>;

template<Integer M>
using dct2_evaler_type   = dct2_impl<dct_config, M,1,basic_codegen,mk::two>;

template<Integer M>
using dct3_evaler_type   = dct3_impl<dct_config, M,1,basic_codegen,mk::one>;

template<Integer M>
using dct4_evaler_type   = dct4_impl<dct_config, M,1,basic_codegen,mk::one>;

template<Integer M, class Config>
using dct1_evaler_conf   = dct1_impl<typename Config::config_dct, M,1,typename Config::config_codegen,mk::two>;

template<Integer M, class Config>
using dct2_evaler_conf   = dct2_impl<typename Config::config_dct, M,1,typename Config::config_codegen,mk::two>;

template<Integer M, class Config>
using dct3_evaler_conf   = dct3_impl<typename Config::config_dct, M,1,typename Config::config_codegen,mk::one>;

template<Integer M, class Config>
using dct4_evaler_conf   = dct4_impl<typename Config::config_dct, M,1,typename Config::config_codegen,mk::one>;

template<Integer M>
class dct_kernels_test
{
    public:
        static void     print();
        static void     check();
        static void     eval_versions();
        static void     test_perf();

    public:
        template<class Evaler>
        static double   eval_ver(Integer T, const double* in, double* out, Integer N, Integer in_ld, Integer out_ld);

        template<template<Integer M, class Config> class Eval_Config, class Configs, class Base_Config>
        static void     eval_versions_1();
};

}};
