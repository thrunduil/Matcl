#include "dct_kernels_test.h"
#include "static_dct.h"
#include "matcl-fft/dct_kernels.h"
#include "matcl-linalg\matcl_linalg.h"

#include <iostream>

#define DO_EVAL     1
#define EVAL_DCT1   0
#define EVAL_DCT2   0
#define EVAL_DCT3   1
#define EVAL_DCT4   0

#define ENABLE_IACA 0

#if ENABLE_IACA == 1
    #include <iacaMarks.h>
#endif

namespace matcl { namespace fft
{

namespace mk = matcl::mkgen;

template<Integer M>
void dct_kernels_test<M>::print()
{
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT1" << "\n";
    std::cout << std::string(80, '-') << "\n";

    dct1_evaler_type<M>::print();

    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT2" << "\n";
    std::cout << std::string(80, '-') << "\n";

    dct2_evaler_type<M>::print();

    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT3" << "\n";
    std::cout << std::string(80, '-') << "\n";

    dct3_evaler_type<M>::print();

    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT4" << "\n";
    std::cout << std::string(80, '-') << "\n";

    dct4_evaler_type<M>::print();
};

template<Integer M>
void dct_kernels_test<M>::check()
{
    Integer N           = 1;
    Matrix x            = matcl::randn(M,N);
    const Real* X       = x.get_array<Real>();

    Real* Y;
    Matrix y            = matcl::make_real_dense_noinit(M,N, Y);

    Integer X_ld        = M;
    Integer Y_ld        = M;

    {
        matcl::fft::dct1_kernel<M>(X,Y, N, X_ld, Y_ld);
        Matrix Z        = matcl::fft::dct1(x,1);

        Real dif        = matcl::norm(y - Z, 1);
        Matrix d        = div(y,Z);

        std::cout << "\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::string(40, ' ') << "DCT1" << "\n";
        std::cout << std::string(80, '-') << "\n";

        disp(d);

        std::cout << "\n";
        std::cout << dif << "\n";
    };

    {
        matcl::fft::dct2_kernel<M>(X,Y, N, X_ld, Y_ld);
        Matrix Z        = matcl::fft::dct2(x,1);

        Real dif        = matcl::norm(y - Z, 1);
        Matrix d        = div(y,Z);

        std::cout << "\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::string(40, ' ') << "DCT2" << "\n";
        std::cout << std::string(80, '-') << "\n";

        disp(d);

        std::cout << "\n";
        std::cout << dif << "\n";
    };

    {
        matcl::fft::dct3_kernel<M>(X,Y, N, X_ld, Y_ld);
        Matrix Z        = matcl::fft::dct3(x,1);

        Real dif        = matcl::norm(y - Z, 1);
        Matrix d        = div(y,Z);

        std::cout << "\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::string(40, ' ') << "DCT3" << "\n";
        std::cout << std::string(80, '-') << "\n";

        disp(d);

        std::cout << "\n";
        std::cout << dif << "\n";
    };
};

template<class DCT_Config>
struct make_config_default
{
    using config_dct        = DCT_Config;
    using config_codegen    = basic_codegen;

    static void print(std::ostream& os)
    {
        Integer Max_Genmat_1    = DCT_Config::max_genmat_dct1;
        Integer Max_Genmat_2    = DCT_Config::max_genmat_dct2;
        Integer Max_Genmat_3    = DCT_Config::max_genmat_dct3;
        Integer Max_Genmat_4    = DCT_Config::max_genmat_dct4;

        Integer Min_Rec_1       = DCT_Config::min_rec_length_dct1;
        Integer Min_Rec_2       = DCT_Config::min_rec_length_dct2;
        Integer Min_Rec_3       = DCT_Config::min_rec_length_dct3;
        Integer Min_Rec_4       = DCT_Config::min_rec_length_dct4;

        os  << "\t\t" 
            << config_codegen::simd_enable << " " << config_codegen::simd_half_allow << " " 
            << "gm: " << Max_Genmat_1 << " " << Max_Genmat_2  << " "  << Max_Genmat_3 << " "  << Max_Genmat_4 << " " 
            << "rec: " << Min_Rec_1    << " " << Min_Rec_2     << " "  << Min_Rec_3    << " "  << Min_Rec_4;
    };
};

template<Integer Seed, bool Simd_Enable, bool Allow_Half, 
    Integer Max_Genmat_1, Integer Max_Genmat_2, Integer Max_Genmat_3, Integer Max_Genmat_4, 
    Integer Min_Rec_1, Integer Min_Rec_2, Integer Min_Rec_3, Integer Min_Rec_4>
struct make_config
{
    struct basic_codegen
    {
        static const bool simd_enable                   = Simd_Enable;
        static const bool simd_allow_negative_step      = true;
        static const bool simd_allow_unaligned          = true;

        static const bool simd_half_allow               = Allow_Half;
        static const bool simd_half_allow_negative_step = true;
        static const bool simd_half_allow_unaligned     = true;

        using simd_type                                 = simd::maximum_tag;
        static const int simd_bits                      = simd::maximum_bits;
    };

    struct dct_config
    {
        static const Integer seed                       = Seed;

        static const Integer min_rec_length_dct1        = Min_Rec_1;
        static const Integer min_rec_length_dct2        = Min_Rec_2;
        static const Integer min_rec_length_dct3        = Min_Rec_3;
        static const Integer min_rec_length_dct4        = Min_Rec_4;

        static const Integer max_genmat_dct1            = Max_Genmat_1;    
        static const Integer max_genmat_dct2            = Max_Genmat_2;    
        static const Integer max_genmat_dct3            = Max_Genmat_3;    
        static const Integer max_genmat_dct4            = Max_Genmat_4;    
    };

    using config_dct        = dct_config;
    using config_codegen    = basic_codegen;

    static void print(std::ostream& os)
    {
        os  << "\t\t"
            << Simd_Enable  << " " << Allow_Half << " " 
            << "gm: "  << Max_Genmat_1 << " " << Max_Genmat_2  << " "  << Max_Genmat_3 << " "  << Max_Genmat_4 << " " 
            << "rec: " << Min_Rec_1    << " " << Min_Rec_2     << " "  << Min_Rec_3    << " "  << Min_Rec_4;
    };
};

template<class Config>
struct print_config
{
    friend std::ostream& operator<<(std::ostream& os, const print_config<Config>& c)
    {
        (void)c;
        Config::print(os);
        return os;
    };
};

template<Integer M>
void dct_kernels_test<M>::eval_versions()
{
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT1" << "\n";
    std::cout << std::string(80, '-') << "\n";

  #if EVAL_DCT1
    {
        using configs = mk::list<   
                                    make_config<-2,true,true,0,0,0,0,2,2,2,2>,
                                    make_config<-3,true,true,0,0,0,0,4,2,2,2>,
                                    make_config<-4,true,true,0,0,0,0,2,4,2,2>,
                                    make_config<-5,true,true,0,0,0,0,2,2,2,4>,
                                    make_config<-6,true,true,0,0,0,0,2,4,2,4>,
                                    make_config<-7,true,true,0,0,0,0,4,2,2,4>,
                                    make_config<-8,true,true,0,0,0,0,4,4,2,2>,
                                    make_config<-9,true,true,0,0,0,0,4,4,2,4>,

                                    make_config<-10,true,true,3,2,2,2,2,2,2,2>,
                                    make_config<-11,true,true,3,2,2,2,4,2,2,2>,
                                    make_config<-12,true,true,3,2,2,2,2,4,2,2>,
                                    make_config<-13,true,true,3,2,2,2,2,2,2,4>,
                                    make_config<-14,true,true,3,2,2,2,2,4,2,4>,
                                    make_config<-15,true,true,3,2,2,2,4,2,2,4>,
                                    make_config<-16,true,true,3,2,2,2,4,4,2,2>,
                                    make_config<-17,true,true,3,2,2,2,4,4,2,4>,

                                    make_config<-18,true,true,5,4,4,4,2,2,2,2>,
                                    make_config<-19,true,true,5,4,4,4,4,2,2,2>,
                                    make_config<-20,true,true,5,4,4,4,2,4,2,2>,
                                    make_config<-21,true,true,5,4,4,4,2,2,2,4>,
                                    make_config<-22,true,true,5,4,4,4,2,4,2,4>,
                                    make_config<-23,true,true,5,4,4,4,4,2,2,4>,
                                    make_config<-24,true,true,5,4,4,4,4,4,2,2>,
                                    make_config<-25,true,true,5,4,4,4,4,4,2,4>,
                                    make_config<-26,true,true,0,0,0,0,0,2,2,2>
                                >;
        eval_versions_1<dct1_evaler_conf,configs,make_config_default<dct_config_1>>();
    };
  #endif

  #if EVAL_DCT2
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT2" << "\n";
    std::cout << std::string(80, '-') << "\n";

    {
        using configs = mk::list<   
                                    make_config<-2,true,true,0,0,0,0,2,2,2,2>,
                                    make_config<-3,true,true,0,0,0,0,2,4,2,2>,
                                    make_config<-4,true,true,0,0,0,0,2,2,2,4>,
                                    make_config<-5,true,true,0,0,0,0,2,4,2,4>,

                                    make_config<-6,true,true,2,2,2,2,2,2,2,2>,
                                    make_config<-7,true,true,2,2,2,2,2,4,2,2>,
                                    make_config<-8,true,true,2,2,2,2,2,2,2,4>,
                                    make_config<-9,true,true,2,2,2,2,2,4,2,4>,

                                    make_config<-10,true,true,4,4,4,4,2,2,2,2>,
                                    make_config<-11,true,true,4,4,4,4,2,4,2,2>,
                                    make_config<-12,true,true,4,4,4,4,2,2,2,4>,
                                    make_config<-13,true,true,4,4,4,4,2,4,2,4>
                                >;
        eval_versions_1<dct2_evaler_conf,configs,make_config_default<dct_config>>();
    };
  #endif

  #if EVAL_DCT3
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT3" << "\n";
    std::cout << std::string(80, '-') << "\n";

    {
        using configs = mk::list<   
                                    make_config<-2,true,true,0,0,0,0,2,2,2,2>,
                                    make_config<-3,true,true,0,0,0,0,2,4,2,2>,
                                    make_config<-4,true,true,0,0,0,0,2,2,4,2>,
                                    make_config<-5,true,true,0,0,0,0,2,2,2,4>,
                                    make_config<-6,true,true,0,0,0,0,2,4,4,2>,
                                    make_config<-7,true,true,0,0,0,0,2,4,2,4>,
                                    make_config<-8,true,true,0,0,0,0,2,2,4,4>,
                                    make_config<-9,true,true,0,0,0,0,2,4,4,4>,

                                    make_config<-10,true,true,2,2,2,2,2,2,2,2>,
                                    make_config<-11,true,true,2,2,2,2,2,4,2,2>,
                                    make_config<-12,true,true,2,2,2,2,2,2,4,2>,
                                    make_config<-13,true,true,2,2,2,2,2,2,2,4>,
                                    make_config<-14,true,true,2,2,2,2,2,4,4,2>,
                                    make_config<-15,true,true,2,2,2,2,2,4,2,4>,
                                    make_config<-16,true,true,2,2,2,2,2,2,4,4>,
                                    make_config<-17,true,true,2,2,2,2,2,4,4,4>,

                                    make_config<-18,true,true,4,4,4,4,2,2,2,2>,
                                    make_config<-19,true,true,4,4,4,4,2,4,2,2>,
                                    make_config<-20,true,true,4,4,4,4,2,2,4,2>,
                                    make_config<-21,true,true,4,4,4,4,2,2,2,4>,
                                    make_config<-22,true,true,4,4,4,4,2,4,4,2>,
                                    make_config<-23,true,true,4,4,4,4,2,4,2,4>,
                                    make_config<-24,true,true,4,4,4,4,2,2,4,4>,
                                    make_config<-25,true,true,4,4,4,4,2,4,4,4>
                                >;
        eval_versions_1<dct3_evaler_conf,configs,make_config_default<dct_config>>();
    };
  #endif

  #if EVAL_DCT4
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::string(40, ' ') << "DCT4" << "\n";
    std::cout << std::string(80, '-') << "\n";

    {
        using configs = mk::list<   
                                    make_config<-2,true,true,0,0,0,0,2,2,2,2>,
                                    make_config<-3,true,true,0,0,0,0,2,4,2,2>,
                                    make_config<-4,true,true,0,0,0,0,2,2,2,4>,
                                    make_config<-5,true,true,0,0,0,0,2,4,2,4>,

                                    make_config<-6,true,true,2,2,2,2,2,2,2,2>,
                                    make_config<-7,true,true,2,2,2,2,2,4,2,2>,
                                    make_config<-8,true,true,2,2,2,2,2,2,2,4>,
                                    make_config<-9,true,true,2,2,2,2,2,4,2,4>,

                                    make_config<-10,true,true,4,4,4,4,2,2,2,2>,
                                    make_config<-11,true,true,4,4,4,4,2,4,2,2>,
                                    make_config<-12,true,true,4,4,4,4,2,2,2,4>,
                                    make_config<-13,true,true,4,4,4,4,2,4,2,4>
                                >;
        eval_versions_1<dct4_evaler_conf, configs,make_config_default<dct_config>>();
    };
  #endif
};        

template<Integer M, template<Integer M, class Config> class Eval_Config, class Configs, Integer Size, Integer Pos>
struct eval_versions_impl
{
    static void eval(Integer T, const double* X, double* Y, Integer N, Integer X_ld, Integer Y_ld, Integer& ver, Real t0)
    {
        using config    = typename mk::get_elem_at_pos<Configs,Pos>::type;
        using evaler    = Eval_Config<M,config>;
        Real t          = dct_kernels_test<M>::eval_ver<evaler>(T, X, Y, N, X_ld, Y_ld);
        std::cout       << "ver " << ver++ << ": " << t << " " << t / t0 << " " 
                        << print_config<config>() << "\n";

        eval_versions_impl<M,Eval_Config,Configs,Size,Pos+1>
            ::eval(T, X, Y, N, X_ld, Y_ld, ver, t0);
    };
};
template<Integer M, template<Integer M, class Config> class Eval_Config, class Configs, Integer Size>
struct eval_versions_impl<M,Eval_Config,Configs,Size,Size>
{
    static void eval(Integer T, const double* X, double* Y, Integer N, Integer X_ld, Integer Y_ld, Integer& ver, Real t0)
    {
        (void)T;
        (void)X;
        (void)Y;
        (void)N;
        (void)X_ld;
        (void)Y_ld;
        (void)ver;
        (void)t0;
    };
};

template<Integer M>
template<template<Integer M, class Config> class Eval_Config, class Configs, class Base_Config>
void dct_kernels_test<M>::eval_versions_1()
{
  #if DO_EVAL == 1
    Integer N           = 50;
    Integer T           = 100000;
    //Integer N           = 10;
    //Integer T           = 1;

    Matrix x            = matcl::randn(M,N);;
    const Real* X       = x.get_array<Real>();

    Real* restricted Y;
    Matrix y            = matcl::make_real_dense_noinit(M,N, Y);

    Integer X_ld        = M;
    Integer Y_ld        = M;

    Real t0;
    {
        using evaler    = Eval_Config<M,Base_Config>;
        t0              = eval_ver<evaler>(T, X, Y, N, X_ld, Y_ld);
    };

    std::cout           << "ver 0: " << t0 << " " << t0 / t0 << " " 
                        << print_config<Base_Config>() << "\n";
    int ver             = 1;
    
    eval_versions_impl<M,Eval_Config,Configs,mk::list_size<Configs>::value,0>
        ::eval(T, X, Y, N, X_ld, Y_ld, ver, t0);
  #endif
};

template<Integer M>
template<class Evaler>
double dct_kernels_test<M>::eval_ver(Integer T, const double* X, double* Y, Integer N, Integer in_ld, Integer out_ld)
{
    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const double* in    = X;
        double* out         = Y;

        #if ENABLE_IACA == 1
            IACA_START
        #endif

        for(Integer j = 0; j < N; ++j)
        {
            Evaler::eval(out,in);

            in  += in_ld;
            out += out_ld;
        };

        #if ENABLE_IACA == 1
            IACA_END
        #endif        
    };
    Real t1 = matcl::toc();
    return t1;
};

template<Integer M>
void dct_kernels_test<M>::test_perf()
{
    Integer N           = 50;
    Integer T           = 100000;

    Integer M2          = ((M - 1)/4 + 1) * 4;

    Matrix x            = matcl::randn(M2,N);
    x                   = x(matcl::colon(1,M),matcl::colon());
    const Real* X       = x.get_array<Real>();

    Real* restricted Y;
    Matrix y            = matcl::make_real_dense_noinit(M2,N, Y);
    y                   = y(matcl::colon(1,M),matcl::colon());

    Integer X_ld        = M2;
    Integer Y_ld        = M2;

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Y;
        matcl::fft::dct1_kernel<M>(X2,Y2, N, X_ld, Y_ld);
    }
    Real t1 = matcl::toc();

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Y;
        dct2_kernel<M>(X2,Y2, N, X_ld, Y_ld);
    }
    Real t2 = matcl::toc();

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Y;
        dct3_kernel<M>(X2,Y2, N, X_ld, Y_ld);
    }
    Real t3 = matcl::toc();

    //TODO
    /*
    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Y;
        dct4_kernel<M>(X2,Y2, N, X_ld, Y_ld);
    }
    Real t4 = matcl::toc();
    */

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        Matrix Y2        = matcl::fft::dct1(x,1);
        (void)Y2;
    }
    Real t11 = matcl::toc();

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        Matrix Y2        = matcl::fft::dct2(x,1);
        (void)Y2;
    }
    Real t12 = matcl::toc();

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        Matrix Y2        = matcl::fft::dct3(x,1);
        (void)Y2;
    }
    Real t13 = matcl::toc();

    std::cout << M << " " << N << " " << T << "\n"
              << "dct1: " << t1 << " " << t11 << " " << t11 / t1 << "\n"
              << "dct2: " << t2 << " " << t12 << " " << t12 / t2 << "\n"
              << "dct3: " << t3 << " " << t13 << " " << t13 / t3 << "\n"
              //<< "dct4: " << t4 << " " << 1.0 << " " << 1.0 / t4 << "\n" TODO
            ;
};

template class dct_kernels_test<DCT_SIZE>;

}};
