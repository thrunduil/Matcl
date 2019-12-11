#include "dct_tests_perf.h"
#include "matcl-matrep\matcl_matrep.h"
#include "matcl-fft/matcl_fft.h"
#include "matcl-fft/matcl_dct.h"
#include "matcl-linalg\matcl_linalg.h"

#include <iostream>


namespace matcl { namespace fft { namespace test
{

void test_dct_perf::test_dct1(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(M+1,N);

    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct1(x,1,1.0,cont);
    };
    Real t1         = toc();

    std::cout << "dct1: "<< t1 << "\n";
};
void test_dct_perf::test_dct1_inpl(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(M+1,N);
    Real scale  = M / 2.0;

    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        x           = dct1(std::move(x),1,1.0,cont);
        x           = dct1(std::move(x),1,scale,cont);
    };
    Real t1         = toc();

    std::cout << "dct1 inpl: "<< t1 << "\n";
};
void test_dct_perf::test_dct1_dim2(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(N,M+1);

    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct1(x,2,1.0,cont);
    };
    Real t1         = toc();

    std::cout   << "dct1: " << t1 
                << "\n";
};

void test_dct_perf::test_dct2(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(M,N);
    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct2(x,1, 1.0, cont);
    };
    Real t1     = toc();

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct2_MH(x,1, cont);
    };
    Real t2     = toc();

    std::cout   << "dct2: " << t1 << ", " << t2 << ", rel: " << t2/t1
                << "\n";
};
void test_dct_perf::test_dct23_inpl(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(M,N);
    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        x           = dct2(std::move(x),1, 1.0, cont);
        x           = dct3(std::move(x),1, 1.0, cont);
    };
    Real t1     = toc();

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        x           = dct2_MH(std::move(x),1.0, cont);
        x           = dct3(std::move(x),1, 1.0, cont);
    };
    Real t2     = toc();

    std::cout   << "dct2 inpl: " 
                << " ver 1: " << t1 
                << ", ver 2: " << t2 << " rel: " << t2/ t1 
                << "\n";
};
void test_dct_perf::test_dct2_dim2(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(N,M);

    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct2(x,2,1.0,cont);
    };
    Real t1         = toc();

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix x2   = trans(x);
        Matrix y    = dct2(x2,1,1.0,cont);
        y           = trans(y);
    };
    Real t2         = toc();

    std::cout   << "dct2: "
                << " ver 1: "   << t1 
                << ", ver 2: "  << t2 << " rel: " << t2 / t1
                << "\n";
};

void test_dct_perf::test_dct3(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(M,N);
    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct3(x,1, 1.0, cont);
    };
    Real t1         = toc();

    std::cout   << "dct3: " << t1 << "\n";
};
void test_dct_perf::test_dct3_dim2(Integer M, Integer N, Integer T)
{
    Matrix x    = randn(N,M);

    fft_context cont;

    tic();
    for (Integer i = 0; i < T; ++i)
    {
        Matrix y    = dct3(x,2,1.0,cont);
    };
    Real t1         = toc();

    std::cout   << "dct3: " << t1 
                << "\n";
};
void test_dct_perf::make()
{
    Integer D   = 4;
    Integer M   = 8*D;
    Integer N   = 16*8*32*4/D;
    Integer T   = 1000;

    Integer dim     = 2;
    bool inpl       = false;

    if (inpl)
    {
        test_dct1_inpl(M,N,T);
        test_dct23_inpl(M,N,T);
    }
    else
    {
        if (dim == 2)
        {        
            test_dct2_dim2(M,N,T);
            test_dct3_dim2(M,N,T);
            test_dct1_dim2(M,N,T);
        }
        else
        {
            test_dct1(M,N,T);
            test_dct2(M,N,T);
            test_dct3(M,N,T);
        };
    };
};

}}};
