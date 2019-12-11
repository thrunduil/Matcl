#include "matcl-fft/matcl_fft.h"
#include "matcl-fft/matcl_dct.h"
#include "fft_tests.h"
#include "dct_tests.h"
#include "dct_tests_perf.h"

#include "matcl-matrep\matcl_matrep.h"
#include "matcl-fft/matcl_fft.h"
#include "matcl-fft/matcl_dct.h"

#include <iostream>

namespace mtf = matcl::fft::test;

namespace matcl { namespace fft { namespace test
{
    void test_dct_II();
}}}

using namespace matcl;

int main()
{
    for (int i = 2; i < 10; ++i)
    {
        Matrix A = randn(i, 1);
        Matrix B = fft::dct1(A, 1);
    };

    for (int i = 2; i < 10; ++i)
    {
        Matrix A = randn(i, 1);
        Matrix B = fft::dct2(A, 1);
    };

    for (int i = 2; i < 10; ++i)
    {
        Matrix A = randn(i, 1);
        Matrix B = fft::dct3(A, 1);
    };


    try
    {
        //odyn::fft::test_dct_II();
        //return 0;

        if (0)
        {
            mtf::test_dct_perf ts;
            ts.make();
            return 0 ;
        }
        {
            mtf::test_dct ts;
            ts.make();
        }
        {
            mtf::test_fft ts;
            ts.make();
        };
    }
    catch(std::exception& ex)
    {
        std::cout << ex.what() << "\n";
    };
};