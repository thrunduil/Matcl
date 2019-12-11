#pragma once

#include "matcl-fft/matcl_fft.h"

namespace matcl { namespace fft { namespace test
{

class test_dct_perf
{
    private:
        void        test_dct1(Integer M, Integer N, Integer T);
        void        test_dct2(Integer M, Integer N, Integer T);
        void        test_dct3(Integer M, Integer N, Integer T);

        void        test_dct1_dim2(Integer M, Integer N, Integer T);
        void        test_dct2_dim2(Integer M, Integer N, Integer T);
        void        test_dct3_dim2(Integer M, Integer N, Integer T);

        void        test_dct1_inpl(Integer M, Integer N, Integer T);
        void        test_dct23_inpl(Integer M, Integer N, Integer T);

    public:
        void        make();
};

}}}
