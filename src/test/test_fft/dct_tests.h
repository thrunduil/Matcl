#pragma once

#include "matcl-fft/matcl_fft.h"

namespace matcl { namespace fft { namespace test
{

class test_dct
{
    private:
        bool        test_dct1();
        bool        test_dct1_inpl();
        bool        test_dct1_1(Integer r, Integer c, Real scal);
        bool        test_dct1_inpl_1(Integer r, Integer c, Real scal);

        bool        test_dct2();
        bool        test_dct2_inpl();
        bool        test_dct2_1(Integer r, Integer c, Real scal);
        bool        test_dct2_inpl_1(Integer r, Integer c, Real scal);

        bool        test_dct2_MH();
        bool        test_dct2_MH_inpl();
        bool        test_dct2_MH_1(Integer r, Integer c, Real scal);
        bool        test_dct2_MH_inpl_1(Integer r, Integer c, Real scal);

        bool        test_dct3();
        bool        test_dct3_1(Integer r, Integer c, Real scal);
        bool        test_dct3_inpl();
        bool        test_dct3_inpl_1(Integer r, Integer c, Real scal);

    public:
        void        make();
};

}}}
