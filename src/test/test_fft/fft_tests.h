#pragma once

#include "matcl-fft/matcl_fft.h"

namespace matcl { namespace fft { namespace test
{

class test_fft
{
    private:
        bool        test_fft_real();
        bool        test_fft_pack();
        bool        test_fft_complex();
        bool        test_fft_type(bool complex);
        bool        test_fft_pack_type(bool complex);
        bool        test_fft_1(Integer r, Integer c, bool complex);
        bool        test_fft_pack_1(Integer r, Integer c, bool complex);

        bool        test_ifft_real();
        bool        test_ifft_complex();
        bool        test_ifft_real_conj_even();
        bool        test_ifft_complex_conj_even();

        bool        test_ifft_type(bool complex, bool conj_even);
        bool        test_ifft_1(Integer r, Integer c, bool complex, bool conj_even);

    public:
        void        make();
};

}}}
