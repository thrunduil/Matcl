#include "fft_tests.h"
#include "matcl-matrep\matcl_matrep.h"
#include "matcl-fft/matcl_fft.h"
#include "matcl-linalg\matcl_linalg.h"

#include <iostream>

namespace matcl { namespace fft { namespace test
{

bool test_fft::test_ifft_1(Integer r, Integer c, bool complex, bool conj_even)
{
    bool ok     = true;
    Real tol    = 1e-13;

    Matrix x;

    if (complex == false)
    {
        x       = randn(r,c);
        if (conj_even == true)
        {
            x   = (mat_col(), x, flipud(x(colon(2,end),colon())));
        }
    }        
    else
    {
        if (conj_even == false)
            x   = crandn(r,c);
        else
        {
            x   = randn(r,c);
            x   = fft(x,1);
        }
    }

    Matrix y    = ifft(x,1,1.0,conj_even);
    Matrix z    = ifft(trans(x),2,1.0,conj_even);
    z           = trans(z);

    Real dif1   = norm(y-z,1);

    Matrix x2   = fft(y,1, 1.0/ y.rows());
    Real dif2   = norm(x2-x,1);

    Real dif3   = 0;
    if (complex == false)
    {
        Matrix y2   = ifft(convert(x,mat_code::complex_dense),1);
        dif3        = norm(y-y2,1);
    };
    if (conj_even == true)
    {
        matcl::value_code vt   = y.get_value_code();
        if (vt != value_code::v_real)
            dif3    += 1.0;
    };

    Real dif    = dif1 + dif2 + dif3;
    if (dif > tol)
    {
        //disp(dif);
        //disp(x);
        //disp(y);
        //disp(x2);
        //disp(div(x2,x));

        if (complex == false)
        {
            Matrix y1   = ifft(x,1,1.0,conj_even);
            Matrix y2   = ifft(convert(x,mat_code::complex_dense),1);

            disp(x);
            disp(y1);
            disp(y2);
        };

        ok      = false;
    };

    return ok;
};
bool test_fft::test_fft_1(Integer r, Integer c, bool complex)
{
    bool ok     = true;
    Real tol    = 1e-14;

    Matrix x;

    if (complex == false)
        x       = randn(r,c);
    else
        x       = crandn(r,c);

    Matrix y    = fft(x,1);
    Matrix z    = fft(trans(x),2);
    z           = trans(z);

    Real dif    = norm(y-z,1);
    Real dif2   = 0;

    if (complex == false)
    {
        Matrix y2   = fft(convert(x,mat_code::complex_dense),1);
        dif2        = norm(y-y2,1);
    };

    if (dif+dif2 > tol)
    {
        disp(dif);
        disp(x);
        disp(y);
        disp(z);
        ok      = false;
    };

    return ok;
};
bool test_fft::test_fft_pack_1(Integer r, Integer c, bool complex)
{
    bool ok     = true;
    Real tol    = 1e-14;

    Matrix x    = randn(r,c);

    if (complex == true)
    {
        x       = convert(x, mat_code::complex_dense);
    };

    Real dif    = 0;
    {
        Matrix y    = fft_pack(x,1);
        Matrix z    = fft_pack(trans(x),2);
        z           = trans(z);

        Real d      = norm(y-z,1);

        d           = dif + d;

        if (d > tol)
        {
            disp(d);
        }
    };
    {
        Matrix y    = ifft_pack(x,1);
        Matrix z    = ifft_pack(trans(x),2);
        z           = trans(z);

        Real d      = norm(y-z,1);

        d           = dif + d;

        if (d > tol)
        {
            disp(d);
        }
    };
    {
        Matrix y    = fft_pack(x,1);
        Matrix z    = ifft_pack(y,1, 1.0 / x.rows());

        Real d      = norm(x-z,1);

        d           = dif + d;

        if (d > tol)
        {
            disp(d);
        }
    };
    {
        Matrix y    = fft_pack(x,2);
        Matrix z    = ifft_pack(y,2, 1.0 / x.cols());

        Real d      = norm(x-z,1);

        d           = dif + d;

        if (d > tol)
        {
            disp(d);
        }
    };

    if (dif > tol)
        ok          = false;

    return ok;
};

bool test_fft::test_fft_type(bool complex)
{
    bool ok     = true;
    ok          &= test_fft_1(0,0,complex);
    ok          &= test_fft_1(0,1,complex);
    ok          &= test_fft_1(1,0,complex);
    
    ok          &= test_fft_1(1,1,complex);
    
    ok          &= test_fft_1(1,5,complex);
    ok          &= test_fft_1(5,1,complex);
    ok          &= test_fft_1(1,6,complex);
    ok          &= test_fft_1(6,1,complex);

    ok          &= test_fft_1(1,5,complex);
    ok          &= test_fft_1(5,1,complex);
    ok          &= test_fft_1(1,6,complex);
    ok          &= test_fft_1(6,1,complex);

    ok          &= test_fft_1(2,5,complex);
    ok          &= test_fft_1(5,2,complex);
    ok          &= test_fft_1(2,6,complex);
    ok          &= test_fft_1(6,2,complex);

    ok          &= test_fft_1(5,5,complex);
    ok          &= test_fft_1(6,6,complex);

    return ok;
};
bool test_fft::test_fft_pack_type(bool complex)
{
    bool ok     = true;
    ok          &= test_fft_pack_1(0,0,complex);
    ok          &= test_fft_pack_1(0,1,complex);
    ok          &= test_fft_pack_1(1,0,complex);
    
    ok          &= test_fft_pack_1(1,1,complex);
    
    ok          &= test_fft_pack_1(1,5,complex);
    ok          &= test_fft_pack_1(5,1,complex);
    ok          &= test_fft_pack_1(1,6,complex);
    ok          &= test_fft_pack_1(6,1,complex);

    ok          &= test_fft_pack_1(1,5,complex);
    ok          &= test_fft_pack_1(5,1,complex);
    ok          &= test_fft_pack_1(1,6,complex);
    ok          &= test_fft_pack_1(6,1,complex);

    ok          &= test_fft_pack_1(2,5,complex);
    ok          &= test_fft_pack_1(5,2,complex);
    ok          &= test_fft_pack_1(2,6,complex);
    ok          &= test_fft_pack_1(6,2,complex);

    ok          &= test_fft_pack_1(5,5,complex);
    ok          &= test_fft_pack_1(6,6,complex);

    return ok;
};
bool test_fft::test_fft_pack()
{
    bool ok = true;
    ok          &= test_fft_pack_type(false);
    ok          &= test_fft_pack_type(true);
    return ok;
};
bool test_fft::test_fft_real()
{
    return test_fft_type(false);
};
bool test_fft::test_fft_complex()
{
    return test_fft_type(true);
};

bool test_fft::test_ifft_type(bool complex, bool conj_even)
{
    bool ok     = true;
    ok          &= test_ifft_1(0,0,complex,conj_even);
    ok          &= test_ifft_1(0,1,complex,conj_even);
    ok          &= test_ifft_1(1,0,complex,conj_even);
    
    ok          &= test_ifft_1(1,1,complex,conj_even);
    
    ok          &= test_ifft_1(1,5,complex,conj_even);
    ok          &= test_ifft_1(5,1,complex,conj_even);
    ok          &= test_ifft_1(1,6,complex,conj_even);
    ok          &= test_ifft_1(6,1,complex,conj_even);

    ok          &= test_ifft_1(1,5,complex,conj_even);
    ok          &= test_ifft_1(5,1,complex,conj_even);
    ok          &= test_ifft_1(1,6,complex,conj_even);
    ok          &= test_ifft_1(6,1,complex,conj_even);

    ok          &= test_ifft_1(2,5,complex,conj_even);
    ok          &= test_ifft_1(5,2,complex,conj_even);
    ok          &= test_ifft_1(2,6,complex,conj_even);
    ok          &= test_ifft_1(6,2,complex,conj_even);

    ok          &= test_ifft_1(5,5,complex,conj_even);
    ok          &= test_ifft_1(6,6,complex,conj_even);

    return ok;
};

bool test_fft::test_ifft_real()
{
    return test_ifft_type(false,false);
};
bool test_fft::test_ifft_complex()
{
    return test_ifft_type(true,false);
};
bool test_fft::test_ifft_real_conj_even()
{
    return test_ifft_type(false, true);
}
bool test_fft::test_ifft_complex_conj_even()
{
    return test_ifft_type(true, true);
};

void test_fft::make()
{
    bool ok;
    
    ok      = test_fft_pack();
    if (ok == true)
        std::cout << "test_fft_pack: passed" << "\n";
    else
        std::cout << "test_fft_pack: FAILED!" << "\n";    

    ok      = test_ifft_real_conj_even();
    if (ok == true)
        std::cout << "test_ifft_real_conj_even: passed" << "\n";
    else
        std::cout << "test_ifft_real_conj_even: FAILED!" << "\n";    

    ok      = test_ifft_complex_conj_even();
    if (ok == true)
        std::cout << "test_ifft_complex_conj_even: passed" << "\n";
    else
        std::cout << "test_ifft_complex_conj_even: FAILED!" << "\n";    

    ok      = test_fft_complex();
    if (ok == true)
        std::cout << "test_fft_complex: passed" << "\n";
    else
        std::cout << "test_fft_complex: FAILED!" << "\n";    

    ok      = test_fft_real();
    if (ok == true)
        std::cout << "test_fft_real: passed" << "\n";
    else
        std::cout << "test_fft_real: FAILED!" << "\n";    

    ok      = test_ifft_real();
    if (ok == true)
        std::cout << "test_ifft_real: passed" << "\n";
    else
        std::cout << "test_ifft_real: FAILED!" << "\n";    

    ok      = test_ifft_complex();
    if (ok == true)
        std::cout << "test_ifft_complex: passed" << "\n";
    else
        std::cout << "test_ifft_complex: FAILED!" << "\n";    
};

}}}
