#include "dct_tests.h"
#include "matcl-matrep\matcl_matrep.h"
#include "matcl-fft/matcl_fft.h"
#include "matcl-fft/matcl_dct.h"
#include "matcl-linalg\matcl_linalg.h"

#include <iostream>

namespace matcl { namespace fft { namespace test
{

bool test_dct::test_dct1_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Real dif    = 0;

    Matrix x    = randn(r,c);
    Matrix y1   = dct1(x,1,scal);

    {        
        Matrix x2   = trans(x);        
        Matrix y2   = dct1(x2,2,scal);
        y2          = trans(y2);

        dif         += norm(y1 - y2, 1);
    };    

    Matrix z;
    if (r > 1 && c > 0)
    {
        Matrix n    = matcl::range(1,r-2);
        Matrix m    = matcl::range(0,r-1);
        Matrix T    = cos(constants::pi() * trans(m)*n/Real(r-1));
        Matrix x2   = x;
        Matrix Z2   = T * x2(colon(2,end-1),colon()) 
                    + 0.5 * (ones(length(m),1) * x2(1,colon()) + pow(-1,trans(m)) * x2(end,colon()));
        z           = 2.0 * scal * Z2;   
    }
    else
    {
        z           = y1;
    };

    dif             += norm(y1 - z, 1);

    if (dif > tol)
        ok      = false;

    return ok;
};

bool test_dct::test_dct2_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct2(x,1,scal);
    Matrix y2   = dct2(x2,2,scal);
    y2          = trans(y2);

    Real dif1   = norm(y1 - y2, 1);

    Matrix z;
    if (r > 1 && c > 0)
    {
        Matrix n    = matcl::range(0,r-1);
        Matrix m    = matcl::range(0,r-1);
        Matrix T    = cos(constants::pi() * trans(m)*(2*n+1)/2.0/Real(r));
        Matrix x3   = x;//flipud(x);
        Matrix Z3   = T * x3 * (2.0*scal); 
        z           = Z3;
    }
    else
    {
        z           = y1;
    };

    Real dif2   = norm(y1 - z, 1);

    Real dif    = dif1 + dif2;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};
bool test_dct::test_dct2_MH_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct2_MH(x,scal);

    Matrix z;
    if (r > 1 && c > 0)
    {
        Matrix n    = matcl::range(0,r-1);
        Matrix m    = matcl::range(0,r-1);
        Matrix T    = cos(constants::pi() * trans(m)*(2*n+1)/2.0/Real(r));
        Matrix x3   = x;//flipud(x);
        Matrix Z3   = T * x3 * (2.0 *scal); 
        z           = Z3;
    }
    else
    {
        z           = y1;
    };

    Real dif2   = norm(y1 - z, 1);

    Real dif    = dif2;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};

bool test_dct::test_dct3_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct3(x,1,scal);
    Matrix y2   = dct3(x2,2,scal);
    y2          = trans(y2);

    Real dif1   = norm(y1 - y2, 1);

    Matrix z;
    if (r > 1 && c > 0)
    {
        Matrix n    = matcl::range(0,r-1);
        Matrix m    = matcl::range(1,r-1);
        Matrix T    = cos(constants::pi() * (2*trans(n)+1) * m /2.0/Real(r));
        Matrix x3   = x;//flipud(x);
        Matrix Z3   = (T * x3(colon(2,end),colon()) + ones(r,1) * x(1, colon()) / 2.0); 
        z           = Z3 * scal;
    }
    else
    {
        z           = y1;
    };

    Real dif2   = norm(y1 - z, 1);

    Matrix w    = dct2(y1, 1);
    Real scale  = 1.0 / r / scal;
    Real dif3   = norm(scale*w-x,1);

    Real dif    = dif1 + dif2 + dif3;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};

bool test_dct::test_dct1_inpl_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct1(x,1,scal);
    Matrix y2   = dct1(x2,2,scal);

    Matrix ix   = x;
    Matrix ix2  = x2;
    
    ix.make_unique();
    ix2.make_unique();

    Matrix z1   = dct1(std::move(ix),1,scal);
    Matrix z2   = dct1(std::move(ix2),2,scal);

    Real dif1   = norm(y1 - z1, 1);
    Real dif2   = norm(y2 - z2, 1);

    Real dif    = dif1 + dif2;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};

bool test_dct::test_dct2_inpl_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct2(x,1,scal);
    Matrix y2   = dct2(x2,2,scal);

    Matrix ix   = x;
    Matrix ix2  = x2;
    
    ix.make_unique();
    ix2.make_unique();

    Matrix z1   = dct2(std::move(ix),1,scal);
    Matrix z2   = dct2(std::move(ix2),2,scal);

    Real dif1   = norm(y1 - z1, 1);
    Real dif2   = norm(y2 - z2, 1);

    Real dif    = dif1 + dif2;

    if (dif > tol)
    {
        disp(x);
        disp(y1);
        disp(y2);
        disp(z1);
        disp(z2);

        ok      = false;
    }

    return ok;
};
bool test_dct::test_dct2_MH_inpl_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);

    Matrix y1   = dct2_MH(x,scal);

    Matrix ix   = x;
    
    ix.make_unique();

    Matrix z1   = dct2_MH(std::move(ix),scal);

    Real dif1   = norm(y1 - z1, 1);

    Real dif    = dif1;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};

bool test_dct::test_dct3_inpl_1(Integer r, Integer c, Real scal)
{
    bool ok     = true;
    Real tol    = 1e-13*r;

    Matrix x    = randn(r,c);
    Matrix x2   = trans(x);

    Matrix y1   = dct3(x,1,scal);
    Matrix y2   = dct3(x2,2,scal);

    Matrix ix   = x;
    Matrix ix2  = x2;
    
    ix.make_unique();
    ix2.make_unique();

    Matrix z1   = dct3(std::move(ix),1,scal);
    Matrix z2   = dct3(std::move(ix2),2,scal);

    Real dif1   = norm(y1 - z1, 1);
    Real dif2   = norm(y2 - z2, 1);

    Real dif    = dif1 + dif2;

    if (dif > tol)
    {
        ok      = false;
    }

    return ok;
};

bool test_dct::test_dct1()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct1_1(r,c, 1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct1_1(r,c, 2.0);

    return ok;
};
bool test_dct::test_dct1_inpl()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct1_inpl_1(r,c, 1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct1_inpl_1(r,c, 2.0);

    return ok;
};

bool test_dct::test_dct2()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_1(r,c, 1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_1(r,c, 2.0);

    return ok;
};
bool test_dct::test_dct2_MH()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_MH_1(r,c,1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_MH_1(r,c,2.0);

    return ok;
};

bool test_dct::test_dct2_inpl()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_inpl_1(r,c,1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_inpl_1(r,c,2.0);

    return ok;
};
bool test_dct::test_dct2_MH_inpl()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_MH_inpl_1(r,c,1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct2_MH_inpl_1(r,c,2.0);

    return ok;
};

bool test_dct::test_dct3()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct3_1(r,c,1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct3_1(r,c,2.0);

    return ok;
};
bool test_dct::test_dct3_inpl()
{
    bool ok     = true;

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct3_inpl_1(r,c,1.0);

    for (Integer c = 0; c < 40; ++c)
        for (Integer r = 0; r < 40; ++r)
            ok  &= test_dct3_inpl_1(r,c,2.0);

    return ok;
};
void test_dct::make()
{
    bool ok;

    ok      = test_dct2_MH();
    if (ok == true)
        std::cout << "test_dct2_MH: passed" << "\n";
    else
        std::cout << "test_dct2_MH: FAILED!" << "\n";  

    ok      = test_dct2_MH_inpl();
    if (ok == true)
        std::cout << "test_dct2_MH_inpl: passed" << "\n";
    else
        std::cout << "test_dct2_MH_inpl: FAILED!" << "\n";  

    ok      = test_dct1();
    if (ok == true)
        std::cout << "test_dct1: passed" << "\n";
    else
        std::cout << "test_dct1: FAILED!" << "\n";    

    ok      = test_dct1_inpl();
    if (ok == true)
        std::cout << "test_dct1_inpl: passed" << "\n";
    else
        std::cout << "test_dct1_inpl: FAILED!" << "\n";    

    ok      = test_dct2();
    if (ok == true)
        std::cout << "test_dct2: passed" << "\n";
    else
        std::cout << "test_dct2: FAILED!" << "\n";  

    ok      = test_dct2_inpl();
    if (ok == true)
        std::cout << "test_dct2_inpl: passed" << "\n";
    else
        std::cout << "test_dct2_inpl: FAILED!" << "\n";  

    ok      = test_dct3();
    if (ok == true)
        std::cout << "test_dct3: passed" << "\n";
    else
        std::cout << "test_dct3: FAILED!" << "\n";  

    ok      = test_dct3_inpl();
    if (ok == true)
        std::cout << "test_dct3_inpl: passed" << "\n";
    else
        std::cout << "test_dct3_inpl: FAILED!" << "\n";  

};

}}};
