
#pragma warning(push)
#pragma warning (disable : 4244)  //possible loss of data

#include "matcl-blas/level1/level1.h"
#include "matcl-scalar/matcl_scalar.h"

#pragma warning(pop)

using compl_t   = matcl::level1::simd_double_complex;
using fcompl_t  = matcl::level1::simd_single_complex;

struct X
{
    operator float() const          { return float(); };
    operator fcompl_t() const       { return fcompl_t(); };
    operator compl_t() const        { return compl_t(); };
};

X abs(const X&)                         { return X(); };
X operator-(const X&)                   { return X(); };
X operator+(const X&, const X&)         { return X(); };
X operator-(const X&, const X&)         { return X(); };
X operator*(const X&, const X&)         { return X(); };

X operator+(const X&, const compl_t&)   { return X(); };
X operator-(const X&, const compl_t&)   { return X(); };
X operator*(const X&, const compl_t&)   { return X(); };
X operator+(const X&, const float&)     { return X(); };
X operator-(const X&, const float&)     { return X(); };
X operator*(const X&, const float&)     { return X(); };

X operator+(const compl_t&, const X&)   { return X(); };
X operator-(const compl_t&, const X&)   { return X(); };
X operator*(const compl_t&, const X&)   { return X(); };
X operator+(const float&, const X&)     { return X(); };
X operator-(const float&, const X&)     { return X(); };
X operator*(const float&, const X&)     { return X(); };

namespace matcl { namespace level1
{

#define INSTANTIATE_MACRO(func)     \
template struct func<int, 0>;       \
template struct func<float, 0>;     \
template struct func<double, 0>;    \
template struct func<compl_t, 0>;   \
template struct func<fcompl_t, 0>;  \
template struct func<X, 0>;

#define INSTANTIATE_TEST_MACRO(func,sel, R, C)\
template struct func<sel, int, R,C,0>;       \
template struct func<sel, float, R,C,0>;     \
template struct func<sel, double, R,C,0>;    \
template struct func<sel, compl_t, R,C,0>;   \
template struct func<sel, fcompl_t, R,C,0>;  \
template struct func<sel, X, R,C,0>;

#define INSTANTIATE_TEST2_MACRO(func)      \
INSTANTIATE_TEST_MACRO(func, true, 0, 0)   \
INSTANTIATE_TEST_MACRO(func, true, 0, 1)   \
INSTANTIATE_TEST_MACRO(func, true, 1, 0)   \
INSTANTIATE_TEST_MACRO(func, true, 1, 1)   \
INSTANTIATE_TEST_MACRO(func, false, 0, 0)  \
INSTANTIATE_TEST_MACRO(func, false, 0, 1)  \
INSTANTIATE_TEST_MACRO(func, false, 1, 0)  \
INSTANTIATE_TEST_MACRO(func, false, 1, 1)

#define INSTANTIATE_TEST_MACRO2(func,sel, R, C)    \
template struct func<sel,int, int, R, C, 0>;  \
template struct func<sel,int, float, R, C, 0>;\
                                               \
template struct func<sel,float, int, R, C, 0>;        \
template struct func<sel,float, float, R, C, 0>;      \
template struct func<sel,float, double, R, C, 0>;     \
                                                       \
template struct func<sel,double, int, R, C, 0>;       \
template struct func<sel,double, float, R, C, 0>;     \
template struct func<sel,double, double, R, C, 0>;    \
                                                       \
template struct func<sel,compl_t, int, R, C, 0>;      \
template struct func<sel,compl_t, float, R, C, 0>;    \
template struct func<sel,compl_t, double, R, C, 0>;   \
template struct func<sel,compl_t, fcompl_t, R, C, 0>; \
template struct func<sel,compl_t, compl_t, R, C, 0>;  \
                                                       \
template struct func<sel,fcompl_t, float, R, C, 0>;   \
template struct func<sel,fcompl_t, fcompl_t, R, C, 0>;

//template struct func<true,fcompl_t, double, 0, 0, 0>;
//template struct func<true,fcompl_t, int, 0, 0, 0>;
//template struct func<true,fcompl_t, compl_t, 0, 0, 0>;

#define INSTANTIATE_TEST2_MACRO2(func)      \
INSTANTIATE_TEST_MACRO2(func, true, 0, 0)   \
INSTANTIATE_TEST_MACRO2(func, true, 0, 1)   \
INSTANTIATE_TEST_MACRO2(func, true, 1, 0)   \
INSTANTIATE_TEST_MACRO2(func, true, 1, 1)   \
INSTANTIATE_TEST_MACRO2(func, false, 0, 0)  \
INSTANTIATE_TEST_MACRO2(func, false, 0, 1)  \
INSTANTIATE_TEST_MACRO2(func, false, 1, 0)  \
INSTANTIATE_TEST_MACRO2(func, false, 1, 1)

#define INSTANTIATE_MACRO2(func)    \
template struct func<int, int, 0>;  \
template struct func<int, float, 0>;\
template struct func<X, X, 0>;      \
                                    \
template struct func<float, int, 0>;        \
template struct func<float, float, 0>;      \
template struct func<float, double, 0>;     \
template struct func<float, X, 0>;          \
                                            \
template struct func<double, int, 0>;       \
template struct func<double, float, 0>;     \
template struct func<double, double, 0>;    \
                                            \
template struct func<compl_t, int, 0>;      \
template struct func<compl_t, float, 0>;    \
template struct func<compl_t, double, 0>;   \
template struct func<compl_t, compl_t, 0>;  \
template struct func<compl_t, fcompl_t, 0>; \
template struct func<compl_t, X, 0>;        \
                                            \
template struct func<fcompl_t, float, 0>;   \
template struct func<fcompl_t, fcompl_t, 0>;

//template struct func<fcompl_t, double, 0>;
//template struct func<fcompl_t, int, 0>;
//template struct func<fcompl_t, compl_t, 0>;

#define INSTANTIATE_MACRO3(func)                    \
template struct func<int, int, int, 0>;             \
template struct func<float, float, float, 0>;       \
template struct func<float, float, double, 0>;      \
template struct func<double, float, float, 0>;      \
template struct func<double, double, double, 0>;    \
template struct func<compl_t, float, float, 0>;     \
template struct func<compl_t, double, compl_t, 0>;  \
template struct func<fcompl_t, float, float, 0>;    \
template struct func<fcompl_t, float, fcompl_t, 0>; \
template struct func<fcompl_t, fcompl_t, float, 0>; \
template struct func<fcompl_t, fcompl_t, fcompl_t, 0>;

#define INSTANTIATE_TEST_MACRO3(func,sel,R,C)       \
template struct func<sel,int, int, int, R,C,0>;             \
template struct func<sel,float, float, float, R,C,0>;       \
template struct func<sel,float, float, double, R,C,0>;      \
template struct func<sel,double, float, float, R,C,0>;      \
template struct func<sel,double, double, double, R,C,0>;    \
template struct func<sel,compl_t, float, float, R,C,0>;     \
template struct func<sel,compl_t, double, compl_t, R,C,0>;  \
template struct func<sel,fcompl_t, float, float, R,C,0>;    \
template struct func<sel,fcompl_t, float, fcompl_t, R,C,0>; \
template struct func<sel,fcompl_t, fcompl_t, float, R,C,0>; \
template struct func<sel,fcompl_t, fcompl_t, fcompl_t, R,C,0>;

#define INSTANTIATE_TEST2_MACRO3(func)      \
INSTANTIATE_TEST_MACRO3(func, true, 0, 0)   \
INSTANTIATE_TEST_MACRO3(func, true, 0, 1)   \
INSTANTIATE_TEST_MACRO3(func, true, 1, 0)   \
INSTANTIATE_TEST_MACRO3(func, true, 1, 1)   \
INSTANTIATE_TEST_MACRO3(func, false, 0, 0)  \
INSTANTIATE_TEST_MACRO3(func, false, 0, 1)  \
INSTANTIATE_TEST_MACRO3(func, false, 1, 0)  \
INSTANTIATE_TEST_MACRO3(func, false, 1, 1)

#define INSTANTIATE_MACRO4(func)                            \
template struct func<int, int, int, int, 0>;                \
template struct func<float, float, float, float, 0>;        \
template struct func<float, float, double, float, 0>;       \
template struct func<float, float, double, double, 0>;      \
template struct func<double, float, float, double, 0>;      \
template struct func<double, float, float, float, 0>;       \
template struct func<double, double, double, float, 0>;     \
template struct func<double, double, double, double, 0>;    \
template struct func<compl_t, float, float, float, 0>;      \
template struct func<compl_t, float, float, double, 0>;     \
template struct func<compl_t, float, float, compl_t, 0>;    \
template struct func<compl_t, double, compl_t, float, 0>;   \
template struct func<compl_t, double, compl_t, double, 0>;  \
template struct func<compl_t, double, compl_t, compl_t, 0>; \
template struct func<fcompl_t, float, float, float, 0>;     \
template struct func<fcompl_t, float, float, fcompl_t, 0>;  \
template struct func<fcompl_t, float, fcompl_t, float, 0>;  \
template struct func<fcompl_t, float, fcompl_t, fcompl_t, 0>;\
template struct func<fcompl_t, fcompl_t, float, float, 0>;  \
template struct func<fcompl_t, fcompl_t, float, fcompl_t, 0>;\
template struct func<fcompl_t, fcompl_t, fcompl_t, float, 0>;\
template struct func<fcompl_t, fcompl_t, fcompl_t, fcompl_t, 0>;

#define INSTANTIATE_TEST_MACRO4(func,sel, R, C)    \
template struct func<sel, int, int, int, int, R, C, 0>;                \
template struct func<sel, float, float, float, float, R, C, 0>;        \
template struct func<sel, float, float, double, float, R, C, 0>;       \
template struct func<sel, float, float, double, double, R, C, 0>;      \
template struct func<sel, double, float, float, double, R, C, 0>;      \
template struct func<sel, double, float, float, float, R, C, 0>;       \
template struct func<sel, double, double, double, float, R, C, 0>;     \
template struct func<sel, double, double, double, double, R, C, 0>;    \
template struct func<sel, compl_t, float, float, float, R, C, 0>;      \
template struct func<sel, compl_t, float, float, double, R, C, 0>;     \
template struct func<sel, compl_t, float, float, compl_t, R, C, 0>;    \
template struct func<sel, compl_t, double, compl_t, float, R, C, 0>;   \
template struct func<sel, compl_t, double, compl_t, double, R, C, 0>;  \
template struct func<sel, compl_t, double, compl_t, compl_t, R, C, 0>; \
template struct func<sel, fcompl_t, float, float, float, R, C, 0>;     \
template struct func<sel, fcompl_t, float, float, fcompl_t, R, C, 0>;  \
template struct func<sel, fcompl_t, float, fcompl_t, float, R, C, 0>;  \
template struct func<sel, fcompl_t, float, fcompl_t, fcompl_t, R, C, 0>;\
template struct func<sel, fcompl_t, fcompl_t, float, float, R, C, 0>;  \
template struct func<sel, fcompl_t, fcompl_t, float, fcompl_t, R, C, 0>;\
template struct func<sel, fcompl_t, fcompl_t, fcompl_t, float, R, C, 0>;\
template struct func<sel, fcompl_t, fcompl_t, fcompl_t, fcompl_t, R, C, 0>;

#define INSTANTIATE_TEST2_MACRO4(func)      \
INSTANTIATE_TEST_MACRO4(func, true, 0, 0)   \
INSTANTIATE_TEST_MACRO4(func, true, 0, 1)   \
INSTANTIATE_TEST_MACRO4(func, true, 1, 0)   \
INSTANTIATE_TEST_MACRO4(func, true, 1, 1)   \
INSTANTIATE_TEST_MACRO4(func, false, 0, 0)  \
INSTANTIATE_TEST_MACRO4(func, false, 0, 1)  \
INSTANTIATE_TEST_MACRO4(func, false, 1, 0)  \
INSTANTIATE_TEST_MACRO4(func, false, 1, 1)

struct Fun1
{
    template<class TY>
    using enable_simd   = details::check_simd<TY>;

    template<class T>
    static T eval(const T& in)
    {
        return -in;
    };

    template<class T>
    static T eval_simd(const T& in)
    {
        return -in;
    };

};

struct Fun2
{
    template<class TY, class TA>
    using enable_simd   = details::check_simd<TY>;

    template<class T, class TA>
    static T eval(const T& in, const TA& a)
    {
        return a * in;
    };

    template<class T, class TA>
    static T eval_simd(const T& in, const TA& a)
    {
        return T(a) * in;
    };
};

template<class TY, Integer Rows>
struct eval_scal_func_Y_impl1
{
    using Func = Fun1;

    static void eval(TY* Y, Integer rows, const Func& f)
    {
        return eval_scal_func_Y<TY, Rows, Func>::eval(Y,rows, f);
    }
    static void eval(TY* Y, Integer Y_step, Integer rows, const Func& f)
    {
        return eval_scal_func_Y<TY, Rows, Func>::eval(Y,Y_step,rows, f);
    }
};

template<class TY, Integer Rows>
struct eval_scal_func_Y_impl2
{
    using TA    = TY;
    using Func  = Fun2;

    static void eval(TY* Y, Integer rows, const TA& a, const Func& f)
    {
        return eval_scal_func_Y<TY, Rows, Func, const TA&>::eval(Y,rows, f, a);
    }
    static void eval(TY* Y, Integer Y_step, Integer rows, const TA& a, const Func& f)
    {
        return eval_scal_func_Y<TY, Rows, Func, const TA& >::eval(Y, Y_step, rows, f, a);
    }
};

// one argument
INSTANTIATE_MACRO(my)
INSTANTIATE_MACRO(swap)
INSTANTIATE_MACRO(dot)
INSTANTIATE_MACRO(apy)
INSTANTIATE_MACRO(amy)
INSTANTIATE_MACRO(set_val)
INSTANTIATE_MACRO(eval_scal_func_Y_impl1)
INSTANTIATE_MACRO(eval_scal_func_Y_impl2)
INSTANTIATE_TEST2_MACRO(set_val_mat)

// two arguments
INSTANTIATE_MACRO2(mx)
INSTANTIATE_MACRO2(ay)
INSTANTIATE_TEST2_MACRO2(ay_test_mat)
INSTANTIATE_MACRO2(xmy)
INSTANTIATE_MACRO2(ymx)
INSTANTIATE_TEST2_MACRO2(ymx_mat)
INSTANTIATE_MACRO2(ypx)
INSTANTIATE_MACRO2(ypxm)
INSTANTIATE_MACRO2(copy)
INSTANTIATE_TEST2_MACRO2(copy_mat)
INSTANTIATE_MACRO2(apby)
INSTANTIATE_TEST2_MACRO2(apby_test_mat)
INSTANTIATE_MACRO2(apx_abs)
INSTANTIATE_MACRO2(x_abs)
INSTANTIATE_TEST2_MACRO2(apx_abs_test_mat)

//three arguments
INSTANTIATE_MACRO3(ax)
INSTANTIATE_MACRO3(axpy)
INSTANTIATE_TEST2_MACRO3(axpy_test_mat)
INSTANTIATE_TEST2_MACRO3(ax_test_mat)
INSTANTIATE_MACRO3(axmy)
INSTANTIATE_MACRO3(xpya)
INSTANTIATE_MACRO3(xmya)
INSTANTIATE_MACRO3(xpby)
INSTANTIATE_MACRO3(mxpby)
INSTANTIATE_MACRO3(xy)

//four arguments
INSTANTIATE_MACRO4(axpby)
INSTANTIATE_MACRO4(axpby_test)
INSTANTIATE_MACRO4(axypz)
INSTANTIATE_TEST2_MACRO4(axpby_test_mat)

}}