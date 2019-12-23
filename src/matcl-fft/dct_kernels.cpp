#include "matcl-fft/dct_kernels.h"
#include "dct_kernels_test.h"
#include "matcl-blas-lapack/level1/level1.h"

namespace matcl { namespace fft
{

//-----------------------------------------------------------------------
//                          dct1_kernel
//-----------------------------------------------------------------------
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel(const double* in, double* out, Integer N, 
                                  Integer in_ld, Integer out_ld)
{
    using evaler    = dct1_evaler_type<M>;

    //std::cout << std::string(50, '-') << "\n";
    //std::cout << "DCT1 " << M << "\n";
    //evaler::print();

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);

        in  += in_ld;
        out += out_ld;
    };
};

template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, double scal)
{
    using evaler    = dct1_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);
        level1::ay<double, double, M, 1>::eval(out, M, scal);

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            Integer in_stride, Integer out_stride)
{
    using evaler    = dct1_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_1[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_1,in2);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_1[j];

        in  += in_ld;
        out += out_ld;
    };
};

template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            double scal, Integer in_stride, Integer out_stride)
{
    using evaler    = dct1_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_2[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_2,in2);
        level1::ay<double, double, M, 1>::eval(out2_2, M, scal);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_2[j];

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_inpl(double* in, Integer N, Integer in_ld)
{
    using evaler    = dct1_evaler_type<M>;

    double tmp_1[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_1,in);

        level1::copy_mat<false,double,double,M,1,1>::eval(in, in_ld, tmp_1, M, M, 1);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_inpl(double* in, Integer N, Integer in_ld, double scal)
{
    using evaler    = dct1_evaler_type<M>;

    double tmp_2[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_2,in);

        level1::ax<double, double, double, M, 1>::eval(in,tmp_2, M, scal);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_stride_inpl(double* in, Integer N, Integer in_ld, Integer in_stride)
{
    using evaler    = dct1_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_3[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_3,in2);

        for (Integer j = 0; j < M; ++j)
        {
            in[j*in_stride]     = out2_3[j];
        };

        in  += in_ld;
    };
};

template<Integer M, class Test>
MATCL_FFT_EXPORT void dct1_kernel_stride_inpl(double* in, Integer N, Integer in_ld, double scal, Integer in_stride)
{
    using evaler    = dct1_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_4[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_4, in2);
        level1::ay<double, double, M, 1>::eval(out2_4, M, scal);

        for (Integer j = 0; j < M; ++j)
            in[j*in_stride]     = out2_4[j];

        in  += in_ld;
    };
};

//-----------------------------------------------------------------------
//                          dct2_kernel
//-----------------------------------------------------------------------
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, double scal)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);
        level1::ay<double, double, M, 1>::eval(out, M, scal);
        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            Integer in_stride, Integer out_stride)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_5[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_5, in2);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_5[j];

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            double scal, Integer in_stride, Integer out_stride)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_6[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_6, in2);
        level1::ay<double, double, M, 1>::eval(out2_6, M, scal);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_6[j];

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_inpl(double* in, Integer N, Integer in_ld)
{
    using evaler    = dct2_evaler_type<M>;

    double tmp_3[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_3,in);

        level1::copy_mat<false, double, double, M, 1, 1>::eval(in, in_ld, tmp_3, M, M, 1);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_inpl(double* in, Integer N, Integer in_ld, double scal)
{
    using evaler    = dct2_evaler_type<M>;

    double tmp_4[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_4,in);

        level1::ax<double, double, double, M, 1>::eval(in,tmp_4, M, scal);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_stride_inpl(double* in, Integer N, Integer in_ld, Integer in_stride)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_7[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_7,in2);

        for (Integer j = 0; j < M; ++j)
            in[j*in_stride]     = out2_7[j];

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct2_kernel_stride_inpl(double* in, Integer N, Integer in_ld, double scal, Integer in_stride)
{
    using evaler    = dct2_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_8[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_8, in2);
        level1::ay<double, double, M, 1>::eval(out2_8, M, scal);

        for (Integer j = 0; j < M; ++j)
            in[j*in_stride]     = out2_8[j];

        in  += in_ld;
    };
};
//-----------------------------------------------------------------------
//                          dct3_kernel
//-----------------------------------------------------------------------
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, double scal)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(out,in);
        level1::ay<double, double, M, 1>::eval(out, M, scal);
        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            Integer in_stride, Integer out_stride)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_9[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_9, in2);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_9[j];

        in  += in_ld;
        out += out_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
                            double scal, Integer in_stride, Integer out_stride)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_10[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_10, in2);
        level1::ay<double, double, M, 1>::eval(out2_10, M, scal);

        for (Integer j = 0; j < M; ++j)
            out[j*out_stride]   = out2_10[j];

        in  += in_ld;
        out += out_ld;
    };
};

template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_inpl(double* in, Integer N, Integer in_ld)
{
    using evaler    = dct3_evaler_type<M>;

    double tmp_5[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_5, in);

        level1::copy_mat<false,double,double,M,1,1>::eval(in, in_ld, tmp_5, M, M, 1);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_inpl(double* in, Integer N, Integer in_ld, double scal)
{
    using evaler    = dct3_evaler_type<M>;

    double tmp_6[M];

    for(Integer i = 0; i < N; ++i)
    {
        evaler::eval(tmp_6,in);

        level1::ax<double, double, double, M, 1>::eval(in,tmp_6, M, scal);

        in  += in_ld;
    };
};
template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_stride_inpl(double* in, Integer N, Integer in_ld, Integer in_stride)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_11[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_11,in2);

        for (Integer j = 0; j < M; ++j)
            in[j*in_stride]     = out2_11[j];

        in  += in_ld;
    };
};

template<Integer M, class Test>
MATCL_FFT_EXPORT void dct3_kernel_stride_inpl(double* in, Integer N, Integer in_ld, double scal, Integer in_stride)
{
    using evaler    = dct3_evaler_type<M>;

    for(Integer i = 0; i < N; ++i)
    {
        double in2[M];
        double out2_12[M];

        for (Integer j = 0; j < M; ++j)
            in2[j]              = in[j*in_stride];

        evaler::eval(out2_12, in2);
        level1::ay<double, double, M, 1>::eval(out2_12, M, scal);

        for (Integer j = 0; j < M; ++j)
            in[j*in_stride]     = out2_12[j];

        in  += in_ld;
    };
};

#define instantiate_dct(func,M)                                                                         \
template MATCL_FFT_EXPORT void func<M>(const double*, double*, Integer, Integer, Integer);                    \
template MATCL_FFT_EXPORT void func<M>(const double*, double*, Integer, Integer, Integer, double);            \
template MATCL_FFT_EXPORT void func##_stride<M>(const double*, double*, Integer, Integer, Integer, Integer, Integer);         \
template MATCL_FFT_EXPORT void func##_stride<M>(const double*, double*, Integer, Integer, Integer, double, Integer, Integer); \
template MATCL_FFT_EXPORT void func##_inpl<M>(double*, Integer, Integer);                                     \
template MATCL_FFT_EXPORT void func##_inpl<M>(double*, Integer, Integer, double);                             \
template MATCL_FFT_EXPORT void func##_stride_inpl<M>(double*, Integer, Integer, Integer);                     \
template MATCL_FFT_EXPORT void func##_stride_inpl<M>(double*, Integer, Integer, double, Integer);                    

instantiate_dct(dct1_kernel,2)
instantiate_dct(dct1_kernel,3)
instantiate_dct(dct1_kernel,4)
instantiate_dct(dct1_kernel,5)
instantiate_dct(dct1_kernel,6)
instantiate_dct(dct1_kernel,7)
instantiate_dct(dct1_kernel,8)
instantiate_dct(dct1_kernel,9)

instantiate_dct(dct2_kernel,2)
instantiate_dct(dct2_kernel,3)
instantiate_dct(dct2_kernel,4)
instantiate_dct(dct2_kernel,5)
instantiate_dct(dct2_kernel,6)
instantiate_dct(dct2_kernel,7)
instantiate_dct(dct2_kernel,8)

instantiate_dct(dct3_kernel,2)
instantiate_dct(dct3_kernel,3)
instantiate_dct(dct3_kernel,4)
instantiate_dct(dct3_kernel,5)
instantiate_dct(dct3_kernel,6)
instantiate_dct(dct3_kernel,7)
instantiate_dct(dct3_kernel,8)

}};
