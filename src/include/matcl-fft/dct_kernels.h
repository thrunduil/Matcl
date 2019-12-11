#pragma once

#include "matcl-fft/matcl_fft_config.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace fft
{

//-----------------------------------------------------------------------
//                          HELPERS
//-----------------------------------------------------------------------
namespace details
{

template<Integer M, Integer Max_M>
struct check_size
{
    static_assert(M <= Max_M && M >= 2, "unsupported size");

    using type = void;
};

};

//-----------------------------------------------------------------------
//                          OPTIMIZED KERNELS
//-----------------------------------------------------------------------

//------------------------- dct1 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void 
dct1_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void 
dct1_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scal);

//------------------------- dct2 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct2_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct2_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scal);

//------------------------- dct3 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scal);

//-----------------------------------------------------------------------
//                  OPTIMIZED KERNELS WITH STRIDES
//-----------------------------------------------------------------------
//These versions reorder data first and then call appropriate kernel.
//These versions can be much slower.

//------------------------- dct1 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, Integer in_stride, Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal, Integer in_stride, 
                            Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, 
                                 Integer inout_stride);

template<Integer M, class Test = typename details::check_size<M,9>::type>
MATCL_FFT_EXPORT void
dct1_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, double scal,
                                 Integer inout_stride);

//------------------------- dct2 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void 
dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, Integer in_stride, Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal, Integer in_stride, 
                            Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void 
dct2_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, 
                                 Integer inout_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void 
dct2_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, double scal,
                                 Integer inout_stride);

//------------------------- dct3 ----------------------------------------
template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, Integer in_stride, Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void 
dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
                            Integer out_ld, double scal, Integer in_stride, 
                            Integer out_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void 
dct3_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, 
                                 Integer inout_stride);

template<Integer M, class Test = typename details::check_size<M,8>::type>
MATCL_FFT_EXPORT void
dct3_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, double scal,
                                 Integer inout_stride);

//-----------------------------------------------------------------------
//                          dct1_kernel
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the first type in the form:
*
*       C_k     = 0.5*[x_0 + (-1)^k*x_{M-1}] + sum_{n=1}^{M-2} x_n*cos(pi*n*k/(M-1))
*       X_k     = 2 * scale * C_k
*
*  for k = 0 .. M-1, {x_k} : vector of orginal values of length M, 
*  {X_k} : vector of transformed values.
*  
*  Template Arguments
*  =======
*
*  M            - order of the transformation, 2 <= M <= 9
*
*  Inputs
*  =======
*
*  in           - array representing a matrix of size MxN to transform
*  out          - array representing resulting matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  in_ld        - leading dimension of in matrix
*  out_ld       - leading dimension of out matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  in_stride    - optional, define layout of in matrix, address of in(i,j) is 
*                   given by in + i * in_stride + j * in_ld, (i,j are 0-based)
*                   default value is 1.
*  out_stride   - optional, define layout of out matrix, address of out(i,j) is 
*                   given by out + i * out_stride + j * out_ld, (i,j are 0-based)
*                   default value is 1.
*       
*  Output
*  =======
*
*           Matrix out containing result of the transform. This transform is its own 
*           inverse up to scalling : dct1 ( dct1(x) ) = 2.0 * (M-1) * x.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
        Integer in_stride, Integer out_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
        double scal, Integer in_stride, Integer out_stride);

//-----------------------------------------------------------------------
//                          dct1_kernel_inpl
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  inplace version of dct1_kernel, see dct1_kernel for details
*
*  Inputs
*  =======
*
*  inout        - array representing a matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  inout_ld     - leading dimension of inout matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  inout_stride - optional, define layout of inout matrix, address of A(i,j) is 
*                   given by A + i * inout_stride + j * inout_ld, (i,j are 0-based)
*                   default value is 1.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, Integer inout_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct1_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, double scal, 
        Integer inout_stride);

//-----------------------------------------------------------------------
//                          dct2_kernel
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the second type in the form:
*
*       X_k     = 2.0 * scale * sum_{n=0}^{M-1} x_n*cos( pi*k*(n+0.5)/M )
*
*  for k = 0 .. M-1, {x_k} : vector of orginal values of length M, 
*  {X_k} : vector of transformed values.
*       
*  Template Arguments
*  =======
*
*  M        - order of the transformation, 2 <= M <= 8
*
*  Inputs
*  =======
*
*  in           - array representing a matrix of size MxN to transform
*  out          - array representing resulting matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  in_ld        - leading dimension of in matrix
*  out_ld       - leading dimension of out matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  in_stride    - optional, define layout of in matrix, address of in(i,j) is 
*                   given by in + i * in_stride + j * in_ld, (i,j are 0-based)
*                   default value is 1.
*  out_stride   - optional, define layout of out matrix, address of out(i,j) is 
*                   given by out + i * out_stride + j * out_ld, (i,j are 0-based)
*                   default value is 1.
*       
*  Output
*  =======
*
*           Matrix out containing result of the transform. The dct2 transform is the
*           inverse of dct3, i.e. (1/M) * dct2 ( dct3(x) ) = x.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void
dct2_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, 
        double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
        Integer out_ld, Integer in_stride, Integer out_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, 
        Integer out_ld, double scal, Integer in_stride, Integer out_stride);

//-----------------------------------------------------------------------
//                          dct2_kernel_inpl
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  inplace version of dct2_kernel, see dct2_kernel for details
*
*  Inputs
*  =======
*
*  inout        - array representing a matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  inout_ld     - leading dimension of inout matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  inout_stride - optional, define layout of inout matrix, address of A(i,j) is 
*                   given by A + i * inout_stride + j * inout_ld, (i,j are 0-based)
*                   default value is 1.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void
dct2_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, 
        Integer inout_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct2_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, 
        double scal, Integer inout_stride);

//-----------------------------------------------------------------------
//                          dct3
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the second type in the form:
*
*       X_k     = scale * [0.5*x_0 + sum_{n=1}^{M-1} x_n*cos( pi*n*(k+0.5)/M )
*
*  for k = 0 .. M-1, {x_k} : vector of orginal values of length M, 
*  {X_k} : vector of transformed values.
*       
*  Template Arguments
*  =======
*
*  M        - order of the transformation, 2 <= M <= 8
*
*  Inputs
*  =======
*
*  in           - array representing a matrix of size MxN to transform
*  out          - array representing resulting matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  in_ld        - leading dimension of in matrix
*  out_ld       - leading dimension of out matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  in_stride    - optional, define layout of in matrix, address of in(i,j) is 
*                   given by in + i * in_stride + j * in_ld, (i,j are 0-based)
*                   default value is 1.
*  out_stride   - optional, define layout of out matrix, address of out(i,j) is 
*                   given by out + i * out_stride + j * out_ld, (i,j are 0-based)
*                   default value is 1.
*       
*  Output
*  =======
*
*           Matrix out containing result of the transform. The dct3 transform is the
*           inverse of dct2, i.e. (1/M) * dct3 ( dct2(x) ) = x.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld,
        double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, 
        Integer in_stride, Integer out_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel_stride(const double* in, double* out, Integer N, Integer in_ld, Integer out_ld, 
        double scal, Integer in_stride, Integer out_stride);

//-----------------------------------------------------------------------
//                          dct3_kernel_inpl
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  inplace version of dct3_kernel, see dct3_kernel for details
*
*  Inputs
*  =======
*
*  inout        - array representing a matrix of size MxN to transform
*  N            - numer of columns of input and output matrix
*  inout_ld     - leading dimension of inout matrix
*  scale        - optional scaling (default = 1.0), no tests for special values
*                   (e.g. 1.0, 0.0) are performed
*  inout_stride - optional, define layout of inout matrix, address of A(i,j) is 
*                   given by A + i * inout_stride + j * inout_ld, (i,j are 0-based)
*                   default value is 1.
*/
template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel_inpl(double* inout, Integer N, Integer inout_ld);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel_inpl(double* inout, Integer N, Integer inout_ld, double scale);

template<Integer M, class Test>
MATCL_FFT_EXPORT void
dct3_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, Integer inout_stride);

template<Integer M, class Test>
MATCL_FFT_EXPORT void 
dct3_kernel_stride_inpl(double* inout, Integer N, Integer inout_ld, double scal, 
        Integer inout_stride);

}};
