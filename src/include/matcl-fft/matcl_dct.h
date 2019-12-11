#pragma once

#include "matcl-fft/matcl_fft_config.h"

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-fft/matcl_fft_exception.h"
#include "matcl-fft/fft_context.h"

#include "matcl-matrep/matrix/matrix.h"

namespace matcl { namespace fft
{

//-----------------------------------------------------------------------
//                          low level dct
//-----------------------------------------------------------------------
MATCL_FFT_EXPORT void
dct1(const Real* X, Integer X_ld, Integer rows, Integer cols,
    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
dct1_inplace(Real* Y, Integer Y_ld, Integer rows, Integer cols,
    Integer dim, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
dct2(const Real* X, Integer X_ld, Integer rows, Integer cols,
    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
dct2_inplace(Real* X, Integer X_ld, Integer rows, Integer cols,
    Integer dim, Real scale, fft_context& cont);

//work is of size at least rows if dim == 1 or cols if dim == 2
MATCL_FFT_EXPORT void
dct2_MH(const Real* X, Integer X_ld, Integer rows, Integer cols,
    Real* Y, Integer Y_ld, Real scale, Real* work, Integer work_size, fft_context& cont);

MATCL_FFT_EXPORT void
dct2_MH_inplace(Real* X, Integer X_ld, Integer rows, Integer cols,
    Real scale, Real* work, Integer work_size, fft_context& cont);

MATCL_FFT_EXPORT void
dct3(const Real* X, Integer X_ld, Integer rows, Integer cols,
    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
dct3_inplace(Real* X, Integer X_ld, Integer rows, Integer cols,
    Integer dim, Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          high level dct
//-----------------------------------------------------------------------
MATCL_FFT_EXPORT Matrix dct1(const matcl::Matrix& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct1(matcl::Matrix&& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct1(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct1(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

MATCL_FFT_EXPORT Matrix dct2(const matcl::Matrix& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct2(matcl::Matrix&& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct2(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct2(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

MATCL_FFT_EXPORT Matrix dct2_MH(const matcl::Matrix& x, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct2_MH(const matcl::Matrix& x, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct2_MH(matcl::Matrix&& x, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct2_MH(matcl::Matrix&& x, Real scale, fft_context& cont);

MATCL_FFT_EXPORT Matrix dct3(const matcl::Matrix& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct3(matcl::Matrix&& x, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix dct3(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct3(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          dct1
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the first type in the form:
*
*       C_k     = 0.5*[x_0 + (-1)^k*x_{N-1}] + sum_{n=1}^{N-2} x_n*cos(pi*n*k/(N-1))
*       X_k     = 2 * scale * C_k
*
*  for k = 0 .. N-1, {x_k} : vector of orginal values of length N, 
*  {X_k} : vector of transformed values.
*       
*  Inputs
*  =======
*
*  mat      - matrix of size MxK to transform, if mat is a temporary and unique
*             then the transform is done in place if possible
*  dim      - dimension along which transormation must be conducted, if dim == 1,
*             then N = M, if dim == 2, then N = K.
*  scal     - scaling factor, default value = 1.0.
*  cont     - optional parameter, stores internal data. If on input cont is empty,
*             i.e. cont = fft_context(), or not initialized for problems of given type,
*             then cont is initialized, and can be reused in next computations of the
*             same type. If given context is used in computations of different type,
*             then data stored in cont are not destroyed and can be used later. A problem
*             is determined by length (number of rows if dim == 1, or number of columns 
*             if dim == 2), and type (dct1, dct_23, or fft). On return cont is a context
*             appended with data for problems of given type.
*       
*  Output
*  =======
*
*           Matrix containing result of the transform. This transform is its own 
*           inverse up to scalling : dct1 ( dct1(x) ) = 2.0 * (N-1) * x.
*/
MATCL_FFT_EXPORT Matrix dct1(const matcl::Matrix& mat, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct1(matcl::Matrix&& x, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct1(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct1(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          dct2
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the second type in the form:
*
*       X_k     = 2.0 * scale * sum_{n=0}^{N-1} x_n*cos( pi*k*(n+0.5)/N )
*
*  for k = 0 .. N-1, {x_k} : vector of orginal values of length N, 
*  {X_k} : vector of transformed values.
*       
*  Inputs
*  =======
*
*  mat      - matrix of size MxK to transform, if mat is a temporary and unique
*             then the transform is done in place if possible
*  dim      - dimension along which transormation must be conducted, if dim == 1,
*             then N = M, if dim == 2, then N = K.
*  dim      - dimension along which transormation must be conducted
*  scal     - scaling factor, default value = 1.0.
*  cont     - optional parameter, stores internal data. If on input cont is empty,
*             i.e. cont = fft_context(), or not initialized for problems of given type,
*             then cont is initialized, and can be reused in next computations of the
*             same type. If given context is used in computations of different type,
*             then data stored in cont are not destroyed and can be used later. A problem
*             is determined by length (number of rows if dim == 1, or number of columns 
*             if dim == 2), and type (dct1, dct_23, or fft). On return cont is a context
*             appended with data for problems of given type.
*       
*  Output
*  =======
*
*           Matrix containing result of the transform. The dct2 transform is the
*           inverse of dct3, i.e. (1/N) * dct2 ( dct3(x) ) = x.
*/
MATCL_FFT_EXPORT Matrix dct2(const matcl::Matrix& mat, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct2(matcl::Matrix&& x, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct2(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct2(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          dct3
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete cosine transform of the second type in the form:
*
*       X_k     = scale * [0.5*x_0 + sum_{n=1}^{N-1} x_n*cos( pi*n*(k+0.5)/N )
*
*  for k = 0 .. N-1, {x_k} : vector of orginal values of length N, 
*  {X_k} : vector of transformed values.
*       
*  Inputs
*  =======
*
*  mat      - matrix of size MxK to transform, if mat is a temporary and unique
*             then the transform is done in place if possible
*  dim      - dimension along which transormation must be conducted, if dim == 1,
*             then N = M, if dim == 2, then N = K.
*  scal     - scaling factor, default value = 1.0.
*  cont     - optional parameter, stores internal data. If on input cont is empty,
*             i.e. cont = fft_context(), or not initialized for problems of given type,
*             then cont is initialized, and can be reused in next computations of the
*             same type. If given context is used in computations of different type,
*             then data stored in cont are not destroyed and can be used later. A problem
*             is determined by length (number of rows if dim == 1, or number of columns 
*             if dim == 2), and type (dct1, dct_23, or fft). On return cont is a context
*             appended with data for problems of given type.
*       
*  Output
*  =======
*
*           Matrix containing result of the transform. The dct3 transform is the
*           inverse of dct2, i.e. (1/N) * dct3 ( dct2(x) ) = x.
*/
MATCL_FFT_EXPORT Matrix dct3(const matcl::Matrix& mat, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct3(matcl::Matrix&& x, Integer dim, Real scale);
MATCL_FFT_EXPORT Matrix dct3(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont);
MATCL_FFT_EXPORT Matrix dct3(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont);

};};
