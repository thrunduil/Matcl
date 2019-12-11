#pragma once

#include "matcl-fft/matcl_fft_config.h"

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-fft/matcl_fft_exception.h"
#include "matcl-fft/fft_context.h"

#include "matcl-matrep/matrix/matrix.h"

namespace matcl { namespace fft
{

//-----------------------------------------------------------------------
//                          low level fft
//-----------------------------------------------------------------------
MATCL_FFT_EXPORT void
fft(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim, 
    Complex* Y, Integer Y_ld, Real scale, bool packed, fft_context& cont);

MATCL_FFT_EXPORT void
fft_pack(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim, 
    Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
fft_pack_inplace(Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
fft(const Complex* X, Integer X_ld, Integer rows, Integer cols, Integer dim, 
    Complex* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
fft_inplace(Complex* X, Integer X_ld, Integer rows, Integer cols,  Integer dim,
    Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
ifft(const Complex* X, Integer X_ld, Integer rows, Integer cols, Integer dim, 
     Complex* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
ifft_inplace(Complex* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
ifft(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim, 
     Complex* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
ifft_conj_even(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void 
ifft_conj_even(const Complex* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void 
ifft_pack(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real* Y, Integer Y_ld, Real scale, fft_context& cont);

MATCL_FFT_EXPORT void
ifft_pack_inplace(Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
    Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          high level fft
//-----------------------------------------------------------------------
MATCL_FFT_EXPORT Matrix fft(const matcl::Matrix& mat, Integer dim, Real scale = 1.0);    
MATCL_FFT_EXPORT Matrix fft(matcl::Matrix&& mat, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix fft(const matcl::Matrix& mat, Integer dim, Real scale, fft_context& cont);    
MATCL_FFT_EXPORT Matrix   fft(matcl::Matrix&& mat, Integer dim, Real scale, fft_context& cont);

MATCL_FFT_EXPORT Matrix ifft(const matcl::Matrix& mat, Integer dim, Real scale = 1.0,
                            bool conjugate_even = false);
MATCL_FFT_EXPORT Matrix ifft(matcl::Matrix&& mat, Integer dim, Real scale = 1.0,
                            bool conjugate_even = false);    
MATCL_FFT_EXPORT Matrix ifft(const matcl::Matrix& mat, Integer dim, Real scale,
                            bool conjugate_even, fft_context& cont);
MATCL_FFT_EXPORT Matrix ifft(matcl::Matrix&& mat, Integer dim, Real scale,
                            bool conjugate_even, fft_context& cont);    

//only for real matrices
MATCL_FFT_EXPORT Matrix fft_pack(const matcl::Matrix& mat, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix fft_pack(matcl::Matrix&& mat, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix fft_pack(const matcl::Matrix& mat, Integer dim, Real scale,
                            fft_context& cont);
MATCL_FFT_EXPORT Matrix fft_pack(matcl::Matrix&& mat, Integer dim, Real scale,
                            fft_context& cont);

MATCL_FFT_EXPORT Matrix ifft_pack(const matcl::Matrix& mat, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix ifft_pack(matcl::Matrix&& mat, Integer dim, Real scale = 1.0);
MATCL_FFT_EXPORT Matrix ifft_pack(const matcl::Matrix& mat, Integer dim, Real scale,
                            fft_context& cont);
MATCL_FFT_EXPORT Matrix ifft_pack(matcl::Matrix&& mat, Integer dim, Real scale,
                            fft_context& cont);

//-----------------------------------------------------------------------
//                          fft
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete fast fourier transform in the form:
*
*       X_k = scal * sum_{j=1}^{K-1} x_j * exp(-2 pi i * k * j / K)
*
*  for k = 0 .. K-1, {x_k} : vector of orginal values, {X_k} : vector of
*  transformed values.
*       
*  Inputs
*  =======
*
*  mat      - matrix of size MxN to transform, if mat is a temporary and unique
*             then the transform is done in place if possible
*  dim      - dimension along which transormation must be conducted. If dim = 1 
*             then K = M, if dim == 2 then K = N.
*  scal     - scaling factor, default value = 1.0.
*  cont     - optional parameter, stores internal data. If on input cont is empty,
*             i.e. cont = fft_context(), or not initialized for problems of given type,
*             then cont is initialized, and can be reused in next computations of the
*             same type. If given context is used in computations of different type,
*             then data stored in cont are not destroyed and can be used later. A problem
*             is determined by length (number of rows if dim == 1, or number of columns 
*             if dim == 2), and type (fft_real, fft_complex, or dct1). On return cont 
*             is a context appended with data for problems of given type.
*       
*  Output
*  =======
*
*           Matrix containing result of the transform. If imaginary part of mat is
*           zero, then resulting matrix A is conjugate even, i.e. for some x and 
*           for v = 0 or v = 1: A = [x; conj( flipud( x(2:end-v,:) ) )] if dim == 1
*           or transposition of such matrix if dim == 2.
*/
MATCL_FFT_EXPORT Matrix fft(const matcl::Matrix& mat, Integer dim, Real scal);    
MATCL_FFT_EXPORT Matrix fft(matcl::Matrix&& mat, Integer dim, Real scal);
MATCL_FFT_EXPORT Matrix fft(const matcl::Matrix& mat, Integer dim, Real scale, fft_context& cont);    
MATCL_FFT_EXPORT Matrix fft(matcl::Matrix&& mat, Integer dim, Real scale, fft_context& cont);

//-----------------------------------------------------------------------
//                          ifft
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*
*  calculate discrete inverse fast fourier transform in the form:
*
*       x_k = scal * sum_{j=1}^{K-1} X_j * exp( 2 pi i * k * j / K)
*
*  for k = 0 .. K-1, {X_k} : vector of orginal values, {x_k} : vector of
*  transformed values.
*       
*  Inputs
*  =======
*
*  mat      - matrix of size MxN to transform, if mat is a temporary and unique
*             then the transform is done in place if possible
*  dim      - dimension along which transormation must be conducted. If dim = 1 
*             then K = M, if dim == 2 then K = N.
*  scal     - scaling factor, if scal = 1/M if dim == 1 or scal = 1/N if dim = 2
*             then this transform is the inverse of fft, i.e. x = ifft(fft(x)).
*             Such scaling is required to obtain results consistent with Matlab's
*             ifft. Default value = 1.0.
*  conj_even- if true then it is assumed (but not checked), that mat is conjugate
*             even along dimension d, i.e. mat = [x; conj( flipud( x(2:end-v,:) ) )]
*             if dim == 1 or mat = [x, conj( fliplr( x(:,2:end-v) ) )] if dim == 2
*             for some matrix x and v = 0 or 1. One should obtaine the same results
*             for conj_even == true and false, but version with conj_even == true 
*             is faster. Additionally only elements in x part of the matrix mat are
*             accessed. Default value = false.
*  cont     - optional parameter, stores internal data. If on input cont is empty,
*             i.e. cont = fft_context(), or not initialized for problems of given type,
*             then cont is initialized, and can be reused in next computations of the
*             same type. If given context is used in computations of different type,
*             then data stored in cont are not destroyed and can be used later. A problem
*             is determined by length (number of rows if dim == 1, or number of columns 
*             if dim == 2), and type (fft_real, fft_complex, or dct1). On return cont 
*             is a context appended with data for problems of given type.
*       
*  Output
*  =======
*
*           Matrix containing result of the transform. If mat is conjugate even,
*           then resulting matrix is real.
*/
MATCL_FFT_EXPORT Matrix ifft(const matcl::Matrix& mat, Integer dim, Real scale, bool conjugate_even);
MATCL_FFT_EXPORT Matrix ifft(matcl::Matrix&& mat, Integer dim, Real scale, bool conj_even);    
MATCL_FFT_EXPORT Matrix ifft(const matcl::Matrix& mat, Integer dim, Real scale, bool conjugate_even, 
                            fft_context& cont);
MATCL_FFT_EXPORT Matrix ifft(matcl::Matrix&& mat, Integer dim, Real scale, bool conjugate_even, 
                            fft_context& cont);    

};};

