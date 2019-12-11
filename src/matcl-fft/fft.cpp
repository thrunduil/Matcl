#include "matcl-fft/matcl_fft.h"
#include "matcl-matrep/details/extract_type_switch.h" 

//#include "matcl-core/details/raw_fwd.h"
#include "matcl-internals/container/sparse_ccs.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

#include "mkl_wrapper/fft_wrapper.h"
#include "mkl_wrapper/fft_utils.h"

#include "matcl-blas-lapack/level1/level1.h"

namespace matcl { namespace fft
{

using details::mkl_fft::m2fft;

//--------------------------------------------------------------------------------
//                                  HELPERS
//--------------------------------------------------------------------------------
static void throw_error_fft(long status)
{
    std::ostringstream msg;

    msg <<"mkl fft error: ";
    details::mkl_fft::fft_wrapper().error_string(status, msg);

    throw matcl::error::matcl_fft_exception(msg.str());
};

//--------------------------------------------------------------------------------
//                                  LOW LEVEL
//--------------------------------------------------------------------------------
void fft::fft(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
              Complex* Y, Integer Y_ld, Real scale, bool packed, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_real(X, X_ld, m2fft(Y), Y_ld, rows, cols, dim, scale, cont);

    if (status != 0)
        throw_error_fft(status);

    if (packed == true)
        return;

    // retrieve remaining elements of conjugate even transform
    // for real inputs mkl does not return conjugate even part.
    if (dim == 1)
    {
        Integer max_r       = rows / 2 + 1;

        for (Integer i = 0; i < cols; i++)
        {
            for (Integer j = rows-1, k = 1; j >= max_r; --j, ++k)
                Y[j]        = conj(Y[k]);

            Y               += Y_ld;
        };
    }
    else
    {
        Integer max_c       = cols / 2 + 1;
        Complex* ptr_dest   = Y + (cols - 1) * Y_ld;
        Complex* ptr_source = Y + 1 * Y_ld;

        for (Integer j = cols-1; j >= max_c; --j)
        {
            for (Integer i = 0; i < rows; i++)
                ptr_dest[i] = conj(ptr_source[i]);

            ptr_dest        -= Y_ld;
            ptr_source      += Y_ld;
        };
    }
};
void fft::fft(const Complex* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Complex* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_comp(m2fft(X), X_ld, m2fft(Y), Y_ld, rows, cols, 
                                    dim, scale, false, false, cont);

    if (status != 0)
        throw_error_fft(status);
};

void fft::fft_inplace(Complex* X, Integer X_ld, Integer rows, Integer cols, 
                    Integer dim, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_comp(m2fft(X), X_ld, m2fft(X), X_ld, rows, cols, 
                                    dim, scale, false, true, cont);

    if (status != 0)
        throw_error_fft(status);
};

void fft::ifft(const Complex* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Complex* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_comp(m2fft(X), X_ld, m2fft(Y), Y_ld, rows, cols, 
                                    dim, scale, true, false, cont);

    if (status != 0)
        throw_error_fft(status);
}

void fft::ifft_inplace(Complex* X, Integer X_ld, Integer rows, Integer cols, 
                    Integer dim, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_comp(m2fft(X), X_ld, m2fft(X), X_ld, rows, cols, 
                                    dim, scale, true, true, cont);

    if (status != 0)
        throw_error_fft(status);
}

void fft::ifft_conj_even(const Complex* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_inv_conj_even(m2fft(X), X_ld, m2fft(Y), Y_ld, rows, 
                                             cols, dim, scale, cont);

    if (status != 0)
        throw_error_fft(status);
};

void fft::ifft(const Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Complex* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    matcl::level1::copy_mat<true, Complex, Real, 0, 0, 0>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols);

    details::mkl_fft::fft_wrapper wrapper;

    long status = wrapper.eval_comp(m2fft(Y), Y_ld, m2fft(Y), Y_ld, rows, cols, dim,
                                    scale, true, true, cont);

    if (status != 0)
        throw_error_fft(status);
};

void fft::ifft_conj_even(const Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    Real* Y_s = Y;

    if (dim == 1)
    {
        //create packed storage
        for (Integer i = 0; i < cols; ++i)
        {
            Y[0]        = X[0];
            for (Integer j = 1, k = 1; j < rows; j += 2, ++k)
            {  
                Y[j]    = X[k];
                Y[j+1]  = 0.0;
            };

            Y           += Y_ld;
            X           += X_ld;
        };
    }
    else
    {
        //create packed storage
        for (Integer j = 0; j < rows; ++j)
        {  
            Y[j]        = X[j];
        };

        Y               += Y_ld;
        X               += X_ld;

        for (Integer i = 1; i < cols; i += 2)
        {
            for (Integer j = 0; j < rows; ++j)
            {  
                Y[j]    = X[j];
            };
            Y           += Y_ld;

            for (Integer j = 0; j < rows; ++j)
            {  
                Y[j]    = 0;
            };

            X           += X_ld;
            Y           += Y_ld;
        };
    }

    details::mkl_fft::fft_wrapper wrapper;

    long status = wrapper.eval_real_inv_conj_even(m2fft(Y_s), Y_ld, rows, cols, dim, scale, cont); 

    if (status != 0)
        throw_error_fft(status);

};

void fft::fft_pack(const Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_real_pack(X, X_ld, Y, Y_ld, rows, cols, dim, scale, 
                                         false,false, cont);

    if (status != 0)
    {
        throw_error_fft(status);
    }
};
void fft::fft_pack_inplace(Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_real_pack(X, X_ld, X, X_ld, rows, cols, dim, scale, 
                                         false,true, cont);

    if (status != 0)
    {
        throw_error_fft(status);
    }
}
void fft::ifft_pack(const Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_real_pack(X, X_ld, Y, Y_ld, rows, cols, dim, scale, 
                                         true,false, cont);

    if (status != 0)
    {
        throw_error_fft(status);
    }
}
void fft::ifft_pack_inplace(Real* X, Integer X_ld, Integer rows, Integer cols,
                    Integer dim, Real scale, fft_context& cont)
{
    details::mkl_fft::fft_wrapper wrapper;
    long status = wrapper.eval_real_pack(X, X_ld, X, X_ld, rows, cols, dim, scale, 
                                         true,true, cont);

    if (status != 0)
    {
        throw_error_fft(status);
    }
};


namespace details
{

static void check_dim(Integer dim)
{
    matcl::error::check_dim(dim);
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V>
struct fft_val{};

template<>
struct fft_val<Complex>
{
    using CMat          = raw::Matrix<Complex, struct_dense>;
    using RMat          = raw::Matrix<Real, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<CMat>;
    
    static Matrix eval(RValue_type&& A0, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {        
        if (conjugate_even == true)
        {
            CMat& A = A0.get();
            RMat res(ti::ti_real(), A.rows(), A.cols());  

            fft::ifft_conj_even(A.ptr(), A.ld(), A.rows(), A.cols(), dim, 
                                res.ptr(), res.ld(), scale, cont);

            return Matrix(std::move(res),false);
        }
        else
        {
            CMat& A = A0.get();

            if (inverse == false)
                fft::fft_inplace(A.ptr(), A.ld(), A.rows(), A.cols(), dim, scale, cont);
            else
                fft::ifft_inplace(A.ptr(), A.ld(), A.rows(), A.cols(), dim, scale, cont);

            return Matrix(std::move(A),false);
        }
    };

    static Matrix eval(const CMat& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {        
        if (conjugate_even == true)
        {
            RMat res(ti::ti_real(), A.rows(), A.cols());  

            fft::ifft_conj_even(A.ptr(), A.ld(), A.rows(), A.cols(), dim, 
                                res.ptr(), res.ld(), scale, cont);

            return Matrix(std::move(res),false);
        }
        else
        {
            CMat res(ti::ti_compl(), A.rows(), A.cols());  

            if (inverse == false)
                fft::fft(A.ptr(), A.ld(), A.rows(), A.cols(), dim, res.ptr(), res.ld(), scale, cont);
            else
                fft::ifft(A.ptr(), A.ld(), A.rows(), A.cols(), dim, res.ptr(), res.ld(), scale, cont);

            return Matrix(std::move(res),false);
        }
    }
};

template<>
struct fft_val<Real>
{
    using CMat          = raw::Matrix<Complex, struct_dense>;
    using Mat           = raw::Matrix<Real, struct_dense>;
    using RValue_Mat    = matcl::details::rvalue_holder<Mat>;
    
    static void eval_inverse_ce(const Mat& A, Integer dim, Real scale, Matrix& ret,
                                fft_context& cont)
    {
        Integer r       = A.rows();
        Integer c       = A.cols();

        Mat res(ti::ti_real(), r, c);  

        fft::ifft_conj_even(A.ptr(), A.ld(), r, c, dim, 
                            res.ptr(), res.ld(), scale, cont);

        ret = Matrix(std::move(res),false);
        return;
    };

    static void eval_inverse_nce(const Mat& A, Integer dim, Real scale, Matrix& ret,
                                 fft_context& cont)
    {
        Integer r       = A.rows();
        Integer c       = A.cols();

        CMat res(ti::ti_compl(), r, c); 
        fft::ifft(A.ptr(), A.ld(), A.rows(), A.cols(), dim, res.ptr(), res.ld(), scale, cont);

        ret = Matrix(std::move(res),false);
        return;
    };

    static Matrix eval(RValue_Mat&& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {
        return eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    }
    static Matrix eval(const Mat& A, bool inverse, bool conjugate_even, Integer dim, Real scale,
                       fft_context& cont)
    {
        Matrix ret;

        if (inverse == true)
        {            
            if (conjugate_even == true)
                eval_inverse_ce(A, dim, scale, ret, cont);
            else
                eval_inverse_nce(A, dim, scale, ret, cont);

            return ret;
        }

        Integer r       = A.rows();
        Integer c       = A.cols();

        CMat res(ti::ti_compl(), r, c);
        fft::fft(A.ptr(), A.ld(), r, c, dim, res.ptr(), res.ld(), scale, false, cont);

        ret = Matrix(std::move(res),false);
        return ret;
    }

    static Matrix eval_pack(RValue_Mat&& A0, Integer dim, Real scale, fft_context& cont)
    {
        Mat& A = A0.get();

        fft::fft_pack_inplace(A.ptr(), A.ld(), A.rows(), A.cols(), dim, scale, cont);
        return Matrix(std::move(A),false);
    }

    static Matrix eval_pack(const Mat& A, Integer dim, Real scale, fft_context& cont)
    {
        Matrix ret;

        Integer r       = A.rows();
        Integer c       = A.cols();

        Mat res(ti::ti_compl(), r, c);
        fft::fft_pack(A.ptr(), A.ld(), r, c, dim, res.ptr(), res.ld(), scale, cont);

        ret = Matrix(std::move(res),false);
        return ret;
    }

    static Matrix eval_pack_inverse(RValue_Mat&& A0, Integer dim, Real scale, fft_context& cont)
    {
        Mat& A = A0.get();

        fft::ifft_pack_inplace(A.ptr(), A.ld(), A.rows(), A.cols(), dim, scale, cont);
        return Matrix(std::move(A),false);
    }

    static Matrix eval_pack_inverse(const Mat& A, Integer dim, Real scale, fft_context& cont)
    {
        Matrix ret;

        Integer r       = A.rows();
        Integer c       = A.cols();

        Mat res(ti::ti_compl(), r, c);
        fft::ifft_pack(A.ptr(), A.ld(), r, c, dim, res.ptr(), res.ld(), scale, cont);

        ret = Matrix(std::move(res),false);
        return ret;
    }
};

template<>
struct fft_val<Object>
{
    using M             = raw::Matrix<Object, struct_dense>;
    using RValue_Mat    = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_Mat&&, const bool, const bool, Integer, Real, fft_context&)
    {
        throw matcl::error::object_value_type_not_allowed("fft");
    };
    static Matrix eval(const M&, const bool, const bool, Integer, Real, fft_context&)
    {
        throw matcl::error::object_value_type_not_allowed("fft");
    }
};

//--------------------------------------------------------------------------------
//                                  STRUCT TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct fft_impl
{
    //TODO
    using M             = raw::Matrix<V, S>;
    using VC            = typename matcl::details::unify_types<V, Real>::type;
    using DM            = raw::Matrix<VC, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        return fft_impl::eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    };

    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {
        DM Ac = raw::converter<DM,M>::eval(A);
        return fft_val<VC>::eval(std::move(Ac), inverse, conjugate_even, dim, scale, cont);
    }
};

template<class S>
struct fft_impl<Integer,S>
{
    using M             = raw::Matrix<Integer, S>;
    using DM            = raw::Matrix<Real, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        return fft_impl::eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    };

    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {
        DM Ac = raw::converter<DM,M>::eval(A);
        return fft_val<Real>::eval(std::move(Ac), inverse, conjugate_even, dim, scale, cont);
    }
};

template<>
struct fft_impl<Integer,struct_dense>
{
    using M             = raw::Matrix<Integer, struct_dense>;
    using DM            = raw::Matrix<Real, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        return fft_impl::eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    };
    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        DM Ac = raw::converter<DM,M>::eval(A);
        return fft_val<Real>::eval(std::move(Ac), inverse, conjugate_even, dim, scale, cont);
    }
};

//TODO
template<>
struct fft_impl<Float,struct_dense>
{
    using M             = raw::Matrix<Float, struct_dense>;
    using DM            = raw::Matrix<Real, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        return fft_impl::eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    };
    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        DM Ac = raw::converter<DM,M>::eval(A);
        return fft_val<Real>::eval(std::move(Ac), inverse, conjugate_even, dim, scale, cont);
    }
};

//TODO
template<>
struct fft_impl<Float_complex,struct_dense>
{
    using M             = raw::Matrix<Float_complex, struct_dense>;
    using DM            = raw::Matrix<Complex, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        return fft_impl::eval(A.get(), inverse, conjugate_even, dim, scale, cont);
    };
    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, Integer dim,
                       Real scale, fft_context& cont)
    {
        DM Ac = raw::converter<DM,M>::eval(A);
        return fft_val<Complex>::eval(std::move(Ac), inverse, conjugate_even, dim, scale, cont);
    }
};

template<class V>
struct fft_impl<V,struct_dense>
{
    using M             = raw::Matrix<V, struct_dense>;
    using DM            = raw::Matrix<V, struct_dense>;
    using RValue_type   = matcl::details::rvalue_holder<M>;

    static Matrix eval(RValue_type&& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {
        return fft_val<V>::eval(std::move(A), inverse, conjugate_even, dim, scale, cont);
    }
    static Matrix eval(const M& A, const bool inverse, const bool conjugate_even, 
                       Integer dim, Real scale, fft_context& cont)
    {
        return fft_val<V>::eval(A, inverse, conjugate_even, dim, scale, cont);
    }
};

//--------------------------------------------------------------------------------
//                                  VISITOR
//--------------------------------------------------------------------------------
struct unary_visitor_fft : public matcl::details::extract_type_switch<Matrix, unary_visitor_fft,true>
{
    template<class T>
    static Matrix eval(const Matrix&, const T& mat, const bool inverse, const bool conjugate_even,
                       Integer dim, Real scale, fft_context& cont)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::fft_impl<V,S>::eval(mat, inverse, conjugate_even, dim, scale, cont);
    }

    template<class T>
    static Matrix eval(T&& mat, const bool inverse, const bool conjugate_even,
                       Integer dim, Real scale, fft_context& cont)
    {
        using M = typename T::matrix_type;
        using V = typename M::value_type;
        using S = typename M::struct_type;

        return details::fft_impl<V,S>::eval(std::move(mat), inverse, conjugate_even, dim, scale, cont);
    }

    template<class T>
    static Matrix eval_scalar(const Matrix& handle, const T& val, bool, bool, Integer, Real scale, 
                              fft_context&)
    {
        if (scale == 1.0)
            return handle;
        else
            return val * scale;
    }

    template<class T>
    static Matrix eval_scalar(T&& val, const bool, const bool, Integer, Real scale, fft_context&)
    {
        if (scale == 1.0)
            return val;
        else
            return val * scale;
    }
};

};

MATCL_FFT_EXPORT Matrix fft::fft(const matcl::Matrix& mat, Integer dim, Real scale)
{
    fft_context cont;
    return fft(mat,dim,scale,cont);
}
MATCL_FFT_EXPORT Matrix fft::fft(const matcl::Matrix& mat, Integer dim, Real scale,
                           fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;    

    return details::unary_visitor_fft
        ::make<const Matrix&>(mat, false, false, dim, scale, cont);
}
MATCL_FFT_EXPORT Matrix fft::fft(matcl::Matrix&& mat, Integer dim, Real scale)
{
    fft_context cont;
    return fft(std::move(mat),dim,scale,cont);
}
MATCL_FFT_EXPORT Matrix fft::fft(matcl::Matrix&& mat, Integer dim, Real scale, fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return std::move(mat);    

    return details::unary_visitor_fft
        ::make<Matrix&&>(std::move(mat), false, false, dim, scale, cont);
};

MATCL_FFT_EXPORT Matrix fft::ifft(const Matrix& mat, Integer dim, Real scale, bool conjugate_even)
{
    fft_context cont;
    return ifft(mat,dim,scale,conjugate_even,cont);
};
MATCL_FFT_EXPORT Matrix fft::ifft(const Matrix& mat, Integer dim, Real scale, bool conjugate_even,
                            fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    return details::unary_visitor_fft
        ::make<const Matrix&>(mat, true, conjugate_even, dim, scale, cont);
}

MATCL_FFT_EXPORT Matrix fft::ifft(Matrix&& mat, Integer dim, Real scale, bool conjugate_even)
{
    fft_context cont;
    return ifft(std::move(mat),dim,scale,conjugate_even,cont);
};

MATCL_FFT_EXPORT Matrix fft::ifft(Matrix&& mat, Integer dim, Real scale, bool conjugate_even,
                            fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    return details::unary_visitor_fft
        ::make<Matrix&&>(std::move(mat), true, conjugate_even, dim, scale, cont);
}

Matrix fft::fft_pack(const matcl::Matrix& mat, Integer dim, Real scale)
{
    fft_context cont;
    return fft_pack(mat,dim,scale,cont);
};
Matrix fft::fft_pack(const matcl::Matrix& mat, Integer dim, Real scale, fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    using Mat   = raw::Matrix<Real, struct_dense>;

    if (mat.get_matrix_code() != mat_code::real_dense)
    {
        Matrix tmp = convert(mat,mat_code::real_dense);
        return fft_pack(std::move(tmp), dim, scale, cont);
    }
    else
    {    
        const Mat& rep = mat.get_impl<Mat>();
        return details::fft_val<Real>::eval_pack(rep, dim, scale, cont);
    };
}

Matrix fft::ifft_pack(const matcl::Matrix& mat, Integer dim, Real scale)
{
    fft_context cont;
    return ifft_pack(mat,dim,scale,cont);
};
Matrix fft::ifft_pack(const matcl::Matrix& mat, Integer dim, Real scale, fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    using Mat   = raw::Matrix<Real, struct_dense>;

    if (mat.get_matrix_code() != mat_code::real_dense)
    {
        Matrix tmp = convert(mat,mat_code::real_dense);
        return ifft_pack(std::move(tmp), dim, scale, cont);
    }
    else
    {    
        const Mat& rep = mat.get_impl<Mat>();
        return details::fft_val<Real>::eval_pack_inverse(rep, dim, scale, cont);
    };    
}

Matrix fft::fft_pack(matcl::Matrix&& mat, Integer dim, Real scale)
{
    fft_context cont;
    return fft_pack(std::move(mat),dim,scale,cont);
};
Matrix fft::fft_pack(matcl::Matrix&& mat, Integer dim, Real scale, fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    using Mat   = raw::Matrix<Real, struct_dense>;

    if (mat.get_matrix_code() != mat_code::real_dense)
    {
        Matrix tmp = convert(std::move(mat),mat_code::real_dense);
        return fft_pack(std::move(tmp), dim, scale, cont);
    }

    if (mat.is_unique() == true)
    {
        matcl::details::rvalue_holder<Mat> rep = mat.move_impl<Mat>();    
        return details::fft_val<Real>::eval_pack(std::move(rep), dim, scale, cont);
    }        
    else
    {
        const Mat& rep = mat.get_impl<Mat>();    
        return details::fft_val<Real>::eval_pack(rep, dim, scale, cont);
    };
}

Matrix fft::ifft_pack(matcl::Matrix&& mat, Integer dim, Real scale)
{
    fft_context cont;
    return ifft_pack(std::move(mat), dim, scale, cont);
};
Matrix fft::ifft_pack(matcl::Matrix&& mat, Integer dim, Real scale, fft_context& cont)
{
    details::check_dim(dim);

    if (mat.numel() == 0)
        return mat;

    using Mat   = raw::Matrix<Real, struct_dense>;

    if (mat.get_matrix_code() != mat_code::real_dense)
    {
        Matrix tmp = convert(std::move(mat),mat_code::real_dense);
        return ifft_pack(std::move(tmp), dim, scale, cont);
    }    

    if (mat.is_unique() == true)
    {
        matcl::details::rvalue_holder<Mat> rep = mat.move_impl<Mat>();    
        return details::fft_val<Real>::eval_pack_inverse(std::move(rep), dim, scale, cont);
    }        
    else
    {
        const Mat& rep = mat.get_impl<Mat>();    
        return details::fft_val<Real>::eval_pack_inverse(rep, dim, scale, cont);
    };
}

}}
