#include "matcl-fft/matcl_dct.h"
#include "matcl-fft/matcl_fft.h"
#include "dct_helpers.h"
#include "matcl-fft/dct_kernels.h"

//#include "mmlib/details/raw_fwd.h"
#include "matcl-internals/container/sparse_ccs.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

#include "matcl-blas-lapack/level1/level1.h"
#include "matcl-core/utils/workspace.h"

#include "matcl-matrep\matrix\matrix_rep_dense.h"
#include "matcl-matrep\lib_functions\matrix_gen.h"

namespace matcl { namespace fft
{

using namespace matcl;

void dct3_1_Rows_dim1(const Real* X, Integer X_ld, Integer rows, Integer cols, 
               Real* Y, Integer Y_ld, Real scale, Real* work, Integer work_size,
               fft_context& cont)
{    
    check_work(work_size, rows);

    Real mult           = scale * 0.5;
    cont.initialize(rows, fft_context::dct_23_1);
    const Real* cs      = cont.get_table(rows, fft_context::dct_23_1);

    make_pack_dct3_dim1::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
    matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 1, mult, cont);    
    make_permute_dct3_dim1::eval(Y, Y_ld, rows, cols, work);
}
void dct3_1_Rows_inpl_dim1(Real* Y, Integer Y_ld, Integer rows, Integer cols, 
                         Real scale, Real* work, Integer work_size, fft_context& cont)
{    
    check_work(work_size, rows);

    Real mult           = scale * 0.5;
    cont.initialize(rows, fft_context::dct_23_1);
    const Real* cs      = cont.get_table(rows, fft_context::dct_23_1);

    make_pack_dct3_inpl_dim1::eval(Y, Y_ld, rows, cols, cs, work);
    matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 1, mult, cont);    
    make_permute_dct3_dim1::eval(Y, Y_ld, rows, cols, work);
}
void dct3_2_Rows_dim1(const Real* X, Integer X_ld, Integer rows, Integer cols, 
               Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{    
    Real mult           = scale * 0.25;
    cont.initialize(rows, fft_context::dct_23_2);
    const Real* cs      = cont.get_table(rows, fft_context::dct_23_2);    

    prepare_dct3_dim1::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
    matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 1, mult, cont);    
    recurrence_dct3_dim1::eval(Y, Y_ld, rows, cols, cs);
}
void dct3_2_Rows_inpl_dim1(Real* Y, Integer Y_ld, Integer rows, Integer cols, 
                         Real scale, fft_context& cont)
{    
    Real mult           = scale * 0.25;
    cont.initialize(rows, fft_context::dct_23_2);
    const Real* cs      = cont.get_table(rows, fft_context::dct_23_2);

    prepare_dct3_inpl_dim1::eval(Y, Y_ld, rows, cols, cs);
    matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 1, mult, cont);    
    recurrence_dct3_dim1::eval(Y, Y_ld, rows, cols, cs);
}
template<Integer Rows>
void dct3_1_Rows_dim2(const Real* X, Integer X_ld, Integer rows0, Integer cols, 
               Real* Y, Integer Y_ld, Real scale, Real* work, Integer work_size,
               fft_context& cont)
{    
    Integer rows        = get_int<Rows>::eval(rows0);    

    check_work(work_size, cols);

    Real mult           = scale * 0.5;
    cont.initialize(cols, fft_context::dct_23_1);
    const Real* cs      = cont.get_table(cols, fft_context::dct_23_1);

    if (Y_ld == Rows)
    {
        if (X_ld == Rows)
            make_pack_dct3_dim2<Rows,Rows,Rows>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
        else
            make_pack_dct3_dim2<Rows,Rows,0>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);

        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        make_permute_dct3_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, work);
    }
    else
    {
        if (X_ld == Rows)
            make_pack_dct3_dim2<Rows,0,Rows>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
        else
            make_pack_dct3_dim2<Rows,0,0>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);

        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        make_permute_dct3_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, work);
    };
}
template<Integer Rows>
void dct3_1_Rows_inpl_dim2(Real* Y, Integer Y_ld, Integer rows0, Integer cols, 
                         Real scale, Real* work, Integer work_size, fft_context& cont)
{    
    Integer rows        = get_int<Rows>::eval(rows0);

    check_work(work_size, cols);

    Real mult           = scale * 0.5;
    cont.initialize(cols, fft_context::dct_23_1);
    const Real* cs      = cont.get_table(cols, fft_context::dct_23_1);

    if (Y_ld == Rows)
    {
        make_pack_dct3_inpl_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, cs, work);
        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        make_permute_dct3_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, work);
    }
    else
    {
        make_pack_dct3_inpl_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, cs, work);
        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        make_permute_dct3_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, work);
    };
}

template<Integer Rows>
void dct3_2_Rows_dim2(const Real* X, Integer X_ld, Integer rows0, Integer cols, 
               Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{    
    Integer rows        = get_int<Rows>::eval(rows0);

    Real mult           = scale * 0.25;
    cont.initialize(cols, fft_context::dct_23_2);
    const Real* cs      = cont.get_table(cols, fft_context::dct_23_2);

    if (Y_ld == Rows)
    {
        if (X_ld == Rows)
            prepare_dct3_dim2<Rows,Rows,Rows>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
        else
            prepare_dct3_dim2<Rows,Rows,0>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);

        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        recurrence_dct3_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, cs);
    }
    else
    {
        if (X_ld == Rows)
            prepare_dct3_dim2<Rows,0,Rows>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);
        else
            prepare_dct3_dim2<Rows,0,0>::eval(Y, Y_ld, X, X_ld, rows, cols, cs);

        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        recurrence_dct3_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, cs);
    };
}
template<Integer Rows>
void dct3_2_Rows_inpl_dim2(Real* Y, Integer Y_ld, Integer rows0, Integer cols, 
                         Real scale, fft_context& cont)
{    
    Integer rows        = get_int<Rows>::eval(rows0);

    Real mult           = scale * 0.25;
    cont.initialize(cols, fft_context::dct_23_2);
    const Real* cs      = cont.get_table(cols, fft_context::dct_23_2);

    if (Y_ld == Rows)
    {
        prepare_dct3_inpl_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, cs);
        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        recurrence_dct3_dim2<Rows,Rows>::eval(Y, Y_ld, rows, cols, cs);
    }
    else
    {
        prepare_dct3_inpl_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, cs);
        matcl::fft::ifft_pack_inplace(Y, Y_ld, rows, cols, 2, mult, cont);    
        recurrence_dct3_dim2<Rows,0>::eval(Y, Y_ld, rows, cols, cs);
    };
}

void fft::dct3(const Real* X, Integer X_ld, Integer rows, Integer cols, Integer dim,
               Real* Y, Integer Y_ld, Real scale, fft_context& cont)
{    
    if (rows == 0 || cols == 0)
        return;

    if (dim == 2)
    {
        switch(cols)
        {
            case 1:
            {
                scale           = scale * 0.5;
                if (scale != 1.0)
                    matcl::level1::ax<Real,Real,Real,0>::eval(Y, X, rows, scale);
                else
                    matcl::level1::copy_mat<true,Real,Real,0,1,1>::eval(Y, Y_ld, X, X_ld, rows, cols);
                return;
            }
            case 2:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<2>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<2>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<2>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<2>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 3:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<3>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<3>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<3>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<3>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 4:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<4>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<4>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<4>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<4>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 5:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<5>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<5>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<5>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<5>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 6:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<6>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<6>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<6>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<6>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 7:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<7>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<7>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<7>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<7>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            case 8:
            {
                if (X_ld == 1 && Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel<8>(X, Y, rows, 1, 1, scale);
                    else
                        return dct3_kernel<8>(X, Y, rows, 1, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride<8>(X, Y, rows, 1, 1, scale, X_ld, Y_ld);
                    else
                        return dct3_kernel_stride<8>(X, Y, rows, 1, 1, X_ld, Y_ld);
                };
            }
            default:
                ;
        };

        switch(rows)
        {
            case 1:
                return dct3_2_Rows_dim2<1>(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
            case 2:
                return dct3_2_Rows_dim2<2>(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
            case 3:
                return dct3_2_Rows_dim2<3>(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
            case 4:
                return dct3_2_Rows_dim2<4>(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
            default:
                return dct3_2_Rows_dim2<0>(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
        };
    }
    else
    {
        switch(rows)
        {
            case 1:
            {
                scale           = scale * 0.5;
                if (scale != 1.0)
                    matcl::level1::ax<Real,Real,Real,0>::eval(Y, Y_ld, X, X_ld, cols, scale);
                else
                    matcl::level1::copy_mat<true,Real,Real,1,0,0>::eval(Y, Y_ld, X, X_ld, rows, cols);

                return;
            }
            case 2:
            {
                if (scale != 1.0)
                    return dct3_kernel<2>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<2>(X,Y, cols, X_ld, Y_ld);                
            }
            case 3:
            {
                if (scale != 1.0)
                    return dct3_kernel<3>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<3>(X,Y, cols, X_ld, Y_ld);                
            }
            case 4:
            {
                if (scale != 1.0)
                    return dct3_kernel<4>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<4>(X,Y, cols, X_ld, Y_ld);                
            }
            case 5:
            {
                if (scale != 1.0)
                    return dct3_kernel<5>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<5>(X,Y, cols, X_ld, Y_ld);                
            }
            case 6:
            {
                if (scale != 1.0)
                    return dct3_kernel<6>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<6>(X,Y, cols, X_ld, Y_ld);                
            }
            case 7:
            {
                if (scale != 1.0)
                    return dct3_kernel<7>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<7>(X,Y, cols, X_ld, Y_ld);                
            }
            case 8:
            {
                if (scale != 1.0)
                    return dct3_kernel<8>(X,Y, cols, X_ld, Y_ld, scale);
                else
                    return dct3_kernel<8>(X,Y, cols, X_ld, Y_ld);                
            }
            default:
            {
                return dct3_2_Rows_dim1(X, X_ld, rows, cols, Y, Y_ld, scale, cont);
            }
        };                
    };
}
void fft::dct3_inplace(Real* Y, Integer Y_ld, Integer rows, Integer cols, Integer dim,
                        Real scale, fft_context& cont)
{    
    if (rows == 0 || cols == 0)
        return;

    if (dim == 2)
    {
        switch(cols)
        {
            case 1:
            {
                scale           = scale * 0.5;
                if (scale != 1.0)
                    matcl::level1::ay<Real,Real,0>::eval(Y, rows, scale);

                return;
            }
            case 2:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<2>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<2>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<2>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<2>(Y, rows, 1, Y_ld);
                };
            }
            case 3:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<3>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<3>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<3>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<3>(Y, rows, 1, Y_ld);
                };
            }
            case 4:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<4>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<4>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<4>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<4>(Y, rows, 1, Y_ld);
                };
            }
            case 5:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<5>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<5>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<5>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<5>(Y, rows, 1, Y_ld);
                };
            }
            case 6:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<6>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<6>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<6>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<6>(Y, rows, 1, Y_ld);
                };
            }
            case 7:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<7>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<7>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<7>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<7>(Y, rows, 1, Y_ld);
                };
            }
            case 8:
            {
                if (Y_ld == 1)
                {
                    if (scale != 1.0)
                        return dct3_kernel_inpl<8>(Y, rows, 1, scale);
                    else
                        return dct3_kernel_inpl<8>(Y, rows, 1);
                }
                else
                {
                    if (scale != 1.0)
                        return dct3_kernel_stride_inpl<8>(Y, rows, 1, scale, Y_ld);
                    else
                        return dct3_kernel_stride_inpl<8>(Y, rows, 1, Y_ld);
                };
            }
            default:
                ;
        };

        switch(rows)
        {
            case 1:
                return dct3_2_Rows_inpl_dim2<1>(Y, Y_ld, rows, cols, scale, cont);
            case 2:
                return dct3_2_Rows_inpl_dim2<2>(Y, Y_ld, rows, cols, scale, cont);
            case 3:
                return dct3_2_Rows_inpl_dim2<3>(Y, Y_ld, rows, cols, scale, cont);
            case 4:
                return dct3_2_Rows_inpl_dim2<4>(Y, Y_ld, rows, cols, scale, cont);
            default:
                return dct3_2_Rows_inpl_dim2<0>(Y, Y_ld, rows, cols, scale, cont);
        };
    }
    else
    {
        switch(rows)
        {
            case 1:
            {
                scale           = scale * 0.5;
                if (scale != 1.0)
                    matcl::level1::ay<Real,Real, 0>::eval(Y, Y_ld, cols, scale);

                return;
            }
            case 2:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<2>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<2>(Y, cols, Y_ld);                
            }
            case 3:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<3>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<3>(Y, cols, Y_ld);                
            }
            case 4:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<4>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<4>(Y, cols, Y_ld);                
            }
            case 5:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<5>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<5>(Y, cols, Y_ld);                   
            }
            case 6:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<6>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<6>(Y, cols, Y_ld);                   
            }
            case 7:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<7>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<7>(Y, cols, Y_ld);                   
            }
            case 8:
            {
                if (scale != 1.0)
                    return dct3_kernel_inpl<8>(Y, cols, Y_ld, scale);
                else
                    return dct3_kernel_inpl<8>(Y, cols, Y_ld);                   
            }
            default:
            {
                return dct3_2_Rows_inpl_dim1(Y, Y_ld, rows, cols, scale, cont);
            }
        };                        
    };
}
Matrix fft::dct3(const Matrix& x, Integer dim, Real scale)
{
    fft_context cont;
    return dct3(x,dim,scale,cont);
};

Matrix fft::dct3(Matrix&& x, Integer dim, Real scale)
{
    fft_context cont;
    return dct3(std::move(x),dim,scale,cont);
};

Matrix fft::dct3(const matcl::Matrix& x, Integer dim, Real scale, fft_context& cont)
{
    Integer rows           = x.rows();
    Integer cols           = x.cols();

    if (rows == 0 || cols == 0)
        return x;

    using const_rep     = matcl::dense_matrix<Real, true>;
    const_rep rep_X     = const_rep(x);
    const Real* X       = rep_X.ptr();
    Integer X_ld        = rep_X.ld();

    Real* Y;
    Matrix mY           = matcl::make_real_dense_noinit(rows,cols, Y);

    dct3(X, X_ld, rows, cols, dim, Y, rows, scale, cont);

    return mY;
};
Matrix fft::dct3(matcl::Matrix&& x, Integer dim, Real scale, fft_context& cont)
{
    if (x.is_unique() == false)
        return dct3(x,dim,scale,cont);

    Integer rows           = x.rows();
    Integer cols           = x.cols();

    if (rows == 0 || cols == 0)
        return x;

    using const_rep     = matcl::dense_matrix<Real, true>;
    const_rep rep_X     = const_rep(x);
    const Real* X       = rep_X.ptr();
    Integer X_ld        = rep_X.ld();

    //x is unique and can be destroyed
    Real* Y             = const_cast<Real*>(X);

    dct3_inplace(Y, X_ld, rows, cols, dim, scale, cont);

    return x;
};

}}
