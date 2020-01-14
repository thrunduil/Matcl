/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "matcl-linalg/decompositions/eig/schur_utils.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/details/matrix.inl"
#include "matcl-matrep/matcl_matrep.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl
{

void details::schur_reorder_check(const Matrix& I, const Integer N, Integer& M)
{  
    if (I.is_vector() == false)
        throw error::error_schur_sel();

    if (I.length() != N)
        throw error::error_schur_sel();

    const Integer* ptr  = I.get_array<Integer>();
    M                   = 0;

    for (Integer i = 0; i < N; ++i)
    {
        Integer val = ptr[i];
        if (val != 1 && val != 0)
            throw error::error_schur_sel();

        if (val == 1)
            ++M;
    };
}

bool details::schur_is_trivial_reorder(const Matrix& I)
{
    Matrix Ic                   = convert(I, mat_code::integer_dense);
    const raw::integer_dense& I_r= Ic.get_impl<raw::integer_dense>();
    const Integer* I_r_ptr      = I_r.ptr();
    bool is_trivial             = true;
    int idx_trivial             = 0;
    Integer I_size              = I_r.size();

    for (idx_trivial = 0; idx_trivial < I_size; idx_trivial++)
    {
        if (I_r_ptr[idx_trivial] == 0) 
            break;
    }
    for (; idx_trivial < I_size; idx_trivial++)
    {
        if (I_r_ptr[idx_trivial] == 1)
        {
            is_trivial = false;
            break;
        }
    }
    return is_trivial;
}

//-----------------------------------------------------------------------------
//                  make_complex_eigenvectors
//-----------------------------------------------------------------------------
template<class V> 
Integer details::make_complex_eigenvectors<V,false>::calc_sel_size(const Mat& TA, const Integer* SELECT, Integer N)
{
    (void)TA;

    Integer M = 0;
    for (Integer j = 0; j < N; ++j)
    {
        if (SELECT[j] != 0)
            ++M;
    };

    return M;
};

template<class V> 
void details::make_complex_eigenvectors<V,false>::eval(const Mat& TA, const Matrix& ind, const Matrix& mat_VL,
                    const Matrix& mat_VR, Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
{
    (void)ind;
    (void)TA;
    (void)comp_left;
    (void)comp_right;

    //nothing to do
    XL = mat_VL;
    XR = mat_VR;

    return;
};

template<class V> 
void details::make_complex_eigenvectors<V,false>::complex_vectors_to_real(const Mat& TA, const Matrix& WL,
                const Matrix& WR, Mat& WL_R, Mat& WR_R, const Matrix& ind)
{
    (void)TA;
    (void)ind;

    //nothing to do

    Matrix WLc  = convert(WL, Mat::matrix_code);
    Matrix WRc  = convert(WR, Mat::matrix_code);

    WL_R.assign_to_fresh(WLc.get_impl<Mat>().make_unique());
    WR_R.assign_to_fresh(WRc.get_impl<Mat>().make_unique());
};

template<class V> 
Integer details::make_complex_eigenvectors<V,true>::calc_sel_size(const Mat& TA, const Integer* SELECT, Integer N)
{
    Integer M = 0;
        
    const V* ptr_T  = TA.ptr();
    Integer LDT     = TA.ld();

    for (Integer i = 0; i < N;)
    {
        V sub               = (i < N-1) ? ptr_T[i+1 + i*LDT] : V(0);
        Integer step        = (sub == V()) ? 1 : 2;
        bool sel            = (SELECT[i] == 1) || (step == 2 && SELECT[i+1] == 1);

        if (sel == true)
            M += step;

        i += step;
    };

    return M;
};

template<class V> 
void details::make_complex_eigenvectors<V,true>::eval(const Mat& TA, const Matrix& ind, const Matrix& mat_VL,
            const Matrix& mat_VR, Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
{
    Integer N               = TA.rows();
    const V* ptr_T          = TA.ptr();
    const Integer* ptr_I    = ind.get_array<Integer>();
    Integer LDT             = TA.ld();

    // real and complex part of eigenvectors are stored separately, we need to combine them
    Integer i = 0;
    for (; i < N-1; ++i)
    {
        V sub               = ptr_T[i+1 + i*LDT];

        if (sub == V(0.0))
            continue;

        bool sel            = (ptr_I[i] == 1 || ptr_I[i+1] == 1);

        if (sel == true)
            break;
    };

    if (i == N-1)
    {
        //no complex eigenvalues
        XL  = mat_VL;
        XR  = mat_VR;
        return;
    };

    using VC    = typename md::complex_type<V>::type;
    using Mat_C = raw::Matrix<VC,struct_dense>;

    Matrix mat_VLc  = convert(mat_VL, Mat::matrix_code);
    Matrix mat_VRc  = convert(mat_VR, Mat::matrix_code);

    const Mat& VL   = mat_VLc.get_impl<Mat>();
    const Mat& VR   = mat_VRc.get_impl<Mat>();

    Mat_C VLC       = Mat_C(VL.get_type(), VL.rows(), VL.cols());
    Mat_C VRC       = Mat_C(VR.get_type(), VR.rows(), VR.cols());

    const V* ptr_VL = VL.ptr();
    const V* ptr_VR = VR.ptr();
    VC* ptr_VLC     = VLC.ptr();
    VC* ptr_VRC     = VRC.ptr();

    Integer VL_ld   = VL.ld();
    Integer VR_ld   = VR.ld();
    Integer VLC_ld  = VLC.ld();
    Integer VRC_ld  = VRC.ld();

    for (i = 0; i < N; ++i)
    {
        V sub               = (i < N-1) ? ptr_T[i+1 + i*LDT] : V(0);
        Integer step        = (sub == V()) ? 1 : 2;
        bool sel            = (ptr_I[i] == 1) || (step == 2 && ptr_I[i+1] == 1);

        if (sel == false)
            continue;

        if (step == 1)
        {
            if (comp_left)
            {
                for (Integer j = 0; j < N; ++j)
                    ptr_VLC[j]  = VC(ptr_VL[j]);
            };

            if (comp_right)
            {
                for (Integer j = 0; j < N; ++j)
                    ptr_VRC[j]  = VC(ptr_VR[j]);
            };
        }
        else
        {
            if (comp_left)
            {
                const V* ptr_VL_re  =  ptr_VL;
                const V* ptr_VL_im  =  ptr_VL + VL_ld;
                VC* ptr_VLC_1       = ptr_VLC;
                VC* ptr_VLC_2       = ptr_VLC + VLC_ld;

                for (Integer j = 0; j < N; ++j)
                {
                    ptr_VLC_1[j]    = VC(ptr_VL_re[j], ptr_VL_im[j]);
                    ptr_VLC_2[j]    = VC(ptr_VL_re[j], -ptr_VL_im[j]);
                }
            };

            if (comp_right)
            {
                const V* ptr_VR_re  =  ptr_VR;
                const V* ptr_VR_im  =  ptr_VR + VR_ld;
                VC* ptr_VRC_1       = ptr_VRC;
                VC* ptr_VRC_2       = ptr_VRC + VRC_ld;

                for (Integer j = 0; j < N; ++j)
                {
                    ptr_VRC_1[j]    = VC(ptr_VR_re[j], ptr_VR_im[j]);
                    ptr_VRC_2[j]    = VC(ptr_VR_re[j], -ptr_VR_im[j]);
                }
            };

            ++i;
        };

        if (comp_left)
        {
            ptr_VL      += step * VL_ld;
            ptr_VLC     += step * VLC_ld;
        };

        if (comp_right)
        {
            ptr_VR      += step * VR_ld;
            ptr_VRC     += step * VRC_ld;
        };
    };

    XL  = Matrix(VLC,true);
    XR  = Matrix(VRC,true);
};

template<class V> 
void details::make_complex_eigenvectors<V,true>::complex_vectors_to_real(const Mat& TA, const Matrix& WL, const Matrix& WR, 
                                    Mat& WL_R, Mat& WR_R, const Matrix& ind)
{
    Integer N               = TA.rows();
    const V* ptr_T          = TA.ptr();
    const Integer* ptr_I    = ind.get_array<Integer>();
    Integer LDT             = TA.ld();

    // real and complex part of eigenvectors are combined, we need to separate them
    Integer i = 0;
    for (; i < N-1; ++i)
    {
        V sub               = ptr_T[i+1 + i*LDT];

        if (sub == V(0.0))
            continue;

        bool sel            = (ptr_I[i] == 1 || ptr_I[i+1] == 1);

        if (sel == true)
            break;
    };

    if (i == N-1)
    {
        Matrix WLc      = convert(WL, Mat::matrix_code);
        Matrix WRc      = convert(WR, Mat::matrix_code);

        //no complex eigenvalues
        WL_R.assign_to_fresh(WLc.get_impl<Mat>().make_unique());
        WR_R.assign_to_fresh(WRc.get_impl<Mat>().make_unique());
        return;
    };

    using VC    = typename md::complex_type<V>::type;
    using Mat_C = raw::Matrix<VC,struct_dense>;

    Matrix WLc  = convert(WL, Mat_C::matrix_code);
    Matrix WRc  = convert(WR, Mat_C::matrix_code);

    Mat_C VLC   = WLc.get_impl<Mat_C>().make_unique();
    Mat_C VRC   = WRc.get_impl<Mat_C>().make_unique();
    
    Mat VLR     = Mat(TA.get_type(), WL.rows(), WL.cols());
    Mat VRR     = Mat(TA.get_type(), WR.rows(), WR.cols());

    const VC* ptr_VLC   = VLC.ptr();
    const VC* ptr_VRC   = VRC.ptr();
    V* ptr_VLR          = VLR.ptr();
    V* ptr_VRR          = VRR.ptr();

    Integer VLR_ld      = VLR.ld();
    Integer VRR_ld      = VRR.ld();
    Integer VLC_ld      = VLC.ld();
    Integer VRC_ld      = VRC.ld();

    for (i = 0; i < N; ++i)
    {
        V sub               = (i < N-1) ? ptr_T[i+1 + i*LDT] : V(0);
        Integer step        = (sub == V()) ? 1 : 2;
        bool sel            = (ptr_I[i] == 1) || (step == 2 && ptr_I[i+1] == 1);

        if (sel == false)
            continue;

        if (step == 1)
        {
            for (Integer j = 0; j < N; ++j)
                ptr_VLR[j]  = real(ptr_VLC[j]);

            for (Integer j = 0; j < N; ++j)
                ptr_VRR[j]  = real(ptr_VRC[j]);
        }
        else
        {
            {
                V* ptr_VLR_re       =  ptr_VLR;
                V* ptr_VLR_im       =  ptr_VLR + VLR_ld;

                for (Integer j = 0; j < N; ++j)
                {
                    ptr_VLR_re[j]   = real(ptr_VLC[j]);
                    ptr_VLR_im[j]   = imag(ptr_VLC[j]);
                }
            };

            {
                V* ptr_VRR_re       =  ptr_VRR;
                V* ptr_VRR_im       =  ptr_VRR + VRR_ld;

                for (Integer j = 0; j < N; ++j)
                {
                    ptr_VRR_re[j]   = real(ptr_VRC[j]);
                    ptr_VRR_im[j]   = imag(ptr_VRC[j]);
                }
            };

            ++i;
        };

        ptr_VLR         += step * VLR_ld;
        ptr_VLC         += step * VLC_ld;
        ptr_VRR         += step * VRR_ld;
        ptr_VRC         += step * VRC_ld;
    };

    WL_R.assign_to_fresh(VLR);
    WR_R.assign_to_fresh(VRR);
};

template struct details::make_complex_eigenvectors<Real>;
template struct details::make_complex_eigenvectors<Float>;
template struct details::make_complex_eigenvectors<Complex>;
template struct details::make_complex_eigenvectors<Float_complex>;

template<class V>
Matrix details::make_tridiag_subdiag_real(const Matrix& D1, V* ptr_U, bool is_D1_compl)
{
    if (is_D1_compl == false)
        return D1;

    using VR            = typename md::real_type<V>::type;

    Integer N           = D1.length()+1;
    Matrix ret          = abs(D1);
    const V* ptr        = D1.get_array<V>();
    const VR* ptr_re    = ret.get_array<VR>();

    //scaling elements u_i solves: u_i'*e_i*u_{i+1} = abs(e_i)
    ptr_U[N-1]          = VR(1.0);
    for (Integer i = N-2; i >= 0; --i)
    {
        V tmp           = ptr_re[i]/(ptr_U[i+1]*ptr[i]);
        ptr_U[i]        = conj(tmp);
    };

    return ret;
};

template
Matrix details::make_tridiag_subdiag_real(const Matrix& D1, Real* ptr_U, bool is_D1_compl);

template
Matrix details::make_tridiag_subdiag_real(const Matrix& D1, Float* ptr_U, bool is_D1_compl);

template
Matrix details::make_tridiag_subdiag_real(const Matrix& D1, Complex* ptr_U, bool is_D1_compl);

template
Matrix details::make_tridiag_subdiag_real(const Matrix& D1, Float_complex* ptr_U, bool is_D1_compl);
}