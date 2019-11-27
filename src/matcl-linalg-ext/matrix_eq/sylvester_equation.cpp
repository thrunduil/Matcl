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

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-linalg/matrix_eq/sylvester_equation.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/linear_eq/linsolve.h"

//TODO
//#include "extern_lib_src/slicot/include/slicot.h"

#if 0
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-core/utils/workspace.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

template<class V>
Matrix solve_sylvester_val(const V& A, const V& B, const V& C)
{
    using constants::nan;
    typedef typename V::value_type  vt;
    typedef typename V::struct_type st;

    //matrices should already be tested

    Integer M   = C.rows();
    Integer N   = C.cols();    
    Matrix Am   = Matrix(A,false);
    Matrix Bm   = Matrix(B,false);

    Matrix Q,T,Z,S;

    // A = QTQ', T schur form
    // B = ZSZ', S schur form
    tie(Q, T)   = schur(Am);
    tie(Z, S)   = schur(Bm);

    Matrix rhs  = chain_mult(ctrans(Q), Matrix(C,true), Z);        // TY + YS = Q'CZ
    V Tc        = T.impl_unique<V>();
    V Sc        = S.impl_unique<V>();
    V rhsc      = rhs.impl_unique<V>();

    rhsc.set_struct(struct_flag());
    
    Integer info    = 0;
    char trana      = 'N';
    char tranb      = 'N';
    Integer isgn    = 1;
    Real scale      = 0;

    lapack::trsyl(&trana, &tranb, isgn, M, N, lap(Tc.ptr()), Tc.ld(), lap(Sc.ptr()), Sc.ld(),
                  lap(rhsc.ptr()), rhsc.ld(), lap(&scale), lap(&info));
    
    if (info != 0 || scale == 0.)
    {
        throw error::error_sylv();
    }

    // X = QYZ'
    return chain_mult(Q, Matrix(rhsc,true), ctrans(Z), 1./ scale);
};

template<class V>
matcl::mat_tup_2 solve_gsylvester_val(const V& A1, const V& A2, const V& B1, const V& B2,
                                      const V& C1, const V& C2)
{
    using constants::nan;
    typedef typename V::value_type  vt;
    typedef typename V::struct_type st;

    //matrices should already be tested

    Integer M   = A1.rows();
    Integer N   = A2.cols();

    Matrix Q1, Z1, TA1, TB1;
    Matrix Q2, Z2, TA2, TB2;

    //A1 = Q1 * TA1 * Z1',  B1 = Q1 * TB1 * Z1'.
    tie(Q1,Z1,TA1,TB1)  = gschur(Matrix(A1,false), Matrix(B1,false));

    //A2 = Q2 * TA2 * Z2',  B2 = Q2 * TB2 * Z2'.
    tie(Q2,Z2,TA2,TB2)  = gschur(Matrix(A2,false), Matrix(B2,false));
    TA2                 = -TA2;
    TB2                 = -TB2;

    //TA1 * (Z1' * X * Z2) - (Q1' * Y * Q2) * TA2 = Q1' * C1 * Z2
    //TB1 * (Z1' * X * Z2) - (Q1' * Y * Q2) * TB2 = Q1' * C2 * Z2

    Matrix C1_tr    = chain_mult(ctrans(Q1), Matrix(C1,false), Z2);
    Matrix C2_tr    = chain_mult(ctrans(Q1), Matrix(C2,false), Z2);
   
    C1_tr.set_struct(struct_flag());
    C2_tr.set_struct(struct_flag());

    V C1_rep        = C1_tr.impl_unique<V>();
    V C2_rep        = C2_tr.impl_unique<V>();

    V TA1_rep       = TA1.impl_unique<V>();
    V TA2_rep       = TA2.impl_unique<V>();
    V TB1_rep       = TB1.impl_unique<V>();
    V TB2_rep       = TB2.impl_unique<V>();
    
    Integer info    = 0;
    char trana      = 'N';
    Integer ijob    = 0;
    Real scale      = 0;
    vt work         = 0;

    typedef typename matcl::details::lapack_value_type<vt>::type VL;

    using iworkspace        = matcl::pod_workspace<Integer>;
    Integer liwork          = M + N + 6;
    iworkspace IWORK        = iworkspace(liwork);
    Integer* ptr_IWORK      = IWORK.ptr();

    lapack::tgsyl<VL>(&trana, ijob, M, N, 
                        lap(TA1_rep.ptr()), TA1_rep.ld(), 
                        lap(TA2_rep.ptr()), TA2_rep.ld(),
                        lap(C1_rep.ptr()), C1_rep.ld(),
                        lap(TB1_rep.ptr()), TB1_rep.ld(),    
                        lap(TB2_rep.ptr()), TB2_rep.ld(),    
                        lap(C2_rep.ptr()), C2_rep.ld(),
                        lap(&scale), nullptr, lap(&work), 1, ptr_IWORK, lap(&info));
    
    if (info != 0 || scale == 0.)
    {
        throw error::error_gsylv();
    }

    //Z1 * C1_rep * Z2' = X
    //Q1 * C2_rep * Q2' = Y

    Matrix X        = chain_mult(Z1, Matrix(C1_rep,false), ctrans(Z2), 1./scale);
    Matrix Y        = chain_mult(Q1, Matrix(C2_rep,false), ctrans(Q2), 1./scale);
    
    return mat_tup_2(X,Y);
};

}

Matrix solve_sylvester(const Matrix& A, const Matrix& B, const Matrix& C,
                       const bool fast_exit_linsolve_allowed)
{
    if (!A.is_square() || !B.is_square() || A.rows() != C.rows() || B.cols() != C.cols())
    {
        throw error::error_size_sylv(A.rows(),A.cols(), B.rows(), B.cols(), C.rows(), C.cols());
    }

    matcl::value_code vA = A.get_value_code();
    matcl::value_code vB = B.get_value_code();
    matcl::value_code vC = C.get_value_code();
    matcl::value_code v  = (matcl::value_code)max(max((Integer)vC, (Integer)vA), (Integer)vB);

    if (v == value_code::v_object)
    {
        throw error::object_value_type_not_allowed("sylvester");
    }

    Integer N           = A.rows();
    Integer M           = B.cols();

    bool isv =    A.all_finite() && B.all_finite() && C.all_finite();

    if (isv == false)
    {
        Matrix X    = repmat(constants::nan(), N, M);

        return X;
    };

    if (C.structural_nnz() == 0)
    {
        return spzeros(N,M);
    }
    
    if (C.numel() == 1)
    {
        if ((A + B == 0))
        {
            throw error::error_sylv();
        }
        return div(C , A + B);
    }
    if (A.structural_nnz() == 0 && B.structural_nnz() == 0)
    {
        throw error::error_sylv();
    }
    
    if (   ((A.structural_nnz() == 0) && B.get_struct().is_id())
        || ((B.structural_nnz() == 0) && A.get_struct().is_id()))
    {
        return C;
    }
    
    if (A.get_struct().is_id() && B.get_struct().is_id())
    {
        return div(C, 2.0f);
    }

    // fast exits below require linsolve to have nan treatment
    // for now - only for finite inputs
    if(fast_exit_linsolve_allowed == true)
    {
        try
        {
            if (B.structural_nnz() == 0)
            {
                return linsolve(A, C);
            }
            if (A.structural_nnz() == 0)
            {
                //return trans(linsolve(trans(B), trans(C)));
                return linsolve_rev(B, C);                
            }
            if (B.get_struct().is_id())
            {
                return linsolve(speye(A.rows(), A.cols()) + A, C);
            }
            if (A.get_struct().is_id())
            {
                //return trans(linsolve(trans(speye(B.rows(), B.cols()) + B), trans(C)));
                return linsolve_rev(speye(B.rows(), B.cols()) + B, C);
            }
            if (B.numel() == 1)
            {
                Matrix newrhs = B * speye(A.rows(), A.cols()) + A;
                return linsolve(newrhs, C);
            }
            if (A.numel() == 1)
            {
                Matrix newrhs = A * speye(B.rows(), B.cols()) + B;
                return linsolve_rev(newrhs, C);

                //Matrix newrhs = trans(A * speye(B.rows(), B.cols()) + B);
                //return trans(linsolve(newrhs, trans(C)));
            }
        }
        catch(error::error_singular&)
        {
            throw error::error_sylv();
        }
    }   
    
    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
        {
            Matrix Ac = convert(A,mat_code::real_dense);
            Matrix Bc = convert(B,mat_code::real_dense);
            Matrix Cc = convert(C,mat_code::real_dense);
            
            typedef matcl::raw::Matrix<Real,struct_dense> DM;
            
            return details::solve_sylvester_val<DM>(Ac.get_impl<DM>(),Bc.get_impl<DM>(),
                                                    Cc.get_impl<DM>());
        }
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {
            Matrix Ac = convert(A,mat_code::complex_dense);
            Matrix Bc = convert(B,mat_code::complex_dense);
            Matrix Cc = convert(C,mat_code::complex_dense);
            
            typedef matcl::raw::Matrix<Complex,struct_dense> DM;
            
            return details::solve_sylvester_val<DM>(Ac.get_impl<DM>(),Bc.get_impl<DM>(),
                                                    Cc.get_impl<DM>());
        }
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("sylvester");
        }
    };
    throw error::error_general("impossible type case in solve_sylvester");
}

Matrix solve_dsylvester(const Matrix& A, const Matrix& B, const Matrix& C)
{
    if (!A.is_square() || !B.is_square() || A.rows() != C.rows() || B.cols() != C.cols())
    {
        throw error::error_size_sylv(A.rows(),A.cols(), B.rows(), B.cols(), C.rows(), C.cols());
    }

    matcl::value_code vA = A.get_value_code();
    matcl::value_code vB = B.get_value_code();
    matcl::value_code vC = C.get_value_code();
    matcl::value_code v  = (matcl::value_code)max(max((Integer)vC, (Integer)vA), (Integer)vB);

    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
        {
            Integer N           = A.rows();
            Integer M           = B.cols();

            bool isv =    A.all_finite() && B.all_finite() && C.all_finite();
            if (isv == false)
            {
                Matrix X    = repmat(constants::nan(), N, M);
                return X;
            }
            if (A.structural_nnz() == 0 || B.structural_nnz() == 0)
            {
                return C;
            }
            if (C.numel() == 0 || C.structural_nnz() == 0)
            {
                return spzeros(N,M);
            }
            if (C.numel() == 1)
            {
                if ((1 + A * B == 0))
                {
                    throw error::error_sylv();
                }
                return div(C, 1.0f + A * B);
            }
            if (A.structural_nnz() == 0 && B.structural_nnz() == 0)
            {
                return C;
            }
            Matrix Ac = convert(A,mat_code::real_dense);
            Matrix Bc = convert(B,mat_code::real_dense);
            Matrix Cc = convert(C,mat_code::real_dense);
            
            using raw::real_dense;
            using raw::integer_dense;
            using details::lap;

            raw::real_dense A_impl = Ac.get_impl_unique<raw::real_dense>();
            raw::real_dense B_impl = Bc.get_impl_unique<raw::real_dense>();
            raw::real_dense C_impl = Cc.get_impl_unique<raw::real_dense>();
            C_impl.set_struct(struct_flag());

            raw::real_dense Z_impl(ti::ti_real(), M, M);
            raw::integer_dense iwork(ti::ti_int(), 4 * N, 1);
            // (For optimum performance LDWORK should be larger.)
            // choosing 2x minimum
            Integer         info        = 0;
            Integer         ldwork      = 2 * std::max(std::max(1, 2*N*N + 9*N), std::max(5*M, N + M));
            raw::real_dense      dwork(ti::ti_real(), ldwork, 1);
            
            Integer lda = A_impl.ld();
            Integer ldb = B_impl.ld();
            Integer ldc = C_impl.ld();
            Integer ldz = Z_impl.ld();

            slicot::sb04qd_(&N, &M, lap(A_impl.ptr()), lap(&lda), lap(B_impl.ptr()), lap(&ldb),
                 lap(C_impl.ptr()), lap(&ldc), lap(Z_impl.ptr()), (&ldz), lap(iwork.ptr()), lap(dwork.ptr()),
                 lap(&ldwork), lap(&info));
    
            if (info != 0)
            {
                throw error::error_sylv();
            }
            return Matrix(C_impl, true);
        }
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {
            throw error::error_complex_value_type_not_allowed();
        }
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("sylvester");
        }
    };
    throw error::error_general("impossible type case in solve_dsylvester");
}

template<class V>
struct eval_22
{
    V   A1, A2, C1;
    V   B1, B2, C2;

    V   X, Y;

    eval_22(const matcl::Matrix& A1_, const matcl::Matrix& A2_, 
            const matcl::Matrix& B1_, const matcl::Matrix& B2_,
            const matcl::Matrix& C1_, const matcl::Matrix& C2_)
        :A1(A1_.get_scalar<V>()), A2(A2_.get_scalar<V>())
        ,B1(B1_.get_scalar<V>()), B2(B2_.get_scalar<V>())
        ,C1(C1_.get_scalar<V>()), C2(C2_.get_scalar<V>())
    {
        make();
    }
    void make()
    {
        V abs_11    = abs(A1);
        V abs_12    = abs(A2);
        
        if (abs_11 >= abs_12)
        {
            if (abs_11 == 0.)
                throw error::error_gsylv();

            //X  = -A2/A1 * Y + C1/A1
            //(B2-B1 * A2/A1) * Y = C2 - B1 * C1/A1

            V D1    = B2 - B1 * A2 / A1;
            V E1    = C2 - B1 * C1 / A1;

            if (D1 == 0.)
                throw error::error_gsylv();

            Y       = E1 / D1;
            X       = -A2/A1 * Y + C1/A1;

            return;
        }
        else
        {
            if (abs_12 == 0.)
                throw error::error_gsylv();

            //Y = -A1/A2 * X + C1/A2
            //(B1 - B2 * A1/A2) * X = C2 - B2 * C1/A2

            V D1    = B1 - B2 * A1 / A2;
            V E1    = C2 - B2 * C1 / A2;

            if (D1 == 0.)
                throw error::error_gsylv();

            X       = E1 / D1;
            Y       = -A1/A2 * X + C1/A2;

            return;
        }
    };
};

mat_tup_2 solve_gsylvester(const Matrix& A1, const Matrix& A2, 
                           const Matrix& B1, const Matrix& B2, 
                           const Matrix& C1, const Matrix& C2)
{
    if (!A1.is_square() || !A2.is_square() || !B1.is_square() || !B2.is_square()
        || A1.rows() != B1.rows() || A2.rows() != B2.rows()
        || C1.rows() != A1.rows() || C1.cols() != A2.cols()
        || C2.rows() != B1.rows() || C2.cols() != B2.cols())
    {
        throw error::error_size_gsylv(A1.rows(), A1.cols(), A2.rows(), A2.cols(),
                                      B1.rows(), B1.cols(), B2.rows(), B2.cols(),
                                      C1.rows(), C1.cols(), C2.rows(), C2.cols());
    }

    Integer N           = A1.rows();
    Integer M           = A2.rows();

    value_code v;
    {
        Integer vA1     = (Integer)A1.get_value_code();
        Integer vA2     = (Integer)A2.get_value_code();
        Integer vB1     = (Integer)B1.get_value_code();
        Integer vB2     = (Integer)B2.get_value_code();
        Integer vC1     = (Integer)C1.get_value_code();
        Integer vC2     = (Integer)C2.get_value_code();

        Integer c1      = max(vA1, vB1);
        Integer c2      = max(vA2, vB2);
        Integer c3      = max(vC1, vC2);

        v               = (value_code)max(max(c1, c2), c3);
    };

    if (v == value_code::v_object)
    {
        throw error::object_value_type_not_allowed("sylvester");
    }

    bool isv    = A1.all_finite() && A2.all_finite() && B1.all_finite() && B2.all_finite()
                && C1.all_finite() && C2.all_finite();

    if (isv == false)
    {
        Matrix X    = repmat(constants::nan(), N, M);
        Matrix Y    = repmat(constants::nan(), N, M);

        return mat_tup_2(X,Y);
    };

    bool is_zero_C1     = C1.numel() == 0 || C1.structural_nnz() == 0;
    bool is_zero_C2     = C2.numel() == 0 || C2.structural_nnz() == 0;

    if (is_zero_C1 && is_zero_C2)
    {
        Matrix X    = spzeros(N,M);
        Matrix Y    = spzeros(N,M);

        return mat_tup_2(X,Y);
    }

    //reduction to the standard sylvester equation
    if (B1.get_struct().is_id() == true && B2.get_struct().is_id() == true)
    {
        //X = -Y + C2;

        Matrix A    = -A1;
        Matrix B    = A2;
        Matrix C    = C1 - A1 * C2;

        Matrix Y;
        
        try
        {
            Y = solve_sylvester(A, B, C);        
        }
        catch(error::error_sylv&)
        {
            throw error::error_gsylv();
        };

        Matrix X    = - Y + C2;
        
        return mat_tup_2(X,Y);
    };
    if (A1.get_struct().is_id() == true && A2.get_struct().is_id() == true)
    {
        // X + Y = -Y + C1
        Matrix A    = -B1;
        Matrix B    = B2;
        Matrix C    = C2 - B1 * C1;

        Matrix Y;

        try
        {
            Y = solve_sylvester(A, B, C);        
        }
        catch(error::error_sylv&)
        {
            throw error::error_gsylv();
        };

        Matrix X    = - Y + C1;
        
        return mat_tup_2(X,Y);
    };

    if (C1.numel() == 1)
    {
        switch (v)
        {
            case value_code::v_integer:
            case value_code::v_real:
            case value_code::v_float:
            {
                //TODO: impl float
            }
            {
                eval_22<Real> obj(A1, A2, B1, B2, C1, C2);

                Real X = obj.X;
                Real Y = obj.Y;

                return mat_tup_2(X,Y);
            }
            case value_code::v_float_complex:
            {
                //TODO: impl float
            }
            case value_code::v_complex:
            {
                eval_22<Complex> obj(A1, A2, B1, B2, C1, C2);

                Complex X = obj.X;
                Complex Y = obj.Y;

                return mat_tup_2(X,Y);
            }
        };

        throw error::error_general("impossible type case in solve_gsylvester");
    };

    bool is_zero_A1     = A1.numel() == 0 || A1.structural_nnz() == 0;
    bool is_zero_A2     = A2.numel() == 0 || A2.structural_nnz() == 0;

    if (is_zero_A1 == true && is_zero_A2 == true)
    {
        throw error::error_gsylv();
    }

    bool is_zero_B1     = B1.numel() == 0 || B1.structural_nnz() == 0;
    bool is_zero_B2     = B2.numel() == 0 || B2.structural_nnz() == 0;

    if (is_zero_B1 == true && is_zero_B2 == true)
    {
        throw error::error_gsylv();
    }

    try
    {
        if (is_zero_A1 == true)
        {
            Matrix Y    = linsolve_rev(A2,C1);
            Matrix X    = linsolve(B1, C2 - Y * B2);
            return mat_tup_2(X,Y);
        }
        else if (is_zero_A2 == true)
        {
            Matrix X    = linsolve(A1, C1);
            Matrix Y    = linsolve_rev(B2,C2 - B1 * X);
            return mat_tup_2(X,Y);
        }
        else if (is_zero_B1 == true)
        {
            Matrix Y    = linsolve_rev(B2,C2);
            Matrix X    = linsolve(A1, C1 - Y * A2);
            return mat_tup_2(X,Y);
        }
        else if (is_zero_B2 == true)
        {
            Matrix X    = linsolve(B1, C2);
            Matrix Y    = linsolve_rev(A2,C1 - A1 * X);
            return mat_tup_2(X,Y);
        };

        if (A1.numel() == 1)
        {
            if ((abs(A1) >= abs(B1)))
            {
                //A1 != 0

                Matrix scal = div(B1,A1);

                Matrix D    = B2 - scal * A2;
                Matrix E    = C2 - scal * C1;
                
                Matrix Y    = linsolve_rev(D,E);

                Matrix sc2  = div(1.0f,A1);

                Matrix X    = sc2 * (C1 - Y * A2);

                return mat_tup_2(X,Y);
            }
            else
            {
                //B1 != 0

                Matrix scal = div(A1,B1);

                Matrix D    = A2 - scal * B2;
                Matrix E    = C1 - scal * C2;
                
                Matrix Y    = linsolve_rev(D,E);

                Matrix sc2  = div(1.0f,B1);

                Matrix X    = sc2 * (C2 - Y * B2);

                return mat_tup_2(X,Y);
            };
        };

        if (A2.numel() == 1)
        {
            if ((abs(A2) >= abs(B2)))
            {
                //A2 != 0
                Matrix scal = div(B2,A2);

                Matrix D    = B1 - scal * A1;
                Matrix E    = C2 - scal * C1;
                
                Matrix X    = linsolve(D,E);

                Matrix sc2  = div(1.0f,A2);

                Matrix Y    = sc2 * (C1 - A1*X);

                return mat_tup_2(X,Y);
            }
            else
            {
                //B2 != 0
                Matrix scal = div(A2,B2);

                Matrix D    = A1 - scal * B1;
                Matrix E    = C1 - scal * C2;
                
                Matrix X    = linsolve(D,E);

                Matrix sc2  = div(1.0f,B2);

                Matrix Y    = sc2 * (C2 - B1*X);

                return mat_tup_2(X,Y);
            };
        };
    }
    catch(error::error_singular&)
    {
        throw error::error_gsylv();
    };    

    //general case
    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
        {
            Matrix A1c  = convert(A1,mat_code::real_dense);
            Matrix B1c  = convert(B1,mat_code::real_dense);
            Matrix C1c  = convert(C1,mat_code::real_dense);

            Matrix A2c  = convert(A2,mat_code::real_dense);
            Matrix B2c  = convert(B2,mat_code::real_dense);
            Matrix C2c  = convert(C2,mat_code::real_dense);
            
            typedef matcl::raw::Matrix<Real,struct_dense> DM;
            
            return details::solve_gsylvester_val<DM>
                        (A1c.get_impl<DM>(),A2c.get_impl<DM>(),
                         B1c.get_impl<DM>(),B2c.get_impl<DM>(),
                         C1c.get_impl<DM>(),C2c.get_impl<DM>());
                                                    
        }
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {                        
            Matrix A1c  = convert(A1,mat_code::complex_dense);
            Matrix B1c  = convert(B1,mat_code::complex_dense);
            Matrix C1c  = convert(C1,mat_code::complex_dense);

            Matrix A2c  = convert(A2,mat_code::complex_dense);
            Matrix B2c  = convert(B2,mat_code::complex_dense);
            Matrix C2c  = convert(C2,mat_code::complex_dense);
            
            typedef matcl::raw::Matrix<Complex,struct_dense> DM;
            
            return details::solve_gsylvester_val<DM>
                        (A1c.get_impl<DM>(),A2c.get_impl<DM>(),
                         B1c.get_impl<DM>(),B2c.get_impl<DM>(),
                         C1c.get_impl<DM>(),C2c.get_impl<DM>());
        }
    };

    throw error::error_general("impossible type case in solve_gsylvester");
};

}

#pragma warning( pop )

#endif