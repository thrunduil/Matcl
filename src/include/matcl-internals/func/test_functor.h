/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#pragma once

#include "matcl-core/matrix/complex_type.h"

namespace matcl { namespace details
{

//=====================================================================
//						TEST MATRIX
//=====================================================================
struct test_geq_zero
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val >= Val(); };
};

struct test_m1_1
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val >= Val(-1) && val <= Val(1);  };
};

struct test_m1_inf
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val >= Val(-1);  };
};

struct test_inf_m1_1_inf
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val <= Val(-1) || val >= Val(1);  };
};

struct test_1_inf
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val >= Val(1);  };
};

struct test_0_1
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return val >= Val(0) && val <= Val(1); };
};

struct test_is_int
{
    static const bool test_complex = false;

    template<class Val>
    bool eval(Val val)      { return (val - (Integer)val) == Val(0); };
    bool eval(Integer)      { return true; };
};

template<bool iso,class vt, class st, class test>
struct test_range_impl
{};

template<class vt, class st, class test>
struct test_range_impl<true,vt,st,test>
{
    using M1 = raw::Matrix<vt,st>;
    static bool eval(const M1&, test& )
    {
        return true;
    };
};

template<class vt, class test>
struct test_range_impl<false,vt,struct_dense,test>
{
    using M1 = raw::Matrix<vt,struct_dense>;
    static bool eval(const M1& A, test& test_object)
    {
        const vt* ptr_A = A.ptr();
        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
            {
                if (!test_object.eval(ptr_A[i]))
                    return false;
            };

            ptr_A += A_ld;
        };

        return true;
    };
};

template<class vt, class test>
struct test_range_impl<false,vt,struct_sparse,test>
{
    using M1 = raw::Matrix<vt,struct_sparse>;
    static bool eval(const M1& A, test& test_object)
    {
        Integer r = A.rows(), c = A.cols(), nz = A.nnz();
        if (Real(nz) < Real(r)*Real(c))
        {
            if (!test_object.eval(0.))
                return false;
        };

        if (nz == 0)
            return true;

        const raw::details::sparse_ccs<vt>& rep = A.rep();
        const vt* d_x = rep.ptr_x() + rep.offset();

        for (Integer i = 0; i < nz; ++i)
        {
            if (!test_object.eval(d_x[i]))
                return false;
        }

        return true;
    };
};

template<class vt, class test>
struct test_range_impl<false,vt,struct_banded,test>
{
    using M1 = raw::Matrix<vt,struct_banded>;
    static bool eval(const M1& A, test& test_object)
    {
        Integer r = A.rows(), c = A.cols(), nz = A.nnz();
        if (Real(nz) < Real(r)*Real(c))
        {
            if (!test_object.eval(0.))
                return false;
        };

        if (nz == 0)
            return true;

        const vt* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j <  c; ++j)
        {
            Integer fr = A.first_row(j);
            Integer lr = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                if (!test_object.eval(ptr_A[pos]))
                    return false;
            };

            ptr_A += A_ld;
        }
        return true;
    };
};

template<class vt, class st, class test>
struct test_range
{
    using M1 = raw::Matrix<vt,st>;
    static const bool iso = std::is_same<vt,Object>::value;	
    static bool eval(const M1& A, test& test_object)
    {
        return test_range_impl<iso,vt,st,test>::eval(A, test_object);
    };
};

template<bool isc, class vt, class st,class test>
struct test_range_compl
{
    using M1 = raw::Matrix<vt,st>;
    static bool eval(const M1& A, test& test_object)
    {
        return test_range_impl<false,vt,st,test>::eval(A, test_object);
    };
};

template<class vt,class st,class test>
struct test_range_compl<false,vt,st,test>
{
    using M1 = raw::Matrix<vt,st>;
    static bool eval(const M1&, test& )
    {
        return true;
    };
};

template<class st, class test>
struct test_range<Complex,st,test>
{
    using M1 = raw::Matrix<Complex,st>;
    static bool eval(const M1& A, test& test_object)
    {
        return test_range_compl<test::test_complex,Complex,st,test>::eval(A, test_object);
    };
};

template<class st, class test>
struct test_range<Float_complex,st,test>
{
    using M1 = raw::Matrix<Float_complex,st>;
    static bool eval(const M1& A, test& test_object)
    {
        return test_range_compl<test::test_complex,Float_complex,st,test>::eval(A, test_object);
    };
};

//=====================================================================
//						TEST MATRIX MATRIX
//=====================================================================
template<class T>	
struct is_conv_to_integer
{};

template<>
struct is_conv_to_integer<Integer>
{
    static const bool eval(Integer )        {	return true;	};
};

template<>
struct is_conv_to_integer<Complex>
{
    static const bool eval(Complex )        {	return false;	};
};

template<>
struct is_conv_to_integer<Float_complex>
{
    static const bool eval(Float_complex )  {	return false;	};
};

template<>
struct is_conv_to_integer<Real>
{
    static const bool eval(Real A)			
    {	
        return (A - (Integer)A) == 0;	
    };
};

template<>
struct is_conv_to_integer<Float>
{
    static const bool eval(Float A)			
    {	
        return (A - (Integer)A) == 0;	
    };
};

template<class test, bool is_inv, class vt1, class vt2>
struct eval_F2
{
    static bool eval(vt1 A, vt2 B)
    {
        return test::eval(A,B);
    };
};

template<class test, class vt1, class vt2>
struct eval_F2<test,true,vt1,vt2>
{
    static bool eval(vt1 A, vt2 B)
    {
        return test::eval(B,A);
    };
};

template<typename value_type_1, typename value_type_2, typename Derived, typename MT1, typename MT2>
struct test_range2_vt_base
{
    static bool base_eval(const MT1& A, const MT2& B)
    {
        return Derived::derived_eval(A, B);
    }
};

template<class vt1, class vt2, class st1, class st2, class test, bool is_inv>
struct test_range2_impl : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1, vt2, 
                    st1, st2, test, is_inv>, raw::Matrix<vt1,st1>, raw::Matrix<vt2,st2>>
{
    using test_range = test_range2_impl<vt1, vt2, st1, st2, test, is_inv>;
    using MT1        = raw::Matrix<vt1,st1> ;
    using MT2        = raw::Matrix<vt2,st2>;
    using base_class = test_range2_vt_base<vt1, vt2, test_range, MT1, MT2> ;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        static_assert(dependent_false<MT1>::value,"this function should not be instantiated");
    }
};

//--------------------------------------------------------------------
//				GM - GM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_dense,struct_dense,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2,test_range2_impl<vt1,vt2,struct_dense,struct_dense,test,is_inv>,
            raw::Matrix<vt1,struct_dense>,raw::Matrix<vt2,struct_dense>>
{
    using MT1           = raw::Matrix<vt1,struct_dense>;
    using MT2           = raw::Matrix<vt2,struct_dense>;
    using test_cl       = test_range2_impl<vt1,vt2,struct_dense,struct_dense,test,is_inv>;
    using base_class    = test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        using eval_func     = eval_F2<test,is_inv,vt1,vt2>;

        const vt1* ptr_A    = A.ptr();
        const vt2* ptr_B    = B.ptr();
        Integer r           = A.rows();
        Integer c           = A.cols();

        Integer A_ld        = A.ld();
        Integer B_ld        = B.ld();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
            {
                if (!eval_func::eval(ptr_A[i],ptr_B[i]))
                    return false;
            };

            ptr_A += A_ld;
            ptr_B += B_ld;
        };
        return true;
    };
};

//--------------------------------------------------------------------
//				GM - SM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_dense,struct_sparse,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_dense,struct_sparse,test,is_inv>,
            raw::Matrix<vt1,struct_dense>,raw::Matrix<vt2,struct_sparse>>
{
    using MT1       = raw::Matrix<vt1,struct_dense>;
    using MT2       = raw::Matrix<vt2,struct_sparse>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_dense,struct_sparse,test,is_inv>;
    using base_class= test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        if (A.rows() == 0 || A.cols() == 0 || B.nnz() == 0)
            return true;

        Integer c = A.cols(); 			
     
        const raw::details::sparse_ccs<vt2>& Bd = B.rep();
        const Integer* Bd_c		= Bd.ptr_c();
        const Integer* Bd_r		= Bd.ptr_r();
        const vt2* Bd_x	        = Bd.ptr_x();
        const vt1* ptr_A        = A.ptr();
        Integer A_ld            = A.ld();

        using eval_func         = eval_F2<test,is_inv,vt1,vt2>;

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer k = Bd_c[j]; k < Bd_c[j + 1] ; ++k)
            {
                Integer p		= Bd_r[k];
                if (!eval_func::eval(ptr_A[p],Bd_x[k]))
                    return false;
            };
            ptr_A += A_ld;
        };
 
        return true;
    };
};

//--------------------------------------------------------------------
//				GM - BM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_dense,struct_banded,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_dense,struct_banded,test,is_inv>,
            raw::Matrix<vt1,struct_dense>, raw::Matrix<vt2,struct_banded>>
{
    using MT1       = raw::Matrix<vt1,struct_dense>;
    using MT2       = raw::Matrix<vt2,struct_banded>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_dense,struct_banded,test,is_inv>;
    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {	
        Integer c = A.cols();

        using eval_func = eval_F2<test,is_inv,vt1,vt2>;

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        if (B.first_diag() == B.last_diag())
        {
            Integer rc          = B.diag_length(B.first_diag());
            const vt1* ptr_A    = A.ptr();
            const vt2* ptr_B    = B.rep_ptr() + B.first_elem_diag(B.first_diag());

            for (Integer j = 0; j < rc; ++j)
            {
                if (!eval_func::eval(ptr_A[0],ptr_B[0]))
                    return false;

                ptr_A += A_ld + 1;
                ptr_B += B_ld;
            };
        }
        else
        {
            const vt1* ptr_A    = A.ptr();
            const vt2* ptr_B    = B.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row = B.first_row(j);
                Integer last_row  = B.last_row(j);
                Integer first_elem= B.first_elem_pos(j);

                for (Integer i = first_row, pos_B = first_elem; i <= last_row; ++i, ++pos_B)
                {
                    if (!eval_func::eval(ptr_A[i],ptr_B[pos_B]))
                        return false;
                };

                ptr_A   += A_ld;
                ptr_B   += B_ld;
            };
        }; 
        return true;
    };
};
//--------------------------------------------------------------------
//				SM - GM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_sparse,struct_dense,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_sparse,struct_dense,test,is_inv>,
                raw::Matrix<vt1,struct_sparse>,raw::Matrix<vt2,struct_dense>>
{
    using MT1       = raw::Matrix<vt1,struct_sparse>;
    using MT2       = raw::Matrix<vt2,struct_dense>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_sparse,struct_dense,test,is_inv>;
    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        return test_range2_impl<vt2,vt1,struct_dense,struct_sparse,test,!is_inv>
            ::eval(B,A);
    };
};

//--------------------------------------------------------------------
//				SM - SM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_sparse,struct_sparse,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_sparse,struct_sparse,test,is_inv>,
        raw::Matrix<vt1,struct_sparse>,raw::Matrix<vt2,struct_sparse>>
{
    using MT1       = raw::Matrix<vt1,struct_sparse>;
    using MT2       = raw::Matrix<vt2,struct_sparse>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_sparse,struct_sparse,test,is_inv>;

    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1, MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        Integer r = A.rows(), c = A.cols();

        if ((r == 0) || (c == 0) )
            return true;
        
        Integer nzret = std::min(A.nnz(),B.nnz());

        if (nzret == 0)
            return true;

        const raw::details::sparse_ccs<vt1>&		Ad = A.rep();
        const raw::details::sparse_ccs<vt2>&		Bd = B.rep();

        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const vt1* Ad_x				= Ad.ptr_x();

        const Integer* Bd_c			= Bd.ptr_c();
        const Integer* Bd_r			= Bd.ptr_r();
        const vt2* Bd_x				= Bd.ptr_x();

        using eval_func = eval_F2<test,is_inv,vt1,vt2>;
         
        for (Integer j = 0; j < c; ++j)
        {
            Integer ka				= Ad_c[j]; 
            Integer kb				= Bd_c[j];

            while (ka < Ad_c[j+1] && kb < Bd_c[j+1])
            {
                if (Ad_r[ka] < Bd_r[kb])
                    ++ka;
                else if (Ad_r[ka] > Bd_r[kb])
                    ++kb; 
                else
                { 
                    if (!eval_func::eval(Ad_x[ka],Bd_x[kb]))
                        return false;

                    ++ka; 
                    ++kb;
                };
            };
        };

        return true;
    };
};

//--------------------------------------------------------------------
//				SM - BM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_sparse,struct_banded,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_sparse,struct_banded,test,is_inv>,
        raw::Matrix<vt1,struct_sparse>,raw::Matrix<vt2,struct_banded>>
{
    using MT1       = raw::Matrix<vt1,struct_sparse>;
    using MT2       = raw::Matrix<vt2,struct_banded>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_sparse,struct_banded,test,is_inv>;
    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        Integer r = A.rows(), c = A.cols();

        if ((r == 0) || (c == 0) || A.nnz() == 0)
            return true;

        const raw::details::sparse_ccs<vt1>& Ad	= A.rep();
        const Integer* Ad_c			= Ad.ptr_c();
        const Integer* Ad_r			= Ad.ptr_r();
        const vt1* Ad_x				= Ad.ptr_x();

        Integer B_ld                = B.ld();

        using eval_func = eval_F2<test,is_inv,vt1,vt2>;

        if (B.first_diag() == B.last_diag())
        {
            Integer rc          = B.diag_length(B.first_diag());
            const vt2* ptr_B    = B.rep_ptr() + B.first_elem_diag(B.first_diag());

            for (Integer j = 0; j < rc; ++j)
            {
                Integer ka;
                bool exist = Ad.has_element(j,j,ka);

                if (exist)
                {
                    if (!eval_func::eval(Ad_x[ka],ptr_B[0]))
                        return false;
                };

                ptr_B += B_ld;
            };
        }
        else
        {
            const vt2* ptr_B = B.rep_ptr();

            for (Integer j = 0; j < c; ++j, ptr_B += B_ld)
            {
                Integer first_row		= B.first_row(j);
                Integer last_row		= B.last_row(j);
                Integer pos_B			= B.first_elem_pos(j);
                Integer kb				= first_row;

                if (first_row >= r)
                    continue;

                Integer ka;
                Ad.has_element(kb,j,ka);

                if (ka == Ad_c[j+1] || Ad_r[ka] > last_row)
                    continue;

                while (ka < Ad_c[j+1] && kb <= last_row)
                {
                    if (Ad_r[ka] > kb)
                    {
                        ++kb; 
                        ++pos_B;
                    }
                    else
                    { 
                        if (!eval_func::eval(Ad_x[ka],ptr_B[pos_B]))
                            return false;

                        ++ka; 
                        ++kb;
                        ++pos_B;
                    };
                };                
            };
        };

        return true;
    };
};
//--------------------------------------------------------------------
//				BM - GM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_banded,struct_dense,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_banded,struct_dense,test,is_inv>,
        raw::Matrix<vt1,struct_banded>,raw::Matrix<vt2,struct_dense>>
{
    using MT1       = raw::Matrix<vt1,struct_banded>;
    using MT2       = raw::Matrix<vt2,struct_dense>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_banded,struct_dense,test,is_inv>;

    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1, MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        return test_range2_impl<vt2,vt1,struct_dense,struct_banded,test,!is_inv>
            ::eval(B,A);
    };
};

//--------------------------------------------------------------------
//				BM - SM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_banded,struct_sparse,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_banded,struct_sparse,test,is_inv>,
        raw::Matrix<vt1,struct_banded>,raw::Matrix<vt2,struct_sparse>>
{
    using MT1       = raw::Matrix<vt1,struct_banded>;
    using MT2       = raw::Matrix<vt2,struct_sparse>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_banded,struct_sparse,test,is_inv>;
    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1,MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        return test_range2_impl<vt2,vt1,struct_sparse,struct_banded,test,!is_inv>
            ::eval(B,A);
    };
};

//--------------------------------------------------------------------
//				BM - BM
//--------------------------------------------------------------------
template<class vt1, class vt2, class test, bool is_inv>
struct test_range2_impl<vt1,vt2,struct_banded,struct_banded,test,is_inv> 
    : public test_range2_vt_base<vt1, vt2, test_range2_impl<vt1,vt2,struct_banded,struct_banded,test,is_inv>,
        raw::Matrix<vt1,struct_banded>,raw::Matrix<vt2,struct_banded>>
{
    using MT1       = raw::Matrix<vt1,struct_banded>;
    using MT2       = raw::Matrix<vt2,struct_banded>;
    using test_cl   = test_range2_impl<vt1,vt2,struct_banded,struct_banded,test,is_inv>;
    using base_class = test_range2_vt_base<vt1, vt2, test_cl, MT1, MT2>;

    static bool eval(const MT1& A, const MT2& B)
    {
        return base_class::base_eval( A, B );
    }

    static bool derived_eval(const MT1& A, const MT2& B)
    {
        Integer fda     = A.first_diag();
        Integer lda     = A.last_diag();
        Integer fdb     = B.first_diag();
        Integer ldb     = B.last_diag();
    
        Integer fd      = std::max(fda, fdb);
        Integer ld      = std::min(lda, ldb);

        Integer A_ld    = A.ld();
        Integer B_ld    = B.ld();

        using eval_func = eval_F2<test,is_inv,vt1,vt2>;

        if (fd == ld)
        {
            Integer rc          = A.diag_length(fd);
            const vt2* ptr_B    = B.rep_ptr() + B.first_elem_diag(fd);
            const vt1* ptr_A    = A.rep_ptr() + A.first_elem_diag(fd);

            for (Integer i = 0; i < rc; ++i)
            {
                if (!eval_func::eval(ptr_A[0], ptr_B[0]))
                    return false;

                ptr_A += A_ld;
                ptr_B += B_ld;
            };

            return true;
        };

        for (Integer d = fd; d <= ld; ++d)
        {
            const vt2* ptr_B    = B.rep_ptr() + B.first_elem_diag(d);
            const vt1* ptr_A    = A.rep_ptr() + A.first_elem_diag(d);
            Integer s           = A.diag_length(d);

            for (Integer j = 0; j < s; ++j)
            {
                if (!eval_func::eval(ptr_A[0], ptr_B[0]))
                    return false;

                ptr_A           += A_ld;
                ptr_B           += B_ld;
            };
        };

        return true;
    };
};

template<class M1, class M2, class test>
struct test_range2
{
    static bool eval(const M1& A, const M2& B)
    {
        using vt1 = typename M1::value_type;
        using vt2 = typename M2::value_type;
        using st1 = typename M1::struct_type;
        using st2 = typename M2::struct_type;
        return test_range2_impl<vt1,vt2,st1,st2,test,false>::eval(A,B);
    };
};

};};