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

#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-core/details/integer.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant

namespace matcl { namespace raw 
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

namespace details
{
    template<class From, class To>
    struct is_precision_increased
    {
        static const bool is_from_single    = md::is_single_precision<From>::value;
        static const bool is_to_double      = md::is_double_precision<To>::value;

        static const bool value             = (is_to_double == true) && (is_from_single == true);
    };

    struct warn_state 
    {
        bool allow_warn_real_compl;
        bool allow_warn_int_real;
        bool allow_warn_int_float;
        bool allow_warn_float_float_compl;
        bool allow_warn_float_real;

        warn_state()
        {
            allow_warn_real_compl           = true;
            allow_warn_int_real             = true;
            allow_warn_int_float            = true;
            allow_warn_float_float_compl    = true;
            allow_warn_float_real           = true;
        };

        bool any_warnings() const
        {
            return allow_warn_real_compl == false || allow_warn_int_real == false
                    || allow_warn_int_float == false || allow_warn_float_float_compl == false
                    || allow_warn_float_real == false;
        };

        bool any_prec_lost() const
        {
            return allow_warn_int_real == false || allow_warn_int_float == false 
                    || allow_warn_float_real == false;
        };

        bool any_compl_to_real() const
        {
            return allow_warn_real_compl == false || allow_warn_float_float_compl == false;
        };
    };

    template<class ret, class in>
    struct conver_scal_warn_impl
    {};

    template<>
    struct conver_scal_warn_impl<Integer,Float>
    {
        static Integer eval(ti::ti_empty, warn_state& w, const Float& val)
        {
            if (w.allow_warn_int_float && val - (Integer)val != 0)
            {
                w.allow_warn_int_float = false;
                error::get_global_messanger()->warning_precision_lost_float_to_int(val);                    
            };
            return (Integer)val;
        };
    };

    template<>
    struct conver_scal_warn_impl<Integer,Real>
    {
        static Integer eval(ti::ti_empty, warn_state& w, const Real& val)
        {
            if (w.allow_warn_int_real && val - (Integer)val != 0)
            {
                w.allow_warn_int_real = false;
                error::get_global_messanger()->warning_precision_lost_real_to_int(val);
            };
            return (Integer)val;
        };
    };

    template<>
    struct conver_scal_warn_impl<Float,Float_complex>
    {
        static Float eval(ti::ti_empty, warn_state& w, const Float_complex& val)
        {
            if (w.allow_warn_float_float_compl && imag(val) != 0)
            {
                w.allow_warn_float_float_compl = false;
                error::get_global_messanger()->warning_precision_lost_float_compl_to_float(val);            
            };
            return real(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Integer,Float_complex>
    {
        static Integer eval(ti::ti_empty ti, warn_state& w, const Float_complex& val)
        {
            Float tmp = conver_scal_warn_impl<Float,Float_complex>::eval(ti, w, val);
            return conver_scal_warn_impl<Integer,Float>::eval(ti, w, tmp);
        };
    };

    template<>
    struct conver_scal_warn_impl<Real,Complex>
    {
        static Real eval(ti::ti_empty, warn_state& w, const Complex& val)
        {
            if (w.allow_warn_real_compl && imag(val) != 0)
            {
                w.allow_warn_real_compl = false;
                error::get_global_messanger()->warning_precision_lost_compl_to_real(val);            
            };
            return real(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Integer,Complex>
    {
        static Integer eval(ti::ti_empty ti, warn_state& w, const Complex& val)
        {
            Real tmp = conver_scal_warn_impl<Real,Complex>::eval(ti, w, val);
            return conver_scal_warn_impl<Integer,Real>::eval(ti, w, tmp);
        };
    };

    template<>
    struct conver_scal_warn_impl<Integer,Object>
    {
        static Integer eval(ti::ti_empty, warn_state& , const Object& val)
        {            
            return cast_integer(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float,Real>
    {
        static Float eval(ti::ti_empty, warn_state& w, const Real& val)
        {
            if (w.allow_warn_float_real && Float(val) != val)
            {
                w.allow_warn_float_real = false;
                error::get_global_messanger()->warning_precision_lost_real_to_float(val);            
            };
            return Float(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float,Complex>
    {
        static Float eval(ti::ti_empty ti, warn_state& w, const Complex& val)
        {
            Real tmp = conver_scal_warn_impl<Real,Complex>::eval(ti, w, val);
            return conver_scal_warn_impl<Float,Real>::eval(ti, w, tmp);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float,Object>
    {
        static Float eval(ti::ti_empty, warn_state&, const Object& val)
        {            
            return cast_float(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Real,Float_complex>
    {
        static Real eval(ti::ti_empty ti, warn_state& w, const Float_complex& val)
        {
            Float tmp = conver_scal_warn_impl<Float,Float_complex>::eval(ti, w, val);
            return tmp;
        };
    };

    template<>
    struct conver_scal_warn_impl<Real,Object>
    {
        static Real eval(ti::ti_empty, warn_state& , const Object& val)
        {            
            return cast_real(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float_complex,Real>
    {
        static Float_complex eval(ti::ti_empty ti, warn_state& w, const Real& val)
        {   
            Float tmp = conver_scal_warn_impl<Float,Real>::eval(ti, w, val);
            return Float_complex(tmp);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float_complex,Complex>
    {
        static Float_complex eval(ti::ti_empty ti, warn_state& w, const Complex& val)
        {   
            Float r = conver_scal_warn_impl<Float,Real>::eval(ti, w, real(val));
            Float i = conver_scal_warn_impl<Float,Real>::eval(ti, w, imag(val));
            return Float_complex(r,i);
        };
    };

    template<>
    struct conver_scal_warn_impl<Float_complex,Object>
    {
        static Float_complex eval(ti::ti_empty, warn_state&, const Object& val)
        {   
            return cast_float_complex(val);
        };
    };

    template<>
    struct conver_scal_warn_impl<Complex,Object>
    {
        static Complex eval(ti::ti_empty, warn_state& , const Object& val)
        {            
            return cast_complex(val);
        };
    };

    template<class T>
    struct conver_scal_warn_impl<Object,T>
    {
        static Object eval(ti::ti_object ti, warn_state& , const T& val)
        {
            return convert_to_object(ti, val);
        }; 
    };

    template<class ret>
    struct conver_scal_warn_impl<ret,ret>
    {
        static ret eval(ti::ti_empty, warn_state& , const ret& val)
        {
            return val;
        }; 
    };

    template<> struct conver_scal_warn_impl<Object,Object>
    {
        static Object eval(ti::ti_object ti, warn_state& warn, const Object& val)  
        { 
            (void)warn;
            return Object(ti,val); 
        };
    };

    template<>
    struct conver_scal_warn_impl<Float,Integer>
    {
        static Float eval(ti::ti_empty, warn_state& , const Integer& val)
        {
            return Float(val);
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Real,Integer>
    {
        static Real eval(ti::ti_empty, warn_state& , const Integer& val)
        {
            return Real(val);
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Float_complex,Integer>
    {
        static Float_complex eval(ti::ti_empty, warn_state& , const Integer& val)
        {
            return Float_complex(Float(val));
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Complex,Integer>
    {
        static Complex eval(ti::ti_empty, warn_state& , const Integer& val)
        {
            return Complex(Real(val));
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Real,Float>
    {
        static Real eval(ti::ti_empty, warn_state& , const Float& val)
        {
            return Real(val);
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Float_complex,Float>
    {
        static Float_complex eval(ti::ti_empty, warn_state& , const Float& val)
        {
            return Float_complex(val);
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Complex,Float>
    {
        static Complex eval(ti::ti_empty, warn_state& , const Float& val)
        {
            return Complex(Real(val));
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Complex,Real>
    {
        static Complex eval(ti::ti_empty, warn_state& , const Real& val)
        {
            return Complex(val);
        }; 
    };

    template<>
    struct conver_scal_warn_impl<Complex,Float_complex>
    {
        static Complex eval(ti::ti_empty, warn_state& , const Float_complex& val)
        {
            return Complex(real(val), imag(val));
        }; 
    };

    template<class ret>
    struct conver_scal_warn
    {
        template<class T>
        static ret eval(typename ti::get_ti_type<ret>::type ti, warn_state& warn, const T& val)
        {
            return conver_scal_warn_impl<ret,T>::eval(ti,warn,val);
        };
    };

    template<class ret, class T>
    struct converter_scalars_impl
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& val)
        {
            warn_state warn;
            return conver_scal_warn_impl<ret, T>::eval(ti, warn, val);
        };
    };

    template<class V, class T>
    struct converter_scalars_impl<raw::Matrix<V,struct_dense>,T>
    {
        using ret = raw::Matrix<V,struct_dense>;
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& val)
        {
            V val2 = converter_scalars_impl<V,T>::eval(ti, val);
            return ret(ti, val2, 1, 1);
        };
    };

    template<class V, class T>
    struct converter_scalars_impl<raw::Matrix<V,struct_banded>,T>
    {
        using ret = raw::Matrix<V,struct_banded>;
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& val)
        {
            V val2 = converter_scalars_impl<V,T>::eval(ti, val);
            return ret(ti, val2, 1, 1, 0, 0);
        };
    };

    template<class V, class T>
    struct converter_scalars_impl<raw::Matrix<V,struct_sparse>,T>
    {
        using ret = raw::Matrix<V,struct_sparse>;
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& val)
        {
            V val2      = converter_scalars_impl<V,T>::eval(ti, val);
            Integer ri  = 0;
            Integer ci  = 0;
            return ret(ti, &ri, &ci, &val2, 1, 1, 1);
        };
    };

    template<class ret, class T>
    struct converter_sparse_full
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& mat)
        {
            using value_type    = typename ret::value_type;
            using in_type       = typename T::value_type;

            static bool prec_increased  = is_precision_increased<in_type, value_type>::value;

            Integer j, k, c = mat.cols();
            value_type Z = matcl::details::default_value<value_type>(ti);
            Matrix<value_type,struct_dense> tmp(ti,Z, mat.rows(), c);
            const sparse_ccs<in_type>& rep = mat.rep();
        
            warn_state warn;

            if (mat.nnz() != 0)
            {
                const Integer* rep_c = rep.ptr_c();
                const Integer* rep_r = rep.ptr_r();
                const in_type* rep_x    = rep.ptr_x();
                value_type* ptr_tmp     = tmp.ptr();
                Integer tmp_ld          = tmp.ld();

                for (j = 0; j < c; ++j)
                {
                    for (k = rep_c[j]; k < rep_c[j+1]; ++k)
                    {
                        value_type val = conver_scal_warn<value_type>::eval(ti,warn,rep_x[k]);
                        mrd::reset_helper(ptr_tmp[rep_r[k]],val);
                    };
                    ptr_tmp += tmp_ld;
                };

                tmp.set_struct(mat.get_struct());
                if (warn.any_warnings())
                {
                    bool prec_lost      = warn.any_prec_lost();
                    bool compl_to_real  = warn.any_compl_to_real();
                    struct_flag sf      = md::predefined_struct::set_warnings(tmp.get_struct(), 
                                                prec_lost, compl_to_real);
                    tmp.set_struct(sf);
                }
                else if (prec_increased)
                {
                    struct_flag sf      = md::predefined_struct::precision_increased(tmp.get_struct());
                    tmp.set_struct(sf);
                }
            }
            else
            {
                tmp.get_struct().set(predefined_struct_type::diag);
            };

            return tmp;
        }
    };

    template<class ret, class T>
    struct converter_sparse_band
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& mat)
        {
            Integer r           = mat.rows();
            Integer c           = mat.cols();

            using val_type      = typename T::value_type;
            using val_type_ret  = typename ret::value_type;

            static bool prec_increased  = is_precision_increased<val_type, val_type_ret>::value;

            val_type_ret Z      = matcl::details::default_value<val_type_ret>(ti);

            if (mat.nnz() == 0)
            {
                ret out(ti,r,c,0,0);
                val_type_ret* ptr_out = out.rep_ptr();

                Integer s = std::min(r,c);
                for (Integer i = 0; i < s; ++i)
                    mrd::reset_helper(ptr_out[i],Z);

                out.get_struct().set(predefined_struct_type::diag);
                return out;
            };            

            Integer ld = raw::get_ld(mat,-1);
            Integer ud = raw::get_ud(mat,-1);

            ret res(ti,r,c,-ld,ud);

            val_type_ret* ptr_res   = res.rep_ptr();
            Integer res_ld          = res.ld();

            Integer res_size = res.impl_size();
            for (Integer i = 0; i < res_size; ++i)
                mrd::reset_helper(ptr_res[i],Z);

            warn_state warn;

            const sparse_ccs<val_type>& Ad = mat.rep();
            const Integer * Ad_c = Ad.ptr_c();
            const Integer * Ad_r = Ad.ptr_r();
            const val_type* Ad_x = Ad.ptr_x();

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row  = res.first_row(j);
                Integer first_pos  = res.first_elem_pos(j) - first_row;

                for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
                {
                    Integer row     = Ad_r[k];

                    //element can be outside the band
                    if (row - j < -ud || row - j > ld)
                        continue;

                    mrd::reset_helper(ptr_res[first_pos + row],
                                       conver_scal_warn<val_type_ret>::eval(ti,warn,Ad_x[k]));
                };

                ptr_res += res_ld;
            };

            res.set_struct(mat.get_struct());
            
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(res.get_struct(), 
                                            prec_lost, compl_to_real);
                res.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(res.get_struct());
                res.set_struct(sf);
            }

            return res;
        };
    };

    //A bit more generic version of fully proxy caller
    template<template <class,class> class caller_, class ret, class T, bool ise>
    struct converter_generic_eval_proxy
    {
        using caller = caller_<ret,T>;

        static ret eval( typename ti::get_ti_type<ret>::type ti, const T& m )
        {
            return caller::eval_f( ti, m );
        }
    };

    template<template <class,class> class caller_, class ret, class T>
    struct converter_generic_eval_proxy<caller_, ret, T, true>
    {
        using caller = caller_<ret,T>;

        static ret eval( typename ti::get_ti_type<ret>::type ti, const T& m )
        {
            return caller::eval_t( ti, m );
        }
    };

    template<class ret, class T>
    struct converter_sparse_sparse
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            static const bool ise = std::is_same<typename ret::value_type,typename T::value_type>::value;            
            return converter_generic_eval_proxy< matcl::raw::details::converter_sparse_sparse, ret, T, ise >::eval( ti, m );
        };

        static ret eval_f(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            Integer r = m.rows(), c = m.cols(), n = m.nnz();
            using val_type_ret  = typename ret::value_type;
            using val_type_in   = typename T::value_type;

            static bool prec_increased  = is_precision_increased<val_type_in, val_type_ret>::value;

            ret res(ti, r,c,n);
            sparse_ccs<val_type_ret>& d = res.rep();
            const sparse_ccs<val_type_in>& Ad = m.rep();

            warn_state warn;

            if (n == 0)
                return res;

            Integer* d_c            = d.ptr_c();
            Integer* d_r            = d.ptr_r();
            val_type_ret* d_x       = d.ptr_x();

            const Integer* Ad_c    = Ad.ptr_c();
            const Integer* Ad_r    = Ad.ptr_r() + Ad.offset();
            const val_type_in* Ad_x= Ad.ptr_x() + Ad.offset();


            for (Integer i = 1; i <= c; ++i)
                d_c[i] = Ad_c[i] - Ad.offset();

            for (Integer i = 0; i < n; ++i)
            {
                d_r[i] = Ad_r[i];
                mrd::reset_helper(d_x[i],conver_scal_warn<val_type_ret>::eval(ti,warn,Ad_x[i]));
            };

            res.set_struct(m.get_struct());
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(res.get_struct(), 
                                            prec_lost, compl_to_real);
                res.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(res.get_struct());
                res.set_struct(sf);
            }

            return res;
        }

        static ret eval_t(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            static const bool iso = std::is_same<typename ret::value_type,Object>::value;
            if (iso)
            {
                if (ti == m.get_type())
                    return m;
                else
                    return eval_f(ti,m);
            }
            else
            {
                return m;
            };
        };
    };

    template<class ret, class T>
    struct converter_dense_dense
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            static const bool ise = std::is_same<typename ret::value_type,typename T::value_type>::value;
            return converter_generic_eval_proxy< matcl::raw::details::converter_dense_dense, ret, T, ise >::eval( ti, m );
        };

        static ret eval_f(typename ti::get_ti_type<ret>::type ti, const T& mat)
        {
            using ret_type  = typename ret::value_type;
            using in_type   = typename T::value_type;

            static bool prec_increased  = is_precision_increased<in_type, ret_type>::value;

            Integer r = mat.rows();
            Integer c = mat.cols();

            Matrix<ret_type,struct_dense> tmp(ti,r,c);

            warn_state warn;

            ret_type* ptr_tmp       = tmp.ptr();
            const in_type* ptr_mat  = mat.ptr();
            Integer tmp_ld          = tmp.ld();
            Integer mat_ld          = mat.ld();

            for(Integer j = 0; j < c; ++j)
            {
                for(Integer i = 0; i < r; ++i)
                    mrd::reset_helper(ptr_tmp[i],conver_scal_warn<ret_type>::eval(ti,warn, ptr_mat[i]));

                ptr_tmp += tmp_ld;
                ptr_mat += mat_ld;
            };

            tmp.set_struct(mat.get_struct());
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(tmp.get_struct(), 
                                            prec_lost, compl_to_real);
                tmp.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(tmp.get_struct());
                tmp.set_struct(sf);
            }

            return tmp;
        };

        static ret eval_t(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            static const bool iso = std::is_same<typename ret::value_type,Object>::value;
            if (iso)
            {
                if (ti == m.get_type())
                    return m;
                else
                    return eval_f(ti,m);
            }
            else
            {
                return m;
            };
        };
    };

    template<class ret, class T>
    struct converter_dense_band
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& mat)
        {
            Integer r = mat.rows(), c = mat.cols();

            using val_type_ret  = typename ret::value_type;
            using val_type      = typename T::value_type;

            static bool prec_increased  = is_precision_increased<val_type, val_type_ret>::value;

            if (r == 0 || c == 0)
                return ret(ti,r,c,0,0);

            Integer ld = raw::get_ld(mat,-1);
            Integer ud = raw::get_ud(mat,-1);

            ret res(ti,r,c,-ld,ud);
            
            const val_type* ptr     = mat.ptr();

            val_type_ret* ptr_ret   = res.rep_ptr();
            Integer mat_ld          = mat.ld();
            Integer res_ld          = res.ld();

            warn_state warn;

            for (Integer j = 0; j < c; ++j)
            {
                Integer first_row   = res.first_row(j);
                Integer first_pos   = res.first_elem_pos(j) - first_row;
                Integer end         = std::min(j+ld,r-1);

                for (Integer k = std::max(j-ud,0); k <= end; ++k)
                {
                    mrd::reset_helper(ptr_ret[first_pos+k],
                                       conver_scal_warn<val_type_ret>::eval(ti,warn,ptr[k]));
                };

                ptr         += mat_ld;
                ptr_ret     += res_ld;
            };

            res.set_struct(mat.get_struct());
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(res.get_struct(), 
                                            prec_lost, compl_to_real);
                res.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(res.get_struct());
                res.set_struct(sf);
            }

            return res;
        };
    };

    template<class ret, class T>
    struct converter_dense_sparse
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& m)
        {
            Integer nnz     = 0;
            Integer r       = m.rows();
            Integer c       = m.cols();
            Integer m_ld    = m.ld();

            if (m.get_struct().is_diag())
            {
                nnz = std::min(r,c);
            }
            else
            {
                const typename T::value_type* ptr_m = m.ptr();

                for (Integer j = 0; j < c; ++j)
                {
                    for (Integer i = 0; i < r; ++i)
                    {
                        if (!mrd::is_zero(ptr_m[i])) 
                            ++nnz;
                    };

                    ptr_m += m_ld;
                };
            };

            ret res(ti,r, c, nnz);
            if (nnz == 0)
                return res;

            using in_val    = typename T::value_type;
            using ret_val   = typename ret::value_type;

            static bool prec_increased  = is_precision_increased<in_val, ret_val>::value;

            warn_state warn;

            sparse_ccs<ret_val>& dat = res.rep();

            if (m.get_struct().is_diag())
            {
                Integer rc = std::min(r,c);

                using V         = typename T::value_type;
                const V* ptr    = m.ptr();

                Integer nz = 0;

                Integer* dat_c  = dat.ptr_c();
                Integer* dat_r  = dat.ptr_r();
                ret_val* dat_x  = dat.ptr_x();

                for (Integer j = 0; j < rc; ++j)
                {
                    dat_c[j] = nz;
                    if (!mrd::is_zero(ptr[0]))
                    {
                        mrd::reset_helper(dat_x[nz],conver_scal_warn<ret_val>::eval(ti,warn,ptr[0]));
                        dat_r[nz] = j;
                        ++nz;
                    };

                    ptr += m_ld + 1;
                };

                for (Integer j = rc; j <= c; ++j)
                    dat_c[j] = nz;
            }
            else
            {
                Integer* dat_c = dat.ptr_c();
                Integer* dat_r = dat.ptr_r();
                ret_val* dat_x = dat.ptr_x();

                using V         = typename T::value_type;
                const V* ptr_m  = m.ptr();

                Integer k = 0;

                for (Integer j = 0; j < c; ++j)
                {
                    dat_c[j] = k;

                    for (Integer i = 0; i < r; ++i)
                    {
                        const V& val = ptr_m[i];
                        if (!mrd::is_zero(val))
                        {
                            mrd::reset_helper(dat_x[k],conver_scal_warn<ret_val>::eval(ti,warn,val));
                            dat_r[k] = i;
                            ++k;
                        };
                    };

                    ptr_m += m_ld;
                };

                dat_c[c] = k;
            };

            res.set_struct(m.get_struct());
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(res.get_struct(), 
                                            prec_lost, compl_to_real);
                res.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(res.get_struct());
                res.set_struct(sf);
            }

            return res;
        };
    };

    template<class ret, class T>
    struct converter_dense
    {
        using converter_type = typename matcl::details::select_if
                <
                    std::is_same<typename ret::struct_type,struct_dense>::value,
                    converter_dense_dense<ret,T>,
                    typename matcl::details::select_if
                        <
                            std::is_same<typename ret::struct_type,struct_banded>::value,
                            converter_dense_band<ret,T>,
                            converter_dense_sparse<ret,T>
                        >::type
                >::type;

        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& mat)
        {
            return converter_type::eval(ti,mat);
        };
    };

    template<class ret_type,class T>
    struct converter_band_band
    {
        static ret_type eval(typename ti::get_ti_type<ret_type>::type ti, const T& m)
        {
            static const bool ise = std::is_same<typename ret_type::value_type,typename T::value_type>::value;
            return converter_generic_eval_proxy< matcl::raw::details::converter_band_band, ret_type, T, ise >
                ::eval( ti, m );
        };

        static ret_type eval_f(typename ti::get_ti_type<ret_type>::type ti, const T& m)
        {
            Integer r   = m.rows();
            Integer c   = m.cols();
            Integer fd  = m.first_diag();
            Integer ld  = m.last_diag();

            ret_type out(ti,r, c, fd, ld);

            using VTR   = typename ret_type::value_type;
            using VT    = typename T::value_type;

            static bool prec_increased  = is_precision_increased<VT, VTR>::value;

            warn_state warn;

            const VT* ptr_m = m.rep_ptr();
            VTR* ptr_out    = out.rep_ptr();
            Integer out_ld  = out.ld();
            Integer m_ld    = m.ld();

            for (Integer j = 0; j < c; ++j)
            {
                Integer fr  = m.first_row(j);
                Integer lr  = m.last_row(j);
                Integer pos = m.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pos)
                    mrd::reset_helper(ptr_out[pos],conver_scal_warn<VTR>::eval(ti, warn, ptr_m[pos]));
                
                ptr_out     += out_ld;
                ptr_m       += m_ld;
            };

            out.set_struct(m.get_struct());

            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(out.get_struct(), 
                                            prec_lost, compl_to_real);
                out.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(out.get_struct());
                out.set_struct(sf);
            }

            return out;
        };

        static ret_type eval_t(typename ti::get_ti_type<ret_type>::type ti, const T& m)
        {
            static const bool iso = std::is_same<typename ret_type::value_type,Object>::value;
            if (iso)
            {
                if (ti == m.get_type())
                    return m;
                else
                    return eval_f(ti,m);
            }
            else
            {
                return m;
            };
        };
    };

    template<class ret_type,class T>
    struct converter_band_dense
    {
        static ret_type eval(typename ti::get_ti_type<ret_type>::type ti, const T& m)
        {
            using VTR   = typename ret_type::value_type;
            using VT    = typename T::value_type;

            static bool prec_increased  = is_precision_increased<VT, VTR>::value;

            Integer r   = m.rows();
            Integer c   = m.cols();

            ret_type out(ti,matcl::details::default_value<VTR>(ti),r, c);

            warn_state warn;

            if (r == 0 || c == 0)
                return out;

            const VT* ptr_m = m.rep_ptr();
            VTR* ptr_out    = out.ptr();
            Integer out_ld  = out.ld();
            Integer m_ld    = m.ld();

            for (Integer j = 0; j < c; ++j)
            {
                Integer fr  = m.first_row(j);
                Integer lr  = m.last_row(j);
                Integer pos = m.first_elem_pos(j);

                for (Integer k = fr; k <= lr; ++k, ++ pos)
                    mrd::reset_helper(ptr_out[k],conver_scal_warn<VTR>::eval(ti,warn,ptr_m[pos]));

                ptr_out     += out_ld;
                ptr_m       += m_ld;
            };

            out.set_struct(m.get_struct());

            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(out.get_struct(), 
                                            prec_lost, compl_to_real);
                out.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(out.get_struct());
                out.set_struct(sf);
            }

            Integer ld  = m.number_subdiagonals();
            Integer ud  = m.number_superdiagonals();

            if (ld == 0)
                out.add_struct(predefined_struct_type::triu);
            else if (ld == 1)
                out.add_struct(predefined_struct_type::hessu);

            if (ud == 0)
                out.add_struct(predefined_struct_type::tril);
            else if (ud == 1)
                out.add_struct(predefined_struct_type::hessl);

            return out;
        };
    };

    template<class ret_type,class T>
    struct converter_band_sparse
    {
        static ret_type eval(typename ti::get_ti_type<ret_type>::type ti, const T& bm)
        {
            Integer r   = bm.rows();
            Integer c   = bm.cols();
            Integer nnz = bm.nnz();

            ret_type res(ti,r, c, nnz);

            if (nnz == 0)
                return res;

            using VTR   = typename ret_type::value_type;
            using VT    = typename T::value_type;

            static bool prec_increased  = is_precision_increased<VT, VTR>::value;

            warn_state warn;            

            sparse_ccs<VTR>& dat = res.rep();

            Integer* d_c        = dat.ptr_c();
            Integer* d_r        = dat.ptr_r();
            VTR* d_x            = dat.ptr_x();
            Integer bm_ld       = bm.ld();

            const VT* ptr_m     = bm.rep_ptr();
            Integer nzr         = 0;

            for (Integer j = 0; j < c; ++j)
            {
                d_c[j]          = nzr;
                Integer fr      = bm.first_row(j);
                Integer lr      = bm.last_row(j);
                Integer pos     = bm.first_elem_pos(j);

                for (Integer i = fr; i <= lr; ++i, ++pos)
                {
                    VTR tmp     = conver_scal_warn<VTR> ::eval(ti,warn,ptr_m[pos]);

                    d_r[nzr]    = i;
                    mrd::reset_helper(d_x[nzr],tmp);
                    ++nzr;
                }

                ptr_m           += bm_ld;
            };

            d_c[c]  = nzr;

            res.set_struct(bm.get_struct());
            
            if (warn.any_warnings())
            {
                bool prec_lost      = warn.any_prec_lost();
                bool compl_to_real  = warn.any_compl_to_real();
                struct_flag sf      = md::predefined_struct::set_warnings(res.get_struct(), 
                                            prec_lost, compl_to_real);
                res.set_struct(sf);
            }
            else if (prec_increased)
            {
                struct_flag sf      = md::predefined_struct::precision_increased(res.get_struct());
                res.set_struct(sf);
            }

            Integer ld  = bm.number_subdiagonals();
            Integer ud  = bm.number_superdiagonals();

            if (ld == 0)
                res.add_struct(predefined_struct_type::triu);
            else if (ld == 1)
                res.add_struct(predefined_struct_type::hessu);

            if (ud == 0)
                res.add_struct(predefined_struct_type::tril);
            else if (ud == 1)
                res.add_struct(predefined_struct_type::hessl);

            return res;
        };
    };	

    template<class ret, class T>
    struct converter_sparse
    {
        static ret eval(typename ti::get_ti_type<ret>::type ti, const T& val)
        {
            using converter_type = typename matcl::details::select_if
                    <
                        !matcl::details::is_sparse<ret>::value,
                        typename matcl::details::select_if
                            <
                                std::is_same<typename ret::struct_type,struct_banded>::value,
                                converter_sparse_band<ret,T>,
                                converter_sparse_full<ret,T>
                            >:: type,
                        converter_sparse_sparse<ret,T>
                    >:: type;

            return converter_type::eval(ti,val);
        };
    };

    template<class ret, class T>
    struct converter_unknown
    {};

    template<class ret, class T>
    struct converter_band
    {
        using type = typename matcl::details::select_if
            <
                std::is_same<typename ret::struct_type,struct_dense>::value,
                details::converter_band_dense<ret,T>,
                typename matcl::details::select_if
                <
                    std::is_same<typename ret::struct_type,struct_sparse>::value,
                    details::converter_band_sparse<ret,T>,
                    details::converter_band_band<ret,T>
                >::type
            >::type;
    };

    template<class ret, class T>
    struct converter_selector_4
    {
        using type = typename matcl::details::lazy_select_if
            <
                std::is_same<typename T::struct_type,struct_banded>::value,
                details::converter_band<ret,T>,
                details::converter_unknown<ret,T>
            >::type;
    };

    template<class ret, class T>
    struct converter_selector_3
    {
        using type = typename matcl::details::lazy_select_if
            <
                std::is_same<typename T::struct_type,struct_dense>::value,
                matcl::details::lazy_type<details::converter_dense<ret,T>>,
                details::converter_selector_4<ret,T>
            >::type;
    };

    template<class ret, class T>
    struct converter_selector_2
    {
        using type = typename matcl::details::lazy_select_if
            <
                matcl::details::is_sparse<T>::value,
                matcl::details::lazy_type<details::converter_sparse<ret,T>>,
                details::converter_selector_3<ret,T>
            >::type;
    };

    template<class ret, class T>
    struct converter_selector
    {
        using type = typename matcl::details::lazy_select_if
            <
                matcl::details::is_scalar<T>::value,
                matcl::details::lazy_type<details::converter_scalars_impl<ret,T>>,
                details::converter_selector_2<ret,T>
            >::type;
    };

    template<class ret, class T>
    ret converter_impl<ret,T>::eval(const T& val, const tinfo_ret& ti)
    {
        using converter_type = typename details::converter_selector<ret,T>::type;
        return converter_type::eval(ti, val);
    };

};

};};


MACRO_INSTANTIATE_GG_2_F(matcl::raw::details::converter_impl)
MACRO_INSTANTIATE_GST_2_F(matcl::raw::details::converter_impl)
MACRO_INSTANTIATE_STG_2_F(matcl::raw::details::converter_impl)
MACRO_INSTANTIATE_STST_2_F(matcl::raw::details::converter_impl)

MACRO_INSTANTIATE_GS_2_F(matcl::raw::details::converter_impl)
MACRO_INSTANTIATE_STS_2_F(matcl::raw::details::converter_impl)

template struct matcl::raw::details::converter_impl<matcl::Integer,matcl::Object>;
template struct matcl::raw::details::converter_impl<matcl::Real,matcl::Object>;
template struct matcl::raw::details::converter_impl<matcl::Float,matcl::Object>;
template struct matcl::raw::details::converter_impl<matcl::Float_complex,matcl::Object>;
template struct matcl::raw::details::converter_impl<matcl::Complex,matcl::Object>;

template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Integer>;
template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Float>;
template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Real>;
template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Float_complex>;
template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Complex>;
template struct matcl::raw::details::converter_impl<matcl::Object,matcl::Object>;

#pragma warning(pop)
