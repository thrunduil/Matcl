/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/func/raw/io.h"
#include "matcl-matrep/base/serialize.h"
#include "matcl-core/details/IO/disp_stream_options.h"
#include "matcl-core/options/options_disp.h"
#include "matcl-core/details/IO/printer.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

#pragma warning (disable: 4127)

namespace matcl
{

namespace details
{
    static inline 
    char tolower_impl(char c)
    {
        return (char)::tolower(c);
    };

    template<class value_type>
    struct matrix_value_string{};

    template<> struct matrix_value_string<Integer>
    {
        static std::string eval()	{	return "Integer";	};
    };

    template<> struct matrix_value_string<Real>
    {
        static std::string eval()	{	return "Real";		};
    };

    template<> struct matrix_value_string<Float>
    {
        static std::string eval()	{	return "Float";		};
    };

    template<> struct matrix_value_string<Complex>
    {
        static std::string eval()	{	return "Complex";	};
    };

    template<> struct matrix_value_string<Float_complex>
    {
        static std::string eval()	{	return "Float_complex";	};
    };

    template<> struct matrix_value_string<Object>
    {
        static std::string eval()	{	return "Object";	};
    };

    template<class struct_type>
    struct matrix_struct_string{};

    template<> struct matrix_struct_string<struct_dense>
    {
        static std::string eval()	{	return "dense";	};
    };

    template<> struct matrix_struct_string<struct_banded>
    {
        static std::string eval()	{	return "banded";	};
    };

    template<> struct matrix_struct_string<struct_sparse>
    {
        static std::string eval()	{	return "sparse";	};
    };

    static std::string get_struct_flag_string(struct_flag sf)
    {
        return sf.save_as_string();
    };

    template<class T>
    struct save_ti
    {
        static void eval(std::ostream&, const matcl::ti::ti_type<T>&){};
    };

    template<>
    struct save_ti<Object>
    {
        static void eval(std::ostream& os, const matcl::ti::ti_type<Object>& ti)
        {
            os << ' ';
            raw::save(os, ti);
        };
    };

    struct save_functor : public details::extract_type_switch<std::ostream&,save_functor,true>
    {
        template<class T>
        static std::ostream& eval(const Matrix&, const T& mat, std::ostream& os, const std::string& comments)
        {
            matcl::struct_flag sf = mat.get_struct();

            os << "\n[\n";
            os << matrix_value_string<typename T::value_type>::eval() << ' ';
            os << matrix_struct_string<typename T::struct_type>::eval();
            os << ' ' << get_struct_flag_string(sf);

            save_ti<typename matcl::details::hide_object_type<typename T::value_type>::type>
                    ::eval(os, mat.get_type());
            
            os << '\n';

            matcl::raw::stream_helpers::save_comments(os, comments);

            raw::save(os,mat);
            os << "\n]\n";
            return os;
        };

        template<class T>
        static std::ostream& eval_scalar(const Matrix&, const T& mat, std::ostream& os, const std::string& comments)
        {
            os << "\n[\n";
            os << matrix_value_string<T>::eval() << ' ';
            os << "scalar";
            os << '\n';

            matcl::raw::stream_helpers::save_comments(os, comments);

            boost::io::ios_flags_saver old_flags(os);
            boost::io::ios_precision_saver old_prec(os);

            int precision = mr::get_stream_precision<T>::eval();

            if (precision > 0)
                os << std::setprecision(precision) << std::scientific;

            matcl::raw::stream_helpers::write(os,mat);
            os << "\n]\n";
            return os;
        };
    };

    struct mm_save_functor : public details::extract_type_switch<std::ostream&,mm_save_functor,true>
    {
        template<class T>
        static std::ostream& eval(const Matrix&, const T& mat, std::ostream& os, const std::string& comments)
        {
            matcl::raw::mm_helper_save<T>::eval(os, mat, comments);
            return os;
        };

        template<class T>
        static std::ostream& eval_scalar(const Matrix& h, const T& mat, std::ostream& os, const std::string& comments)
        {
            using Mat = raw::Matrix<T,struct_dense>;
            Mat A       = raw::converter<Mat,T>::eval(mat);
            return eval<Mat>(h, A, os, comments);
        };

        template<class S>
        static std::ostream& eval(const Matrix&, const raw::Matrix<Object,S>&, std::ostream&, const std::string&)
        {
            throw error::object_value_type_not_allowed("mm_save");
        };

        static std::ostream& eval_scalar(const Matrix&, const Object&, std::ostream&, const std::string&)
        {
            throw error::object_value_type_not_allowed("mm_save");
        };
    };

    template<>
    std::ostream& save_functor::eval_scalar<Object>(const Matrix&, const Object& mat, std::ostream& os, 
                                                    const std::string& comments)
    {
        os << "\n[\n";
        os << matrix_value_string<Object>::eval() << ' ';
        os << "scalar";
        save_ti<typename matcl::details::hide_object_type<Object>::type>::eval(os, mat.get_type());
        os << '\n';

        if (comments.size() != 0)
        {
            matcl::raw::stream_helpers::save_comments(os, comments);
        };

        //TODO: formatting?
        matcl::raw::stream_helpers::write(os,mat);
        os << "\n]\n";
        return os;
    };

    template<class T>
    struct load_functor_impl
    {
        static Matrix eval(std::istream& is)
        {
            using value_type    = typename T::value_type;
            using matrix_type   = raw::Matrix<value_type,typename T::struct_type>;

            ti::ti_empty ti;
            matrix_type tmp(ti);
            raw::load(is,tmp);
            return Matrix(tmp,false);
        };
    };

    template<class s_type>
    struct load_functor_impl<raw::Matrix<Object,s_type>>
    {
        static Matrix eval(std::istream& is)
        {
            using matrix_type = raw::Matrix<Object,s_type>;

            ti::ti_object ti;
            raw::load(is,ti);

            matrix_type tmp(ti);
            raw::load(is,tmp);
            return Matrix(tmp,false);
        };
    };

    struct load_functor : public details::extract_type_from_code<Matrix,load_functor>
    {
        template<class T>
        static Matrix eval(std::istream& is)
        {
            return load_functor_impl<T>::eval(is);
        };

        template<class T>
        static Matrix eval_scalar(std::istream& is)
        {
            T tmp;
            matcl::raw::stream_helpers::skip_all_blank(is);
            matcl::raw::stream_helpers::read(is,tmp);
            matcl::raw::stream_helpers::skip_all_blank(is);
            return Matrix(tmp,false);
        };
    };

    struct mm_load_functor
    {
        static void eval(std::istream& is, Matrix& m, std::string& comments)
        {
            return matcl::raw::mm_helper_load::eval(is, m, comments);
        };
    };

    template<>
    Matrix load_functor::eval_scalar<matcl::Object>(std::istream& is)
    {
        ti::ti_object tio;
        raw::load(is,tio);

        Object tmp(tio);
        matcl::raw::stream_helpers::skip_all_blank(is);
        matcl::raw::stream_helpers::read(is,tmp);
        matcl::raw::stream_helpers::skip_all_blank(is);
        return Matrix(tmp,false);
    };

    static void load_begin(std::istream& is)
    {
        matcl::raw::stream_helpers::skip_all_blank(is);

        if (!is.good())
            throw error::unable_to_read_matrix();

        char c;
        is.get(c);

        if (!is.good() || c != '[')
            throw error::unable_to_read_matrix();

        matcl::raw::stream_helpers::skip_all_blank(is);

        if (!is.good())
            throw error::unable_to_read_matrix();
    };

    static void load_end(std::istream& is)
    {
        char c;
        is.get(c);
    
        if (c != ']')
            throw error::unable_to_read_matrix();

        is.get(c);

        if (c != '\n')
            throw error::unable_to_read_matrix();

    };

    static void load_matrix_type(std::istream& is,matcl::value_code& val_type,
                            matcl::struct_code& struct_type, struct_flag& sf,
                                 std::string& comments)
    {
        std::string tmp;

        is >> tmp;

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), &tolower_impl);

        if (!is.good())
            throw error::unable_to_read_matrix();

        if (tmp == "integer" || tmp == "int")
        {
            val_type = value_code::v_integer;
        }
        else if (tmp == "float" || tmp == "single")
        {
            val_type = value_code::v_float;
        }
        else if (tmp == "real" || tmp == "double")
        {
            val_type = value_code::v_real;
        }
        else if (tmp == "float_complex")
        {
            val_type = value_code::v_float_complex;
        }
        else if (tmp == "complex")
        {
            val_type = value_code::v_complex;
        }
        else if (tmp == "object")
        {
            val_type = value_code::v_object;
        }
        else
        {
            throw error::unable_to_read_matrix();
        };

        is >> tmp;
        if (!is.good())
            throw error::unable_to_read_matrix();

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), &tolower_impl);

        if (tmp == "dense")
        {
            struct_type = struct_code::struct_dense;
        }
        else if (tmp == "banded" || tmp == "band" )
        {
            struct_type = struct_code::struct_banded;
        }
        else if (tmp == "sparse")
        {
            struct_type = struct_code::struct_sparse;
        }
        else if (tmp == "scalar")
        {
            struct_type = struct_code::struct_scalar;
        }
        else
        {
            throw error::unable_to_read_matrix();
        };

        if (struct_type != struct_code::struct_scalar)
        {
            is >> tmp;
            sf.load_from_string(tmp);
        }
        else
        {
            sf  = struct_flag();
        };        

        if (matcl::raw::stream_helpers::check_nl(is) == false)
            throw error::unable_to_read_matrix();

        matcl::raw::stream_helpers::load_comments(is, comments);

        matcl::raw::stream_helpers::skip_all_blank(is);

        if (!is.good())
            throw error::unable_to_read_matrix();
    };
};

std::ostream& matcl::operator<<(std::ostream& os, const Matrix& A)
{
    std::string comments;
    return details::save_functor::make<const Matrix&,std::ostream&>(A,os,comments);
};

std::ostream& matcl::save(std::ostream& os, const Matrix& A)
{
    std::string comments;
    return details::save_functor::make<const Matrix&,std::ostream&>(A,os, comments);
};

std::ostream& matcl::save(std::ostream& os, const Matrix& A, const std::string& comments)
{
    return details::save_functor::make<const Matrix&,std::ostream&>(A,os, comments);
}

std::ostream& matcl::mm_save(std::ostream& os, const Matrix& A)
{
    std::string comments;
    return details::mm_save_functor::make<const Matrix&,std::ostream&>(A,os, comments);
};

std::ostream& matcl::mm_save(std::ostream& os, const Matrix& A, const std::string& comments)
{
    return details::mm_save_functor::make<const Matrix&,std::ostream&>(A,os, comments);
};

std::istream& matcl::operator>>(std::istream& is, Matrix& A)
{
    std::string comments;
    return load(is, A, comments);
};

std::istream& matcl::load(std::istream& is, Matrix& m)
{
    std::string comments;
    return load(is, m, comments);
};

std::istream& matcl::load(std::istream& is, Matrix& A, std::string& comments)
{
    matcl::value_code   val_type;
    matcl::struct_code  struct_type;
    struct_flag sf;

    details::load_begin(is);

    details::load_matrix_type(is,val_type,struct_type,sf, comments);

    matcl::mat_code mat_type = matrix_traits::get_matrix_type(val_type,struct_type);

    A = details::load_functor::make<std::istream&>(mat_type,is);
    details::load_end(is);
    A.set_struct(sf);
    return is;
};

std::istream& matcl::mm_load(std::istream& is, Matrix& m)
{
    std::string comments;
    details::mm_load_functor::eval(is, m, comments);
    return is;
};

std::istream& matcl::mm_load(std::istream& is, Matrix& m, std::string& comments)
{
    details::mm_load_functor::eval(is, m, comments);
    return is;
};

void matcl::disp_header(const Matrix& m, const disp_stream_ptr& os, const options& opts)
{
    options opts2 = opts;
    opts2.set(matcl::opt::disp::header_only(true));
    return disp(m, os, opts2);
};

void matcl::disp(const Matrix& m, const disp_stream_ptr& os, const options& opts)
{
    switch(m.get_matrix_code())
    {
        case mat_code::integer_scalar:
            return md::disp_impl(m.get_scalar<Integer>(), os, opts);
        case mat_code::real_scalar:
            return md::disp_impl(m.get_scalar<Real>(), os, opts);
        case mat_code::float_scalar:
            return md::disp_impl(m.get_scalar<Float>(), os, opts);
        case mat_code::complex_scalar:
            return md::disp_impl(m.get_scalar<Complex>(), os, opts);
        case mat_code::float_complex_scalar:
            return md::disp_impl(m.get_scalar<Float_complex>(), os, opts);
        case mat_code::object_scalar:
            return md::disp_impl(m.get_scalar<Object>(), os, opts);
    };

    if (opts.size() == 0)
    {
        return details::matrix_data_accesser::get_ptr(m)->disp(os);
    }
    else
    {
        disp_stream_ptr os2 = std::make_shared<disp_stream_options>(os, opts);
        return details::matrix_data_accesser::get_ptr(m)->disp(os2);
    }
};

namespace details
{
    struct serialize_save_functor : public details::extract_type_switch<void,serialize_save_functor,true>
    {
        template<class T>
        static void eval(const Matrix& h, const T& , oarchive_impl& os)
        {
            using cont_type     = matrix_container<typename T::value_type,typename T::struct_type>;
            const cont_type* c  = static_cast<const cont_type*>(matrix_data_accesser::get_ptr(h));
            os << c;
        };
        template<class T>
        static void eval_scalar(const Matrix& , const T& mat, oarchive_impl& os)
        {
            matcl::details::serialize_save(os,mat,0);
        };
    };

    template<>
    void serialize_save_functor::eval_scalar<Object>(const Matrix& , const Object& mat,
                                                     oarchive_impl& os)
    {
        ti::ti_object ti = mat.get_type();
        matcl::details::serialize_save(os,ti,0);
        matcl::details::serialize_save(os,mat,0);
    };

    template<class T, class S, class Mat>
    static void init(T& owner,S* cont, Mat& mat)
    {
        owner.m_refcount    = mat.get_refstr();
        owner.m_mat_ptr     = cont; 

        if (cont->get_count() > 0)
            owner.m_refcount->increase();

        cont->increase();
    };

    struct serialize_load_functor : public details::extract_type_from_code<Matrix,serialize_load_functor>
    {
        template<class T>
        static Matrix eval(iarchive_impl& is)
        {
            using cont_type = matrix_container<typename T::value_type,typename T::struct_type>;
            cont_type* tmp;
            is >> tmp;

            Matrix out;
            init(matrix_data_accesser::get_base(out).m_value.m_mat,tmp,tmp->get());
            matcl::mat_code mt = matrix_traits::mat_type_info_type<T>::matrix_code;
            matrix_data_accesser::get_base(out).m_type = mt;
            return out;
        };
        template<class T>
        static Matrix eval_scalar(iarchive_impl& is)
        {
            T val;
            matcl::details::serialize_load(is,val,0);
            return val;
        };
    };

    template<>
    Matrix serialize_load_functor::eval_scalar<Object>(iarchive_impl& is)
    {
        ti::ti_object tio;
        matcl::details::serialize_load(is,tio,0);

        Object val(tio);
        matcl::details::serialize_load(is,val,0);
        return val;
    };
};

void matcl::load(iarchive & ar, Matrix& m)
{
    int code;
    ar.get() >> code;

    m = details::serialize_load_functor
            ::make<iarchive_impl&>((matcl::mat_code)code,ar.get());
};

void matcl::save(oarchive & ar, const Matrix& m)
{
    int code = (int)m.get_matrix_code();
    ar.get() << code;
    return details::serialize_save_functor
            ::make<const Matrix&,oarchive_impl&>(m,ar.get());
};

void matcl::convert_mmlib_to_mm(std::istream& mmlib_format, std::ostream& mm_format)
{
    Matrix A;
    std::string comments;
    load(mmlib_format, A, comments);

    mm_save(mm_format, A, comments);
};

void matcl::convert_mm_to_mmlib(std::istream& mm_format, std::ostream& mmlib_format)
{
    Matrix A;
    std::string comments;
    mm_load(mm_format, A, comments);

    save(mmlib_format, A, comments);
};

};
