/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/IO/scalar_io.h"
#include "matcl-core/details/IO/printer.h"
#include "matcl-core/details/IO/disp_stream_impl.h"
#include "matcl-core/details/IO/io_impl.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

namespace matcl { namespace details
{
    namespace mr = matcl :: raw;

    template<class T>
    struct to_string_scalar_helper_impl
    {
        static std::string eval(const T& v)
        {
            md::bufor_info buf;
            md::printer pr(&buf);

            pr.set_precision(options().get_option<Integer>(opt::disp::precision()));
            pr.disp_elem(0, v, align_type::left, 0);

            return buf.get_stream().str();
        };
    };

    template<class T, bool Is_compl = is_complex<T>::value>
    struct saveload_scalar_helper_impl
    {
        static std::ostream& eval_save(std::ostream& os, const T& A)
        {
            boost::io::ios_flags_saver old_flags(os);
            boost::io::ios_precision_saver old_prec(os);

            int precision = mr::get_stream_precision<T>::eval();

            if (precision > 0)
                os << std::setprecision(precision) << std::scientific;

            os << " ";
            mr::stream_helpers::write(os, A);
            os << " ";
            return os;
        }
        static std::istream& eval_load(std::istream& is, T& A)
        {
            mr::stream_helpers::read(is, A);
            return is;
        }
    };

    template<class T>
    struct saveload_scalar_helper_impl<T, true>
    {
        static std::ostream& eval_save(std::ostream& os, const T& A)
        {
            boost::io::ios_flags_saver old_flags(os);
            boost::io::ios_precision_saver old_prec(os);

            int precision = mr::get_stream_precision<T>::eval();

            if (precision > 0)
                os << std::setprecision(precision) << std::scientific;

            os << " ";
            mr::stream_helpers::write(os, A);
            os << " ";
            return os;
        }
        static std::istream& eval_load(std::istream& is, T& A)
        {
            mr::stream_helpers::read(is, A);
            return is;
        }
    };

    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, Integer A)
    {
        return saveload_scalar_helper_impl<Integer>::eval_save(os, A);
    };
    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, Float A)
    {
        return saveload_scalar_helper_impl<Float>::eval_save(os, A);
    };
    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, Real A)
    {
        return saveload_scalar_helper_impl<Real>::eval_save(os, A);
    };
    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, const Complex& A)
    {
        return saveload_scalar_helper_impl<Complex>::eval_save(os, A);
    };
    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, const Float_complex& A)
    {
        return saveload_scalar_helper_impl<Float_complex>::eval_save(os, A);
    };
    /*
    std::ostream& saveload_scalar_helper::eval_save(std::ostream& os, const Object& A)
    {
        return saveload_scalar_helper_impl<Object>::eval_save(os, A);
    };
    */

    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Integer& A)
    {
        return saveload_scalar_helper_impl<Integer>::eval_load(is, A);
    }
    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Real& A)
    {
        return saveload_scalar_helper_impl<Real>::eval_load(is, A);
    }
    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Float& A)
    {
        return saveload_scalar_helper_impl<Float>::eval_load(is, A);
    }
    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Complex& A)
    {
        return saveload_scalar_helper_impl<Complex>::eval_load(is, A);
    }
    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Float_complex& A)
    {
        return saveload_scalar_helper_impl<Float_complex>::eval_load(is, A);
    }
    /*
    std::istream& saveload_scalar_helper::eval_load(std::istream& is, Object& A)
    {
        return saveload_scalar_helper_impl<Object>::eval_load(is, A);
    }
    */

    std::string details::to_string_scalar_helper::eval(Integer v)
    {
        return to_string_scalar_helper_impl<Integer>::eval(v);
    };
    std::string details::to_string_scalar_helper::eval(Float v)
    {
        return to_string_scalar_helper_impl<Float>::eval(v);
    };
    std::string details::to_string_scalar_helper::eval(Real v)
    {
        return to_string_scalar_helper_impl<Real>::eval(v);
    };
    std::string details::to_string_scalar_helper::eval(const Complex& v)
    {
        return to_string_scalar_helper_impl<Complex>::eval(v);
    };
    std::string details::to_string_scalar_helper::eval(const Float_complex& v)
    {
        return to_string_scalar_helper_impl<Float_complex>::eval(v);
    };
    std::string details::to_string_scalar_helper::eval(const Object& v)
    {
        return to_string_scalar_helper_impl<Object>::eval(v);
    };

}};
