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
#include "matcl-core/details/printer.h"
#include "matcl-core/details/disp_stream_impl.h"
#include "matcl-core/details/io_impl.h"

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

    template<class T, bool Is_compl = is_complex<T>::value, bool Is_obj = is_object<T>::value>
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
            mr::stream_helpers::write<T>(os, A);
            os << " ";
            return os;
        }
        static std::istream& eval_load(std::istream& is, T& A)
        {
            mr::stream_helpers::read<T>(is, A);
            return is;
        }
    };

    template<class T>
    struct saveload_scalar_helper_impl<T, true, false>
    {
        static std::ostream& eval_save(std::ostream& os, const T& A)
        {
            boost::io::ios_flags_saver old_flags(os);
            boost::io::ios_precision_saver old_prec(os);

            int precision = mr::get_stream_precision<T>::eval();

            if (precision > 0)
                os << std::setprecision(precision) << std::scientific;

            os << " ";
            mr::stream_helpers::write<T>(os, A);
            os << " ";
            return os;
        }
        static std::istream& eval_load(std::istream& is, T& A)
        {
            mr::stream_helpers::read<T>(is, A);
            return is;
        }
    };

    template<class T>
    struct saveload_scalar_helper_impl<T, false, true>
    {
        static std::ostream& eval_save(std::ostream& os, const T& A)
        {
            os << " ";
            mr::stream_helpers::write<T>(os, A);
            os << " ";
            return os;
        }
        static std::istream& eval_load(std::istream& is, T& A)
        {
            mr::stream_helpers::read<T>(is, A);
            return is;
        }
    };

    template<class T>
    std::ostream& saveload_scalar_helper<T>::eval_save(std::ostream& os, const T& A)
    {
        return saveload_scalar_helper_impl<T>::eval_save(os, A);
    };
    template<class T>
    std::istream& saveload_scalar_helper<T>::eval_load(std::istream& is, T& A)
    {
        return saveload_scalar_helper_impl<T>::eval_load(is, A);
    }

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

namespace matcl
{

template struct details::saveload_scalar_helper<Integer>;
template struct details::saveload_scalar_helper<Real>;
template struct details::saveload_scalar_helper<Float>;
template struct details::saveload_scalar_helper<Complex>;
template struct details::saveload_scalar_helper<Float_complex>;
template struct details::saveload_scalar_helper<Object>;

};