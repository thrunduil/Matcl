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

#include "matcl-matrep/func/raw/io.h"
#include <vector>
#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/func/raw/mvgen.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
//#include "matcl-matrep/lib_functions/func_unary.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

namespace matcl { namespace raw 
{
 
namespace mrd = matcl::raw::details;

static bool checkeof(std::istream &is)
{
    char c;
    is.get(c);
    is.putback(c);

    if (c == ']')   return true;
    else            return false;
}

// check if we have empty object mark '[]'
static bool checkempty(std::istream &is)
{
    char c = 0;

    if (is.eof())
        return true;

    while (is)
    {
        is.get(c);
        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }

    if (c == ']')
    {
        is.putback(c);
        return true;
    }
    else
    {
        is.putback(c);
    };

    return false;
}

template<class V, class S>
struct writemat_format
{};

template<class V, class S>
struct readmat_format
{};

template <class val_type>
struct writemat_format<val_type,struct_dense>
{
    static std::ostream& eval(std::ostream &os, const Matrix<val_type,struct_dense> &m)
    {
        static const Integer NB = 5;

        Integer r       = m.rows();
        Integer c       = m.cols();
        Integer m_ld    = m.ld();

        os << r << ' ' << c;

        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        int precision = get_stream_precision<val_type>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        const val_type* ptr_m = m.ptr();
        Integer n_add   = 0;

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0, pos = 0; i < r; i += NB)
            {
                Integer NBJ = std::min(NB, r - i);

                for (Integer k = 1; k <= NBJ; ++k, ++pos)
                {
                    if (n_add == 0)
                        os << "\n";
                    else
                        os << ' ';

                    stream_helpers::write(os,ptr_m[pos]);

                    ++n_add;
                    if (n_add == NB)
                        n_add = 0;
                };                    
            }

            ptr_m += m_ld;
        }

        return os;
    };
};

template <class val_type>
struct readmat_format<val_type,struct_dense>
{
    static std::istream& eval(std::istream &is, Matrix<val_type,struct_dense> &mat)
    {
        val_type buf = matcl::details::default_value<val_type>(mat.get_type());

        stream_helpers::skip_all_blank(is);

        Integer mr, mc;

        is >> mr; 
        is >> mc;

        stream_helpers::skip_all_blank(is);

        if (!is.good())
            goto err;

        mat.reset_unique(mr, mc);

        if (checkempty(is))
        {
            if (mc == 0 || mr == 0)
                return is;
            else
                goto err;
        };

        val_type* ptr   = mat.ptr();
        Integer pos     = 0;
        Integer size    = mr * mc;

        while (is && pos < size)
        {
            if (stream_helpers::read(is, buf))
                ptr[pos++]  = buf;                
            else
                goto err;

            if (stream_helpers::check_nl_comment(is))
            {
                if (checkeof(is))
                    break;
            };
        };

        if (pos != mr * mc)
            goto err;

        if (checkeof(is) == false)
            goto err;

        return is;

    err:
        throw error::unable_to_read_matrix();
    };
};

template <class val_type>
struct writemat_format<val_type,struct_sparse>
{
    static std::ostream& eval(std::ostream &os, const Matrix<val_type,struct_sparse> &m)
    {
        Integer c   = m.cols();
        Integer n   = m.nnz();

        os << m.rows() << ' ' << m.cols() << ' ' << m.nnz();

        if (n == 0)
            return os;

        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        int precision = get_stream_precision<val_type>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        const details::sparse_ccs<val_type>& tmp = m.rep();

        const Integer* tmp_c    = tmp.ptr_c();
        const Integer* tmp_r    = tmp.ptr_r();
        const val_type* tmp_x   = tmp.ptr_x();

        for (Integer i = 0; i < c; ++i)
        {
            for (Integer k = tmp_c[i]; k < tmp_c[i + 1]; ++k)
            {
                os << '\n';
                os << (tmp_r[k] + 1) << ' ' << (i + 1) << ' ';
                stream_helpers::write(os,tmp_x[k]);				
            };
        };

        return os;
    };
};

template <class val_type>
struct readmat_format<val_type,struct_sparse>
{
    static std::istream& eval(std::istream &is, Matrix<val_type,struct_sparse> &m)
    {
        //Integer i, k, nnz;
        val_type bufx = matcl::details::default_value<val_type>(m.get_type());        

        m.assign_to_fresh(Matrix<val_type,struct_sparse>(m.get_type()));
        stream_helpers::skip_all_blank(is);

        Integer mr, mc, nz;
        is >> mr;
        is >> mc;
        is >> nz;

        std::vector<Integer> ri(nz);
        std::vector<Integer> ci(nz);
        std::vector<val_type> vx;
        vx.reserve(nz);

        stream_helpers::skip_all_blank(is);

        if (!is.good())
            goto err;

        if (checkempty(is))
        {
            if (nz != 0)
                goto err;

            m.assign_to_fresh(Matrix<val_type,struct_sparse>(m.get_type(),mr,mc));
            return is;
        };

        Integer bufi    = 0;
        Integer bufj    = 0;

        Integer pos     = 0;

        while (is && pos < nz)
        {
            if (stream_helpers::read(is, bufi))
                ri[pos] = bufi;
            else
                goto err;
        
            if (stream_helpers::check_nl_comment(is))
                goto err;

            if (stream_helpers::read(is, bufj))
                ci[pos] = bufj;
            else
                goto err;

            if (stream_helpers::check_nl_comment(is)) 
                goto err;

            if (stream_helpers::read(is, bufx))
                vx.push_back(bufx);
            else
                goto err;

            ++pos;

            if (stream_helpers::check_nl_comment(is))
            {
                if (checkeof(is))
                    break;
            }
            else
            {
                goto err;
            };            
        };

        if (pos != nz)
            goto err;

        if (checkeof(is) == false)
            goto err;

        {
            using Mat   = Matrix<val_type,struct_sparse>;
            Mat ret(m.get_type(), ri.data(), ci.data(), vx.data(), mr, mc, nz);

            m.assign_to_fresh(ret);
        };

        return is;

    err:
        throw error::unable_to_read_matrix();
    };
};

template <class val_type>
struct writemat_format<val_type,struct_banded>
{
    static std::ostream& eval(std::ostream &os, const Matrix<val_type,struct_banded> &m)
    {
        static const Integer NB = 5;

        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        Integer r   = m.rows();
        Integer c   = m.cols();
        Integer fd  = m.first_diag();
        Integer ld  = m.last_diag();

        os << r << ' ' << c << ' '<< fd << ' ' << ld;

        if (r == 0 || c == 0)
            return os;
        
        Integer m_ld            = m.ld();
        Integer n_add           = 0;

        int precision = get_stream_precision<val_type>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        for (Integer d = ld; d >= fd; --d)
        {
            const val_type* ptr_m   = m.rep_ptr() + m.first_elem_diag(d);
            Integer s               = m.diag_length(d);

            for (Integer k = 0; k < s; k += NB)
            {
                Integer NBJ = std::min(NB, s - k);

                for (Integer i = 1; i <= NBJ; ++i)
                {
                    if (n_add == 0)
                        os << "\n";
                    else
                        os << ' ';

                    stream_helpers::write(os,ptr_m[0]);

                    ++n_add;

                    if (n_add == NB)
                        n_add = 0;

                    ptr_m   += m_ld;
                };
            };
        };

        return os;
    };
};

template <class val_type>
struct readmat_format<val_type,struct_banded>
{
    static std::istream& eval(std::istream &is, Matrix<val_type,struct_banded> &mat)
    {		
        stream_helpers::skip_all_blank(is);

        Integer mr, mc, mfd, mld;

        is >> mr; 
        is >> mc;
        is >> mfd;
        is >> mld;

        stream_helpers::skip_all_blank(is);

        val_type Z = matcl::details::default_value<val_type>(mat.get_type());
        val_type buf(Z);        

        using Mat = Matrix<val_type,struct_banded>;
        mat.assign_to_fresh(Mat(mat.get_type(),mr, mc, mfd, mld));

        if (!is.good())
            goto err;
        
        if (checkempty(is))
        {
            if (mc == 0 || mr == 0)
                return is;
            else
                goto err;
        };		

        Integer fd      = mat.first_diag();
        Integer ld      = mat.last_diag();
        Integer mat_ld  = mat.ld();

        for (Integer d = ld; d >= fd; --d)
        {
            Integer s       = mat.diag_length(d);
            val_type* ptr   = mat.rep_ptr() + mat.first_elem_diag(d);

            for (Integer i = 0; i < s; ++i)
            {
                if (stream_helpers::read(is, buf))
                    mrd::reset_helper(ptr[0],buf);
                else
                    goto err;

                stream_helpers::check_nl_comment(is);

                ptr += mat_ld;
            };
        };

        if (checkeof(is) == false)
            goto err;

        return is;

    err:
        throw error::unable_to_read_matrix();
    };
};

template<class val_type,class struct_type>
std::istream& load(std::istream &is, Matrix<val_type,struct_type> &mat)
{
    return readmat_format<val_type, struct_type>::eval(is, mat);
};

std::ostream& save(std::ostream& os, ti::ti_object ti)
{
    os << ti;
    return os;
};

std::istream& load(std::istream& is, ti::ti_object& ti)
{
    is >> ti;
    return is;
};

template<class val_type,class struct_type>
std::ostream& save(std::ostream &os, const Matrix<val_type,struct_type> &m)
{
    return writemat_format<val_type,struct_type>::eval(os, m);
};

//--------------------------------------------------------------------
//                      INSTANTIATIONS
//--------------------------------------------------------------------
template std::ostream& save(std::ostream &os, const Matrix<Integer,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Integer,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Integer,struct_banded> &m);

template std::ostream& save(std::ostream &os, const Matrix<Real,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Real,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Real,struct_banded> &m);

template std::ostream& save(std::ostream &os, const Matrix<Float,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Float,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Float,struct_banded> &m);

template std::ostream& save(std::ostream &os, const Matrix<Complex,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Complex,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Complex,struct_banded> &m);

template std::ostream& save(std::ostream &os, const Matrix<Float_complex,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Float_complex,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Float_complex,struct_banded> &m);

template std::ostream& save(std::ostream &os, const Matrix<Object,struct_dense> &m);
template std::ostream& save(std::ostream &os, const Matrix<Object,struct_sparse> &m);
template std::ostream& save(std::ostream &os, const Matrix<Object,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Integer,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Integer,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Integer,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Real,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Real,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Real,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Float,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Float,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Float,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Complex,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Complex,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Complex,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Float_complex,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Float_complex,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Float_complex,struct_banded> &m);

template std::istream& load(std::istream &is, Matrix<Object,struct_dense> &m);
template std::istream& load(std::istream &is, Matrix<Object,struct_sparse> &m);
template std::istream& load(std::istream &is, Matrix<Object,struct_banded> &m);

};};
