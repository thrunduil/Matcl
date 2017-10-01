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

#pragma once

#include "matcl-scalar/config.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/IO/scalar_io.h"
#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-scalar/details/enablers.h"

namespace matcl
{

namespace md = matcl::details;

namespace details
{
    MATCL_SCALAR_EXPORT 
    void            disp_impl(Integer, const disp_stream_ptr& os, const options& opts);

    MATCL_SCALAR_EXPORT 
    void            disp_impl(const Real&, const disp_stream_ptr& os, const options& opts);

    MATCL_SCALAR_EXPORT 
    void            disp_impl(const Complex&, const disp_stream_ptr& os, const options& opts);

    MATCL_SCALAR_EXPORT 
    void            disp_impl(const Object&, const disp_stream_ptr& os, const options& opts);
};

//--------------------------------------------------------------------
//                      PRETTY PRINTING
//--------------------------------------------------------------------

// Display matrix, string, character, or scalars (Integer, Real, Complex,
// Object) using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
MATCL_SCALAR_EXPORT
void                disp(const std::string& , const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options());

MATCL_SCALAR_EXPORT
void                disp(const char*, const disp_stream_ptr& os,
                         const options& opts);

template<class S>
MATCL_SCALAR_EXPORT
typename md::enable_if_scalar<S, void>::type
                    disp(const S& A, const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options());

// Display matrix, string, character, or scalars (Integer, Real, Complex,
// Object) header only; using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
MATCL_SCALAR_EXPORT
void                disp_header(const std::string& , const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options());
MATCL_SCALAR_EXPORT
void                disp_header(const char*, const disp_stream_ptr& os = default_disp_stream(),
                        const options& opts = options());

template<class S>
MATCL_SCALAR_EXPORT
inline typename md::enable_if_scalar<S, void>::type
                    disp_header(const S& A, const disp_stream_ptr& os = default_disp_stream(),
                            const options& opts = options());

// convert a scalar to a string
template<class S, class Enable = typename md::enable_if_scalar<S, void>::type> 
std::string         to_string(const S& A);

//--------------------------------------------------------------------
//                      SERIALIZATION
//--------------------------------------------------------------------

// Serialize matrices and scalars using boost::serialization library.
template<class T> 
MATCL_SCALAR_EXPORT
typename details::enable_if_scalar<T,void>::type
                    save(oarchive& ar,const T& A);

// Deserialize matrices and scalars.
template<class T> 
MATCL_SCALAR_EXPORT
typename details::enable_if_scalar<T,void>::type
                    load(iarchive& ar, T& A);

//--------------------------------------------------------------------
//                  SAVE AND LOAD IN TEXT FORMAT
//--------------------------------------------------------------------

// Output stream operator for a matrix and scalars. This operator is intendent to
// save matrices in text format, not for pretty printing. However output looks
// reasonably good. If you want to get nice looking output, use disp function
template<class T>
MATCL_SCALAR_EXPORT
typename details::enable_if_matcl_scalar<T,std::ostream&>::type
                    operator<<(std::ostream& os, const T& A);

// Input stream operator for a matrix and scalars. This operator is intendent to
// load matrices in text format, saved previously with operator<<..
template<class T>
MATCL_SCALAR_EXPORT
typename details::enable_if_matcl_scalar<T,std::istream&>::type
                    operator>>(std::istream& is, T& A);

};

#include "matcl-scalar/details/scalar_io.inl"