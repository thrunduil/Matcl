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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-core/IO/scalar_io.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/options/matcl_options.h"

namespace matcl
{

namespace md = matcl::details;

//--------------------------------------------------------------------
//                      PRETTY PRINTING
//--------------------------------------------------------------------

// Display matrix, string, character, or scalars (Integer, Real, Complex,
// Object) using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
MATCL_MATREP_EXPORT
void            disp(const Matrix& m, const disp_stream_ptr& os = default_disp_stream(),
                    const options& opts = options());

MATCL_SCALAR_EXPORT
void            disp(const std::string& , const disp_stream_ptr& os,
                     const options& opts);

MATCL_SCALAR_EXPORT
void            disp(const char*, const disp_stream_ptr& os,
                     const options& opts);

template<class S>
MATCL_SCALAR_EXPORT
typename md::enable_if_scalar<S, void>::type
                disp(const S& A, const disp_stream_ptr& os, const options& opts);

// Display any data using the same formatting rules as for matrices. Display
// data using global disp stream or local disp stream. Data are represented by
// abstract class disp_data_provider defined elsewhere. Options controls how printing
// is performed, see options_disp for details
// redefinition of function defined in matcl-core
MATCL_CORE_EXPORT
void            disp(disp_data_provider& data, const disp_stream_ptr& os, 
                    const options& opts);

// Display matrix, string, character, or scalars (Integer, Real, Complex,
// Object) header only; using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
MATCL_MATREP_EXPORT 
void            disp_header(const Matrix& m, const disp_stream_ptr& os = default_disp_stream(),
                    const options& opts = options());

MATCL_SCALAR_EXPORT
void            disp_header(const std::string& , const disp_stream_ptr& os,
                    const options& opts);

MATCL_SCALAR_EXPORT
void            disp_header(const char*, const disp_stream_ptr& os, const options& opts);

template<class S>
MATCL_SCALAR_EXPORT
typename md::enable_if_scalar<S, void>::type
                disp_header(const S& A, const disp_stream_ptr& os, const options& opts);

// Display heder of any data using the same formatting rules as for matrices. Display
// data using global disp stream or local disp stream. Data are represented by
// abstract class disp_data_provider defined elsewhere. Options controls how printing
// is performed, see options_disp for details
// redefinition of function defined in matcl-core
MATCL_CORE_EXPORT
void            disp_header(disp_data_provider& data, const disp_stream_ptr& os,
                    const options& opts);

// convert a scalar to a string
// redefinition of function defined in matcl-scalar
template<class S, class Enable> 
std::string     to_string(const S& A);

//--------------------------------------------------------------------
//                      SERIALIZATION
//--------------------------------------------------------------------

// Serialize matrices and scalars using boost::serialization library.
MATCL_MATREP_EXPORT 
void            save(oarchive& ar,const Matrix& mat);

template<class T>
MATCL_SCALAR_EXPORT
typename details::enable_if_scalar<T,void>::type
                save(oarchive& ar,const T& A);

// Deserialize matrices and scalars.
MATCL_MATREP_EXPORT 
void            load(iarchive& ar,Matrix& mat);

template<class T> 
MATCL_SCALAR_EXPORT
typename details::enable_if_scalar<T,void>::type
                load(iarchive& ar,T& A);

//--------------------------------------------------------------------
//                  SAVE AND LOAD IN TEXT FORMAT
//--------------------------------------------------------------------

// Output stream operator for a matrix and scalars. This operator is intendent to
// save matrices in text format, not for pretty printing. However output looks
// reasonably good. If you want to get nice looking output, use disp function
MATCL_MATREP_EXPORT 
std::ostream&   operator<<(std::ostream&, const Matrix&);

template<class T>
MATCL_SCALAR_EXPORT
typename details::enable_if_matcl_scalar<T,std::ostream&>::type
                operator<<(std::ostream& os, const T& A);

// Input stream operator for a matrix and scalars. This operator is intendent to
// load matrices in text format, saved previously with operator<<..
MATCL_MATREP_EXPORT 
std::istream&   operator>>(std::istream&, Matrix&);

template<class T>
MATCL_SCALAR_EXPORT
typename details::enable_if_matcl_scalar<T,std::istream&>::type
                operator>>(std::istream& is, T& A);

//--------------------------------------------------------------------
//                  SAVE AND LOAD IN TEXT FORMAT
//--------------------------------------------------------------------

// save matrix in matcl text format
MATCL_MATREP_EXPORT 
std::ostream&   save(std::ostream&, const Matrix&);

// save matrix in matcl text format with additional comments
MATCL_MATREP_EXPORT 
std::ostream&   save(std::ostream&, const Matrix&, const std::string& comments);

// load matrix in matcl text format
MATCL_MATREP_EXPORT 
std::istream&   load(std::istream& is, Matrix& m);

// load matrix in matcl text format; additionally return associated
// comments
MATCL_MATREP_EXPORT 
std::istream&   load(std::istream& is, Matrix& m, std::string& comments);

//--------------------------------------------------------------------
//                  EXTERNAL TEXT FORMATS
//--------------------------------------------------------------------

// save matrix m in Matrix Market format
MATCL_MATREP_EXPORT 
std::ostream&   mm_save(std::ostream& os, const Matrix& m);

// save matrix m in Matrix Market format; additionally put comments
MATCL_MATREP_EXPORT 
std::ostream&   mm_save(std::ostream& os, const Matrix& m, 
                    const std::string& comments);

// load matrix in Matrix Market format
MATCL_MATREP_EXPORT 
std::istream&   mm_load(std::istream& is, Matrix& m);

// load matrix in Matrix Market format; additionally return associated
// comments
MATCL_MATREP_EXPORT 
std::istream&   mm_load(std::istream& is, Matrix& m, std::string& comments);

// convert matrix stored in matcl text format to Matrix Market format
// only one matrix is converted; if given stream constains many matrices,
// then this function must be called many times to convert all matrices
MATCL_MATREP_EXPORT 
void            convert_matcl_to_mm(std::istream& matcl_format, std::ostream& mm_format);

// convert matrix stored in Matrix Market format to matcl format
// only one matrix is converted; if given stream constains many matrices,
// then this function must be called many times to convert all matrices
MATCL_MATREP_EXPORT 
void            convert_mm_to_matcl(std::istream& mm_format, std::ostream& matcl_format);

};
