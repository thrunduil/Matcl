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
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-dynamic/type.h"

#include <iostream>

namespace matcl { namespace raw
{

namespace md = matcl::details;

//precision required to represent given type; 
//return nonzero value for floating points only
template<class T>
struct get_stream_precision;

template<>
struct get_stream_precision<Integer>
{
    static Integer eval()   { return 0; };
};
template<>
struct get_stream_precision<Object>
{
    static Integer eval()   { return 0; };
};
template<>
struct get_stream_precision<Real>
{
    static Integer eval()   { return std::numeric_limits<Real>::max_digits10; };
};
template<>
struct get_stream_precision<Float>
{
    static Integer eval()   { return std::numeric_limits<Float>::max_digits10; };
};
template<>
struct get_stream_precision<Complex>
{
    static Integer eval()   { return std::numeric_limits<Real>::max_digits10; };
};
template<>
struct get_stream_precision<Float_complex>
{
    static Integer eval()   { return std::numeric_limits<Float>::max_digits10; };
};

struct MATCL_SCALAR_EXPORT stream_helpers
{
    static void     save_comments(std::ostream& os, const std::string& comments);
    static void     load_comments(std::istream& is, std::string& comments);

    // go to the end of line (or eof)
    static void     skiptoeol(std::istream &is);

    // consume all white characters except new line
    static void     skip_white(std::istream &is);

    //consume all white charecters and comments beginning with % #
    //up to first new line; return true if last consumed charecter is end of line or eof
    static bool     check_nl_comment(std::istream &is);

    //consume all white charecters up to first new line
    //return true if last consumed charecter is end of line of eof
    static bool     check_nl(std::istream &is);

    //consume all white spaces and lines beginning with % #
	static void    skip_all_blank(std::istream&);

    //if line begins with comment_char than all characters up to end of line are 
    //appended to comment argument; if add_newline is true then newline charater 
    //is added at the beginning of comment line; add_newline is set to true
    //otherwise nothing is added to comment and no characters are consumed and
    //add_newline is not changed
    //function return true if line is a comment line
    static bool    read_comment_line(std::istream&, char comment_char, std::ostringstream& comment, bool& add_newline);

    //consume all characted up to end of line and put them to line stream
    //(except end of line)
    static void     read_toeol(std::istream&, std::ostringstream& line);

	static Integer maxwidth(const integer_dense &);

	template<class value_type>
	static bool    read(std::istream&, value_type& val);

    //correct precision must already be set
	template<class value_type>
	static void    write(std::ostream&, const value_type& val);
};

std::ostream& save(std::ostream&, const dynamic::Type& ti);
std::istream& load(std::istream&, dynamic::Type& ti);

};};
