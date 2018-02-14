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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/general/fwd_decls.h"

#include <iostream>

namespace matcl { namespace raw
{

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
struct get_stream_precision<Integer_64>
{
    static Integer eval()   { return 0; };
};
template<>
struct get_stream_precision<dynamic::object>
{
    static Integer eval()   { return 0; };
};

// we need two extra decimal digits to read value back without error
template<>
struct get_stream_precision<Real>
{
    static Integer eval()   { return std::numeric_limits<Real>::max_digits10 + 2; };
};
template<>
struct get_stream_precision<Float>
{
    static Integer eval()   { return std::numeric_limits<Float>::max_digits10 + 2; };
};
template<>
struct get_stream_precision<Complex>
{
    static Integer eval()   { return std::numeric_limits<Real>::max_digits10 + 2; };
};
template<>
struct get_stream_precision<Float_complex>
{
    static Integer eval()   { return std::numeric_limits<Float>::max_digits10 + 2; };
};

struct MATCL_CORE_EXPORT stream_helpers
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
	static void     skip_all_blank(std::istream&);

    //if line begins with comment_char than all characters up to end of line are 
    //appended to comment argument; if add_newline is true then newline charater 
    //is added at the beginning of comment line; add_newline is set to true
    //otherwise nothing is added to comment and no characters are consumed and
    //add_newline is not changed
    //function return true if line is a comment line
    static bool     read_comment_line(std::istream&, char comment_char,
                        std::ostringstream& comment, bool& add_newline);

    //consume all characted up to end of line and put them to line stream
    //(except end of line)
    static void     read_toeol(std::istream&, std::ostringstream& line);

	static Integer  maxwidth(const integer_dense &);

    // if this read returns false, then return value 'val' need not be
    // initialized, additionally failbit is set
	static bool     read(std::istream&, Integer& val);
    static bool     read(std::istream&, Integer_64& val);
    static bool     read(std::istream&, Float& val);
    static bool     read(std::istream&, Real& val);
    static bool     read(std::istream&, Complex& val);
    static bool     read(std::istream&, Float_complex& val);
    static bool     read(std::istream&, std::string& val);
    static bool     read(std::istream&, Object& val);

    //correct precision must be already set in the stream os
	static void     write(std::ostream& os, Integer val);
    static void     write(std::ostream& os, Integer_64 val);
    static void     write(std::ostream& os, Float val);
    static void     write(std::ostream& os, Real val);
    static void     write(std::ostream& os, const Complex& val);
    static void     write(std::ostream& os, const Float_complex& val);
    static void     write(std::ostream& os, const std::string& val);
    static void     write(std::ostream& os, const Object& val);
};

};};
