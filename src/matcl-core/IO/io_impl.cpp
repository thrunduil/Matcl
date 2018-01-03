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

#include "matcl-core/details/IO/io_impl.h"
#include "matcl-core/memory/global_objects.h"
#include <vector>
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-core/details/integer.h"
#include "matcl-core/IO/scalar_io.h"
#include "matcl-core/details/object_interface.h"
#include "matcl-core/details/scalfunc_real.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

namespace matcl { namespace raw 
{
 
namespace mrd = matcl::raw::details;

void stream_helpers::load_comments(std::istream& is, std::string& comments)
{
    std::ostringstream str;
    bool add_newline = false;

    while (stream_helpers::read_comment_line(is, '%', str, add_newline) == true)
        ;

    comments = str.str();
};

void stream_helpers::save_comments(std::ostream& os, const std::string& comments)
{
    if (comments.empty())
        return;

    std::istringstream iss(comments);

    std::string line;
    while (std::getline(iss, line))
        os << '%' << line << "\n";
};

// go to the end of line (or eof)
void stream_helpers::skiptoeol(std::istream &is)
{
    char c;

    while (is)
    {
        is.get(c);
        if (c == '\n') break;
    }
};

// skip all lines, in which the first character that is not whitespace
// is '#' or '%'
void stream_helpers::skip_all_blank(std::istream &is)
{
    char c;

    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t' && c != '\n')
        {
            if (c == '#' || c == '%') 
                skiptoeol(is);
            else 
            {
                is.putback(c);
                break;
            }
        }
    }
};

void stream_helpers::skip_white(std::istream &is)
{
    char c;

    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t')
        {
            is.putback(c);
            break;
        }
    }
};

void stream_helpers::read_toeol(std::istream &is, std::ostringstream& line)
{
    char c;

    while (is)
    {
        is.get(c);
        if (c == '\n')
            break;
        else
            line << c;
    }
};

bool stream_helpers::read_comment_line(std::istream& is, char comment_char, std::ostringstream& comments, 
                                       bool& add_newline)
{
    char c;

    while (is)
    {
        is.get(c);

        if (c == comment_char) 
        {
            if (add_newline)
                comments << "\n";

            read_toeol(is,comments);
            add_newline = true;
            return true;
        }
        else 
        {
            is.putback(c);
            break;
        }
    }
    return false;
};

bool stream_helpers::check_nl_comment(std::istream &is)
{
    char c;

    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t')
        {
            if (c == '\n')
                return true;

            if (c == '#' || c == '%' )
            {
                skiptoeol(is);
                return true;
            }

            is.putback(c);
            return false;
        }
    }

    if (is.eof()) 
        return true;
    else
        return false;
}

bool stream_helpers::check_nl(std::istream &is)
{
    char c;

    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t')
        {
            if (c == '\n')
                return true;

            is.putback(c);
            return false;
        }
    }

    if (is.eof()) 
        return true;
    else
        return false;
}

// set failbit and return false
static bool set_failbit(std::istream &is)
{
    is.setstate(std::ios::failbit);
    return false;
};

static bool Iread(std::istream &is, Integer &i)
{
    int c;

    is >> i;

    if (is.fail() || is.bad()) 
        return false;

    if (is.eof()) 
        return true;
 
    c = is.peek();
        
    if (c == ' ' || c == '\t' || c == '\n') 
        return true;
    
    return set_failbit(is);
}

static bool Iread64(std::istream &is, Integer_64 &i)
{
    int c;

    is >> i;

    if (is.fail() || is.bad()) 
        return false;

    if (is.eof()) 
        return true;
 
    c = is.peek();
        
    if (c == ' ' || c == '\t' || c == '\n') 
        return true;
    
    return set_failbit(is);
}

static bool Rread_helper(std::istream &is, Real &r)
{
    char c = 0;
    bool neg = false;

    while (is)
    {
        is.get(c);
        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }

    if (is.eof())
        return set_failbit(is);

    if (c == '-') 
        neg = true; 
    else 
        is.putback(c);

    if (is)
    {
        is.get(c);
        if (is && (c == 'i' || c == 'I'))
        {
            is.get(c);
            
            if (is && (c == 'n' || c == 'N'))
            {
                is.get(c);

                if ((c == 'f') ||(c == 'F'))
                {
                    r = (neg) ? -constants::inf() : constants::inf();
                    return true;
                }
                else 
                {
                    return set_failbit(is);
                }
            }
            else
            {
                return set_failbit(is);
            }
        }
        if (is && (c == 'n' || c == 'N'))
        {
            is.get(c);
            if (is && (c == 'a' || c == 'A'))
            {
                is.get(c);
                if ((c == 'n' || c == 'N') && !neg)
                {
                    r = constants::nan();
                    return true;
                }
                else 
                {
                    return set_failbit(is);
                };
            }
            else
            {
                return set_failbit(is);
            };
        }

        if (c == '-')
            return set_failbit(is); // double -
    }
    else
    {
        return set_failbit(is);
    }

    is.putback(c);
    is >> r;

    if (is.fail() || is.bad())
        return false;
    
    if (neg) 
        r = -r;

    return true;
};

static bool Rread(std::istream &is, Real &r)
{
    Real tmp;
    int c;

    if (Rread_helper(is, tmp))
    {
        if (is.eof())
        {
            r = tmp;
            return true;
        }
        c = is.peek();
        
        if (c == ' ' || c == '\t' || c == '\n')
        {
            r = tmp;
            return true;
        }
    }

    return set_failbit(is);
}

static bool Cread(std::istream &is, Complex &cm)
{
    Real r, i;
    char c = 0;

    while (is)
    {
        is.get(c);
        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }

    if (is.eof())
        return set_failbit(is);

    if (c == '(') // format '(re, im)'
    {
        while (is)
        {
            is.get(c);
            if (c != ' ' && c != '\t')
                break;
        }

        if (!is)
            return set_failbit(is);
        
        if (c == '\n')
            return set_failbit(is); // everything should be in one line

        is.putback(c);
        
        if (!Rread_helper(is, r)) 
            return set_failbit(is);

        while (is)
        {
            is.get(c);
            if (c != ' ' && c != '\t')
                break;
        }

        if (!is)
            return set_failbit(is);
        
        if (c != ',')
            return set_failbit(is);

        while (is)
        {
            is.get(c);
            if (c != ' ' && c != '\t')
                break;
        }

        if (!is)
            return set_failbit(is);

        if (c == '\n')
            return set_failbit(is); // everything should be in one line

        is.putback(c);
        
        if (!Rread_helper(is, i))
            return set_failbit(is);

        while (is)
        {
            is.get(c);
            if (c != ' ' && c != '\t')
                break;
        }

        if (!is)
            return set_failbit(is);

        if (c != ')')
            return set_failbit(is);

        if (is.eof())
        {
            cm = Complex(r, i);
            return true;
        }
        {
            int c2 = is.peek();
            
            if (c2 != ' ' && c2 != '\t' && c2 != '\n')
                return set_failbit(is);

            cm = Complex(r, i);
        };

        return true;
    }

    is.putback(c);

    // format 're im' (complex number as a pair of reals 
    // separated by whitespace or tab)
    if (!Rread_helper(is, r))
        return set_failbit(is);

    while (is)
    {
        is.get(c);
        if (c != ' ' && c != '\t')
            break;
    }

    if (!is)
        return set_failbit(is);
    
    if (c == '\n')
        return set_failbit(is);; // everything should be in one line

    is.putback(c);
    
    if (!Rread_helper(is, i)) 
        return set_failbit(is);

    {
        int c2 = is.peek();
        if (c2 != ' ' && c2 != '\t' && c2 != '\n')
            return set_failbit(is);

        cm = Complex(r, i);
    };

    return true;
}

bool stream_helpers::read(std::istream& is,Integer& x)
{
    return Iread(is,x);
}
bool stream_helpers::read(std::istream& is,Integer_64& x)
{
    return Iread64(is,x);
}

bool stream_helpers::read(std::istream& is,Real& x)
{
    return Rread(is,x);
}

bool stream_helpers::read(std::istream& is,Float& x)
{
    Real tmp;
    bool ret = Rread(is,tmp);
    x = static_cast<Float>(tmp);
    return ret;
}

bool stream_helpers::read(std::istream& is,Complex& x)
{
    return Cread(is,x);
}

bool stream_helpers::read(std::istream& is,Float_complex& x)
{
    Complex tmp;
    bool ret = Cread(is,tmp);
    x = Float_complex(static_cast<Float>(real(tmp)), static_cast<Float>(imag(tmp)));
    return ret;
}

bool stream_helpers::read(std::istream& is,std::string& str)
{
    is >> str;

    if (is.fail() || is.bad())
        return false;
    else
        return true;
}

bool stream_helpers::read(std::istream& is, Object& x)
{
    md::nonconst_object_interface oi(&x);
    oi.read(is);

    if (is.fail() || is.bad())
        return false;
    else
        return true;
}

void stream_helpers::write(std::ostream& os, Integer x)
{
    os << x;
}

void stream_helpers::write(std::ostream& os, Integer_64 x)
{
    os << x;
}

void stream_helpers::write(std::ostream& os, Real val)
{
    namespace mrds = matcl::raw::details::scal_func;

    if (mrds::finite(val) == false)
    {
        if (mrds::isinf(val))
        {
            if (val < 0)
                os << "-Inf";
            else
                os << "Inf";
        }
        else
            os << "NaN";
    }
    else
    {
        os << val;
    };
}

void stream_helpers::write(std::ostream& os, Float val)
{
    namespace mrds = matcl::raw::details::scal_func;

    if (mrds::finite(val) == false)
    {
        if (mrds::isinf(val))
        {
            if (val < 0)
                os << "-Inf";
            else
                os << "Inf";
        }
        else
            os << "NaN";
    }
    else
    {
        os << val;
    };
};

void stream_helpers::write(std::ostream& os,const Complex& x)
{
    write(os,matcl::real(x));
    os << ' ';
    write(os,matcl::imag(x));
}

void stream_helpers::write(std::ostream& os,const Float_complex& x)
{
    write(os,matcl::real(x));
    os << ' ';
    write(os,matcl::imag(x));
};

void stream_helpers::write(std::ostream& os,const std::string& str)
{
    os << str;
}

void stream_helpers::write(std::ostream& os,const Object& x)
{
    md::const_object_interface oi(&x);
    oi.write(os);
}

};};
