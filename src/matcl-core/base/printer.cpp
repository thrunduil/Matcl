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

#include "matcl-core/details/printer.h"
#include "matcl-core/IO/disp_stream.h"
#include "matcl-core/general/exception.h"
#include "matcl-core/details/integer.h"
#include "matcl-core/details/object_interface.h"

#include <algorithm>

namespace matcl { namespace details
{

namespace mrd = matcl::raw::details;

bufor_info& bufor_info::operator<<(const std::string& str)
{
    m_bufor<<str;
    Integer s       = (Integer)str.size();
    m_line_pos      += s;
    m_size          += s;
    return *this;
};
void bufor_info::new_line()
{
    m_bufor << "\n";
    m_line_pos  = 0;
    m_size      += 1;
};
bufor_info& bufor_info::disp(Integer w, const std::string& s)
{
    if (w < 0)
    {
        w           = -w;
        w           = std::min(w, (Integer)s.size());
        m_bufor     << std::setw(w) << s;
        m_line_pos  += w;
        m_size      += w;
    }
    else
    {
        m_bufor     << std::setw(w) << s;
        m_line_pos  += w;
        m_size      += w;
    };
    return *this;
};

std::string printer::disp_string(Integer w, const std::string& s, align_type at)
{
    std::ostringstream os;
    if (w == 0)
    {
        os << std::setw(s.size()) << s;
        return os.str();
    };

    Integer size	= (Integer)s.size();
    bool fixed		= (w > 0);
    w				= (w > 0)? w : -w;

    if (size > w)
    {
        std::string s2 = s;
        s2.resize(w);
        if (w >= 2)
        {
            s2[w-1] = '.';
            s2[w-2] = '.';
        }
        else
        {
            s2[0]	= '.';
        };

        os << s2;
        return os.str();
    };	

    if (fixed == true)
    {
        disp_string_align(os, s, w, at);
    }
    else
    {
        os << s;
    }

    return os.str();
};
void printer::disp_string_align(std::ostream& os, const std::string& s, Integer w, align_type at)
{
    if ((Integer)s.size() == w)
    {
        os << s;
        return;
    };

    switch(at)
    {
        case align_type::left:
        {
            os << std::setw(w) << std::left << s;
            break;
        }
        case align_type::right:
        {
            os << std::setw(w) << std::right << s;
            break;
        }
        case align_type::center:
        default:
        {
            Integer padding     = std::max(0,w - (Integer)s.size());
            Integer pl          = padding/2;
            Integer pr          = padding - pl;
            os  << std::string(pl, ' ' )
                << s
                << std::string(pr, ' ' );
            break;
        }        
    };

}
void printer::disp_elem(Integer w, const std::string& s, align_type at, Integer)
{
    get_stream() << disp_string(w,s,at);
};
void printer::disp_elem(Integer w, Integer v, align_type at, Integer vp)
{
    std::ostringstream os;

    if (v == 0 && m_disp_zero == false)
    {}
    else
        os << v;

    disp_elem(w,os.str(),at,vp);
};
Integer printer::get_precision(Integer elem_width, Integer off, Integer max_w)
{
    if (elem_width == 0)
    {
        return max_w;
    };

    elem_width	= (elem_width > 0)? elem_width : -elem_width;
    Integer p	= std::min(std::max(elem_width-off,1),max_w);
    return p;
};

std::string printer::disp_real(Integer w, Real v, bool disp_zero, align_type at)
{
    if (std::isnan(v))
        return disp_string(w,"NaN",at);

    if (std::isinf(v))
    {
        if (v > 0)  return disp_string(w,"Inf",at);
        else        return disp_string(w,"-Inf",at);
    };

    //we calculate std::precision required to print the first digit
    //and req_prec following; let m_precision mean total number of 
    //digits
    Integer req_prec    = m_precision - 1;
    Integer min_prec    = 2;
    Integer dprec       = std::max(0,req_prec - min_prec);

    std::ostringstream str;
    if ( v >= 1e100 )
    {
        //[1e100, inf]
        Integer p = get_precision(w,7,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v <= -1e100 )
    {
        //[-inf, -1e100]
        Integer p = get_precision(w,8,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v >= 1e4 )
    {
        //[1e4, 1e100)
        Integer p = get_precision(w,7,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v <= -1e4)
    {
        //(-1e100, -1e4)
        Integer p = get_precision(w,8,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v >= 1e3 )
    {
        //[1e3, 1e4)
        Integer p = get_precision(w,5,std::max(0,-1+dprec));
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v <= -1e3)
    {
        //(-1e4, -1e3]
        Integer p = get_precision(w,6,std::max(0,-1+dprec));
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v >= 1e2 )
    {
        //[1e2, 1e3)
        Integer p = get_precision(w,4,0+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v <= -1e2)
    {
        //(-1e3, -1e2]
        Integer p = get_precision(w,5,0+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v >= 1e1 )
    {
        //[10, 100)
        Integer p = get_precision(w,3,1+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v <= -1e1)
    {
        //(-100, -10]
        Integer p = get_precision(w,4,1+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v > 0  && v <= 1e-3)
    {
        //(0, 0.001]
        Integer p = get_precision(w,7,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v < 0  && v >= -1e-3)
    {
        //[-0.001, 0)
        Integer p = get_precision(w,8,2+dprec);
        str << std::scientific << std::setprecision(p) << v;
    }
    else if ( v > 0   && v <= 1e-2)
    {
        //(0.001, 0.01]
        Integer p = get_precision(w,2,5+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v < 0 && v >= -1e-2)
    {
        //[-0.01, -0.001)
        Integer p = get_precision(w,3,5+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v > 0   && v <= 1e-1)
    {
        //(0.01, 0.1]
        Integer p = get_precision(w,2,4+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v < 0 && v >= -1e-1)
    {
        //[-0.1, -0.01)
        Integer p = get_precision(w,3,4+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v > 0 && v <= 1)
    {
        //(0.1, 1]
        Integer p = get_precision(w,2,3+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v < 0 && v >= -1)
    {
        //[-1, -0.1)
        Integer p = get_precision(w,3,3+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v > 0 )
    {
        //(1, 10]
        Integer p = get_precision(w,2,2+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else if ( v < 0 )
    {
        //[-10, -1)
        Integer p = get_precision(w,3,2+dprec);
        str << std::fixed << std::setprecision(p) << v;
    }
    else
    {
        if (disp_zero == true)
            str << 0;
    };

    //maximum size = 10 + dprec
    return str.str();
};
Integer printer::get_min_width(Real v, Integer precition)
{
    if (std::isnan(v))
        return 3;

    if (std::isinf(v))
    {
        if (v > 0)  return 3;
        else        return 4;
    };

    //we calculate std::precision required to print the first digit
    //and precition following; let precition mean total number of 
    //digits
    precition           = precition - 1;

    Integer min_prec    = 2;
    Integer dprec       = std::max(0,precition - min_prec);

    std::ostringstream str;
    if ( v >= 1e100 )
        return 7 + 2 + dprec;
    else if ( v <= -1e100 )
        return 8 + 2 + dprec;
    else if ( v >= 1e4 )
        return 7 + 2+dprec;
    else if ( v <= -1e4)
        return 8 + 2 + dprec;
    else if ( v >= 1e3 )
        return 5 + std::max(0,-1+dprec);
    else if ( v <= -1e3)
        return 6 + std::max(0,-1+dprec);
    else if ( v >= 1e2 )
        return 4 + 0+dprec;
    else if ( v <= -1e2)
        return 5 + 0 + dprec;
    else if ( v >= 1e1 )
        return 3 + 1+dprec;
    else if ( v <= -1e1)
        return 4+1+dprec;
    else if ( v > 0  && v <= 1e-3)
        return 7 + 2 + dprec;
    else if ( v < 0  && v >= -1e-3)
        return 8 + 2 + dprec;
    else if ( v > 0   && v <= 1e-2)
        return 2 + 5+dprec;
    else if ( v < 0 && v >= -1e-2)
        return 3 + 5+dprec;
    else if ( v > 0   && v <= 1e-1)
        return 2 + 4+dprec;
    else if ( v < 0 && v >= -1e-1)
        return 3 + 4+dprec;
    else if ( v > 0 && v <= 1)
        return 2 + 3 + dprec;
    else if ( v < 0 && v >= -1)
        return 3 + 3 + dprec;
    else if ( v > 0 )
        return 2 + 2 + dprec;
    else if ( v < 0 )
        return 3 + 2 + dprec;
    else
        return 1;
};

void printer::disp_elem(Integer w, Real v, align_type at, Integer vp)
{
    disp_elem(w,disp_real(w,v,m_disp_zero,at),at,vp);
};
void printer::disp_elem(Integer w, Float v, align_type at, Integer vp)
{
    disp_elem(w, Real(v), at, vp);
};

void printer::disp_elem(Integer w, const dynamic::object& v, align_type at, 
                        Integer value_pos)
{
    bool display = true;

    const_object_interface oi = const_object_interface(&v);

    if (m_disp_zero == false)
    {
        if (oi.is_zero() == true)
            display = false;
    };

    if (display == false)
        this->disp_elem(w, " ", at, value_pos);
    else
        oi.disp_elem(*this, w, at, value_pos);
};

void printer::disp_elem(Integer w, const Float_complex& v, align_type at, Integer vp)
{
    Complex tmp(real(v),imag(v));
    disp_elem(w, tmp, at, vp);
};
void printer::disp_elem(Integer w, const Complex& v, align_type at, Integer vp)
{
    Real r = real(v);
    Real i = imag(v);

    if (w == 0)
    {
        if (r == 0. && i == 0. && m_disp_zero == false)
        {
            disp_elem(0,r,at,vp);
            return;
        };

        disp_elem(0,r,at,vp);

        if (i < 0)  get_stream() << " - ";
        else        get_stream() << " + ";

        disp_elem(0, std::abs(i),at,vp);
        get_stream() << "i";
        return;
    };

    if (r == 0. && i == 0. && m_disp_zero == false)
    {
        disp_elem(w,r,at,vp);
        return;
    };

    Integer wa = (w > 0)? w : -w;

    Integer w1 = (wa - 2)/2;
    Integer w0 = wa - w1 - 2;

    std::string str_r	= disp_real(0,r,true,at);
    std::string str_i	= disp_real(0, std::abs(i),true,at);

    Integer tot_size	= (Integer)str_r.size() + (Integer)str_i.size();

    Integer ds          = 2;

    if (w < 0 || tot_size + ds > wa || (Integer)str_r.size() > w0 || (Integer)str_i.size() > w1)
    {
        std::ostringstream os;
        os << str_r;

        if (tot_size + 4 > wa)
        {
            if (i < 0)  os << "-";
            else        os << "+";
        }
        else
        {
            if (i < 0)  os << " - ";
            else        os << " + ";
        };

        os << str_i;
        os << "i";

        disp_elem(w,os.str(),at,vp);
        return;
    };    

    if ((Integer)str_r.size() < w0)
    {
        disp_elem(w0-1,str_r,at,vp);
        get_stream() << " ";
    }
    else
    {
        disp_elem(w0,str_r,at,vp);
    }

    if (i < 0)  get_stream() << "-";
    else        get_stream() << "+";

    if ((Integer)str_i.size() < w1)
    {
        get_stream() << " ";
        str_i = str_i + "i";
        disp_elem(w1-1+1,str_i,at,vp);
    }
    else
    {
        str_i = str_i + "i";
        disp_elem(w1+1,str_i,at,vp);
    };
};

}}
