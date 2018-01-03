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
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/details/exception_details.h"
#include "matcl-core/IO/output_stream.h"
#include "matcl-core/IO/disp_stream.h"

#include <string>
#include <iosfwd>
#include <vector>
#include <iomanip>
#include <sstream>

namespace matcl  { namespace details
{

class bufor_info
{
    private:
        std::stringstream   m_bufor;
        Integer             m_line_pos;
        Integer             m_size;

    public:
        bufor_info()
            :m_line_pos(0), m_size(0)
        {};

        bufor_info& operator<<(const std::string& str);

        std::stringstream&  get_stream()        { return m_bufor; };
        Integer             get_line_pos() const{ return m_line_pos; };
        Integer             get_size() const    { return m_size; };
        void                clear()             { m_bufor = std::stringstream(); m_line_pos = 0; m_size = 0;};
        bufor_info&         disp(Integer w, const std::string& s);
        void                new_line();
};

class MATCL_CORE_EXPORT printer
{
    private:
        bufor_info*         m_buf;
        bool                m_disp_zero;
        Integer             m_precision;

    public:
        printer(bufor_info* buf)
            : m_buf(buf), m_disp_zero(true), m_precision(5)
        {};

        void                set_disp_zero(bool disp_zero)   { m_disp_zero = disp_zero; }
        void                set_precision(Integer prec)     { m_precision = prec; }
        Integer             get_precision() const           { return m_precision; };

                            //if elem_width > 0, then width is fixed; 
                            //if elem_width < 0, then width is not fixed, and width gives maximal width
                            //if elem_width == 0, then width is not specified
        void                disp_elem(Integer elem_width, const std::string& s, 
                                align_type at, Integer value_pos);
        void                disp_elem(Integer elem_width, Integer r, align_type at, 
                                Integer value_pos);
        void                disp_elem(Integer elem_width, Integer_64 r, align_type at, 
                                Integer value_pos);
        void                disp_elem(Integer elem_width, Real r, align_type at, 
                                Integer value_pos);
        void                disp_elem(Integer elem_width, Float r, align_type at,
                                Integer value_pos);
        void                disp_elem(Integer elem_width, const Complex& r, align_type at,
                                Integer value_pos);
        void                disp_elem(Integer elem_width, const Float_complex& r, 
                                align_type at, Integer value_pos);
        void                disp_elem(Integer elem_width, const dynamic::object& r, 
                                align_type at, Integer value_pos);
        
        static Integer      get_min_width(Real r, Integer precision);
        static std::string  disp_string(Integer elem_width, const std::string& s, align_type at);

    private:
        bufor_info&         get_stream()    { return *m_buf; };

        Integer             get_precision(Integer elem_width, Integer required_symbols, Integer max_precision);
        std::string         disp_real(Integer elem_width, Real r, bool disp_zero, align_type at);
        static void         disp_string_align(std::ostream& os, const std::string& s, Integer w, align_type at);
};

}};