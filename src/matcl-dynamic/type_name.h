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

#pragma warning(push)
#pragma warning(disable:4127)       //conditional expression is constant

#include <boost/tokenizer.hpp>

#pragma warning(pop)

#include <map>

namespace matcl { namespace dynamic { namespace details
{

class class_name_parser
{
    private:
        using string_map    = std::map<std::string,std::string>;
        using tokenizer     = boost::tokenizer<boost:: char_separator<char> >;
        using tok_iterator  = tokenizer::iterator;

        std::string			parse_templ(tok_iterator& beg,tok_iterator end, bool& simpl, 
                                bool ignore_templates);
        std::string			parse_type(tok_iterator& beg,tok_iterator end, bool& simpl, 
                                bool ignore_templates);		
        std::string			convert_known_types(const std::string& s, bool& simpl);

        static string_map	init_table();

        void				normalize(string_map& smap);

    public:
        std::string			make(const std::string& s, bool& simpl, bool ignore_templates);
};

};};};
