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

#include "type_name.h"
#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/special_types.h"
#include "matcl-dynamic/exception.h"

namespace matcl { namespace dynamic { namespace details
{

std::string class_name_parser::parse_templ(tok_iterator& beg,tok_iterator end, bool& simpl, 
        bool ignore_templates)
{
    std::string templ_args = "<";

    std::string cur = *beg;
    if (cur != "<")
        throw error::matcl_dynamic_exception("invalid type name format");

    ++beg;
    for(; beg != end; ++beg)
    {        
        std::string name = parse_type(beg,end, simpl, ignore_templates);
        templ_args += name;

        cur     = *beg;
        if  (cur == ">")
        {
            templ_args += ">";
            ++beg;
            break;
        }
        else if (cur == ",")
        {
            templ_args += ",";
            continue;
        }
    };

    if (ignore_templates == true)
        return "";

    return templ_args;
};

std::string class_name_parser::parse_type(tok_iterator& beg,tok_iterator end, 
                                         bool& simpl, bool ignore_templates)
{
    std::string name;
    std::string short_name;
    for(; beg != end; )
    {        
        std::string cur = *beg;
        if (cur == "<")
        {
            std::string templ = parse_templ(beg,end, simpl,ignore_templates);
            name        += templ;
            short_name  += templ;
            name        = convert_known_types(name, simpl);
        }
        else if (cur == ">")
        {
            break;
        }
        else if (cur == ",")
        {
            break;
        }
        else if (cur == " ")
        {
            name += " ";
            ++beg;
        }
        else if (cur == ":")
        {
            name += ":";
            ++beg;
        }
        else
        {
            name        += cur;
            short_name  = cur;
            ++beg;
        };
    };

    bool simpl2 = false;
    name = convert_known_types(name, simpl2);
    if (simpl2)
    {
        simpl = simpl2;
        return name;
    }
    else
    {
        short_name = convert_known_types(short_name, simpl2);
        return short_name;
    };
};

std::string class_name_parser::make(const std::string& s, bool& simpl, bool ignore_templates)
{
    boost::char_separator<char> sep("", " :<>,");        
                
    tokenizer tok(s,sep);

    tok_iterator beg = tok.begin(); 
    tok_iterator end = tok.end();

    std::string name = parse_type(beg,end, simpl, ignore_templates);
    return name;
};

std::string class_name_parser::convert_known_types(const std::string& s, bool& simpl)
{   
    static string_map g_table = init_table();
    static bool initialized = false;

    if (initialized == false)
    {
        initialized = true;
        normalize(g_table);
    };

    using iterator = string_map::const_iterator;
    iterator pos = g_table.find(s);
    if (pos != g_table.end())
    {
        simpl = true;
        return pos->second;
    };

    return s;
};

class_name_parser::string_map class_name_parser::init_table()
{
    string_map ret;

    ret[typeid(bool).name()]                = "bool";
    ret[typeid(unit_type).name()]           = "unit";
    ret[typeid(any_type).name()]            = "any";
    ret[typeid(std::string).name()]         = "string";
    ret[typeid(Integer).name()]             = "Integer";    
    ret[typeid(Real).name()]                = "Real";
    ret[typeid(Float).name()]               = "Float";
    ret[typeid(Complex).name()]             = "Complex";
    ret[typeid(Float_complex).name()]       = "Float_complex";
    ret[typeid(object).name()]              = "Template";

    //TODO:
    //ret[typeid(Matrix).name()]              = "Matrix";

    return ret;
};

void class_name_parser::normalize(string_map& smap)
{
    using iterator      = string_map::iterator;
    using value_type    = string_map::value_type;

    iterator pos = smap.begin();
    while (pos != smap.end())
    {
        const std::string& key = pos->first;
        const std::string& value = pos->second;

        bool simpl = false;
        std::string key_simpl = make(key,simpl,false);

        if (key_simpl != value)
            smap[key_simpl] = value;

        ++pos;
    };
};

};};};