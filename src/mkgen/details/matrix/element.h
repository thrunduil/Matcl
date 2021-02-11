/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

namespace matcl { namespace mkgen { namespace details
{

//TODO:

//-----------------------------------------------------------------------
//                      forward declarations
//-----------------------------------------------------------------------
template<class Elem>
struct get_offset_elem_step;

//-----------------------------------------------------------------------
//                      element
//-----------------------------------------------------------------------
// This class represents a generic element of a matrix tagged with Tag at row
// row, column col.
template<class Tag, Integer Row, Integer Col>
struct element : public mkd::scalar_data<element<Tag, Row, Col>>
{
    //TODO: add checks
    static const Integer    get_offset      = Tag::get_offset(Row,Col);
    using                   root_align_type = typename Tag::root_align_type;

    using tag       = Tag;
    using this_type = element<Tag, Row, Col>;

    // print element
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        (void)prior;

        Tag::print(os, details::prior_start);
        os << "[" << Row << "," << Col << "]";
    };

    //store result of type Val of computation in Array supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    static void assign(const Val& v, const Local_Storage& ls)
    {
        const_cast<Val&>(ls.get_extern<Tag,Row,Col>()) = v;
    };    

    template<class Val, class Data_Provider, class Temp_Storage>
    inline_lev_1
    static const Val* get_array(const Data_Provider& dp, const Temp_Storage* ev)
    {        
        return Tag::get_data_ptr<Val,0>(dp);
    };

    //get data of type Val associated with this element supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    //static const Val& eval(const Local_Storage& ls) TODO
    static Val eval(const Local_Storage& ls)
    {
        //return ls.get_extern<Tag,Row,Col>();
        return Val(ls.get_extern<Tag,Row,Col>());
    };

    template<class Loop_Storage, class Val, class Local_Storage>
    inline_lev_1
    static void eval_loop(Val& ret, Integer off, const Local_Storage& cont)
    {
        return Loop_Storage::get_value<Val, this_type>(cont, ret, off);
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename list::push_back<Arr_List, details::array_item<this_type,Step * Tag::step,
                            details::array_item_extern>> :: type;

    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };

    template<class Visitor>
    static void accept_assign(Visitor& vis)
    {
        vis.visit_store();
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        vis.visit_load();
    };
};

//-----------------------------------------------------------------------
//                      element_step
//-----------------------------------------------------------------------
template<class Elem, Integer Step>
struct element_step : public mkd::scalar_data<element_step<Elem, Step>>
{
    using this_type = element_step<Elem, Step>;

    // print element
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Elem::print<Subs_Context>(os,prior);
    };

    //store result of type Val of computation in Array supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    static void assign(const Val& v, const Local_Storage& ls)
    {
        Elem::assign<Val>(v,ls);
    };

    static const Integer    get_offset = get_offset_elem_step<Elem>::value;

    template<class Val, class Data_Provider, class Temp_Storage>
    inline_lev_1
    static const Val* get_array(const Data_Provider& dp, const Temp_Storage* ev)
    {     
        return Elem::get_array<Val>(dp,ev);
    };

    //get data of type Val associated with this element supplied by data_provided
    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        return Elem::eval<Val>(ls);
    };

    template<class Loop_Storage, class Val, class Local_Storage>
    inline_expr
    static void eval_loop(Val& ret, Integer off, const Local_Storage& cont)
    {
        return Elem::eval_loop<Loop_Storage,Val>(ret,off,cont);
    };

    template<Integer Step_loc, class Arr_List>
    using get_arrays    = typename Elem::template get_arrays<Step*Step_loc,Arr_List>;

    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };

    template<class Visitor>
    static void accept_assign(Visitor& vis)
    {
        Elem::accept_assign<Visitor>(vis);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        Elem::accept<Visitor>(vis);
    };
};

//-----------------------------------------------------------------------
//                      get_offset_elem_step
//-----------------------------------------------------------------------
template<class Elem>
struct get_offset_elem_step
{};

template<class Tag, Integer R, Integer C>
struct get_offset_elem_step<mkd::element<Tag,R,C>>
{
    static const Integer value = mkd::element<Tag,R,C>::get_offset;
};

template<class Tag, Integer MR, Integer MC, Integer R, Integer C>
struct get_offset_elem_step<get_temporary<Tag,MR,MC,R,C>>
{
    using elem = get_temporary<Tag,MR,MC,R,C>;
    static const Integer value = elem::get_offset;
};

template<class Elem, Integer Step>
struct get_offset_elem_step<mkd::element_step<Elem,Step>>
{
    static const Integer value = Elem::get_offset;
};

//TODO: other?
}}}
