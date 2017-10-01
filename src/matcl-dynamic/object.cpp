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

#include "matcl-dynamic/object.h"
#include "matcl-dynamic/details/object_data.inl"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/function.h"
#include "matcl-dynamic/details/type_impl.h"
#include "matcl-dynamic/predefined_type_traits.h"
#include "type_table.h"
#include "matcl-core/details/object_interface.h"

namespace matcl { namespace dynamic
{

//------------------------------------------------------------
//                      object_initializer
//------------------------------------------------------------

struct object_interface : matcl::details::object_interface_impl
{
    bool is_zero(const dynamic::object& v) const override
    {
        return v.is_zero();
    };

    void disp_elem(const dynamic::object& v, matcl::details::printer& p, Integer w, 
                        align_type at, Integer value_pos) const override
    {
        return v.disp(p, w, at, value_pos);
    };

    void read(dynamic::object& v, std::istream& is) const
    {
        is >> v;
    }

    void read_type(dynamic::Type& ty, std::istream& is) const
    {
        is >> ty;
    }
    void write_type(const dynamic::Type& ty, std::ostream& os) const
    {
        os << ty;
    }
};

// nifty counter
static int g_counter = 0;

static void open()
{
    static object_interface oi; 
    set_object_intrface(&oi);
};

static void close()
{
};

object_initializer::object_initializer()
{
    if (g_counter == 0)
        open();

    ++g_counter;
}

object_initializer::~object_initializer()
{
    --g_counter;

    if (g_counter == 0)   
        close();
};

//------------------------------------------------------------
//                      object
//------------------------------------------------------------
object::object()
    :m_data(nullptr)
{};

object::object(const object& other)    
    :m_type(other.m_type),m_data(other.m_data)
{
    if (m_data)
        m_data->increase_refcount();
};

object::object(object&& other)
    :m_type(std::move(other.m_type)), m_data(other.m_data)
{
    other.m_data = nullptr;
};

object::object(Type t, const object& other)
    :m_type(t), m_data(nullptr)
{
    object tmp = convert(t, other);
    swap(*this, tmp);
}

object::object(Type t, object&& other)
    :m_type(t), m_data(nullptr)
{
    object tmp = convert(t, other);
    swap(*this, tmp);
};

object::object(Type t)
    :m_type(t), m_data(nullptr)
{
    if (m_type != Type())
    {
        m_data = details::type_impl::get(m_type)->create();
        m_data->increase_refcount();
    };
};

static void error_no_one(Type t)
{
    details::error_handler eh;
    eh.error_one_not_defined(t);
    eh.report();
};

object object::make_one(Type t)
{
    object ret;
    ret.m_type = t;

    if (t == Type())
        error_no_one(t);

    ret.m_data = details::type_impl::get(t)->create_one();

    if (ret.m_data == nullptr)
        error_no_one(t);

    ret.m_data->increase_refcount();
    return ret;
};

object::object(Type ti, data_type* data)
    :m_type(ti),m_data(data)
{
    if (m_type == Type())
        return;

    if(data == nullptr)
        m_data = details::type_impl::get(m_type)->create();

    m_data->increase_refcount();
};

object::object(bool v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OBool(v));
    swap(*this, tmp);
};

object::object(const char* str)
    :object(str ? std::string(str) : std::string())
{};

object::object(const std::string& s)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OString(s));
    swap(*this, tmp);
};

object::object(std::string&& s)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OString(std::move(s)));
    swap(*this, tmp);
};

object::object(Integer v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OInteger(v));
    swap(*this, tmp);
};

object::object(Float v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OFloat(v));
    swap(*this, tmp);
};

object::object(Real v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OReal(v));
    swap(*this, tmp);
};

object::object(const Float_complex& v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OFloat_complex(v));
    swap(*this, tmp);
};

object::object(const Complex& v)
    :m_data(nullptr)
{
    object tmp = object(dynamic::OComplex(v));
    swap(*this, tmp);
};

object::~object()
{
    if (m_data)
        m_data->decrease_refcount();
};

object& object::operator=(const object& other) &
{
    make_unique();

    Type in1                = this->get_type();
    Type in2                = other.get_type();
    function f              = details::type_table::get()->get_assigner(in1, in2);
    const object* args[]    = {this, &other};
    f.make(2,args);

    return *this;
};

object& object::operator=(object&& other) &
{
    make_unique();

    Type in1                = this->get_type();
    Type in2                = other.get_type();
    function f              = details::type_table::get()->get_assigner(in1, in2);

    object tmp(std::move(other));
    const object* args[]    = {this, &tmp};
    f.make(2,args);

    return *this;
};

object::operator bool() const
{
    return cast_bool(*this);
};

void object::reset(const object& other) &
{
    data_type* ptr = other.m_data;

    if (ptr)
        ptr->increase_refcount();

    if(m_data != nullptr)
        m_data->decrease_refcount();

    m_data  = ptr;
    m_type  = other.m_type;
    return;
};

void object::reset(object&& other) &
{
    swap(*this, other);
    return;
};

object object::clone() const
{
    if (is_null() == true)
        return *this;

    data_type* ptr = details::type_impl::get(m_type)->clone(m_data);
    return object(m_type,ptr);
};

bool object::is_unique() const
{
    if (m_data == nullptr || m_data->is_unique())
        return true;
    else 
        return false;
}

void object::make_unique()
{
    if (m_data == nullptr || m_data->is_unique())
        return;

    data_type* ptr = details::type_impl::get(m_type)->copy(m_data);
    ptr->increase_refcount();
    m_data->decrease_refcount();
    m_data = ptr;
};

void swap(object& o1, object& o2)
{
    std::swap(o1.m_type,o2.m_type);
    std::swap(o1.m_data,o2.m_data);
};

bool object::is_zero() const
{
    if (is_null() == true)
        return true;

    return details::type_impl::get(m_type)->is_zero(m_data);
};

bool object::is_one() const
{
    if (is_null() == true)
        return false;

    return details::type_impl::get(m_type)->is_one(m_data);
};

void object::disp(matcl::details::printer& pr, Integer elem_width, 
                  matcl::align_type at, Integer value_pos) const
{
    return details::type_impl::get(m_type)->disp(m_data,pr,elem_width,at,value_pos);
};

std::string object::to_string(matcl::details::printer& pr) const
{
    return details::type_impl::get(m_type)->to_string(m_data,pr);
};

void object::serialize(oarchive_impl & ar, const unsigned int version) const
{
    m_type.serialize(ar,version);

    if (is_null() == true)
        return;

    if(is_zero())
        details::type_impl::get(m_type)->save(ar, version, nullptr);
    else
        details::type_impl::get(m_type)->save(ar, version, m_data);
};

void object::serialize(iarchive_impl & ar, const unsigned int version)
{    
    m_type.serialize(ar,version);

    if (m_type == Type())
    {
        *this = object(m_type, nullptr);
        return;
    };

    data_type* ptr = details::type_impl::get(m_type)->load(ar,version);

    *this = object(m_type, ptr);
};

std::ostream& dynamic::save_data(std::ostream& os, const object& A)
{
    if (A.is_null() == true)
        return os;

    details::type_impl::get(A.m_type)->save(os, A.m_data);
    return os;
};

std::ostream& dynamic::operator<<(std::ostream& os, const object& A)
{
    os << A.get_type();

    if (A.is_null() == true)
        return os;

    os << " ";
    details::type_impl::get(A.m_type)->save(os, A.m_data);
    return os;
};

std::istream& dynamic::load_data(std::istream& is, object& A)
{
    if (A.is_null() == true)
        return is;

    object::data_type* ptr = details::type_impl::get(A.m_type)->load(is);
    A = object(A.get_type(),ptr);
    return is;
};

std::istream& dynamic::operator>>(std::istream& is, object& A)
{
    is >> A.m_type;

    object::data_type* ptr = nullptr;

    if (A.is_null() == false)
        ptr = details::type_impl::get(A.m_type)->load(is);

    A = object(A.get_type(),ptr);
    return is;
};

MATCL_DYN_EXPORT object dynamic::convert(Type new_type, const object& obj)
{
    Type in1                = obj.get_type();

    if (in1 == new_type)
        return obj;

    function f              = details::type_table::get()->get_converter(new_type, in1, false);
    const object* args[]    = {&obj};

    return f.make(1,args);
};

MATCL_DYN_EXPORT object dynamic::cast(Type new_type, const object& obj)
{
    Type in1                = obj.get_type();

    if (in1 == new_type)
        return obj;

    function f              = details::type_table::get()->get_cast_function(new_type, in1);
    const object* args[]    = {&obj};

    return f.make(1,args);
};

};};
