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

#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/function.h"
#include "matcl-core/details/hash_table/object_table.inl"
#include "matcl-core/details/hash_table/hash_equal.inl"
#include "matcl-core/general/thread.h"
#include "matcl-core/memory/alloc.h"
#include "matcl-dynamic/initialization.h"

namespace matcl { namespace dynamic
{

//-----------------------------------------------------------
//                  identifier_impl
//-----------------------------------------------------------
namespace details
{

class identifier_impl 
{
    private:
        using string_t  = std::string;
        using allocator = default_allocator_simple<true, char>;

        size_t          m_size;
        size_t          m_hash;
        size_t          m_code;
        const char*     m_data;

    public:
        explicit identifier_impl(string_t val);
        identifier_impl(const char* data, size_t size);

        ~identifier_impl();

        const char*     get_data() const    { return m_data; };
        size_t          get_size() const    { return m_size; };
        size_t          get_code() const    { return m_code; };

        bool            is_special() const;

        static size_t   eval_hash(const char* data, size_t size);
        static size_t   eval_hash(const std::string& data);
        std::size_t     hash_value() const;  
        bool            equal(const std::string& data) const;

    private:
        static size_t   get_next_code();
};

inline identifier_impl::identifier_impl(string_t val)
{
    m_hash      = eval_hash(val.data(), val.size());
    size_t N    = val.size();    
    m_size      = N;

    char* tmp   = allocator::malloc((m_size+1)*sizeof(char));

    strcpy_s(tmp, N+1, val.data());
    tmp[N]      = 0;

    m_data      = tmp;
    m_code      = get_next_code();
};

inline identifier_impl::identifier_impl(const char* data, size_t size)
{
    m_hash      = eval_hash(data, size);
    size_t N    = size;    
    m_size      = N;
    
    char* tmp   = allocator::malloc((m_size+1)*sizeof(char));    

    strcpy_s(tmp, N+1, data);
    tmp[N]      = 0;

    m_data      = tmp;
    m_code      = get_next_code();
};

identifier_impl::~identifier_impl()
{
    allocator::free(const_cast<char*>(m_data));
};

size_t identifier_impl::get_next_code()
{
    // only one thread can call this function at this moment
    static size_t code  = 0;
    return code++;
};

inline size_t identifier_impl::eval_hash(const char* data, size_t size)
{
    size_t seed = 0;

    for (size_t i = 0; i < size; ++i)
        boost::hash_combine(seed,data[i]);        

    return seed;
};

inline size_t identifier_impl::eval_hash(const std::string& val)
{
    return eval_hash(val.data(), val.size());
};

inline bool identifier_impl::equal(const std::string& val) const
{
    if (m_size != val.size())
        return false;

    return strcmp(m_data, val.data()) == 0;
};

inline bool identifier_impl::is_special() const
{
    return m_size >= 1 && m_data[0] == '$';
};

inline std::size_t identifier_impl::hash_value() const
{
    return m_hash;
};

//---------------------------------------------------------------
//                  identifier_table
//---------------------------------------------------------------

template<class T>
struct ptr_traits
{
    static void check_free(const T*) {};
    static void copy(const T*) {};
};

struct identifier_ptr
{
    using value_type    = identifier_impl*;
    using type_traits   = ptr_traits<identifier_impl>;

    identifier_impl*    m_ptr;

    static identifier_ptr   make(identifier_impl* p) { return identifier_ptr{p}; };
};

class identifier_table : public matcl_new_delete, global_object
{
    private:        
        using string_hash   = matcl::details::obj_hasher<identifier_impl>;
        using string_equal  = matcl::details::obj_equaler<identifier_impl>;
        using allocator     = matcl::default_allocator_simple<true, char>;
        using string_table  = matcl::details::object_table<identifier_ptr,string_hash,string_equal,
                                allocator>;
        using mutex_type    = matcl::default_spinlock_mutex;
        using lock_type     = std::unique_lock<mutex_type>;

    private:
        string_table        m_table;
        mutex_type          m_mutex;

    public:
        identifier_table();
        ~identifier_table();
        
        void                clear_global() override;
        void                close_global() override;
        identifier_ptr      get(const std::string& val);
};

identifier_ptr identifier_table::get(const std::string& val)
{ 
    lock_type lock(m_mutex);

    return m_table.get(val); 
};

identifier_table::identifier_table()
{
};

void identifier_table::clear_global()
{
    m_table.clear();
};

void identifier_table::close_global()
{
    delete this;
}

identifier_table::~identifier_table()
{}

identifier_table* g_id_table = nullptr;

void details::initialize_identifier_table()
{
    g_id_table = new identifier_table();
};

};

//---------------------------------------------------------------
//                  identifier
//---------------------------------------------------------------
identifier::identifier(const std::string& name)
{
    details::identifier_ptr ptr = details::g_id_table->get(name);
    m_impl  = ptr.m_ptr;    
};

std::string identifier::to_string() const
{
    return m_impl->get_data();
};

size_t identifier::hash_value() const
{ 
    return m_impl->hash_value();
};

size_t identifier::get_unique_code() const
{
    return m_impl->get_code();
};


bool identifier::operator==(identifier other) const
{
    return m_impl == other.m_impl;
};

bool identifier::operator!=(identifier other) const
{
    return m_impl != other.m_impl;
};

bool identifier::operator<(identifier other) const
{
    return m_impl < other.m_impl;
};
bool identifier::operator>(identifier other) const
{
    return m_impl > other.m_impl;
};

//---------------------------------------------------------------
//                  function_name
//---------------------------------------------------------------

bool function_validator::operator==(function_validator o) const
{
    return m_validator == o.m_validator && m_printer == o.m_printer;
};

bool function_validator::operator!=(function_validator o) const
{
    return m_validator != o.m_validator || m_printer != o.m_printer;
};

function_name::function_name(const std::string& name, validator valid)
    :identifier(name), m_validator(valid)
{};

void function_name::disp_requirements(std::ostream& os) const
{
    if (m_validator.get_validator() == nullptr)
    {
        os << "none";
        return;
    }

    m_validator.get_printer()(os);
};

bool function_name::validate_function(function f) const
{
    if (f.is_null() == true)
        return false;

    if (m_validator.get_validator() == nullptr)
        return true;

    return m_validator.get_validator()(f);
};

bool function_name::is_special() const
{
    return m_impl != nullptr && m_impl->is_special();
};

function_name::validator function_name::get_validator() const
{
    return m_validator;
};

};};

