#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen
{

template<class Derived>
class visitor
{
    public:
        void    visit_plus()        { return static_cast<Derived*>(this)->visit_plus(); };
        void    visit_minus()       { return static_cast<Derived*>(this)->visit_minus(); };
        void    visit_mult()        { return static_cast<Derived*>(this)->visit_mult(); };
        void    visit_store()       { return static_cast<Derived*>(this)->visit_store(); };
        void    visit_load()        { return static_cast<Derived*>(this)->visit_load(); };
        void    visit_call_inplace(){ return static_cast<Derived*>(this)->visit_call_inplace(); };
};

class op_count : public visitor<op_count>
{
    private:
        Integer m_plus;
        Integer m_mult;
        Integer m_load;
        Integer m_store;

    public:
        op_count()
            :m_plus(0), m_mult(0), m_load(0), m_store(0)
        {};

        void    visit_plus()        { ++m_plus; };
        void    visit_minus()       { ++m_plus; };
        void    visit_mult()        { ++m_mult; };
        void    visit_store()       { ++m_store; };
        void    visit_load()        { ++m_load; };
        void    visit_call_inplace(){ ; };

        void print(std::ostream& os)
        {
            os << "plus: "      << m_plus
               << ", mult: "    << m_mult
               << ", load: "    << m_load
               << ", store: "   << m_store
               << "\n";
        };
};

}}