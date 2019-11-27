/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/general/config.h"

#include <type_traits>

#include <stdexcept>

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients
namespace matcl
{

struct FunEvalBase;
struct FunEvalBaseWithStructure;
namespace details
{

struct FunEvalBase_impl
{
    friend struct FunEvalBase;
    friend struct FunEvalBaseWithStructure;

    virtual ~FunEvalBase_impl() {};
    virtual Matrix eval(const Matrix&) const = 0;
    virtual Matrix get_structure() const = 0;

  private:
    FunEvalBase_impl() {};
};

} // end details

/**
 * Base class for functor without sparsity structure
 * Implement eval method in the derived class
 */    
struct FunEvalBase : public details::FunEvalBase_impl 
{
    // Dummy implementation, should never be called:
    Matrix get_structure() const 
    { 
        assert(0 && "FunEvalBase::get_structure() was called on unstructured MmlibF"); 
        return zeros(0,0);
    }
};
/**
 * Base class for functor with sparsity structure
 * Implement eval and get_structure methods in the derived class
 */    
struct FunEvalBaseWithStructure : public details::FunEvalBase_impl {};
 
/**
 * Class to manage a functor
 */    
class MmlibF
{
    public: 
/**
 * Creates empty MmlibF
 */    
        MmlibF(): my_has_structure(false) {};
         
/**
 * Creates MmlibF to manage a FunEvalBase or FunEvalBaseWithStructure functor
 *
 * @param   _fe     FunEvalBase or FunEvalBaseWithStructure representing the functor
 */    
#if defined(COMPILE_FOR_SWIG_WRAPPER)
        // special constructor for SWIG wrapper
        // NOTE: MmlibF takes over ownership of fe_in!
        MmlibF(FunEvalBase* fe_in)
            : fe(dynamic_cast<details::FunEvalBase_impl*>(fe_in))
            , my_has_structure(false) {}
        MmlibF(FunEvalBaseWithStructure* fe_in)
            : fe(dynamic_cast<details::FunEvalBase_impl*>(fe_in))
            , my_has_structure(true) {}
#else
        template <typename T>
        MmlibF(const T& _fe)
            : fe(new T(_fe))
            , my_has_structure(std::is_base_of<FunEvalBaseWithStructure,T>::value) {}

#endif
        // we want implicit special member functions, comment this out (all resources are managed)
        // MmlibF(const MmlibF& other)...
        // MmlibF& operator=(const MmlibF& other)...
/**
 * Evaluate the underlying functor (function f)
 *
 * @param   x       Column vector, point x to evaluate f at
 * 							  
 * @return  Matrix with the evaluation result of f(x)
 */    
        Matrix operator()(const Matrix& x) const
        {
            if (fe)
            {
                return fe->eval(x);
            }
            else
            {
                // TODO dorobic
                return zeros(0,0);
            }
        }

        bool isEmpty() const
        {
            return fe ? false : true;
        }
         
/**
 * Return the sparsity structure of the functor result
 * 							  
 * @return  (0-1) Matrix representing the sparsity structure
 */    
        Matrix get_structure() const
        {
            assert(has_structure() && "MmlibF::get_structure() called on unstructured callback");

            if(fe) return fe->get_structure();
            else return zeros(0,0); // TODO: how about some exception?
        }

        bool has_structure() const { return my_has_structure; }

    private:
        std::shared_ptr<details::FunEvalBase_impl> fe;
        bool my_has_structure;
};

}

#pragma warning( pop )
