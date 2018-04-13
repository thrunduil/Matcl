/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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
#pragma warning(disable:4100)
#pragma warning(disable:4127)

#include "morfa_petsc/option_setting_utility.h"
#include "petsc.h"

#pragma warning(pop)

#include "petsc_library.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-linalg/linear_eq/petsc_enums.h"

//-------------------------------------------------------------
//                  petsc_options_handler
//-------------------------------------------------------------

namespace matcl { namespace petsc
{

/// convert typed values to string recognized by PETSC
template<class VT>
struct string_converter
{
    virtual ~string_converter(){};
    virtual std::string to_string(const VT& val) const = 0;
};

template<class VT>
struct default_converter : string_converter<VT>
{
    virtual std::string to_string(const VT& val) const override;
};

struct converter_int_bool : string_converter<Integer>
{
    virtual std::string to_string(const Integer& val) const override;
};

struct converter_solver : public string_converter<Integer>
{
    private:
        using code_map  = std::map<std::string, petsc_solver>;

    public:
        virtual std::string to_string(const Integer& val) const override;
        petsc_solver        to_code(const std::string& name) const;

    private:
        const code_map&     get_map() const;
        static code_map     init_map();
};

class converter_precond : public string_converter<Integer>
{
    private:
        using code_map  = std::map<std::string, petsc_precond>;

    public:
        virtual std::string to_string(const Integer& val) const override;
        petsc_precond       to_code(const std::string& name) const;

    private:
        const code_map&     get_map() const;
        static code_map     init_map();
};

class converter_order : public string_converter<Integer>
{
    public:
        virtual std::string to_string(const Integer& val) const override;
};

class converter_shift_type : public string_converter<Integer>
{
    public:
        virtual std::string to_string(const Integer& val) const override;
};

class converter_precond_side : public string_converter<Integer>
{
    public:
        virtual std::string to_string(const Integer& val) const override;
};
class converter_jacobi_type : public string_converter<Integer>
{
    public:
        virtual std::string to_string(const Integer& val) const override;
};

struct converter_refine : public string_converter<Integer>
{
    virtual std::string to_string(const Integer& val) const override;
};

struct converter_fcg_orth : public string_converter<Integer>
{
    virtual std::string to_string(const Integer& val) const override;
};

struct converter_norm_type : public string_converter<Integer>
{
    virtual std::string to_string(const Integer& val) const override;
};

struct converter_optional : public string_converter<Integer>
{
    Integer m_val;
    converter_optional(Integer v)   :m_val(v){};

    virtual std::string to_string(const Integer& val) const override;
};

class option_stack
{
    private:
        _OptionUtility*     m_old;
        bool                m_valid;

    public:
        option_stack()                      :m_old(nullptr), m_valid(false){};
        option_stack(_OptionUtility* old)   :m_old(old), m_valid(true){};
        option_stack(option_stack&& o)      :m_old(o.m_old), m_valid(o.m_valid){o.m_valid = false;};
        ~option_stack();

        void reset_fresh(option_stack&& o)
        {
            m_old   = o.m_old;
            m_valid = o.m_valid;
            o.m_valid = false;
        };

        option_stack(const option_stack&) = delete;
        option_stack& operator=(const option_stack&) = delete;
};

/// interface which can be bound to Petsc world to feed it with options
class petsc_options_handler
{
    public:
        using ptr_type              = petsc_options_handler*;
        using sptr_type             = std::shared_ptr<petsc_options_handler>;
        using options_utility_ptr   = std::shared_ptr<_OptionUtility>;
        using util_type             = _OptionUtility*;

    private:
        const options*      m_opts;
        options*            m_opts_missing;
        options_utility_ptr m_util;

    public:
        petsc_options_handler(const options* opts, options* missing);
        ~petsc_options_handler();

        /// returns true if a particular options is known to the handler, i.e. is defined
        /// name: name of option, val: value of option if present
        /// name_only: if true, then petsc is testing only for existence, value type may be
        ///            incorrect
        template<class T, class PT>
        bool            get_option(const std::string& prefix, const std::string& name, 
                                const std::string& full_name, PT& val) const;
        template<class T>
        bool            has_option(const std::string& prefix, const std::string& name, 
                                const std::string& full_name) const;

        //set as active option handler in petsc
        option_stack    install();
        const options*  get_current_options() const;

        /// implements the petsc_options_handler interface
        static sptr_type make(const options* opts, options* missing);

    private:
        void                make_util();
        static util_type    null_utility_constr();
        static util_type    thread_utility_constr();
        
        void            reset_options(const options* opts);

        friend class scoped_options;
};

class scoped_options
{
    private:
        petsc_options_handler*  m_handler;
        const options*          m_old_options;
        option_stack            m_stack_opt;

    public:
        scoped_options(petsc_options_handler& handler, const options* opts);
        scoped_options::~scoped_options();

        scoped_options(const scoped_options&) = delete;
        scoped_options& operator=(const scoped_options&) = delete;
};

}};

//-------------------------------------------------------------
//                  C defines injection
//-------------------------------------------------------------
extern "C"
{

struct _PetscOptionsListenerC
{
    using ptr_type  = matcl::options*;
    ptr_type    data;
};

struct _PetscOptionsHandlerC
{
    using ptr_type = matcl::petsc::petsc_options_handler::ptr_type;

    ptr_type    data;
};

struct _std_stringC
{
    std::string data;
};

// C-wrapping of option infrastructure
using PetscOptionsListenerC = _PetscOptionsListenerC*;
using PetscOptionsHandlerC  = _PetscOptionsHandlerC*;
using std_stringC           = _std_stringC*;

}

//-------------------------------------------------------------
//                  Petsc options callbacks
//-------------------------------------------------------------
namespace matcl { namespace petsc
{

// Callback functions passed within OptionUtility C-structure to Petsc
// These callbacks get (set) options from (to) mmlib option infrastructure

PetscErrorCode PetscOptionsInt_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscInt,PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsReal_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscReal,PetscReal*,PetscBool *);
PetscErrorCode PetscOptionsScalar_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscScalar,PetscScalar*,PetscBool *);
PetscErrorCode PetscOptionsName_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscBool *);
PetscErrorCode PetscOptionsEnum_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],const char *const*,PetscEnum,PetscEnum*,PetscBool *);
PetscErrorCode PetscOptionsBool_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscBool ,PetscBool *,PetscBool *);
PetscErrorCode PetscOptionsBool2_impl(OptionUtility, const char* pref, const char *name, PetscBool* ret);

PetscErrorCode PetscOptionsString_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],const char[],char*,size_t,PetscBool *);
PetscErrorCode PetscOptionsList_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscFunctionList,const char[],char[],size_t,PetscBool *);
PetscErrorCode PetscOptionsEList_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],const char*const*,PetscInt,const char[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsRealArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscReal[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsScalarArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscScalar[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsEnumArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],const char *const *list,PetscEnum[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsIntArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscInt[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsStringArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],char*[],PetscInt*,PetscBool *);
PetscErrorCode PetscOptionsBoolArray_impl(OptionUtility, const char[], const char[],const char[],
                    const char[],PetscBool [],PetscInt*,PetscBool *);

}}
