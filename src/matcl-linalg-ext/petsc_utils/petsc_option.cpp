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

#if 0
TODO

#include "petsc_option.h"
//#include "petscviewer.h"
#include "matcl-linalg/general/mmlib_petsc_exception.h"
#include "matcl-linalg/options/options_petsc.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/lib_functions/func_unary.h"

#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant

//-------------------------------------------------------
//              options handling
//-------------------------------------------------------
namespace matcl { namespace petsc
{

template<class T> struct make_option_type{};
template<>        struct make_option_type<const char*>  { using type = std::string; };
template<>        struct make_option_type<std::string>  { using type = std::string; };
template<>        struct make_option_type<Integer>      { using type = Integer; };
template<>        struct make_option_type<PetscReal>    { using type = Real; };
template<>        struct make_option_type<bool>         { using type = bool; };
template<>        struct make_option_type<PetscBool>    { using type = bool; };

template<class T> 
struct make_option_type<std::vector<T>>
{ 
    using TR    = typename make_option_type<T>::type;
    using type  = std::vector<TR>; 
};

template<class T>
struct convert_opt_value
{
    template<class S>
    static T eval(const S& val)
    {
        return T(val);
    };
};

template<>
struct convert_opt_value<std::string>
{
    template<class S>
    static std::string eval(const S& val)
    {
        return std::string(val);
    };
    static std::string eval(const char* val)
    {
        return (val == nullptr) ? "" : std::string(val);
    };
    static std::string eval(char* val)
    {
        return (val == nullptr) ? "" : std::string(val);
    };
};

template<>
struct convert_opt_value<bool>
{
    static bool eval(bool val)
    {
        return val;
    };
    static bool eval(PetscBool val)
    {
        return (val == PETSC_TRUE) ? true : false;
    };
};

template<class T>
struct convert_opt_string
{
    static std::string eval(const T& val)
    {
        std::ostringstream str;
        str << val;
        return str.str();
    }
};
template<>
struct convert_opt_string<std::string>
{
    static std::string eval(const std::string& val)
    {
        return val;
    }
};
template<>
struct convert_opt_string<bool>
{
    static std::string eval(bool val)
    {
        return val ? "true" : "false";
    }
};
template<class T>
struct convert_opt_string<std::vector<T>>
{
    static std::string eval(const std::vector<T>& val)
    {
        size_t N = val.size();

        if (N == 0)
            return "";

        std::ostringstream str;

        str << val[0];

        for (size_t i = 1; i < N; ++i)
            str << "," << val[i];

        return str.str();
    }
};

template<class T>
struct copy_to_petsc_array
{
    static void eval(T* ret, PetscInt* n, const std::vector<T>& ret_val)
    {        
        Integer s   = std::min<Integer>(*n, ret_val.size());

        for (Integer i = 0; i < s; ++i)
            ret[i] = ret_val[i];

        *n = ret_val.size();
    };
};

PetscEnum convert_string_enum(const std::string& str_val, const char *const* list, bool& found)
{
    Integer pos = 0;
    found       = false;

    while(*list)
    {
        if (strcmp(str_val.c_str(), *list) == 0)
        {
            found = true;
            return (PetscEnum)pos;
        }

        ++pos;
        ++list;
    };

    return (PetscEnum)0;
};
PetscInt convert_string_enum(const std::string& str_val, const char *const* list, Integer size, bool& found)
{
    Integer pos = 0;
    found       = false;

    while (pos < size)
    {
        if (strcmp(str_val.c_str(), *list) == 0)
        {
            found = true;
            return (PetscInt)pos;
        }

        ++pos;
        ++list;
    };

    return (PetscInt)0;
};

PetscErrorCode copy_to_petsc_array_enum(PetscEnum values[], PetscInt* n, const std::vector<std::string>& ret, 
                                        const char *const* list, const char* pref, const char* opt)
{
    Integer s = std::min<Integer>(*n, ret.size());

    for (Integer i = 0; i < s; ++i)
    {
        bool found;
        PetscEnum val = convert_string_enum(ret[i], list, found);

        if (!found)
            SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Unknown option %s for -%s%s",ret[i].c_str(),pref ? pref : "",opt+1);
        
        values[i] = val;
    };

    return 0;
};

void copy_to_petsc_array_string(char values[], Integer len, const std::string& str)
{
    Integer s = std::min<Integer>(len - 1, str.size());
    
    for (Integer i = 0; i < s; ++i)
        values[i] = str[i];

    values[s] = 0;
};

// Generic function to handle options getting-setting between Morfa and Petsc
PetscErrorCode petsc::PetscOptionsBool2_impl(OptionUtility util, const char pref[], const char opt[], PetscBool *set)
{
    std::string str_pref        = pref ? std::string(pref) : "";
    std::string option_name     = std::string(opt);

    std::string full_opt_name   = std::string("-") + str_pref + std::string(opt + 1);
    
    using opt_type  = make_option_type<PetscBool>::type;

    *set = PETSC_FALSE;

    // setting options to Petsc
    if (util->handler->data)
    {
        if ( util->handler->data->has_option<opt_type>(str_pref, option_name, full_opt_name) == false )
        {
            if (util->listener->data)
            {                
                string_option<opt_type> opt_1(full_opt_name, false, opt);
                util->listener->data->set(opt_1);
            }
        }
        else
        {            
            *set = PETSC_TRUE;
        }
    }

    return 0;
}   

template <class T>
PetscErrorCode petsc_options_impl(OptionUtility util, const char pref[], const char opt[], const char descr[], T dfault,
                                  T* ret, PetscBool* set)
{
    std::string str_pref        = pref ? std::string(pref) : "";
    std::string option_name     = std::string(opt);

    std::string full_opt_name   = std::string("-") + str_pref + std::string(opt + 1);
    
    using opt_type  = typename make_option_type<T>::type;

    *set = PETSC_FALSE;

    // setting options to Petsc
    if (util->handler->data)
    {
        std::string option_value;

        T ret_val;
        if ( util->handler->data->get_option<opt_type>(str_pref, option_name, full_opt_name, ret_val) == false )
        {
            if (util->listener->data)
            {                
                opt_type val    = convert_opt_value<opt_type>::eval(dfault);

                string_option<opt_type> opt_1(full_opt_name, val, descr);
                util->listener->data->set(opt_1);
            }
        }
        else
        {            
            *set = PETSC_TRUE;
            *ret = ret_val;
        }
    }

    return 0;
}

template <class T>
PetscErrorCode petsc_options_test_impl(OptionUtility util, const char pref[], const char opt[], const char descr[], PetscBool* set)
{
    std::string str_pref        = pref ? std::string(pref) : "";
    std::string option_name     = std::string(opt);

    std::string full_opt_name   = std::string("-") + str_pref + std::string(opt + 1);
    
    using opt_type  = typename make_option_type<T>::type;

    *set = PETSC_FALSE;

    // setting options to Petsc
    if (util->handler->data)
    {
        if ( util->handler->data->has_option<opt_type>(str_pref, option_name, full_opt_name) == false )
        {
            if (util->listener->data)
            {                
                opt_type val    = convert_opt_value<opt_type>::eval(T());

                string_option<opt_type> opt_1(full_opt_name, val, descr);
                util->listener->data->set(opt_1);
            }
        }
        else
        {            
            *set = PETSC_TRUE;
        }
    }

    return 0;
}

template <class T, class TR = typename make_option_type<T>::type>
PetscErrorCode petsc_options_array_impl(OptionUtility util, const char pref[], const char opt[], const char descr[], 
                    const std::vector<TR>& dfault, T* ret, PetscInt* n, PetscBool* set)
{
    std::string str_pref        = pref ? std::string(pref) : "";
    std::string option_name     = std::string(opt);

    std::string full_opt_name   = std::string("-") + str_pref + std::string(opt + 1);
    
    using opt_type  = typename make_option_type<std::vector<T>>::type;

    *set = PETSC_FALSE;

    // setting options to Petsc
    if (util->handler->data)
    {
        std::string option_value;

        std::vector<T> ret_val;
        if ( util->handler->data->get_option<opt_type>(str_pref, option_name, full_opt_name, ret_val) == false )
        {
            if (util->listener->data)
            {                
                opt_type val    = convert_opt_value<opt_type>::eval(dfault);

                string_option<opt_type> opt_1(full_opt_name, val, descr);
                util->listener->data->set(opt_1);
            }
        }
        else
        {            
            *set = PETSC_TRUE;
            copy_to_petsc_array<T>::eval(ret, n, ret_val);
        }
    }

    return 0;
}

template <class T, class TR = typename make_option_type<T>::type>
PetscErrorCode petsc_options_array_str_impl(OptionUtility util, const char pref[], const char opt[], const char descr[], 
                    const std::vector<TR>& dfault, std::vector<std::string>& ret, PetscBool* set)
{
    std::string str_pref        = pref ? std::string(pref) : "";
    std::string option_name     = std::string(opt);

    std::string full_opt_name   = std::string("-") + str_pref + std::string(opt + 1);
    
    using opt_type  = typename make_option_type<std::vector<T>>::type;

    *set = PETSC_FALSE;

    // setting options to Petsc
    if (util->handler->data)
    {
        std::vector<T> ret_val;
        if ( util->handler->data->get_option<opt_type>(str_pref, option_name, full_opt_name, ret_val) == false )
        {
            if (util->listener->data)
            {                
                opt_type val    = convert_opt_value<opt_type>::eval(dfault);

                string_option<opt_type> opt_1(full_opt_name, val, descr);
                util->listener->data->set(opt_1);
            }
        }
        else
        {            
            *set = PETSC_TRUE;
            ret = ret_val;
        }
    }

    return 0;
}

template <class T>
struct default_list
{
    using TR        = typename make_option_type<T>::type;
    using vec_type  = std::vector<TR>;

    static vec_type eval(T values[], int n)
    {
        vec_type def_values;

        for (int i = 0; i < n; i++)
            def_values.push_back(convert_opt_value<TR>::eval(values[i]));

        return def_values;
    }
};

template <>
struct default_list<const char*>
{
    using vec_type  = std::vector<std::string>;

    static vec_type eval(char* values[], const int n)
    {
        vec_type def_values;

        for (int i = 0; i < n; i++)
            def_values.push_back(values[i] != nullptr ? values[i] : "");

        return def_values;
    }
};

template <>
struct default_list<PetscEnum>
{
    using TR        = std::string;
    using vec_type  = std::vector<TR>;

    static vec_type eval(PetscEnum values[], int n, const char *const *list)
    {
        vec_type def_values;

        for (int i = 0; i < n; i++)
        {
            PetscEnum val       = values[i];
            const char* name    = list[val];
            def_values.push_back(name? std::string(name) : "");
        };

        return def_values;
    }
};

PetscErrorCode petsc::PetscOptionsEnum_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], const char *const* list, PetscEnum dfault, PetscEnum* value, 
                        PetscBool *set)
{
    (void)set;
    (void)man;
    (void)value;

    const char * dfault_element = list[dfault];

    std::string def_val = dfault_element? std::string(dfault_element) : "";
    std::string str_val;
    petsc_options_impl<std::string>(util, pref, opt, descr, def_val, &str_val, set);
    
    if (*set == PETSC_FALSE)
        return 0;

    bool found;
    *value = convert_string_enum(str_val, list, found);

    if (!found)
        SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Unknown option %s for -%s%s",str_val.c_str(),pref ? pref : "",opt+1);

    return 0;
}

PetscErrorCode petsc::PetscOptionsEList_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], const char*const* list, PetscInt ntext, const char dfault[],
                        PetscInt* value, PetscBool *set)
{
    (void)man;

    std::string str_val;
    petsc_options_impl<std::string>(util, pref, opt, descr, dfault, &str_val, set);

    if (*set == PETSC_FALSE)
        return 0;

    bool found;
    *value = convert_string_enum(str_val, list, ntext, found);

    if (!found)
        SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Unknown option %s for -%s%s",str_val.c_str(),pref ? pref : "",opt+1);

    return 0;
}

PetscErrorCode petsc::PetscOptionsInt_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscInt dfault, PetscInt* value, PetscBool *set)
{
    (void)man;

    return petsc_options_impl<PetscInt>(util, pref, opt, descr, dfault, value, set);
}
PetscErrorCode petsc::PetscOptionsReal_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscReal dfault, PetscReal* value, PetscBool *set)
{
    (void)man;

    return petsc_options_impl(util, pref, opt, descr, dfault, value, set);
}
PetscErrorCode petsc::PetscOptionsScalar_impl(OptionUtility util, const char pref[], const char opt[], const char descr[],
                        const char arg2[], PetscScalar dfault, PetscScalar* value, PetscBool *set)
{
    (void)arg2;

    return petsc_options_impl(util, pref, opt, descr, dfault, value, set);
}
PetscErrorCode petsc::PetscOptionsName_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscBool *set)
{
    (void)set;
    (void)man;

    // NOTE: for this kind of option, when option absent is "", when present may hold some value
    return petsc_options_test_impl<const char*>(util, pref, opt, descr, set);
}

PetscErrorCode petsc::PetscOptionsBool_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscBool dfault, PetscBool *flg, PetscBool *set)
{
    (void)man;

     // petsc understands strings like "true", "false"
    return petsc_options_impl<PetscBool>(util, pref, opt, descr, dfault, flg, set);
}

PetscErrorCode petsc::PetscOptionsString_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], const char dfault[], char* value, size_t len, PetscBool *set)
{
    (void)man;

    std::string ret;    
    petsc_options_impl<std::string>(util, pref, opt, descr, dfault ? dfault : "", &ret, set);

    if (*set == PETSC_FALSE)
        return 0;

    copy_to_petsc_array_string(value, len, ret);
    return 0;
}

PetscErrorCode petsc::PetscOptionsList_impl(OptionUtility util, const char pref[], const char opt[], const char descr[],
                        const char man[], PetscFunctionList list, const char dfault[], char value[],
                        size_t len, PetscBool *set)
{
    (void)man;
    (void)list;

    std::string ret;    
    petsc_options_impl<std::string>(util, pref, opt, descr, dfault ? dfault : "", &ret, set);

    if (*set == PETSC_FALSE)
        return 0;

    copy_to_petsc_array_string(value, len, ret);
    return 0;
}

PetscErrorCode petsc::PetscOptionsIntArray_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscInt values[], PetscInt* n, PetscBool *set)
{
    (void)set;
    (void)man;

    return petsc_options_array_impl<PetscInt>(util, pref, opt, descr, default_list<PetscInt>::eval(values, *n), values, n, set);
}

PetscErrorCode petsc::PetscOptionsRealArray_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscReal values[], PetscInt *n, PetscBool *set)
{
    (void)set;
    (void)man;

    return petsc_options_array_impl<PetscReal>(util, pref, opt, descr, default_list<PetscReal>::eval(values, *n), values, n, set);
}

PetscErrorCode petsc::PetscOptionsScalarArray_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscScalar values[], PetscInt *n, PetscBool *set)
{
    (void)set;
    (void)man;

    return petsc_options_array_impl<PetscScalar>(util, pref, opt, descr, 
                default_list<PetscScalar>::eval(values, *n), values, n, set);
}

PetscErrorCode petsc::PetscOptionsBoolArray_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[], PetscBool values [], PetscInt* n, PetscBool *set)
{
    (void)set;
    (void)man;

    return petsc_options_array_impl<PetscBool>(util, pref, opt, descr, 
                default_list<PetscBool>::eval(values, *n), values, n, set);
}
PetscErrorCode petsc::PetscOptionsEnumArray_impl(OptionUtility util, const char pref[], const char opt[] ,const char descr[],
                        const char man[],const char *const *list, PetscEnum values[], PetscInt *n, PetscBool *set)
{
    (void)man;

    std::vector<std::string> ret;
    petsc_options_array_str_impl<std::string>(util, pref, opt, descr, 
                default_list<PetscEnum>::eval(values, *n, list), ret, set);

    if (*set == PETSC_FALSE)
        return 0;

    copy_to_petsc_array_enum(values, n, ret, list, pref, opt);
    return 0;
}

PetscErrorCode petsc::PetscOptionsStringArray_impl(OptionUtility util, const char pref[], const char opt[], const char descr[],
                        const char man[], char* values[], PetscInt* n, PetscBool *set)
{
    (void)man;

    std::vector<std::string> ret;
    petsc_options_array_str_impl<std::string>(util, pref, opt, descr, 
                default_list<const char*>::eval(values, *n), ret, set);

    if (*set == PETSC_FALSE)
        return 0;

    Integer s   = std::min<Integer>(*n, ret.size());

    for (Integer i = 0; i < s; ++i)
    {
        PetscErrorCode ierr = PetscStrallocpy(ret[i].c_str(), &values[i]);CHKERRQ(ierr);
    };

    return 0;
}

//-------------------------------------------------------
//              petsc_options_handler
//-------------------------------------------------------
template<class VT>
std::string default_converter<VT>::to_string(const VT& val) const
{
    std::string val_str = convert_opt_string<VT>::eval(val);
    return val_str;
};

std::string converter_int_bool::to_string(const Integer& val) const
{
    std::string val_str = convert_opt_string<bool>::eval(val != 0? true : false);
    return val_str;
};

std::string converter_solver::to_string(const Integer& val) const
{
    petsc_solver solver  = (petsc_solver)val;

    switch(solver)
    {            
        case petsc_solver::preonly:      return "preonly";  
        case petsc_solver::cg:           return "cg";
        case petsc_solver::fcg:          return "fcg";
        case petsc_solver::cgne:         return "cgne";
        case petsc_solver::gmres:        return "gmres";
        case petsc_solver::richardson:   return "richardson";
        case petsc_solver::chebyshev:    return "chebyshev";
        case petsc_solver::tfqmr:        return "tfqmr";
        case petsc_solver::tcqmr:        return "tcqmr";
        case petsc_solver::bcgs:         return "bcgs";
        case petsc_solver::bcgsl:        return "bcgsl";
        case petsc_solver::bicg:         return "bicg";
        case petsc_solver::cgs:          return "cgs";
        case petsc_solver::fgmres:       return "fgmres";
        case petsc_solver::lgmres:       return "lgmres";
        case petsc_solver::gcr:          return "gcr";
        case petsc_solver::cr:           return "cr";
        case petsc_solver::minres:       return "minres";   
        case petsc_solver::symmlq:       return "symmlq";
        case petsc_solver::lsqr:         return "lsqr";
        case petsc_solver::dgmres:       return "dgmres";
        case petsc_solver::ibcgs:        return "ibcgs";
        case petsc_solver::fbcgs:        return "fbcgs";
        case petsc_solver::fbcgsr:       return "fbcgsr";
        case petsc_solver::lcd:          return "lcd";
        case petsc_solver::nash:         return "nash";
        case petsc_solver::stcg:         return "stcg";
        case petsc_solver::gltr:         return "gltr";
        case petsc_solver::qcg:          return "qcg";

        // this case should already be removed
        case petsc_solver::default_val:
            return "gmres";

        default:
        {
            matcl_assert(0,"invalid case");
            return "gmres";
        }
    };
};

petsc_solver converter_solver::to_code(const std::string& name) const
{
    code_map map = get_map();
    
    std::string name_conv   = name;
    std::transform(name_conv.begin(), name_conv.end(), name_conv.begin(), ::tolower);

    auto pos = map.find(name);
    if (pos == map.end())
        throw error::unrecognized_petsc_solver(name_conv);

    return pos->second;
};
const converter_solver::code_map& converter_solver::get_map() const
{
    static code_map map = init_map();
    return map;
}
converter_solver::code_map converter_solver::init_map()
{
    code_map map;

    map["preonly"]      = petsc_solver::preonly;
    map["cg"]           = petsc_solver::cg;
    map["fcg"]          = petsc_solver::fcg;
    map["cgne"]         = petsc_solver::cgne;
    map["gmres"]        = petsc_solver::gmres;
    map["richardson"]   = petsc_solver::richardson;
    map["chebyshev"]    = petsc_solver::chebyshev;
    map["tfqmr"]        = petsc_solver::tfqmr;
    map["tcqmr"]        = petsc_solver::tcqmr;
    map["bcgs"]         = petsc_solver::bcgs;
    map["bcgsl"]        = petsc_solver::bcgsl;
    map["bicg"]         = petsc_solver::bicg;
    map["cgs"]          = petsc_solver::cgs;
    map["fgmres"]       = petsc_solver::fgmres;
    map["lgmres"]       = petsc_solver::lgmres;
    map["gcr"]          = petsc_solver::gcr;
    map["cr"]           = petsc_solver::cr;
    map["minres"]       = petsc_solver::minres;
    map["symmlq"]       = petsc_solver::symmlq;
    map["lsqr"]         = petsc_solver::lsqr;
    map["dgmres"]       = petsc_solver::dgmres;
    map["ibcgs"]        = petsc_solver::ibcgs;
    map["fbcgs"]        = petsc_solver::fbcgs;
    map["fbcgsr"]       = petsc_solver::fbcgsr;
    map["lcd"]          = petsc_solver::lcd;
    map["nash"]         = petsc_solver::nash;
    map["stcg"]         = petsc_solver::stcg;
    map["gltr"]         = petsc_solver::gltr;
    map["qcg"]          = petsc_solver::qcg;

    return map;
};

std::string converter_precond::to_string(const Integer& val) const
{
    petsc_precond precond    = (petsc_precond)val;

    switch(precond)
    {            
        case petsc_precond::none:       return "none";
        case petsc_precond::ilu:        return "ilu";
        case petsc_precond::icc:        return "icc";
        case petsc_precond::lu:         return "lu";
        case petsc_precond::svd:        return "svd";
        case petsc_precond::cholesky:   return "cholesky";
        case petsc_precond::sor:        return "sor";
        case petsc_precond::eisenstat:  return "eisenstat";
        case petsc_precond::jacobi:     return "jacobi";
        case petsc_precond::pbjacobi:   return "pbjacobi";
        case petsc_precond::kaczmarz:   return "kaczmarz";

        // this case should already be removed
        case petsc_precond::default_val:
            return "ilu";

        default:
        {
            matcl_assert(0,"invalid case");
            return "ilu";
        }
    }
};

const converter_precond::code_map& converter_precond::get_map() const
{
    static code_map map = init_map();
    return map;
}
converter_precond::code_map converter_precond::init_map()
{
    code_map map;

    map["none"]     = petsc_precond::none;
    map["ilu"]      = petsc_precond::ilu;
    map["icc"]      = petsc_precond::icc;
    map["lu"]       = petsc_precond::lu;
    map["svd"]      = petsc_precond::svd;
    map["cholesky"] = petsc_precond::cholesky;
    map["sor"]      = petsc_precond::sor;
    map["eisenstat"]= petsc_precond::eisenstat;
    map["jacobi"]   = petsc_precond::jacobi;
    map["pbjacobi"] = petsc_precond::pbjacobi;
    map["kaczmarz"] = petsc_precond::kaczmarz;
    
    return map;
};

petsc_precond converter_precond::to_code(const std::string& name) const
{
    code_map map = get_map();
    
    std::string name_conv   = name;
    std::transform(name_conv.begin(), name_conv.end(), name_conv.begin(), ::tolower);

    auto pos = map.find(name);
    if (pos == map.end())
        throw error::unrecognized_petsc_preconditioner(name_conv);

    return pos->second;
};


std::string converter_refine::to_string(const Integer& val) const
{
    orth_refine_type ref    = (orth_refine_type)val;

    switch(ref)
    {            
        case orth_refine_type::never:       return "REFINE_NEVER";
        case orth_refine_type::if_needed:   return "REFINE_IFNEEDED";
        case orth_refine_type::always:      return "REFINE_ALWAYS";

        // this case should already be removed
        case orth_refine_type::default_val:
            return "REFINE_NEVER";

        default:
        {
            matcl_assert(0,"invalid case");
            return "REFINE_NEVER";
        }
    }
};

std::string converter_fcg_orth::to_string(const Integer& val) const
{
    fcg_orth_type ref    = (fcg_orth_type)val;

    switch(ref)
    {            
        case fcg_orth_type::full:       return "STANDARD";
        case fcg_orth_type::partial:    return "NOTAY";

        // this case should already be removed
        case fcg_orth_type::default_val:
            return "notay";

        default:
        {
            matcl_assert(0,"invalid case");
            return "notay";
        }
    }
};

std::string converter_norm_type::to_string(const Integer& val) const
{
    petsc_conv_norm ref    = (petsc_conv_norm)val;

    switch(ref)
    {            
        case petsc_conv_norm::none:             return "NONE";
        case petsc_conv_norm::preconditioned:   return "PRECONDITIONED";
        case petsc_conv_norm::true_resid:       return "UNPRECONDITIONED";
        case petsc_conv_norm::normal:           return "NATURAL";

        // this case should already be removed
        case petsc_conv_norm::default_val:
            return "DEFAULT";

        default:
        {
            matcl_assert(0,"invalid case");
            return "UNPRECONDITIONED";
        }
    }
};

std::string converter_order::to_string(const Integer& val) const
{
    petsc_ordering order    = (petsc_ordering)val;

    switch(order)
    {            
        case petsc_ordering::natural:   return "natural";
        case petsc_ordering::nd:        return "nd";
        case petsc_ordering::w1d:       return "1wd";
        case petsc_ordering::rcm:       return "rcm";
        case petsc_ordering::qmd:       return "qmd";
        case petsc_ordering::rowlength: return "rowlength";
            
        // this case should already be removed
        case petsc_ordering::default_val:
            return "natural";

        default:
        {
            matcl_assert(0,"invalid case");
            return "natural";
        }
    }
};

std::string converter_shift_type::to_string(const Integer& val) const
{
    petsc_shift_type st    = (petsc_shift_type)val;

    switch(st)
    {            
        case petsc_shift_type::none:                return "NONE";
        case petsc_shift_type::nonzero:             return "NONZERO";
        case petsc_shift_type::positive_definite:   return "POSITIVE_DEFINITE";
        case petsc_shift_type::inblocks:            return "INBLOCKS";
            
        // this case should already be removed
        case petsc_shift_type::default_val:
            return "NONZERO";

        default:
        {
            matcl_assert(0,"invalid case");
            return "natural";
        }
    }
};

std::string converter_jacobi_type::to_string(const Integer& val) const
{
    petsc_jacobi_type st    = (petsc_jacobi_type)val;

    switch(st)
    {            
        case petsc_jacobi_type::diag:       return "DIAGONAL";
        case petsc_jacobi_type::max_abs_row:return "ROWMAX";
        case petsc_jacobi_type::sumrows:    return "SUMROWS";
            
        // this case should already be removed
        case petsc_jacobi_type::default_val:
            return "DIAGONAL";

        default:
        {
            matcl_assert(0,"invalid case");
            return "natural";
        }
    }
};

std::string converter_precond_side::to_string(const Integer& val) const
{
    petsc_precond_side st    = (petsc_precond_side)val;

    switch(st)
    {            
        case petsc_precond_side::left:      return "LEFT";
        case petsc_precond_side::right:     return "RIGHT";
        case petsc_precond_side::symmetric: return "SYMMETRIC";
            
        // this case should already be removed
        case petsc_precond_side::default_val:
            return "LEFT";

        default:
        {
            matcl_assert(0,"invalid case");
            return "natural";
        }
    }
};

std::string converter_optional::to_string(const Integer& val) const
{
    if (val == m_val)
        return "ON";
    else
        return "";
};

template<class VT>
using string_converter_ptr = std::shared_ptr<string_converter<VT>>;

class option_access
{
    public:
        virtual ~option_access(){};
        virtual bool    get_value(const options& opts, std::string& val) const = 0;
        virtual bool    has_value(const options& opts) const = 0;
};

template<class Option>
class option_access_impl : public option_access
{
    private:
        using VT            = typename Option::value_type;
        using str_conv      = string_converter<VT>;
        using str_conv_ptr  = std::shared_ptr<str_conv>;

    private:
        Option          m_opt;
        str_conv_ptr    m_conv;

    public:
        option_access_impl(const Option& opt, const str_conv_ptr& conv)
            :m_opt(opt), m_conv(conv)
        {};

        virtual bool get_value(const options& opts, std::string& val_ret) const
        {
            VT val = opts.get_option<VT>(m_opt);
            if (val == Option::m_default_value)
                return false;

            val_ret = m_conv->to_string(val);
            return true;
        };
        virtual bool has_value(const options& opts) const override
        {
            VT val = opts.get_option<VT>(m_opt);
            if (val == Option::m_default_value)
                return false;
            else 
                return true;
        };
};

template<class T, class VT = typename T::value_type>
std::shared_ptr<option_access> get_option_accesser(const T& opt, const string_converter_ptr<VT>& conv)
{
    using ptr_type  = std::shared_ptr<option_access>;
    return ptr_type(new option_access_impl<T>(opt, conv));
};

class registered_petsc_options
{
    private:
        using option_data   = std::shared_ptr<option_access>;
        using option_map    = std::map<std::string, option_data>;

    public:
        template<class T>
        static bool get_option(const options& opts, const std::string& prefix, const std::string& name, std::string& val)
        {
            if (allowed_prefix(prefix) == false)
                return false;

            const option_map& map   = get_map();
            
            auto pos = map.find(name);

            if (pos == map.end())
                return false;

            return pos->second->get_value(opts, val);
        };
        template<class T>
        static bool has_option(const options& opts, const std::string& prefix, const std::string& name)
        {
            if (allowed_prefix(prefix) == false)
                return false;

            const option_map& map   = get_map();
            
            auto pos = map.find(name);

            if (pos == map.end())
                return false;

            return pos->second->has_value(opts);
        };
    private:
        static bool allowed_prefix(const std::string& str)
        {
            if (str == "")
                return true;
            if (str == "sub_")
                return true;

            if (str.find("fieldsplit",0) != std::string::npos)
                return true;

            if (str.find("mg_lev",0) != std::string::npos)
                return true;

            if (str.find("mg_coarse",0) != std::string::npos)
                return true;

            return false;
        };

        static const option_map& get_map()
        {
            //in C++ static variables are thread safe
            static option_map map = init_map();
            return map;
        };

        static option_map init_map()
        {
            option_map map;
            map["-ksp_max_it"]          = get_option_accesser(opt::petsc::max_it(), conv_int());
            map["-ksp_type"]            = get_option_accesser(opt::petsc::solver(), conv_solver());
            map["-pc_type"]             = get_option_accesser(opt::petsc::preconditioner(), conv_precond());
            map["-ksp_atol"]            = get_option_accesser(opt::petsc::atol(), conv_real());
            map["-ksp_rtol"]            = get_option_accesser(opt::petsc::rtol(), conv_real());
            map["-ksp_divtol"]          = get_option_accesser(opt::petsc::divtol(), conv_real());
            map["-ksp_norm_type"]       = get_option_accesser(opt::petsc::norm_type(), conv_norm());
            map["-ksp_chebyshev_eigenvalues"]
                                        = get_option_accesser(opt::petsc::cheb_eig(), conv_vec());
            map["-ksp_chebyshev_esteig"]= get_option_accesser(opt::petsc::cheb_eig_est(), conv_vec());
            map["-ksp_gmres_restart"]   = get_option_accesser(opt::petsc::gmres_restart(), conv_int());
            map["-ksp_gcr_restart"]     = get_option_accesser(opt::petsc::gmres_restart(), conv_int());
            map["-ksp_lcd_restart"]     = get_option_accesser(opt::petsc::gmres_restart(), conv_int());
            map["-ksp_fcg_mmax"]        = get_option_accesser(opt::petsc::gmres_restart(), conv_int());
            map["-ksp_gmres_haptol"]    = get_option_accesser(opt::petsc::gmres_haptol(), conv_real());
            map["-ksp_lcd_haptol"]      = get_option_accesser(opt::petsc::gmres_haptol(), conv_real());
            map["-ksp_fcg_truncation_type"]
                                        = get_option_accesser(opt::petsc::fcg_orthog(), conv_fcg_orth());
            map["-ksp_gmres_cgs_refinement_type"]
                                        = get_option_accesser(opt::petsc::gmres_refine(), conv_refine());
            map["-ksp_gmres_classicalgramschmidt"]
                                        = get_option_accesser(opt::petsc::gmres_gs(), conv_optional(1));
            map["-ksp_gmres_modifiedgramschmidt"]
                                        = get_option_accesser(opt::petsc::gmres_gs(), conv_optional(0));
            map["-ksp_dgmres_max_eigen"]= get_option_accesser(opt::petsc::dgmres_max_eigen(), conv_int());
            map["-ksp_dgmres_eigen"]    = get_option_accesser(opt::petsc::dgmres_eigen(), conv_int());
            map["-ksp_lgmres_augment"]  = get_option_accesser(opt::petsc::lgmres_augment(), conv_int());
            map["-ksp_gltr_radius"]     = get_option_accesser(opt::petsc::trust_region_radius(), conv_real());
            map["-ksp_nash_radius"]     = get_option_accesser(opt::petsc::trust_region_radius(), conv_real());
            map["-ksp_stcg_radius"]     = get_option_accesser(opt::petsc::trust_region_radius(), conv_real());
            map["-ksp_richardson_scale"]= get_option_accesser(opt::petsc::richardson_damping(), conv_real());
            map["-ksp_richardson_self_scale"]
                                        = get_option_accesser(opt::petsc::richardson_auto(), conv_int_bool());            
            map["-ksp_bcgsl_ell"]       = get_option_accesser(opt::petsc::bcgsl_num_dir(), conv_int());        
            map["-ksp_lsqr_set_standard_error"]
                                        = get_option_accesser(opt::petsc::lsqr_std_est(), conv_optional(1));

            map["-pc_factor_levels"]    = get_option_accesser(opt::petsc::ilu_k(), conv_int());
            map["-pc_factor_diagonal_fill"]
                                        = get_option_accesser(opt::petsc::allow_diagonal_fill(), conv_int_bool());
            map["-pc_factor_mat_ordering_type"]
                                        = get_option_accesser(opt::petsc::ordering(), conv_order());
            map["-pc_factor_nonzeros_along_diagonal"]
                                        = get_option_accesser(opt::petsc::reororder_zero_diagonal(), conv_real());
            map["-pc_factor_zeropivot"]
                                        = get_option_accesser(opt::petsc::zero_pivot(), conv_real());
            map["-pc_factor_shift_amount"]
                                        = get_option_accesser(opt::petsc::shift_amount(), conv_real());
            map["-pc_factor_shift_type"]
                                        = get_option_accesser(opt::petsc::shift_type(), conv_shift_type());            
            map["-pc_sor_diagonal_shift"]
                                        = get_option_accesser(opt::petsc::shift_amount(), conv_real());
            map["-pc_sor_omega"]        = get_option_accesser(opt::petsc::sor_omega(), conv_real());
            map["-pc_eisenstat_omega"]  = get_option_accesser(opt::petsc::sor_omega(), conv_real());
            map["-pc_sor_backward"]     = get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::backward));
            map["-pc_sor_forward"]      = get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::forward));
            map["-pc_sor_symmetric"]    = get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::symmetric));
            map["-pc_sor_local_forward"]= get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::last));
            map["-pc_sor_local_symmetric"]
                                        = get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::last));
            map["-pc_sor_local_backward"]
                                        = get_option_accesser(opt::petsc::sor_type(), 
                                                conv_optional((Integer)petsc_sor_type::last));
            map["-pc_sor_its"]          = get_option_accesser(opt::petsc::sor_iter(), conv_int());
            map["-pc_sor_lits"]         = get_option_accesser(opt::petsc::sor_iter(), conv_int());
            map["-ksp_pc_side"]         = get_option_accesser(opt::petsc::precond_side(), conv_side());
            map["-pc_jacobi_type"]      = get_option_accesser(opt::petsc::jacobi_type(), conv_jacobi_type());
            map["-ksp_chebyshev_esteig_steps"]
                                        = get_option_accesser(opt::petsc::cheb_eig_steps(), conv_int());
            map["-pc_jacobi_abs"]       = get_option_accesser(opt::petsc::jacobi_abs(), conv_int_bool());
            map["-pc_kaczmarz_lambda"]  = get_option_accesser(opt::petsc::kaczmarz_omega(), conv_real());
            map["-pc_kaczmarz_symmetric"]
                                        = get_option_accesser(opt::petsc::kaczmarz_sym(), conv_int_bool());
            map["-pc_fieldsplit_detect_saddle_point"]
                                        = get_option_accesser(opt::petsc::detect_saddle(), conv_int_bool());

            return map;
        };

        static string_converter_ptr<Integer> conv_int()
        {
            return string_converter_ptr<Integer>(new default_converter<Integer>());
        };
        static string_converter_ptr<Integer> conv_int_bool()
        {
            return string_converter_ptr<Integer>(new converter_int_bool());
        };
        static string_converter_ptr<bool> conv_bool()
        {
            return string_converter_ptr<bool>(new default_converter<bool>());
        };

        static string_converter_ptr<Integer> conv_optional(Integer val)
        {
            return string_converter_ptr<Integer>(new converter_optional(val));
        };
        static string_converter_ptr<Real> conv_real()
        {
            return string_converter_ptr<Real>(new default_converter<Real>());
        };

        static string_converter_ptr<Integer> conv_solver()
        {
            return string_converter_ptr<Integer>(new converter_solver());
        }
        
        static string_converter_ptr<Integer> conv_refine()
        {
            return string_converter_ptr<Integer>(new converter_refine());
        }
        static string_converter_ptr<Integer> conv_fcg_orth()
        {
            return string_converter_ptr<Integer>(new converter_fcg_orth());
        }

        static string_converter_ptr<Integer> conv_norm()
        {
            return string_converter_ptr<Integer>(new converter_norm_type());
        }
        
        static string_converter_ptr<Integer> conv_precond()
        {
            return string_converter_ptr<Integer>(new converter_precond());
        }
        static string_converter_ptr<Integer> conv_order()
        {
            return string_converter_ptr<Integer>(new converter_order());
        }
        static string_converter_ptr<Integer> conv_shift_type()
        {
            return string_converter_ptr<Integer>(new converter_shift_type());
        }
        static string_converter_ptr<Integer> conv_side()
        {
            return string_converter_ptr<Integer>(new converter_precond_side());
        }        
        static string_converter_ptr<Integer> conv_jacobi_type()
        {
            return string_converter_ptr<Integer>(new converter_jacobi_type());
        }        

        static string_converter_ptr<std::vector<Real>> conv_vec()
        {
            return string_converter_ptr<std::vector<Real>>(new default_converter<std::vector<Real>>());
        }
        
};

template<class T, class PT>
struct convert_val
{
    static PT eval(const T& val)
    { 
        return val; 
    };
};
template<>
struct convert_val<bool, PetscBool>
{
    static PetscBool eval(const bool val)
    { 
        return val ? PETSC_TRUE : PETSC_FALSE; 
    };
};
template<>
struct convert_val<PetscBool, bool>
{
    static bool eval(PetscBool val)
    { 
        return val == PETSC_TRUE ? true : false; 
    };
};

template<class T, class PT>
struct convert_val<std::vector<T>, std::vector<PT>>
{
    static std::vector<PT> eval(const std::vector<T>& val)
    { 
        std::vector<PT> ret(val.size());

        for (size_t i = 0; i < val.size(); ++i)
            ret[i]  = convert_val<T,PT>::eval(val[i]);

        return ret; 
    };
};
template<class T>
struct convert_val<std::vector<T>, std::vector<T>>
{
    static std::vector<T> eval(const std::vector<T>& val)
    { 
        return val;
    };
};

template<class PT>
struct convert_string_val
{
    static PT eval(const std::string& val);
};

template<>
struct convert_string_val<Integer>
{
    static Integer eval(const std::string& val)
    {
        PetscInt ret;
        PetscOptionsStringToInt(val.c_str(), &ret);
        return ret;
    }
};
template<>
struct convert_string_val<Real>
{
    static Real eval(const std::string& val)
    {
        PetscReal ret;
        PetscOptionsStringToReal(val.c_str(), &ret);
        return ret;
    }
};
template<>
struct convert_string_val<bool>
{
    static bool eval(const std::string& val)
    {
        PetscBool ret;
        PetscOptionsStringToBool(val.c_str(), &ret);
        return convert_val<PetscBool, bool>::eval(ret);
    }
};
template<>
struct convert_string_val<PetscBool>
{
    static PetscBool eval(const std::string& val)
    {
        PetscBool ret;
        PetscOptionsStringToBool(val.c_str(), &ret);
        return ret;
    }
};
template<>
struct convert_string_val<std::string>
{
    static std::string eval(const std::string& val)
    {
        return val;
    }
};
template<class T>
struct convert_string_val<std::vector<T>>
{
    static std::vector<T> eval(const std::string& val)
    {
        std::vector<T> ret;
        if (val.size() == 0)
            return ret;

        std::istringstream is(val);
        std::string elem;

        while (std::getline(is, elem, ','))
        {
            T loc   = convert_string_val<T>::eval(elem);
            ret.push_back(loc);
        }

        return ret;
    }
};
//-------------------------------------------------------
//              petsc_options_handler
//-------------------------------------------------------

petsc_options_handler::petsc_options_handler(const options* opts, options* missing)
    :m_opts(opts), m_opts_missing(missing)
{
    make_util();
}

template<class T, class PT>
bool petsc_options_handler::get_option(const std::string& prefix, const std::string& name, 
                                       const std::string& full_name, PT& val) const
{
    std::string string_val;
    if (registered_petsc_options::get_option<T>(*m_opts, prefix, name, string_val) == true)
    {
        val = convert_string_val<PT>::eval(string_val);
        return true;
    };

    option opt = string_option<T>(full_name);

    if (m_opts->has_option(opt) == false)
        return false;

    optional<T> opt_val = m_opts->get<T>(opt);

    if (!opt_val)
        return false;

    val = convert_val<T,PT>::eval(opt_val.value());
    return true;
}
template<class T>
bool petsc_options_handler::has_option(const std::string& prefix, const std::string& name, 
                                       const std::string& full_name) const
{
    if (registered_petsc_options::has_option<T>(*m_opts, prefix, name) == true)
        return true;

    option opt = string_option<T>(full_name);

    if (m_opts->has_option(opt) == false)
        return false;
    else
        return true;
}

petsc_options_handler::sptr_type
petsc_options_handler::make(const options* opt, options* missing)
{
    return sptr_type(new petsc_options_handler(opt, missing));
}

const options* petsc_options_handler::get_current_options() const
{
    return m_opts;
};
void petsc_options_handler::reset_options(const options* opts)
{
    m_opts = opts;
};

void petsc_options_handler::make_util()
{
    m_util              = options_utility_ptr(new _OptionUtility);

    m_util->listener    = nullptr;
    m_util->handler     = nullptr;
    m_util->prefix      = nullptr;

    try
    {
        m_util->listener    = new _PetscOptionsListenerC;
        m_util->handler     = new _PetscOptionsHandlerC;
        m_util->handler->data = this;

        m_util->listener->data  = m_opts_missing;

        m_util->PetscOptionsBool        = PetscOptionsBool_impl;
        m_util->PetscOptionsBool2       = PetscOptionsBool2_impl;
        m_util->PetscOptionsBoolArray   = PetscOptionsBoolArray_impl;
        m_util->PetscOptionsEList       = PetscOptionsEList_impl;
        m_util->PetscOptionsInt         = PetscOptionsInt_impl;
        m_util->PetscOptionsIntArray    = PetscOptionsIntArray_impl;
        m_util->PetscOptionsList        = PetscOptionsList_impl;
        m_util->PetscOptionsName        = PetscOptionsName_impl;
        m_util->PetscOptionsReal        = PetscOptionsReal_impl;
        m_util->PetscOptionsRealArray   = PetscOptionsRealArray_impl;
        m_util->PetscOptionsScalar      = PetscOptionsScalar_impl;
        m_util->PetscOptionsString      = PetscOptionsString_impl;
        m_util->PetscOptionsStringArray = PetscOptionsStringArray_impl;
        m_util->PetscOptionsScalarArray = PetscOptionsScalarArray_impl;
        m_util->PetscOptionsEnumArray   = PetscOptionsEnumArray_impl;
    }
    catch(...)
    {
        if(m_util->handler != nullptr)
            delete m_util->handler;
        if(m_util->listener != nullptr)
            delete m_util->listener;
        if(m_util->prefix != nullptr) 
            delete m_util->prefix;

        throw;
    };
};

thread_local _OptionUtility* thread_util = nullptr;

option_stack petsc_options_handler::install()
{
    _OptionUtility* old = thread_util;
    thread_util = this->m_util.get();
    ::setOptionUtilityConstructor(&thread_utility_constr);

    return std::move(option_stack(old));
};

option_stack::~option_stack()
{
    if (m_valid)
        thread_util = m_old;
};

petsc_options_handler::~petsc_options_handler()
{
    ::setOptionUtilityConstructor(&null_utility_constr);

    delete m_util->listener;
    delete m_util->handler;
    delete m_util->prefix; 
};

petsc_options_handler::util_type petsc_options_handler::null_utility_constr()
{
    return nullptr;
};

petsc_options_handler::util_type petsc_options_handler::thread_utility_constr()
{
    return thread_util;
};

scoped_options::scoped_options(petsc_options_handler& handler, const options* opts)
    :m_handler(&handler), m_old_options(handler.get_current_options())
{
    m_handler->reset_options(opts);
    m_stack_opt.reset_fresh(handler.install());
}
scoped_options::~scoped_options()
{
    m_handler->reset_options(m_old_options);
};

}}

#pragma warning(pop)

#endif