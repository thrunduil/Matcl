#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen { namespace details
{

//allow use always false conditions in static_assert 
//template<class T>
//struct dependent_false          { static const bool value = false; };

template<class... T>
struct dependent_false_var      { static const bool value = false; };

template<class T, T Val>
struct dependent_value_false    { static const bool value = false; };

// hide type
template<class T>
struct lazy_type
{
    using type = T;
};

}}}

