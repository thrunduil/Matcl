#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              Process dependencies
//----------------------------------------------------------------------------------


template<class New_Elem, class List, bool Is_Member>
struct insert_new_dep
{
    using type = List;
};
template<class New_Elem, class ... Deps>
struct insert_new_dep<New_Elem, dps<Deps...>, false>
{
    using type = dps<New_Elem, Deps...>;
};

template<class Deps1, class Deps2>
struct get_unique
{};
template<class Dep1, class ...Deps1, class Deps2>
struct get_unique<dps<Dep1, Deps1...>, Deps2>
{    
    using type_1                = typename get_unique<dps<Deps1...>, Deps2>::type;
    static const bool is_member = is_member<Dep1, Deps2>::value;
    using type                  = typename insert_new_dep<Dep1, type_1,is_member>::type;    
};
template<class Dep1, class Deps2>
struct get_unique<dps<Dep1>, Deps2>
{
    static const bool is_member = is_member<Dep1, Deps2>::value;
    using type                  = typename std::conditional<is_member, dps<>, dps<Dep1>> :: type;
};

template<class Dep1, class Dep2>
struct merge_deps
{};
template<class ... Dep1, class ... Dep2>
struct merge_deps<dps<Dep1...>, dps<Dep2...>>
{
    using type = dps<Dep1..., Dep2...>;
};
template<class ... Dep1>
struct merge_deps<dps<Dep1...>, void>
{
    using type = dps<Dep1...>;
};
template<class Dep1, class Dep2>
struct link_deps
{
    static_assert(details::dependent_false<Dep1>::value, 
                "this type should not be instantiated");
};
template<class Dep1>
struct link_deps<Dep1,Dep1>
{
    using type = Dep1;
};
template<class Dep1, class... Deps1, class Dep2>
struct link_deps<dps<Dep1, Deps1 ... >, Dep2>
{
    using unique    = typename get_unique<dps<Dep1, Deps1 ... >, Dep2>::type;
    using type      = typename merge_deps<unique, Dep2>::type;
};

template<class Dep1, class... Deps1>
struct link_deps<dps<Dep1, Deps1 ... >, dps<>>
{
    using type  = dps<Dep1, Deps1 ... >;
};
template<class Dep1, class... Deps1>
struct link_deps<dps<Dep1, Deps1 ... >, dps<Dep1, Deps1 ... >>
{
    using type  = dps<Dep1, Deps1 ... >;
};
template<class Dep2>
struct link_deps<dps<>, Dep2>
{
    using type  = Dep2;
};
template<>
struct link_deps<dps<>, dps<>>
{
    using type  = dps<>;
};


template<class Deps, class New_Dep>
struct add_dep
{    
    using temp_type = dps<New_Dep>;
    using type      = typename link_deps<temp_type, Deps>::type;
};

template<class Dep, class Unique>
struct is_member_dps
{
    static_assert(details::dependent_false<Dep>::value,
                  "this type should not be instantiated");
};
template<class Dep>
struct is_member_dps<Dep,dps<>>
{
    static const bool value = false;
};
template<class Dep, class ...Args>
struct is_member_dps<Dep, dps<Dep,Args...>>
{
    static const bool value = true;
};
template<class Dep, class Arg, class ...Args>
struct is_member_dps<Dep, dps<Arg, Args...>>
{
    static const bool value = is_member_dps<Dep,dps<Args...>>::value;
};

}}