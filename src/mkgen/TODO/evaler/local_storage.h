#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              local_storage
//----------------------------------------------------------------------------------
template<class Val, Integer Size>
struct local_arrays
{
    const Val* restricted m_array[Size];
};
template<class Val>
struct local_arrays<Val,0>
{};

template<class Val, Integer Size>
struct local_values
{
    Val m_array[Size];
};
template<class Val>
struct local_values<Val,0>
{};

template<class Val, Integer Pos, class List_Tags>
struct init_storage_arrays
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};

template<class Val, Integer Pos, class Dep>
struct init_storage_arrays_dep
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};
template<class Val, Integer Pos, class Tag,Integer Size>
struct init_storage_arrays_dep<Val,Pos,dep<Tag,Size,dep_temp>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    inline_lev_1
    static void eval(Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        const Val* arr  = &ev->get<1,Data_Provider,Tag>(dp);
        ls->set_array<Pos>(arr);
    };
};
template<class Val, Integer Pos, class Tag>
struct init_storage_arrays_dep<Val,Pos,dep<Tag,0,dep_extern>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    inline_lev_1
    static void eval(Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        (void)ev;

        const Val* arr = Tag::get_data_ptr<Val,0>(dp);
        ls->set_array<Pos>(arr);
    };
};

template<class Val, Integer Pos, class Tag, class ... Tags>
struct init_storage_arrays<Val,Pos,list::list<Tag,Tags...>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    inline_initializer
    static void eval(Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        using dep       = typename make_return_dep_from_tag<Tag>::type;
        init_storage_arrays_dep<Val,Pos,dep>::eval(dp,ev,ls);

        init_storage_arrays<Val,Pos+1,list::list<Tags...>>::eval(dp,ev,ls);
    };
};

template<class Val, Integer Pos>
struct init_storage_arrays<Val,Pos,list::list<>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    static void eval(const Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        (void)dp;
        (void)ev;
        (void)ls;
    };
};

template<class Val, Integer Pos, class List_Deps>
struct init_storage_scalars
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};
template<class Val, Integer Pos, class Tag, class ... Deps>
struct init_storage_scalars<Val,Pos,list::list<Tag,Deps...>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    inline_initializer
    static void eval(Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        Val v = eval_scalar<Val>(Tag(), ls);
        ls->set_scalar<Tag>(v);

        init_storage_scalars<Val,Pos+1,list::list<Deps...>>::eval(dp,ev,ls);
    };
};
template<class Val, Integer Pos>
struct init_storage_scalars<Val,Pos,list::list<>>
{
    template<class Data_Provider, class Temp_Storage, class Local_Storage>
    static void eval(Data_Provider& dp, Temp_Storage* ev, Local_Storage* ls)
    {
        (void)dp;
        (void)ev;
        (void)ls;
    };
};

template<class Deps>
struct make_array_deps
{
    static_assert(md::dependent_false<Deps>::value,
                  "this type should not be instantiated");
};
template<class Dep, class ... Deps>
struct make_array_deps<list::list<Dep, Deps...>>
{
    using list_1    = typename make_array_deps<list::list<Deps...>>::type;
    using type      = typename list::push_front<list_1, Dep>::type;
};
template<class Tag, class ... Deps>
struct make_array_deps<list::list<dep<Tag,0,dep_scalar>, Deps...>>
{
    using type      = typename make_array_deps<list::list<Deps...>>::type;
};
template<class Tag, class ... Deps>
struct make_array_deps<list::list<dep<Tag,0,dep_computation>, Deps...>>
{
    using type      = typename make_array_deps<list::list<Deps...>>::type;
};
template<class ... Deps>
struct make_array_deps<list::list<void, Deps...>>
{
    using type      = typename make_array_deps<list::list<Deps...>>::type;
};
template<>
struct make_array_deps<list::list<>>
{
    using type = list::list<>;
};

template<class Deps>
struct make_scalar_deps
{
    static_assert(md::dependent_false<Deps>::value,
                  "this type should not be instantiated");
};
template<class Dep, class ... Deps>
struct make_scalar_deps<list::list<Dep, Deps...>>
{
    using type      = typename make_scalar_deps<list::list<Deps...>>::type;
};
template<class Tag, class ... Deps>
struct make_scalar_deps<list::list<dep<Tag,0,dep_scalar>, Deps...>>
{
    using list_1    = typename make_scalar_deps<list::list<Deps...>>::type;
    using type      = typename list::push_front<list_1, Tag>::type;
};
template<>
struct make_scalar_deps<list::list<>>
{
    using type = list::list<>;
};

template<class Tag, Integer Offset, class Colon, Integer Rows>
struct array_info
{
    using tag                   = Tag;
    using colon                 = Colon;
    static const Integer offset = Offset;
    static const Integer rows   = Rows;
};

template<class Subs, class Dep>
struct make_array_info_temp
{
    static_assert(md::dependent_false<Subs>::value,
                  "this type should not be instantiated");
};
template<class Tag, class Dep>
struct make_array_info_temp<modif2<Tag,colon_all, -1, -1>, Dep>
{
    using offset_info           = decltype(get_temp_offset_info(Tag()));

    static const Integer off    = offset_info::offset;
    using base_tag              = typename offset_info::base_tag;
    using type                  = array_info<base_tag,off,colon_all, Dep::size>;
};
template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols, class Dep>
struct make_array_info_temp<modif2<Tag,Colon, Mat_Rows, Mat_Cols>, Dep>
{
    static_assert(std::is_base_of<temp_tag_base,Tag>::value == false, "not implemented");

    static const Integer off    = 0;
    using type                  = array_info<Tag,off,Colon, Mat_Rows>;
};

template<class Subs_Context, class Dep>
struct make_array_info_1
{
    static_assert(md::dependent_false<Dep>::value,
                  "this type should not be instantiated");
};
template<class Subs_Context, class Tag>
struct make_array_info_1<Subs_Context, dep<Tag,0,dep_extern>>
{
    using type = array_info<Tag,0,colon_all,0>;
};
template<class Subs_Context, class Tag, Integer Size>
struct make_array_info_1<Subs_Context, dep<Tag,Size,dep_temp>>
{
    using subs  = decltype(get_substitution(Subs_Context(),Tag()));
    using dep   = dep<Tag,Size,dep_temp>;
    using type  = typename make_array_info_temp<subs, dep>::type;
};

template<class Subs_Context, class List_Deps>
struct make_array_info
{
    static_assert(md::dependent_false<List_Deps>::value,
                  "this type should not be instantiated");
};
template<class Subs_Context, class... Deps>
struct make_array_info<Subs_Context, list::list<Deps...>>
{
    using type = list::list<typename make_array_info_1<Subs_Context, Deps>::type...>;
};

template<class List, class Elem, bool Static>
struct insert_elem{};

template<class ... List, class Elem>
struct insert_elem<list::list<List...>,Elem,true>
{
    using type = list::list<List...>;
};
template<class ... List, class Elem>
struct insert_elem<list::list<List...>,Elem,false>
{
    using type = list::list<typename Elem::tag, List...>;
};

template<class List_Deps>
struct make_array_tags
{
    static_assert(md::dependent_false<List_Deps>::value,
                  "this type should not be instantiated");
};
template<class... Info>
struct make_array_tags<list::list<Info...>>
{
    using list_tags = list::list<typename Info::tag ... >;
    using type      = typename list::unique_list<list_tags>::type;
};

template<class Subs_Type, class Dep, class Subs_Dep> 
struct merge_return_subs
{
    using type = Subs_Type;
};
template<class Dep1, class Dep2, class Tag1, class Tag2>
struct merge_return_subs_impl
{
    using type = Tag2;
};
template<class Dep1, class Dep2, class Tag1>
struct merge_return_subs_impl<Dep1,Dep2,Tag1,Tag1>
{
    using type = void;
};
template<class Dep, class Subs_Dep> 
struct merge_return_subs<void,Dep,Subs_Dep>
{
    using type = typename merge_return_subs_impl<Dep,Subs_Dep,typename Dep::tag, typename Subs_Dep::tag>::type;
};


template<class Subs_Context, class ... Deps>
struct get_return_subs
{
    using type = void;
};
template<class Subs_Context, class Dep, class ... Deps>
struct get_return_subs<Subs_Context,Dep,Deps...>
{
    using next      = typename get_return_subs<Subs_Context,Deps...>::type;
    using tag       = typename Dep::tag;
    using subs      = decltype(get_substitution(Subs_Context(),tag()));
    using type      = typename merge_return_subs<next,Dep,subs>::type;
};

template<class Tag, bool Is_Temp>
struct make_return_dep_from_tag_impl
{
    using type = dep<Tag,0,dep_extern>;
};
template<class Tag>
struct make_return_dep_from_tag_impl<Tag,true>
{
    using type = dep<Tag,0,dep_temp>;
};

template<class Tag>
struct is_temporary_tag
{
    static const bool value = std::is_base_of<temp_tag_base,Tag>::value;
};

template<class Tag>
struct make_return_dep_from_tag : make_return_dep_from_tag_impl<Tag, is_temporary_tag<Tag>::value>
{};

template<>
struct make_return_dep_from_tag<void>
{
    using type = void;
};

template<class Subs_Context, class Deps>
struct make_deps_list
{
    static_assert(md::dependent_false<Subs_Context>::value,
                  "this type should not be instantiated");
};

template<class Subs_Context, class ... Deps>
struct make_deps_list<Subs_Context, dps<Deps...>>
{
    using ret_tag   = typename get_return_subs<Subs_Context, Deps...>::type;
    using ret_dep   = typename make_return_dep_from_tag<ret_tag>::type;
    using type      = list::list<ret_dep, Deps...>;
};


template<class Val, class Data_Provider, class Subs_Context>
struct local_storage
{
    using subs_context  = Subs_Context;
    using data_provider = Data_Provider;
    using value_type    = Val;

    using deps_all      = typename subs_context::deps_all;
    using deps_list     = typename make_deps_list<subs_context,deps_all>::type;
    using array_deps    = typename make_array_deps<deps_list>::type;
    using scalar_deps   = typename make_scalar_deps<deps_list>::type;
    using array_info    = typename make_array_info<subs_context, array_deps>::type;
    using array_tags    = typename make_array_tags<array_info>::type;

    using storage       = local_arrays<Val,list::size<array_tags>::value>;        
    using storage_scal  = local_values<Val,list::size<scalar_deps>::value>;
    using dep_return    = typename list::elem_at_pos<array_deps,0>::type;
    using return_tag    = typename dep_return::tag;

    storage         m_data;
    storage_scal    m_data_scal;

    template<class Data_Provider, class Temp_Storage>
    inline_initializer
    void init(Data_Provider& dp, Temp_Storage* ev)
    {
        init_storage_arrays<Val, 0, array_tags>::eval(dp,ev,this);

        //TODO: scalars are initialized too early, must be initialized during temp storage 
        //initialization
        init_storage_scalars<Val, 0, scalar_deps>::eval(dp,ev,this);
    };

    template<Integer Pos>
    inline_lev_1
    void set_array(const Val* arr)
    {
        m_data.m_array[Pos] = arr;
    };

    template<class Dep, Integer Pos>
    inline_lev_1
    const Val& get_temp() const
    {
        static const Integer pos    = list::elem_pos<array_deps,Dep>::value;
        using info                  = typename list::elem_at_pos<array_info,pos>::type;

        using tag                   = typename info::tag;
        using colon                 = typename info::colon;
        static const Integer offb   = info::offset;    
        static const Integer rows   = info::rows;

        static const Integer pos_colon  = colon_func::index<Pos,colon>::value;
        static const Integer col        = (pos_colon-1) / rows + 1;
        static const Integer row        = (pos_colon-1) % rows + 1;
        static const Integer off        = tag::get_offset(row,col);        

        static const Integer pos2   = list::elem_pos<array_tags,tag>::value;

        return m_data.m_array[pos2][off + offb];
    };

    template<class Tag, Integer Row, Integer Col>
    inline_lev_1
    const Val& get_extern() const
    {
        using dep_type              = extern_dep<Tag>;
        static const Integer off    = Tag::get_offset(Row,Col);        
        static const Integer pos    = list::elem_pos<array_deps,dep_type>::value;
        using info                  = typename list::elem_at_pos<array_info,pos>::type;

        using tag                   = typename info::tag;
        static const Integer offb   = info::offset;    
        static const Integer pos2   = list::elem_pos<array_tags,tag>::value;

        return m_data.m_array[pos2][off + offb];
    }

    template<class Tag, Integer Row, Integer Col>
    inline_lev_1
    const Val& get_return() const
    {
        using dep_type              = typename list::elem_at_pos<array_deps,0>::type;
        static_assert(std::is_same<Tag,typename dep_type::tag>::value, "invalid tag");

        return get_extern<Tag,Row,Col>();
    };

    template<class Dep>
    inline_lev_1
    const Val* get_array() const
    {
        static const Integer pos    = list::elem_pos<array_deps,Dep>::value;
        using info                  = typename list::elem_at_pos<array_info,pos>::type;

        using tag                   = typename info::tag;
        static const Integer pos2   = list::elem_pos<array_tags,tag>::value;

        return m_data.m_array[pos2];
    };

    template<class Tag>
    inline_lev_1
    void set_scalar(const Val& v)
    {
        static const Integer pos    = list::elem_pos<scalar_deps,Tag>::value;
        m_data_scal.m_array[pos]    = v;
    };

    template<class Tag>
    inline_lev_1
    Val get_scalar() const
    {
        static const Integer pos    = list::elem_pos<scalar_deps,Tag>::value;
        return m_data_scal.m_array[pos];
    };

    local_storage() = default;
    local_storage(const local_storage&) = delete;
    local_storage& operator=(const local_storage&) = delete;
};

//----------------------------------------------------------------------------------
//                              make_local_storage
//----------------------------------------------------------------------------------


}}