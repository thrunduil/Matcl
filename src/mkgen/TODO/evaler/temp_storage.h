#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"
#include "mkgen/TODO/evaler/storage_initializer.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              DSP_Modifier
//----------------------------------------------------------------------------------
template<class DSP_Modifier, class Dep>
struct is_modified
{
    static_assert(details::dependent_false<DSP_Modifier>::value, 
                "this type should not be instantiated");
};

template<class DSP_Modifier_List, class Dep>
struct get_modif_tag
{
    static_assert(details::dependent_false<DSP_Modifier_List>::value,
                "this type should not be instantiated");
};
template<class ... Modif, class Temp_Tag2, Integer Size, class Type>
struct is_modified<list<Modif...>, dep<Temp_Tag2,Size,Type>>
{
    static_assert(details::dependent_false<Temp_Tag2>::value, 
                "this type should not be instantiated");
};
template<class Temp_Tag2, Integer Size, class Type>
struct is_modified<list<>, dep<Temp_Tag2,Size,Type>>
{
    static const bool value = false;
};
template<class Temp_Tag, class Colon, class ... Args, class Ret_Tag, Integer R, Integer C, bool Init,
         Integer Size, class Type>
struct is_modified<list<dps_modif<Temp_Tag,Ret_Tag,Colon, R, C, Init>, Args...>, dep<Temp_Tag,Size,Type>>
{
    static const bool value = true;
};
template<class Temp_Tag, class Colon, class ... Args, class Ret_Tag, Integer R, Integer C, bool Init,
        class Temp_Tag2, Integer Size, class Type>
struct is_modified<list<dps_modif<Temp_Tag, Ret_Tag, Colon, R, C, Init>, Args...>, dep<Temp_Tag2, Size, Type>>
{
    static const bool value = is_modified<list<Args...>, dep<Temp_Tag2, Size, Type>>::value;
};

template<class Modif, class Dep>
struct is_modif_one
{
    static_assert(details::dependent_false<Modif>::value,
                  "this type should not be instantiated");
};

template<class Temp_Tag, class Ret_Tag, class Colon, Integer R, Integer C, bool Init, 
         class Temp_Tag2, Integer Size, class Type>
struct is_modif_one<dps_modif<Temp_Tag, Ret_Tag, Colon, R, C, Init>, dep<Temp_Tag2, Size, Type>>
{
    static const bool value = std::is_same<Temp_Tag, Temp_Tag2>::value;
};

template<class DPS_Modif_1, class DPS_Modif_2, class Dep, bool Is_Equal>
struct get_modif_tag_impl
{
    using type = DPS_Modif_1;
};
template<class DPS_Modif_1, class DPS_Modif_2, class ... Args, class Dep>
struct get_modif_tag_impl<DPS_Modif_1, list<DPS_Modif_2, Args...>, Dep, false>
{
    using type = typename get_modif_tag_impl<DPS_Modif_2, list<Args...>, Dep, 
                            is_modif_one<DPS_Modif_2, Dep>::value>::type;
};

template<class Arg, class ... Args, class Dep>
struct get_modif_tag<list<Arg, Args...>, Dep>
    :get_modif_tag_impl<Arg, list<Args...>, Dep, is_modif_one<Arg,Dep>::value>
{};

template<class Val, class Ret_Dep, Integer Pos>
struct get_ret_array
{};
template<class Val, class Tag, Integer Pos>
struct get_ret_array<Val, dep<Tag,0,dep_extern>,Pos>
{
    template<class Data_Provider, class Temp_Storage>
    inline_lev_1
    static Val* eval(Data_Provider& dp, Temp_Storage* ts)
    {   
        return Tag::get_data_ptr<Val, Pos>(dp);
    };
};
template<class Val, class Tag, Integer Pos>
struct get_ret_array<Val, dep<Tag,0,dep_temp>,Pos>
{
    template<class Data_Provider, class Temp_Storage>
    inline_lev_1
    static Val* eval(Data_Provider& dp, Temp_Storage* ts)
    {   
        const Val* ptr = &ts->get<Pos,Data_Provider,Tag>(dp);
        return const_cast<Val*>(ptr);
    };
};

//----------------------------------------------------------------------------------
//                              temp_storage_elem
//----------------------------------------------------------------------------------
//all elements are aligned according to VEC_ALIGN
template<class Subs_Context, class Val, class ... Deps>
struct temp_storage_elem;

template<class Val, Integer Size>
struct aligned_size
{
    static const Integer val_size   = sizeof(Val);
    static const Integer value0     = ((Size - 1) * val_size) / VEC_ALIGN;
    static const Integer value1     = (value0 + 1) * VEC_ALIGN;
    static const Integer value      = value1 / val_size;
};

template<Integer Offset, class Base_Tag>
struct offset_info
{
    static const Integer offset     = Offset;
    using base_tag                  = Base_Tag;
};

template<class Tag, class Base>
struct make_base_tag;

template<class Subs_Context, class Val, class Dep, bool Is_Modif, class ... Deps>
struct temp_storage_elem_impl : temp_storage_elem<Subs_Context, Val, Deps...>
{
    using tag       = typename Dep::tag;
    using base      = temp_storage_elem<Subs_Context, Val, Deps...>;
    using base_tag  = typename make_base_tag<tag,base>::type;

    static const Integer size0          = Dep::size;
    static const Integer size           = aligned_size<Val,size0>::value;
    static const Integer total_size     = size + base::total_size;
    static const Integer offset         = base :: next_offset;
    static const Integer next_offset    = offset + size;

    using tag_offset_info               = offset_info<offset, base_tag>;

    friend tag_offset_info get_temp_offset_info(tag)
    {
        return tag_offset_info();
    };

    template<Integer Pos, class Val, class Data_Provider, class Temp_Storage>
    inline_lev_1
    friend const Val& get_impl(Subs_Context, tag,Data_Provider& dp, const Temp_Storage* ts)
    {
        (void)dp;
        return ts->storage.elem<Pos - 1 + offset>();
    };
    template<Integer Pos, class Data_Provider, class Temp_Storage>
    inline_lev_1
    friend void set_impl(Subs_Context, tag, Data_Provider& dp, Temp_Storage* ts, Val&& val)
    {
        ts->storage[Pos - 1 + offset] = std::forward<Val>(val);
    };

    template<class Stored_Matrix, class Local_Storage>
    inline_initializer
    friend void init_temporary_impl(Subs_Context, Dep, Local_Storage& ls)
    {
        storage_initializer<Stored_Matrix, Val, colon_all, Dep, offset,1>
            ::eval<align_full>(ls);
    };
    template<class Stored_Matrix, class Visitor>
    friend void init_temporary_impl_accept(Subs_Context, Dep, Visitor& vis)
    {
        storage_initializer<Stored_Matrix, Val, colon_all, Dep, offset,1>::accept(vis);
    };
};

template<class Subs_Context, class Val, class Dep, class ...Deps>
struct temp_storage_elem_impl<Subs_Context, Val, Dep, true, Deps...>
    : temp_storage_elem<Subs_Context, Val, Deps...>
{
    public:
        using tag                   = typename Dep::tag;

        using dps_modifier          = typename Subs_Context::dps_modifier;
        using modif_type            = typename get_modif_tag<dps_modifier,Dep>::type;
        using ret_tag               = typename modif_type::modified_tag;
        using ret_dep               = typename make_return_dep_from_tag<ret_tag>::type;
        using colon_type            = typename modif_type::colon_type;
        using base                  = temp_storage_elem<Subs_Context, Val, Deps...>;

        static const Integer rows       = modif_type::rows;
        static const Integer cols       = modif_type::cols;
        static const Integer size0      = rows * cols;
        static const Integer size       = aligned_size<Val,size0>::value;
        static const Integer total_size = 0 + base::total_size;
        static const Integer offset     = base::next_offset;
        static const Integer next_offset= offset + 0;
        static const bool init          = modif_type::initialize;

        template<Integer Pos0, class Val, class Data_Provider, class Temp_Storage>
        inline_lev_1
        friend const Val& get_impl(Subs_Context, tag, Data_Provider& dp, const Temp_Storage* ts)
        {
            static const Integer pos    = get_pos_colon<Pos0, colon_type>::value;
            static const Integer pos0   = typename ret_tag::template get_offset<1,1>::value;
            static_assert(pos <= size, "invalid access");
            
            Val* storage                = get_ret_array<Val,ret_dep,pos0+pos-1>::eval(dp,ts);
            return storage[0];
        };

        template<Integer Pos0, class Data_Provider, class Temp_Storage>
        inline_lev_1
        friend void set_impl(Subs_Context, tag, Data_Provider& dp, Temp_Storage* ts, Val&& v)
        {
            static const Integer pos    = get_pos_colon<Pos0, colon_type>::value;
            static const Integer pos0   = typename ret_tag::template get_offset<1,1>::value;
            static_assert(pos <= size, "invalid access");

            using ret_dep               = typename make_return_dep_from_tag<ret_tag>::type;
            Val* storage                = get_ret_array<Val,ret_dep,pos0+pos-1>::eval(dp,ts);

            storage[0]                  = std::forward<Val>(v);
        };

        template<class Stored_Matrix, class Local_Storage>
        inline_initializer
        static void init_temporary(std::true_type, Local_Storage& ls)
        {
            static const Integer pos        = ret_tag::template get_offset<1, 1>::value;
            using align_r                   = ret_tag::root_align_type;
            static const Integer step       = ret_tag::step;
            using data_provider             = typename Local_Storage::data_provider;

            storage_initializer<Stored_Matrix,Val,colon_type,ret_dep,pos,step>::eval<align_r>(ls);
        };

        template<class Stored_Matrix, class Data_Provider>
        static void init_temporary(std::false_type, Data_Provider& dp)
        {};

        template<class Stored_Matrix, class Local_Storage>
        inline_initializer
        friend void init_temporary_impl(Subs_Context, Dep, Local_Storage& ls)
        {  
            using ver           = std::integral_constant<bool, init>;

            temp_storage_elem_impl::init_temporary<Stored_Matrix,Local_Storage>(ver(), ls);
        };

        template<class Stored_Matrix, class Visitor>
        static void init_temporary_accept(std::true_type, Visitor& vis)
        {
            storage_initializer<Stored_Matrix,Val,colon_type,ret_dep, 0,1>::accept(vis);
        };
        template<class Stored_Matrix, class Visitor>
        static void init_temporary_accept(std::false_type, Visitor& vis)
        {};

        template<class Stored_Matrix, class Visitor>
        friend void init_temporary_impl_accept(Subs_Context, Dep, Visitor& vis)
        {   
            using ver           = std::integral_constant<bool, init>;
            temp_storage_elem_impl::init_temporary_accept<Stored_Matrix,Visitor>(ver(),vis);
        };
};

template<class Subs_Context, class Val, class ... Deps>
struct temp_storage_elem
{
    static_assert(details::dependent_false<Val>::value,
                  "this type should not be instantiated");
};
template<class Subs_Context, class Val, class Dep, class ... Deps>
struct temp_storage_elem<Subs_Context, Val, Dep,Deps...> : 
        temp_storage_elem_impl<Subs_Context, Val, Dep, 
            is_modified<typename Subs_Context::dps_modifier, Dep>::value, Deps...>
{};
template<class Subs_Context, class Val>
struct temp_storage_elem<Subs_Context, Val>
{
    static const Integer total_size     = 0;
    static const Integer offset         = 0;
    static const Integer next_offset    = 0;
};

template<class Subs_Context, class Val, class Tag, class ... Deps>
struct temp_storage_elem<Subs_Context, Val, dep<Tag,0,dep_computation>, Deps...>
    :temp_storage_elem<Subs_Context, Val, Deps...>
{
    using base = temp_storage_elem<Subs_Context, Val, Deps...>;

    static const Integer total_size     = 0 + base::total_size;
    static const Integer offset         = base:: next_offset;
    static const Integer next_offset    = offset;
};

template<class Tag, class Base>
struct make_base_tag
{
    using type  = typename Base :: base_tag;
};

template<class Tag, class Subs_Context, class Val>
struct make_base_tag<Tag, temp_storage_elem<Subs_Context,Val>>
{
    using type  = Tag;
};

//----------------------------------------------------------------------------------
//                              temp_storage_initializer
//----------------------------------------------------------------------------------
template<class Val, Integer Length, class Deps>
struct temp_storage_initializer
{
    using first_half        = temp_storage_initializer<Val, Length / 2, Deps>;
    using rem_args_1        = typename first_half::remaining_args;
    
    using second_half       = temp_storage_initializer<Val, Length - Length / 2, rem_args_1>;
    using rem_args_2        = typename second_half::remaining_args;
    
    using remaining_args    = rem_args_2;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>    
    inline_no
    static void init_tags(Local_Storage& ls, Data_Provider& dp, Temp_Storage* owner)
    {
        //initialization in reverse order
        second_half::init_tags(ls, dp, owner);
        first_half::init_tags(ls, dp, owner);
    };

    template<class Visitor, class Temp_Storage>
    static void accept(Visitor& vis, Temp_Storage* owner)
    {
        second_half::accept<Visitor>(vis, owner);
        first_half::accept<Visitor>(vis, owner);
    };
};

template<class Val, class Deps>
struct temp_storage_initializer<Val,0, Deps>
{
    using remaining_args    = Deps;
    
    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    static void init_tags(Local_Storage& ls, Data_Provider& dp, Temp_Storage* owner)
    {
        (void)ls;
        (void)dp;
        (void)owner;
    };

    template<class Visitor, class Temp_Storage>
    static void accept(Visitor& vis, Temp_Storage* owner)
    {
        (void)vis;
        (void)owner;
    };
};

template<class Val, class Dep1, class ...Deps>
struct temp_storage_initializer<Val, 1, dps<Dep1,Deps...>>
{    
    using remaining_args    = dps<Deps...>;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    static void init_tags(Local_Storage& ls, Data_Provider& dp, Temp_Storage *ts)
    {
        using tag1          = typename Dep1::tag;
        tag_initializer(tag1(), ls, dp, ts);
    };

    template<class Visitor, class Temp_Storage>
    static void accept(Visitor& vis, Temp_Storage* ts)
    {
        using tag1          = typename Dep1::tag;
        tag_initializer_accept(tag1(), vis, ts);
    };
};
template<class Val, class Dep1, class Dep2, class ...Deps>
struct temp_storage_initializer<Val, 2, dps<Dep1, Dep2, Deps...>>
{    
    using remaining_args    = dps<Deps...>;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    static void init_tags(Local_Storage& ls, Data_Provider& dp, Temp_Storage* ts)
    {
        using tag1  = typename Dep1::tag;
        using tag2  = typename Dep2::tag;

        tag_initializer(tag2(), ls, dp, ts);
        tag_initializer(tag1(), ls, dp, ts);
    };

    template<class Visitor, class Temp_Storage>
    static void accept(Visitor& vis, Temp_Storage* ts)
    {
        using tag1  = typename Dep1::tag;
        using tag2  = typename Dep2::tag;

        tag_initializer_accept(tag2(), vis, ts);
        tag_initializer_accept(tag1(), vis, ts);
    };
};

template<class Val, class DPS_Modifier, class Deps>
struct temp_storage_members
{
    static_assert(details::dependent_false<Val>::value,
                  "this type should not be instantiated");
};
template<class Subs_Context, class Val, class Dep, class ... Deps>
struct temp_storage_members<Subs_Context, Val, dps<Dep,Deps...>> 
    : temp_storage_elem<Subs_Context, Val, Dep, Deps...>
{
    using base = temp_storage_elem<Subs_Context, Val, Dep, Deps...>;
    static const Integer total_size = base::total_size;
};
template<class Val, class DPS_Modifier>
struct temp_storage_members<Val, DPS_Modifier, dps<>>
{
    static const Integer total_size = 0;
};

template<class D>
struct get_dps_length
{
    static_assert(details::dependent_false<D>::value,
                  "this type should not be instantiated");
};
template<class ... Args>
struct get_dps_length<dps<Args...>>
{
    static const Integer value = sizeof...(Args);
};

//selected if set_impl for given dep is not defined in given temporary storage, try in parent storage
template<Integer Pos, class Data_Provider, class Temp_Storage, class Val, class Parent_Subs_Context, class Tag>
inline_lev_1
void set_impl(Parent_Subs_Context, Tag, Data_Provider& dp, Temp_Storage* ts, Val&& val)
{   
    ts->set_in_parent<Pos,Data_Provider,Tag>(dp,std::forward<Val>(val));
};
//selected if get_impl for given Dep is not defined in given temporary storage, try in parent storage
template<Integer Pos, class Val, class Data_Provider, class Temp_Storage, class Parent_Subs_Context, class Tag>
inline_lev_1
const Val& get_impl(Parent_Subs_Context, Tag, Data_Provider& dp, const Temp_Storage* ts)
{   
    return ts->get_in_parent<Pos,Data_Provider,Tag>(dp);
};

template<class Val>
struct empty_storage
{
    template<Integer Pos, class Data_Provider, class Dep, class Val>
    void set(Data_Provider& dp, Val&& val)
    {               
        static_assert(details::dependent_false<Dep>::value, 
                "temporary given by dep not found");
    };
    template<Integer Pos, class Data_Provider, class Dep>
    const Val& get(Data_Provider& dp) const
    {
        static_assert(details::dependent_false<Dep>::value, 
                "temporary given by dep not found");
    };
};

//----------------------------------------------------------------------------------
//                              temporary_storage
//----------------------------------------------------------------------------------
template<class Parent_Storage>
struct temp_storage_parent
{
    Parent_Storage*     m_parent;

    void init_parent(Parent_Storage* p)
    {
        m_parent = p;
    };
};
template<class Val>
struct temp_storage_parent<empty_storage<Val>>
{
    void init_parent(empty_storage<Val>*){};
};

template<class Val, class Subs_Context0, class Parent_Storage>
class temporary_storage : temp_storage_parent<Parent_Storage>
{
    public: 
        using val_type              = Val;
        using subs_context          = Subs_Context0;
        using base_type             = temp_storage_parent<Parent_Storage>;
        using expr_deps             = typename subs_context::expr_deps;
        using dps_modifier_maker    = typename subs_context::dps_mod_maker;        
        using deps_all              = typename subs_context::deps_temp;

        static const Integer deps_all_length
                                    = get_dps_length<deps_all>::value;

        using dps_modifier          = typename dps_modifier_maker::type;        
        using storage_member        = temp_storage_members<subs_context, Val, deps_all>;
        static const Integer size   = storage_member::total_size;

        stack_array<Val, size>      storage;        

    public:
        temporary_storage(Parent_Storage* parent)   {base_type::init_parent(parent);};

        template<class Data_Provider, class Local_Storage>
        inline_initializer
        void init_storage(Local_Storage& ls, Data_Provider& dp)
        {
            temp_storage_initializer<Val, deps_all_length, deps_all>::init_tags(ls, dp, this);
        };
        template<class Visitor>
        void accept(Visitor& vis)
        {
            temp_storage_initializer<Val, deps_all_length, deps_all>::accept<Visitor>(vis, this);
        };

        template<Integer Pos, class Data_Provider, class Tag>
        inline_lev_1
        const Val& get(Data_Provider& dp) const
        {
            return get_impl<Pos, Val, Data_Provider>(subs_context(), Tag(),dp,this);
        };
        template<Integer Pos, class Data_Provider, class Tag>
        inline_lev_1
        const Val& get_in_parent(Data_Provider& dp) const
        {   
            return m_parent->get<Pos, Data_Provider, Tag>(dp);
        };

        template<Integer Pos, class Data_Provider, class Tag>
        inline_lev_1
        void set(Data_Provider& dp, Val&& v)
        {     
            return set_impl<Pos,Data_Provider>(subs_context(), Tag(), dp, this, 
                                                std::forward<Val>(v));
        };
        template<Integer Pos, class Data_Provider, class Tag>
        inline_lev_1
        void set_in_parent(Data_Provider& dp, Val&& v)
        {   
            return m_parent->set<Pos, Data_Provider, Tag>(dp, std::forward<Val>(v));
        };

        template<class Tag, class Computation, class Local_Storage>
        inline_initializer
        void init_computation(Local_Storage& ls)
        {
            using assignments   = typename Computation::assignments;
            using subject       = typename Computation::subject;

            comp_initializer<subject, assignments, list_size<assignments>::value, 
                    temporary_storage> :: eval_comp<Val,Local_Storage>(ls, this);
        };

        template<class Tag, class Computation, class Visitor>
        void init_computation_accept(Visitor& vis)
        {
            using assignments   = typename Computation::assignments;
            using subject       = typename Computation::subject;

            comp_initializer<subject, assignments, list_size<assignments>::value, 
                    temporary_storage> :: accept<Visitor>(vis);
        };
};

}}