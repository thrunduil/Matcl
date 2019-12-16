#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/TODO/utils/op_count.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              expr_evaler
//----------------------------------------------------------------------------------
template<class Array>
struct is_output_matrix
{
    static const bool value = false;
};

template<class Output_Tag, Integer M, Integer N, class Deps>
struct is_output_matrix<ct_matrix<M,N, mkd::output_array<Output_Tag>, Deps>>
{
    static const bool value = true;
};

template<class Output_Tag, Integer M, Integer N, Integer MR, Integer MC, class Deps>
struct is_output_matrix<ct_matrix<M,N, mkd::temp_output_array<Output_Tag, MR, MC>, Deps>>
{
    static const bool value = true;
};

template<class Mat>
struct is_temporary
{
    static const bool value = false;
};
template<Integer M, Integer N, class Tag, Integer Mat_Rows, Integer Mat_Cols, bool Force, class Deps>
struct is_temporary<ct_matrix<M,N,mkd::mat_temp_array<Tag, Mat_Rows, Mat_Cols, Force>,Deps>>
{
    //for now only temporary matrix are market as temporary, virtual matrices of temporaries
    //are ignored
    static const bool value = (Force == false);
};

template<class Mat, bool With_Forced>
struct is_virtual_temporary
{
    static const bool value = false;
};
template<Integer M, Integer N, class Array, class Deps, bool With_Forced>
struct is_virtual_temporary<ct_matrix<M, N, Array, Deps>,With_Forced>
{
    static const bool value = is_temporary_mat_array<Array, With_Forced>::value;
};

template<bool All_Temps, class Subs_Context, class ... Assign_List>
struct get_temp_tags_virtual
{
    //output is a temporary, which stores assignments, assignments from temporaries
    //can be removed, currently virtual array is removed only if all assigned arrays
    //are temporary
    static_assert(details::dependent_false_var<Assign_List...>::value,
                  "this case is not implemented");
};

template<class Array, class Subs_Context>
struct get_temp_tags;

template<class Tag, class Colon, bool Init>
struct modif{};

template<class Tag_, class Colon_, Integer Mat_Rows, Integer Mat_Cols>
struct modif2
{
    using tag   = Tag_;
    using colon = Colon_;

    static const Integer    rows = Mat_Rows;
    static const Integer    cols = Mat_Cols;
};

template<class Modif, Integer Old>
struct get_mat_rows
{
    static_assert(md::dependent_false<Modif>::value, 
                "this type should not be instantiated");
};
template<class Modif, Integer Old>
struct get_mat_cols
{
    static_assert(md::dependent_false<Modif>::value, 
                "this type should not be instantiated");
};
template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols, Integer Old>
struct get_mat_rows<modif2<Tag,Colon,Mat_Rows,Mat_Cols>,Old>
{
    static const Integer value = Mat_Rows;
};
template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols, Integer Old>
struct get_mat_cols<modif2<Tag,Colon,Mat_Rows,Mat_Cols>,Old>
{
    static const Integer value = Mat_Cols;
};
template<class Tag, class Colon, Integer Old>
struct get_mat_rows<modif2<Tag,Colon,-1,-1>,Old>
{
    static const Integer value = Old;
};
template<class Tag, class Colon, Integer Old>
struct get_mat_cols<modif2<Tag,Colon,-1,-1>,Old>
{
    static const Integer value = Old;
};

template<class Item1, class Colon>
struct make_tag_colon
{
    static_assert(md::dependent_false<Item1>::value,
                  "this type should not be instantiated");
};
template<class Colon1, class Colon2>
struct link_colons
{
    static_assert(md::dependent_false<Colon1>::value,
                  "this type should not be instantiated");
};
template<Integer Start, Integer Step, Integer End>
struct link_colons<colon3<Start, Step, End>, colon_all>
{
    using type = colon3<Start, Step, End>;
};

template<Integer Start1, Integer Step1, Integer End1, 
         Integer Start2, Integer Step2, Integer End2>
struct link_colons<colon3<Start1, Step1, End1>, colon3<Start2, Step2, End2>>
{
    using type = colon3<(Start2-1)*Step1 + Start1, Step1*Step2, (End2 - 1)*Step1 + Start1>;
};

template<Integer Start, Integer End>
struct link_colons<colon2<Start, End>, colon_all>
{
    using type = colon2<Start, End>;
};
template<Integer Pos>
struct link_colons<colon<Pos>, colon_all>
{
    using type = colon<Pos>;
};
template<class Colon>
struct link_colons<colon_all, Colon>
{
    using type = Colon;
};

template<class Tag1, class Colon1, bool Init, class Colon2>
struct make_tag_colon<modif<Tag1,Colon1,Init>,Colon2>
{
    using colon = typename link_colons<Colon2,Colon1>::type;
    using type  = modif<Tag1,colon,Init>;
};

template<class Tags1, class Colon_1>
struct make_tags_colon
{
    static_assert(md::dependent_false<Tags1>::value,
                  "this type should not be instantiated");
};
template<class ... Tags1, class Colon_1>
struct make_tags_colon<list::list<Tags1...>,Colon_1>
{
    using type = list::list<typename make_tag_colon<Tags1,Colon_1>::type ...>;
};
template<class Colon_1>
struct make_tags_colon<list::list<>, Colon_1>
{
    using type = list::list<>;
};

template<class Tags1, class Tags2>
struct link_tags
{
    static_assert(md::dependent_false<Tags1>::value,
                  "this type should not be instantiated");
};
template<class ... Tags1, class ... Tags2>
struct link_tags<list::list<Tags1...>,list::list<Tags2...>>
{
    using type = list::list<Tags1...,Tags2...>;
};

template<class Subs_Context, Integer M1, Integer N1, Integer M2, Integer N2, class Array_T, 
        class Colon_1, class ... Assign_List>
struct get_temp_tags_virtual<true, Subs_Context, mkd::virtual_assign_item<M1, N1, M2, N2, Array_T, Colon_1>,Assign_List...>
{
    using tags1     = typename get_temp_tags<Array_T,Subs_Context>::type;
    using tags_col  = typename make_tags_colon<tags1,Colon_1>::type;
    using tags2     = typename get_temp_tags_virtual<true,Subs_Context, Assign_List...>::type;
    using type      = typename link_tags<tags_col, tags2>::type;
};
template<class Subs_Context, class ... Assign_List>
struct get_temp_tags_virtual<false, Subs_Context, Assign_List...>
{
    using type      = list::list<>;
};
template<class Subs_Context>
struct get_temp_tags_virtual<true, Subs_Context>
{
    using type      = list::list<>;
};

template<class Array, class Subs_Context>
struct get_temp_tags
{
    using type = list::list<>;
};

template<class Virt_Tag, class ... Assign_List, class Subs_Context>
struct get_temp_tags<mkd::virtual_array<Virt_Tag, Assign_List...>, Subs_Context> 
    : get_temp_tags_virtual<is_temporary_mat_array<mkd::virtual_array<Virt_Tag, Assign_List...>, false>::value,  
                            Subs_Context, Assign_List...>
{};

template<class Tag, class Dep_Tags>
struct make_final_tags
{
    static_assert(md::dependent_false<Tag>::value,
                  "this type should not be instantiated");
};
template<class Tag>
struct make_final_tags<Tag,list::list<>>
{
    using type = list::list<modif<Tag,colon_all, true>>;
};
template<class Tag, class Item, class ... Items>
struct make_final_tags<Tag, list::list<Item, Items...>>
{
    using type = list::list<modif<Tag, colon_all, false>, Item, Items...>;
};

template<class Tag, Integer Mat_Rows, Integer Mat_Cols, class Subs_Context>
struct get_temp_tags<mkd::mat_temp_array<Tag, Mat_Rows, Mat_Cols, false>, Subs_Context>
{
    using mat       = decltype(get_stored_matrix<Subs_Context>(Tag()));
    using array     = typename mat::array_type;
    using dep_tags  = typename get_temp_tags<array, Subs_Context>::type;
    using type      = typename make_final_tags<Tag, dep_tags>::type;
};
template<class Tag, Integer Mat_Rows, Integer Mat_Cols, class Subs_Context>
struct get_temp_tags<mkd::mat_temp_array<Tag, Mat_Rows, Mat_Cols, true>, Subs_Context>
{
    using type = list::list<>;
};

template<class Temp_Tag, class Ret_Tag, class Colon, Integer Rows0, Integer Cols0, bool Init>
struct dps_modif
{
    using modified_tag  = Ret_Tag;
    using colon_type    = Colon;

    static const Integer rows   = Rows0;
    static const Integer cols   = Cols0;
    static const bool initialize= Init;
};

template <class Subs_Context, class Tag>
modif2<Tag,colon_all, -1, -1> get_substitution(Subs_Context sc, Tag tag)
{
    static_assert(std::is_base_of<temp_tag_base,Tag>::value, "tag is not derived from tag_base");
    return modif2<Tag, colon_all, -1, -1>();
};
template<class Subs_Context, class Tag>
bool tag_need_initialization(Subs_Context sc, Tag tag)
{
    (void)sc;
    (void)tag;
    return true;
};


template<class Mod, class Ret_Tag, Integer Rows, Integer Cols, class Colon, class Subs_Context>
struct make_dps_modif_1
{
    static_assert(md::dependent_false<Mod>::value,
                  "this case is not implemented");
};
template<class Tag, class Colon_in, class Ret_Tag, Integer Rows, Integer Cols, 
    bool Init, class Colon_Ret, class Subs_Context>
struct make_dps_modif_1<modif<Tag,Colon_in,Init>,Ret_Tag,Rows,Cols,Colon_Ret, Subs_Context>
{
    using colon = typename link_colons<Colon_Ret, Colon_in>::type;
    using type  = dps_modif<Tag,Ret_Tag,colon,Rows,Cols,Init>;

    friend modif2<Ret_Tag,colon,Rows,Cols> get_substitution(Subs_Context sc, Tag t)
    {
        static_assert(std::is_base_of<temp_tag_base,Tag>::value, "tag is not derived from tag_base");
        return modif2<Ret_Tag, colon, Rows,Cols>();
    };
    friend bool tag_need_initialization(Subs_Context sc, Tag t)
    {
        (void)sc;
        (void)t;
        return Init;
    };
};

template<class Temp_Tags, class Ret_Tag, Integer R, Integer C, class Colon, class Subs_Context>
struct make_dps_modifier
{
    static_assert(md::dependent_false<Temp_Tags>::value,
                  "this case is not implemented");
};
template<class ...Items, class Ret_Tag, Integer R, Integer C, class Colon, class Subs_Context>
struct make_dps_modifier<list::list<Items...>, Ret_Tag, R, C, Colon, Subs_Context>
{
    using type = list::list<typename make_dps_modif_1<Items,Ret_Tag, R, C, Colon, Subs_Context>::type ...>;
};

template<class Array, class Ret_Tag0, Integer R0, Integer C0, class Return_Subs, class Subs_Context>
struct dps_modifier_maker
{
    static_assert(md::dependent_false<Array>::value, 
                "this type should not be instantiated");
};
template<class Array, class Ret_Tag0, Integer R0, Integer C0, 
        class Ret_Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols, class Subs_Context>
struct dps_modifier_maker<Array,Ret_Tag0,R0,C0,modif2<Ret_Tag,Colon,Mat_Rows,Mat_Cols>,Subs_Context>
{
    using temp_tags = typename get_temp_tags<Array, Subs_Context>::type;
    using type      = typename make_dps_modifier<temp_tags, Ret_Tag, Mat_Rows, Mat_Cols, Colon, 
                                Subs_Context>::type;
};
template<class Array, class Ret_Tag, Integer R, Integer C, class Ret_Tag0, class Colon, class Subs_Context>
struct dps_modifier_maker<Array,Ret_Tag,R,C,modif2<Ret_Tag0,Colon,-1,-1>,Subs_Context>
{
    using temp_tags = typename get_temp_tags<Array, Subs_Context>::type;
    using type      = typename make_dps_modifier<temp_tags, Ret_Tag, R, C, colon_all,Subs_Context>::type;
};
template<class Array, class Ret_Tag, Integer R, Integer C, class Subs_Context>
struct dps_modifier_maker<Array, Ret_Tag,R,C,void,Subs_Context>
{
    using temp_tags = typename get_temp_tags<Array, Subs_Context>::type;
    using type      = typename make_dps_modifier<temp_tags, Ret_Tag, R, C, colon_all,Subs_Context>::type;
};

template<class Ret_Matrix, class Expr_Matrix, class Ret_Mat_Subs, class Subs_Context, 
        bool Is_Temp = is_temporary<Expr_Matrix>::value || is_virtual_temporary<Expr_Matrix, false>::value>
struct final_expr
{
    static_assert(md::dependent_false<Ret_Matrix>::value, 
                "this type should not be instantiated");
};

template<class Matrix_Type, class Rem_Dps>
struct add_removed_deps
{
    static_assert(md::dependent_false<Matrix_Type>::value, 
                "this type should not be instantiated");
};

template<class Tag, class Return_Subs, class Subs_Context>
struct make_empty_dps_maker
{
    static_assert(md::dependent_false<Return_Subs>::value, 
                "this type should not be instantiated");
};
template<class Tag, class Subs_Context>
struct make_empty_dps_maker<Tag,void,Subs_Context>
{
    using type          = list::list<>;
};
template<class Tag, class Colon, class Subs_Context>
struct make_empty_dps_maker<Tag, modif2<Tag, Colon, -1, -1>,Subs_Context>
{
    using type          = list::list<>;
};
template<class Tag, class Ret_Tag, class Colon, Integer Rows, Integer Cols, class Subs_Context>
struct make_empty_dps_maker<Tag, modif2<Ret_Tag, Colon, Rows, Cols>,Subs_Context>
{
    using type          = dps_modif<Tag,Ret_Tag,Colon,Rows,Cols,false>;

    friend modif2<Ret_Tag,Colon,Rows,Cols> get_substitution(Subs_Context sc, Tag t)
    {
        static_assert(std::is_base_of<temp_tag_base,Tag>::value, "tag is not derived from tag_base");
        return modif2<Ret_Tag, Colon, Rows,Cols>();
    };
    friend bool tag_need_initialization(Subs_Context sc, Tag t)
    {
        return false;
    };
};

template<Integer Rows, Integer Cols, class Array, class Deps, class Rem_Dps>
struct add_removed_deps<ct_matrix<Rows,Cols,Array,Deps>,Rem_Dps>
{
    using all_deps  = typename link_deps<Deps,Rem_Dps>::type;
    using type      = ct_matrix<Rows,Cols,Array,all_deps>;
};

template<Integer Rows, Integer Cols, class Array, class Deps>
struct add_removed_deps<ct_matrix<Rows,Cols,Array,Deps>,void>
{
    using type      = ct_matrix<Rows,Cols,Array,Deps>;
};

template<class Tag>
struct make_hidden_matrix_type
{
    using tag = Tag;
};

template<class Tag, class Ret, class Mat>
struct hide_matrix_type
{
    static_assert(is_output_matrix<Ret>::value == true, "invalid output matrix");

    using type  = make_hidden_matrix_type<Tag>;

    friend Ret get_matrix_return(type)  { return Ret(); };
    friend Mat get_matrix_input(type)   { return Mat(); };
};

template<class T>
struct get_return_matrix
{
    using type  = decltype( get_matrix_return(T()) );
};

template<class T>
struct get_input_matrix
{
    using type  = decltype( get_matrix_input(T()) );
};

template<class Deps_or_void>
struct make_dps_from_subs{};

template<>
struct make_dps_from_subs<void>
{
    using type = dps<>;
};

template<class ... Args>
struct make_dps_from_subs<dps<Args...>>
{
    using type = dps<Args...>;
};

template<Integer M, Integer N, class Ret_Tag, class Ret_DPS, class Expr_Matrix, class Ret_Mat_Subs, 
        class Subs_Context>
struct final_expr<ct_matrix<M,N, mkd::output_array<Ret_Tag>,Ret_DPS>, Expr_Matrix, Ret_Mat_Subs, Subs_Context,
        false>
{
    static_assert(M == Expr_Matrix::rows && N == Expr_Matrix::cols, "invalid assignment");

    using ret_mat           = ct_matrix<M,N, mkd::output_array<Ret_Tag>,Ret_DPS>;
    using matrix_type       = typename mat_assign<ret_mat, Expr_Matrix>::type;
    using return_subs       = typename list::elem_at_pos<Ret_Mat_Subs,0>::type;
    using removed_deps      = typename list::elem_at_pos<Ret_Mat_Subs,1>::type;
    using dps_modif_maker   = make_empty_dps_maker<Ret_Tag,return_subs,Subs_Context>;
    using matrix_type_deps  = typename add_removed_deps<matrix_type,removed_deps>::type;

    using removed_deps_dps  = typename make_dps_from_subs<removed_deps>::type;
    using additional_deps   = typename link_deps<Ret_DPS, removed_deps_dps>::type;
};

template<Integer M, Integer N, class Ret_Tag, Integer MR, Integer MC, class Ret_DPS, class Expr_Matrix, 
        class Ret_Mat_Subs, class Subs_Context>
struct final_expr<ct_matrix<M,N, mkd::temp_output_array<Ret_Tag,MR,MC>,Ret_DPS>, Expr_Matrix, Ret_Mat_Subs, 
        Subs_Context, false>
{
    static_assert(M == Expr_Matrix::rows && N == Expr_Matrix::cols, "invalid assignment");

    using ret_mat           = ct_matrix<M,N, mkd::temp_output_array<Ret_Tag,MR,MC>,Ret_DPS>;
    using matrix_type       = typename mat_assign<ret_mat, Expr_Matrix>::type;
    using return_subs       = typename list::elem_at_pos<Ret_Mat_Subs,0>::type;
    using removed_deps      = typename list::elem_at_pos<Ret_Mat_Subs,1>::type;
    using dps_modif_maker   = make_empty_dps_maker<Ret_Tag,return_subs,Subs_Context>;
    using matrix_type_deps  = typename add_removed_deps<matrix_type,removed_deps>::type;

    using removed_deps_dps  = typename make_dps_from_subs<removed_deps>::type;
    using additional_deps   = typename link_deps<Ret_DPS, removed_deps_dps>::type;
};

template<Integer M, Integer N, class Ret_Tag, class Ret_DPS, 
        Integer M2, Integer N2, class Expr_Array, class Expr_DPS, class Ret_Mat_Subs, class Subs_Context>
struct final_expr<ct_matrix<M,N, mkd::output_array<Ret_Tag>,Ret_DPS>,
                  ct_matrix<M2,N2,Expr_Array,Expr_DPS>,Ret_Mat_Subs, Subs_Context, true>
{
    //final expression is empty but temporary storage is in output Array
    static_assert(M == M2 && N == N2, "invalid assignment");
    
    using array_type        = mkd::empty_array<Ret_Tag>;
    using dps_all           = typename link_deps<Ret_DPS, Expr_DPS>::type;
    using matrix_type       = ct_matrix<M, N, array_type, dps_all>;
    using return_subs       = typename list::elem_at_pos<Ret_Mat_Subs,0>::type;
    using removed_deps      = typename list::elem_at_pos<Ret_Mat_Subs,1>::type;

    using dps_modif_maker   = dps_modifier_maker<Expr_Array, Ret_Tag, M, N, return_subs, Subs_Context>;
    using matrix_type_deps  = typename add_removed_deps<matrix_type,removed_deps>::type;

    using removed_deps_dps  = typename make_dps_from_subs<removed_deps>::type;
    using additional_deps   = typename link_deps<Ret_DPS, removed_deps_dps>::type;
};

template<Integer M, Integer N, class Ret_Tag, Integer MR, Integer MC, class Ret_DPS, 
        Integer M2, Integer N2, class Expr_Array, class Expr_DPS, class Ret_Mat_Subs, class Subs_Context>
struct final_expr<ct_matrix<M,N, mkd::temp_output_array<Ret_Tag, MR, MC>,Ret_DPS>,
                  ct_matrix<M2,N2,Expr_Array,Expr_DPS>,Ret_Mat_Subs, Subs_Context, true>
{
    //final expression is empty but temporary storage is in output Array
    static_assert(M == M2 && N == N2, "invalid assignment");
    
    using array_type        = mkd::empty_array<Ret_Tag>;
    using matrix_type       = ct_matrix<M, N, array_type, Expr_DPS>;
    using return_subs       = typename list::elem_at_pos<Ret_Mat_Subs,0>::type;
    using removed_deps      = typename list::elem_at_pos<Ret_Mat_Subs,1>::type;

    using dps_modif_maker   = dps_modifier_maker<Expr_Array, Ret_Tag, M, N, return_subs, Subs_Context>;
    using matrix_type_deps  = typename add_removed_deps<matrix_type,removed_deps>::type;
    using additional_deps   = removed_deps;
};

template<class Code_Gen, class Hidden_Matrix, class Ret_Mat_Subs0 = list::list<void,void>>
struct expr_evaler_base
{
    public:
        using ret_matrix        = typename get_return_matrix<Hidden_Matrix>::type;
        using in_matrix         = typename get_input_matrix<Hidden_Matrix>::type;
        using in_deps           = typename in_matrix::dps_type;
        using code_gen          = Code_Gen;
        using return_subs       = typename list::elem_at_pos<Ret_Mat_Subs0,0>::type;
        using removed_subs      = typename list::elem_at_pos<Ret_Mat_Subs0,1>::type;
        using expr_deps         = in_deps;

        using expr_matrix       = in_matrix;        
        using final_expr_t      = final_expr<ret_matrix, expr_matrix, Ret_Mat_Subs0, expr_evaler_base>;
        using final_matrix      = typename final_expr_t::matrix_type;                
        using add_deps          = typename final_expr_t::additional_deps;        
        using deps_collected    = typename collect_deps<expr_evaler_base,expr_deps,dps<>>::type;        

        using final_matrix_deps = typename final_expr_t::matrix_type_deps;
        using dps_mod_maker     = typename final_expr_t::dps_modif_maker;
        using dps_modifier      = typename dps_mod_maker::type;
        using deps_temp         = typename collect_deps_temp<deps_collected,dps<>>::type;
        using deps_all          = typename link_deps<deps_collected, add_deps>::type;        

    public:
        static void print(std::ostream& os, int tabs)
        {
            final_matrix::print<expr_evaler_base>(os, tabs);
        };
};

template<class Evaler_Base, class Val, class Parent_Storage = empty_storage<Val>>
struct expr_evaler_impl
{    
    private:        
        using base          = Evaler_Base;
        using parent        = Parent_Storage;
        using temp          = temporary_storage<Val, base, parent>;

        using final_matrix  = typename base::final_matrix_deps;

        temp  m_temp;

    public:
        expr_evaler_impl()                      : m_temp(nullptr) {};
        expr_evaler_impl(Parent_Storage* st)    : m_temp(st) {};

        template<class Data_Provider>
        inline_lev_root
        void eval(Data_Provider& dp)
        {
            using subs_context  = typename temp::subs_context;            
            using loc_storage   = local_storage<Val, Data_Provider, subs_context>;

            loc_storage ls;
            ls.init(dp,&m_temp);

            m_temp.init_storage(ls, dp);            
            expr_evaler_elems<Val, final_matrix>::eval(ls);
        };

        template<class Visitor>
        void accept(Visitor& vis)
        {
            m_temp.accept<Visitor>(vis);
            expr_evaler_elems<Val, final_matrix>::accept<Visitor>(vis);
        };

        static void print(std::ostream& os, int tabs)
        {
            Evaler_Base::print(os, tabs);
        };

        op_count get_op_count()
        {
            op_count vis;
            accept<op_count>(vis);
            return vis;
        };
};

template<class Tag, class Code_Gen, class Ret_Matrix, class Mat, class Val>
struct expr_evaler: expr_evaler_impl<expr_evaler_base<Code_Gen, 
           typename hide_matrix_type<Tag, Ret_Matrix, Mat>::type >, Val>
{};

template<class Val, class Mat>
struct expr_evaler_elems_expand
{};
template<class Val, class Mat>
struct expr_evaler_elems
{};

template<class Val, Integer M, Integer N, class Array, class Deps, class ... Elems>
struct expr_evaler_elems_expand<Val, list::list<ct_matrix<M,N,Array,Deps>, Elems...>>
{
    using matrix_type   = ct_matrix<M,N,Array,Deps>;

    template<class Local_Storage>
    inline_initializer
    static void eval(Local_Storage& ls)
    {        
        using subs_context          = typename Local_Storage::subs_context;
        static const bool is_cont   =  simd_enable<subs_context,matrix_type>::value
                                    && (N == 1);
        using is_cont_t             = std::integral_constant<bool, is_cont>;

        eval_impl(ls,is_cont_t());

        expr_evaler_elems_expand<Val,list::list<Elems...>>::eval(ls);
    };
    template<class Local_Storage>
    inline_initializer
    static void eval_impl(Local_Storage& ls, std::true_type)
    {
        using assign_info   = typename make_assign_info<matrix_type>::type;
        using rhs           = typename assign_info::rhs;
        using rhs_array     = typename rhs::array_type;
        using elem          = typename rhs_array::template get_element<1,1>::type;
        using subs_context  = typename Local_Storage::subs_context;
        using loop_context  = typename make_loop_context<elem>::type;
        using ret_colon     = typename assign_info::colon_type;
        using ret_dep       = typename assign_info::return_dep;
        using data_provider = typename Local_Storage::data_provider;

        static const Integer ret_step   = assign_info::return_step;
        static const Integer ret_offset = assign_info::return_offset;
        using root_align                = assign_info::ret_root_align;

        using align_type                = typename link_alignment<root_align,
                                            get_offset_alignment<Val,ret_offset,ret_step>::type>::type;

        loop_evaler<loop_context, M, elem, ret_dep, align_type, ret_step, ret_offset>::eval<Val>(ls);

        return;
    };

    template<class Local_Storage>
    inline_initializer
    static void eval_impl(Local_Storage& ls, std::false_type)
    {
        return expr_evaler_elems_cols<M,1,N,N,Array,Val>::eval(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        expr_evaler_elems_cols<M,1,N,N,Array,Val>::accept<Visitor>(vis);
        expr_evaler_elems_expand<Val,list::list<Elems...>>::accept<Visitor>(vis);
    }
};
template<class Val, Integer M, Integer N, class Ret_Tag, class Deps, class ... Elems>
struct expr_evaler_elems_expand<Val, list::list<ct_matrix<M,N, mkd::empty_array<Ret_Tag>,Deps>, Elems...>>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(Local_Storage& ls)
    {
        //nothing to do
        expr_evaler_elems_expand<Val,list::list<Elems...>>::eval(ls);
    }
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        expr_evaler_elems_expand<Val,list::list<Elems...>>::accept<Visitor>(vis);
    }
};
template<class Val>
struct expr_evaler_elems_expand<Val, list::list<>>
{
    template<class Local_Storage>
    static void eval(Local_Storage& ls)
    {
        (void)ls;
        //nothing to do
    }
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        (void)vis;
    }
};

template<class Val, Integer M, Integer N, class Array, class Deps>
struct expr_evaler_elems<Val, ct_matrix<M,N,Array,Deps>>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(Local_Storage& ls)
    {
        using matrix_type   = ct_matrix<M,N,Array,Deps>;
        using expand_vm     = typename mkd::expand_virtual_matrix<matrix_type>::type;

        return expr_evaler_elems_expand<Val,expand_vm>::eval(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using matrix_type   = ct_matrix<M,N,Array,Deps>;
        using expand_vm     = typename mkd::expand_virtual_matrix<matrix_type>::type;

        return expr_evaler_elems_expand<Val,expand_vm>::accept(vis);
    }
};

template<Integer First_Row, Integer Last_Row, Integer Length, class Array, class Val, Integer Col>
struct expr_evaler_elems_row
{
    template<class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        static const Integer last_row_1 = First_Row + Length/2 - 1;

        expr_evaler_elems_row<First_Row, last_row_1, Length/2, Array, Val, Col>
            ::eval(ls);

        expr_evaler_elems_row<last_row_1+1, Last_Row, Length - Length/2, Array, Val, Col>
            ::eval(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer last_row_1 = First_Row + Length/2 - 1;

        expr_evaler_elems_row<First_Row, last_row_1, Length/2, Array, Val, Col>
            ::accept<Visitor>(vis);

        expr_evaler_elems_row<last_row_1+1, Last_Row, Length - Length/2, Array, Val, Col>
            ::accept<Visitor>(vis);
    };
};
template<Integer First_Row, Integer Last_Row, class Array, class Val, Integer Col>
struct expr_evaler_elems_row<First_Row, Last_Row, 0, Array,Val,Col>
{
    template<class Local_Storage>
    static void eval(const Local_Storage& ls)
    {};
    template<class Visitor>
    static void accept(Visitor& vis)
    {};
};
template<Integer First_Row, Integer Last_Row, class Array, class Val, Integer Col>
struct expr_evaler_elems_row<First_Row, Last_Row, 1, Array,Val,Col>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        using elem1 = typename Array::template get_element<First_Row,Col>::type;
        elem1::eval<Val>(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using elem1 = typename Array :: template get_element<First_Row,Col>::type;
        elem1::accept<Visitor>(vis);
    };
};
template<Integer First_Row, Integer Last_Row, class Array, class Val, Integer Col>
struct expr_evaler_elems_row<First_Row, Last_Row, 2, Array,Val,Col>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        using elem1     = typename Array :: template get_element<First_Row, Col>::type;
        elem1::eval<Val>(ls);

        using elem2     = Array :: template get_element<First_Row+1, Col>::type;
        elem2::eval<Val>(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using elem1     = typename Array :: template get_element<First_Row, Col>::type;
        elem1::accept<Visitor>(vis);

        using elem2     = typename Array :: template get_element<First_Row+1,Col>::type;
        elem2::accept<Visitor>(vis);
    };
};

template<Integer M, Integer First_Col, Integer Last_Col, Integer Length, class Array, class Val>
struct expr_evaler_elems_cols
{
    template<class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        static const Integer last_col_1 = First_Col + Length/2 - 1;

        expr_evaler_elems_cols<M, First_Col, last_col_1, Length/2, Array, Val>
            ::eval(ls);

        expr_evaler_elems_cols<M,last_col_1+1, Last_Col, Length - Length/2, Array, Val>
            ::eval(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer last_col_1 = First_Col + Length/2 - 1;

        expr_evaler_elems_cols<M, First_Col, last_col_1, Length/2, Array, Val>
            ::accept<Visitor>(vis);

        expr_evaler_elems_cols<M,last_col_1+1, Last_Col, Length - Length/2, Array, Val>
            ::accept<Visitor>(vis);
    }
};
template<Integer M, Integer First_Col, Integer Last_Col, class Array, class Val>
struct expr_evaler_elems_cols<M, First_Col, Last_Col, 0, Array, Val>
{
    template<class Local_Storage>
    static void eval(const Local_Storage&)
    {};
    template<class Visitor>
    static void accept(Visitor& vis)
    {};
};
template<Integer M, Integer First_Col, Integer Last_Col, class Array, class Val>
struct expr_evaler_elems_cols<M, First_Col, Last_Col, 1, Array, Val>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col>::eval(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col>::accept<Visitor>(vis);
    };
};
template<Integer M, Integer First_Col, Integer Last_Col, class Array, class Val>
struct expr_evaler_elems_cols<M, First_Col, Last_Col, 2, Array, Val>
{
    template<class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col>::eval(ls);
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col+1>::eval(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col>::accept<Visitor>(vis);
        expr_evaler_elems_row<1,M,M,Array,Val,First_Col+1>::accept<Visitor>(vis);
    };
};

}}