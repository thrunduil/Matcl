#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              make_call_inline
//----------------------------------------------------------------------------------

struct call_array_type{};

template<class Tag, template<class Arg> class Func, 
         Integer M_arg, Integer N_arg, Mat_array Array_Arg, class Deps0>
struct make_call_inline<Tag, Func, ct_matrix<M_arg,N_arg,Array_Arg,Deps0>>
{
    //remove dependencies
    using arg_mat               = ct_matrix<M_arg,N_arg,Array_Arg, dps<>>;

    using result                = typename Func<arg_mat>::expression;
    static const Integer rows   = result::rows;
    static const Integer cols   = result::cols;    

    using new_dep               = dep<Tag, rows * cols, dep_temp>;
    using temp_deps             = dps<new_dep>;

    using array_type            = mkd::mat_temp_array<Tag, rows, cols, false>;
    using type                  = ct_matrix<rows, cols, array_type, temp_deps>;

    using child_deps            = Deps0;
    using call_array            = call_array_type;
    using call_matrix           = ct_matrix<M_arg, N_arg, call_array, child_deps>;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    friend void tag_initializer(Tag, Local_Storage& ls, Data_Provider& dp, Temp_Storage* ts)
    {
        using ret               = temp_output_mat<rows, cols, Tag>;
        using subs_context      = typename Temp_Storage::subs_context;
        using code_gen          = typename subs_context::code_gen;
        using ret_subs          = decltype(get_substitution(subs_context(), Tag()));
        using rem_subs          = list::list<ret_subs,Deps0>;
        using evaler_type       = expr_evaler_base<code_gen, ret, result, rem_subs>;       
        using val               = typename Temp_Storage::val_type;
        using evaler            = expr_evaler_impl<evaler_type, val, Temp_Storage>;

        evaler ev(ts);
        ev.eval(dp);
    };    
    template<class Visitor, class Temp_Storage>
    friend void tag_initializer_accept(Tag, Visitor& vis, Temp_Storage* ts)
    {
        using ret               = temp_output_mat<rows, cols, Tag>;
        using subs_context      = typename Temp_Storage::subs_context;
        using code_gen          = typename subs_context::code_gen;
        using ret_subs          = decltype(get_substitution(subs_context(), Tag()));
        using rem_subs          = list::list<ret_subs,Deps0>;
        using evaler_type       = expr_evaler_base<code_gen, ret, result, rem_subs>;       
        using val               = typename Temp_Storage::val_type;
        using evaler            = expr_evaler_impl<evaler_type, val, Temp_Storage>;

        vis.visit_call_inplace();

        evaler ev(ts);
        return ev.accept(vis);
    };

    template<class Subs_Context>
    friend void tag_printer(Tag, std::ostream& os, int tabs)
    {
        print_whitespace(os, tabs);
        os << "inline function with tag: ";

        Tag::print(os,details::prior_start);
        os << "; code:" << "\n";

        print_whitespace(os, tabs);
        os << std::string(80,'-') << "\n";

        using ret               = temp_output_mat<rows, cols, Tag>;
        using ret_subs          = decltype(get_substitution(Subs_Context(), Tag()));
        using code_gen          = typename Subs_Context::code_gen;
        using rem_subs          = list::list<ret_subs,Deps0>;
        using evaler_type       = expr_evaler_base<code_gen, ret, result, rem_subs>;

        evaler_type::print(os, tabs+4);

        print_whitespace(os, tabs);
        os << std::string(80,'-') << "\n";
    };

    template<class Subs_Context>
    friend child_deps get_child_deps(new_dep)
    {
        return child_deps();
    };

    template<class Subs_Context>
    friend call_matrix get_stored_matrix(Tag)
    {
        return call_matrix();
    };
};

//----------------------------------------------------------------------------------
//                              make_call_external
//----------------------------------------------------------------------------------
template<class Colon>
struct print_colon
{
    static_assert(md::dependent_false<Colon>::value, 
                "this type should not be instantiated");
};
template<>
struct print_colon<colon_all>
{
    static void eval(std::ostream& os)
    {
        os << "(:)";
    };
};
template<Integer Pos>
struct print_colon<colon<Pos>>
{
    static void eval(std::ostream& os)
    {
        os << "(" << Pos << ")";
    };
};
template<Integer Start, Integer End>
struct print_colon<colon2<Start,End>>
{
    static void eval(std::ostream& os)
    {
        os << "(" << Start << ":" << End << ")";
    };
};
template<Integer Start, Integer Step, Integer End>
struct print_colon<colon3<Start,Step,End>>
{
    static void eval(std::ostream& os)
    {
        os << "(" << Start << ":" << Step << ":" << End << ")";
    };
};

template<class Mat>
struct is_value_matrix
{
    static_assert(md::dependent_false<Mat>::value, 
                "this type should not be instantiated");
};
template<class Arr>
struct get_arg_tag
{
    static_assert(md::dependent_false<Arr>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, class Tag, class Deps>
struct is_value_matrix<ct_matrix<M,N, mkd::gen_array<Tag>,Deps>>
{
    static const Integer value = true;
};
template<class Tag>
struct get_arg_tag<mkd::gen_array<Tag>>
{
    using type = Tag;
};

template<Integer M, Integer N, class Tag, class Deps>
struct is_value_matrix<ct_matrix<M,N, mkd::output_array<Tag>,Deps>>
{
    static const Integer value = true;
};
template<class Tag>
struct get_arg_tag<mkd::output_array<Tag>>
{
    using type = Tag;
};

template<Integer M, Integer N, class Tag, Integer MR, Integer MC, class Deps>
struct is_value_matrix<ct_matrix<M,N, mkd::temp_output_array<Tag,MR,MC>,Deps>>
{
    static const Integer value = true;
};
template<class Tag, Integer MR, Integer MC>
struct get_arg_tag< mkd::temp_output_array<Tag,MR,MC>>
{
    using type = Tag;
};

template<Integer M, Integer N, class Tag, Integer Rows, Integer Cols, bool Force, class Deps>
struct is_value_matrix<ct_matrix<M,N, mkd::mat_temp_array<Tag, Rows, Cols, Force>,Deps>>
{
    static const Integer value = true;
};
template<class Tag, Integer MR, Integer MC, bool Force>
struct get_arg_tag<mkd::mat_temp_array<Tag,MR,MC,Force>>
{
    using type = Tag;
};

template<Integer M, Integer N, Mat_array Array, class Deps>
struct is_value_matrix<ct_matrix<M,N,Array,Deps>>
{
    static const Integer value = false;
};

template<class Tag, class Func, class Mat>
struct make_call_external<Tag, Func, Mat, false>
{
    static_assert(md::dependent_false<Tag>::value, 
                "argument must be value matrix");
};

template<Integer M, Integer N>
struct make_func_input
{
    static const Integer rows   = M;
    static const Integer cols   = N;

    using type = make_func_input;
};

template<Integer Rows0, Integer Cols0, Integer Mat_Rows, Integer Mat_Cols, class Colon, Integer Base_Step>
struct make_func_output
{
    static const Integer rows   = Rows0;
    static const Integer cols   = Cols0;
    static const Integer start  = colon_func::offset<Colon>::value + 1;
    static const Integer step   = colon_func::step<Colon>::value * Base_Step;

    using type = make_func_output;
};

template<class Tag, class Func, Integer M_arg, Integer N_arg, Mat_array Array_Arg, class Deps0>
struct make_call_external<Tag, Func, ct_matrix<M_arg,N_arg,Array_Arg,Deps0>, true>
{
    using return_size           = typename Func::template return_size<M_arg, N_arg>;
    static const Integer rows   = return_size::rows;
    static const Integer cols   = return_size::cols;    

    using new_dep               = dep<Tag, rows * cols, dep_temp>;
    using temp_deps             = dps<new_dep>;

    using array_type            = mkd::mat_temp_array<Tag, rows, cols,false>;
    using type                  = ct_matrix<rows, cols, array_type, temp_deps>;

    using child_deps            = Deps0;
    using call_array            = call_array_type;
    using call_matrix           = ct_matrix<M_arg, N_arg, call_array, child_deps>;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    friend void tag_initializer(Tag, Local_Storage& ls, Data_Provider& dp, Temp_Storage* ts)
    {
        using subs_context      = typename Temp_Storage::subs_context;
        using subs              = typename get_temporary_elem<Tag,subs_context>::type;
        using ret_tag           = typename subs::tag;
        using colon             = typename subs::colon;

        using arg_tag           = typename get_arg_tag<Array_Arg>::type;
        using in_dep            = dep<arg_tag,M_arg * N_arg, dep_temp>;
        using out_dep           = dep<Tag,rows*cols,dep_temp>;

        static const Integer mat_rows   = get_mat_rows<subs, rows>::value;
        static const Integer mat_cols   = get_mat_cols<subs, cols>::value;
        static const Integer out_step   = ret_tag::step;
        static const Integer out_start  = colon_func::offset<colon>::value;

        using input             = typename make_func_input<M_arg, N_arg>::type;
        using output            = typename make_func_output<rows, cols, mat_rows, 
                                        mat_cols, colon, out_step>::type;

        using evaler            = typename Func::template evaler<output, input>;
        using val               = typename Temp_Storage::val_type;

        const val* in           = &ts->get<1, data_provider, arg_tag>(dp);
        val* out                = const_cast<val*>(&ts->get<1, data_provider, Tag>(dp));

        evaler::eval(out, in);
    };    
    
    template<class Visitor, class Temp_Storage>
    friend void tag_initializer_accept(Tag, Visitor& vis, Temp_Storage* ts)
    {
        //TODO
    };

    template<class Subs_Context>
    friend void tag_printer(Tag, std::ostream& os, int tabs)
    {
        print_whitespace(os, tabs);
        os << "external call:" << "\n";

        using subs      = typename get_temporary_elem<Tag,Subs_Context>::type;
        using ret_tag   = typename subs::tag;
        using colon     = typename subs::colon;

        using arg_tag   = typename get_arg_tag<Array_Arg>::type;

        print_whitespace(os, tabs);
        os << "[" << rows << "x" << cols << "] : ";
        ret_tag::print(os, details::prior_start);
        print_colon<colon>::eval(os);

        os << " = " << Func::name() << "(";
        os << "[" << M_arg << "x" << N_arg << "] : ";
        arg_tag::print(os, details::prior_start);
        os << ")";
        os << "\n";
    };

    template<class Subs_Context>
    friend child_deps get_child_deps(new_dep)
    {
        return child_deps();
    };

    template<class Subs_Context>
    friend call_matrix get_stored_matrix(Tag)
    {
        return call_matrix();
    };

};

}}