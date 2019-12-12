#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"
#include "mkgen/TODO/evaler/loop_promotion.h"
#include "matcl-simd/simd.h"
#include "mkgen/TODO/evaler/local_storage.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              storage_initializer
//----------------------------------------------------------------------------------
template<class Val, Integer First_Row, Integer Last_Row, Integer Length, Integer Col>
struct storage_initializer_rows
{
    template<class Storage, class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        static const Integer last_row_1 = First_Row + Length/2 - 1;

        storage_initializer_rows<Val, First_Row, last_row_1, Length/2, Col>
            ::eval<Storage>(ls);

        storage_initializer_rows<Val, last_row_1+1, Last_Row, Length - Length/2, Col>
            ::eval<Storage>(ls);
    };
    
    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer last_row_1 = First_Row + Length/2 - 1;

        storage_initializer_rows<Val, First_Row, last_row_1, Length/2, Col>
            ::accept<Storage>(vis);

        storage_initializer_rows<Val, last_row_1+1, Last_Row, Length - Length/2, Col>
            ::accept<Storage>(vis);
    };
};

template<class Val, Integer First_Row, Integer Last_Row, Integer Col>
struct storage_initializer_rows<Val, First_Row,Last_Row,0,Col>
{
    template<class Storage, class Local_Storage>
    static void eval(const Local_Storage& )
    {};
    
    template<class Visitor>
    static void accept(Visitor& vis)
    {};
};
template<class Val, Integer First_Row, Integer Last_Row, Integer Col>
struct storage_initializer_rows<Val, First_Row,Last_Row,1,Col>
{
    template<class Storage, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        Storage::init_elem<First_Row,Col>(ls);
    };

    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        Storage::init_elem_accept<First_Row,Col>(vis);
    };
};
template<class Val, Integer First_Row, Integer Last_Row, Integer Col>
struct storage_initializer_rows<Val, First_Row,Last_Row,2,Col>
{
    template<class Storage, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        Storage::init_elem<First_Row,Col>(ls);
        Storage::init_elem<First_Row+1,Col>(ls);
    };

    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        Storage::init_elem_accept<First_Row,Col>(vis);
        Storage::init_elem_accept<First_Row+1,Col>(vis);
    };
};

template<class Val, Integer First_Row, Integer Last_Row, Integer First_Col, Integer Last_Col, Integer Length>
struct storage_initializer_cols
{
    template<class Storage, class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        static const Integer last_col_1 = First_Col + Length/2 - 1;

        storage_initializer_cols<Val, First_Row, Last_Row, First_Col, last_col_1, Length/2>
            ::eval<Storage>(ls);

        storage_initializer_cols<Val, First_Row,Last_Row, last_col_1+1, Last_Col, Length - Length/2>
            ::eval<Storage>(ls);
    };

    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer last_col_1 = First_Col + Length/2 - 1;

        storage_initializer_cols<Val, First_Row, Last_Row, First_Col, last_col_1, Length/2>
            ::accept<Storage>(vis);

        storage_initializer_cols<Val, First_Row,Last_Row, last_col_1+1, Last_Col, Length - Length/2>
            ::accept<Storage>(vis);
    };
};
template<class Val, Integer First_Row, Integer Last_Row, Integer First_Col, Integer Last_Col>
struct storage_initializer_cols<Val, First_Row,Last_Row,First_Col,Last_Col,0>
{
    template<class Storage, class Local_Storage>
    static void eval(const Local_Storage& ls)
    {};
    
    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {};
};
template<class Val, Integer First_Row, Integer Last_Row, Integer First_Col, Integer Last_Col>
struct storage_initializer_cols<Val,First_Row,Last_Row,First_Col,Last_Col,1>
{
    template<class Storage, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        static const Integer length = Last_Row - First_Row + 1;
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col>::eval<Storage>(ls);
    };

    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer length = Last_Row - First_Row + 1;
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col>::accept<Storage>(vis);
    };
};
template<class Val, Integer First_Row, Integer Last_Row, Integer First_Col, Integer Last_Col>
struct storage_initializer_cols<Val,First_Row,Last_Row,First_Col,Last_Col,2>
{
    template<class Storage, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        static const Integer length = Last_Row - First_Row + 1;
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col>::eval<Storage>(ls);
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col+1>::eval<Storage>(ls);
    };

    template<class Storage, class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer length = Last_Row - First_Row + 1;
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col>
            ::accept<Storage>(vis);
        storage_initializer_rows<Val, First_Row,Last_Row,length,First_Col+1>
            ::accept<Storage>(vis);
    };
};

template<class Val, class expand_vm, class Colon_Type, class Ret_Dep, Integer Offset, Integer Step>
struct storage_initializer_expand
{};

template<class Val, class Colon_Type, class Ret_Dep, Integer Offset, Integer Step>
struct storage_initializer_expand<Val,list::list<>,Colon_Type,Ret_Dep, Offset,Step>
{
    template<class Array_Aligned, class Local_Storage>
    static void eval(Local_Storage& ls)
    {
        (void)ls;
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        (void)vis;
    };
};

template<class Val, class Mat_Colon, class ... Mats, class Colon_Type, class Ret_Dep, 
    Integer Offset, Integer Step>
struct storage_initializer_expand<Val, list::list<Mat_Colon, Mats...>,Colon_Type,Ret_Dep,Offset,Step>
{
    using colon_assign              = typename list::elem_at_pos<Mat_Colon,0>::type;
    using matrix_type               = typename list::elem_at_pos<Mat_Colon,1>::type;
    static const Integer rows       = matrix_type::rows;
    static const Integer cols       = matrix_type::cols;

    using array_type                = typename matrix_type::array_type;

    template<class Array_Aligned, class Local_Storage>
    inline_initializer
    static void eval(Local_Storage& ls)
    {
        using subs_context          = typename Local_Storage::subs_context;        
        static const bool is_cont   = (cols == 1)
                                    && simd_enable<subs_context,matrix_type>::value;
        using is_cont_t             = std::integral_constant<bool,is_cont>;

        eval_impl<Array_Aligned>(is_cont_t(), ls);

        using list_type             = list::list<Mats...>;


        using storate_initializer   = storage_initializer_expand
                                        <Val, list_type, Colon_Type,Ret_Dep, Offset,Step>;

        
        storate_initializer::eval<Array_Aligned>(ls);
    }

    template<class Array_Aligned, class Local_Storage>
    inline_initializer
    static void eval_impl(std::false_type, Local_Storage& ls)
    {
        storage_initializer_cols<Val, 1,rows,1,cols, cols>
            ::eval<storage_initializer_expand>(ls);
    };

    template<class Array_Aligned, class Local_Storage>
    inline_initializer
    static void eval_impl(std::true_type, Local_Storage& ls)
    {
        using elem          = typename get_array_elem<array_type,1,1>::type;
        using subs_context  = typename Local_Storage::subs_context;
        using loop_context  = typename make_loop_context<elem>::type;

        static const Integer off_c  = get_offset_colon<colon_assign>::value;
        static const Integer step_c = get_step_colon<colon_assign>::value;
        static const Integer pos    = get_pos_colon<1,Colon_Type>::value;        
        static const Integer offset = Step * (off_c + (pos-1)*step_c) + Offset; 
        static const Integer step   = step_c * Step * get_step_colon<Colon_Type>::value;

        using align_t               = typename link_alignment<Array_Aligned, 
                                        typename get_offset_alignment<Val,offset,step>::type> :: type;

        loop_evaler<loop_context, rows, elem, Ret_Dep,align_t, step, offset>::eval<Val>(ls);

        return;
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        storage_initializer_cols<Val, 1,rows,1,cols, cols>
            ::accept<storage_initializer_expand>(vis);

        storage_initializer_expand<Val,list::list<Mats...>,Colon_Type,Ret_Dep,Offset,Step>::accept(vis);
    };

    template<Integer Row, Integer Col, class Local_Storate>
    inline_initializer
    static void init_elem(const Local_Storate& ls)
    {        
        using elem = typename get_array_elem<array_type,Row,Col>::type;
        static const Integer pos0   = (Col - 1)*rows + Row;
        static const Integer pos    = get_pos_colon<pos0,Colon_Type>::value;

        static const Integer off_c  = get_offset_colon<colon_assign>::value;
        static const Integer step_c = get_step_colon<colon_assign>::value;

        static const Integer offset = Step * (off_c + (pos-1) * step_c) + Offset;

        Val* arr                    = const_cast<Val*>(ls.get_array<Ret_Dep>());
        arr[offset]                 = elem::eval<Val>(ls);
    };

    template<Integer Row, Integer Col, class Visitor>
    static void init_elem_accept(Visitor& vis)
    {        
        using elem = typename get_array_elem<array_type,Row,Col>::type;
        elem::accept<Visitor>(vis);
        vis.visit_store();
    };
};

template<class Stored_Matrix, class Val, class Colon_Type, class Ret_Dep, Integer Offset, Integer Step>
struct storage_initializer
{
    template<class Array_Align, class Local_Storate>
    inline_initializer
    static void eval(Local_Storate& ls)
    {
        using expand_vm             = typename expand_virtual_matrix2<Stored_Matrix>::type;

        return storage_initializer_expand<Val,expand_vm, Colon_Type, Ret_Dep, Offset, Step>
            ::eval<Array_Align>(ls);
    }

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using expand_vm             = typename expand_virtual_matrix2<Stored_Matrix>::type;

        return storage_initializer_expand<Val,expand_vm, Colon_Type, Ret_Dep, Offset, Step>::accept(vis);
    };
};

//----------------------------------------------------------------------------------
//                              comp_initializer
//----------------------------------------------------------------------------------
template<class Subject, class Assignments, Integer Length>
struct comp_initializer_impl
{
    using first_half    = comp_initializer_impl<Subject, Assignments, Length / 2>;
    using rem_args_1    = typename first_half::remaining_args;

    using second_half   = comp_initializer_impl<Subject, rem_args_1, Length - Length / 2>;
    using rem_args_2    = typename second_half::remaining_args;

    using remaining_args = rem_args_2;

    template<class Val, class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        first_half::eval<Val>(ls);
        second_half::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        first_half::accept(vis);
        second_half::accept(vis);
    };
};
template<class Subject, class Assignments>
struct comp_initializer_impl<Subject, Assignments,0>
{
    using remaining_args    = Assignments;

    template<class Val, class Local_Storage>
    static void eval(const Local_Storage& ls)
    {};

    template<class Visitor>
    static void accept(Visitor& vis)
    {};
};

template<class Subject, class Assign_Type>
struct comp_initializer_1
{};
template<class Subject, Integer Pos1, class Scal1>
struct comp_initializer_1<Subject,assign_colon_scal<Pos1,Scal1>>
{
    template<class Val, class Local_Storage>
    inline_lev_1
    static void eval(const Local_Storage& ls)
    {
        static const Integer rows   = Subject::rows;
        static const Integer col    = (Pos1 - 1) / rows + 1;
        static const Integer row    = (Pos1 - 1) % rows + 1;

        Val val     = Scal1::eval<Val>(ls);
        assign_elem<Subject, row, col>::eval<Val>(val, ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer rows   = Subject::rows;
        static const Integer col    = (Pos1 - 1) / rows + 1;
        static const Integer row    = (Pos1 - 1) % rows + 1;

        Scal1::accept<Visitor>(vis);
        assign_elem<Subject, row, col>::accept<Visitor>(vis);
    };
};

template<class Subject, class Colon, class Mat, Integer Start, Integer Size>
struct comp_initializer_1_impl
{
    template<class Val, class Local_Storage>
    inline_expr_split
    static void eval(const Local_Storage& ls)
    {
        static const Integer size_half  = Size / 2;

        comp_initializer_1_impl<Subject,Colon,Mat,Start,size_half>::eval<Val>(ls);
        comp_initializer_1_impl<Subject,Colon,Mat,Start+size_half,Size - size_half>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer size_half  = Size / 2;

        comp_initializer_1_impl<Subject,Colon,Mat,Start,size_half>::accept(vis);
        comp_initializer_1_impl<Subject,Colon,Mat,Start+size_half,Size - size_half>::accept(vis);
    };
};
template<class Subject, class Colon, class Mat, Integer Start>
struct comp_initializer_1_impl<Subject,Colon,Mat,Start,0>
{
    template<class Val, class Local_Storage>
    static void eval(const Local_Storage& ls)
    {};

    template<class Visitor>
    static void accept(Visitor& vis)
    {};
};
template<class Subject, class Colon, class Mat, Integer Start>
struct comp_initializer_1_impl<Subject,Colon,Mat,Start,1>
{
    template<class Val, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        static const Integer pos    = get_pos_colon<Start,Colon>::value;
        using mat_array             = typename Mat::array_type;
        using elem                  = typename get_array_elem<mat_array,Start,1>::type;

        comp_initializer_1<Subject,assign_colon_scal<pos,elem>>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer pos    = get_pos_colon<Start,Colon>::value;
        using mat_array             = typename Mat::array_type;
        using elem                  = typename get_array_elem<mat_array,Start,1>::type;

        comp_initializer_1<Subject,assign_colon_scal<pos,elem>>::accept(vis);
    };
};
template<class Subject, class Colon, class Mat, Integer Start>
struct comp_initializer_1_impl<Subject,Colon,Mat,Start,2>
{
    template<class Val, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        using mat_array             = typename Mat::array_type;

        static const Integer pos1   = get_pos_colon<Start,Colon>::value;
        using elem1                 = typename get_array_elem<mat_array,Start,1>::type;
        comp_initializer_1<Subject,assign_colon_scal<pos1,elem1>>::eval<Val>(ls);

        static const Integer pos2   = get_pos_colon<Start+1,Colon>::value;
        using elem2                 = typename get_array_elem<mat_array,Start+1,1>::type;
        comp_initializer_1<Subject,assign_colon_scal<pos2,elem2>>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using mat_array             = typename Mat::array_type;

        static const Integer pos1   = get_pos_colon<Start,Colon>::value;
        using elem1                 = typename get_array_elem<mat_array,Start,1>::type;
        comp_initializer_1<Subject,assign_colon_scal<pos1,elem1>>::accept(vis);

        static const Integer pos2   = get_pos_colon<Start+1,Colon>::value;
        using elem2                 = typename get_array_elem<mat_array,Start+1,1>::type;
        comp_initializer_1<Subject,assign_colon_scal<pos2,elem2>>::accept(vis);
    };
};

template<class Subject, class Colon, class Mat>
struct comp_initializer_1<Subject,assign_colon<Colon,Mat>>
{
    template<class Val, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        static const Integer size   = get_size_colon<Colon,Subject::size>::value;
        comp_initializer_1_impl<Subject,Colon,Mat,1,size>::eval<Val>(ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        static const Integer size   = get_size_colon<Colon,Subject::size>::value;
        comp_initializer_1_impl<Subject,Colon,Mat,1,size>::accept(vis);
    };
};

template<class Subject, class Assign_Type, class ... Items>
struct comp_initializer_impl<Subject, list::list<Assign_Type, Items...>, 1>
{
    using remaining_args    = list::list<Items...>;

    template<class Val, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        comp_initializer_1<Subject,Assign_Type>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        comp_initializer_1<Subject,Assign_Type>::accept<Visitor>(vis);
    };
};
template<class Subject, class Assign_Type_1, class Assign_Type_2, class ... Items>
struct comp_initializer_impl<Subject, list::list<Assign_Type_1, Assign_Type_2, Items...>, 2>
{
    using remaining_args = list::list<Items...>;

    template<class Val, class Local_Storage>
    inline_initializer
    static void eval(const Local_Storage& ls)
    {
        comp_initializer_1<Subject,Assign_Type_1>::eval<Val>(ls);
        comp_initializer_1<Subject,Assign_Type_2>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        comp_initializer_1<Subject,Assign_Type_1>::accept<Visitor>(vis);
        comp_initializer_1<Subject,Assign_Type_2>::accept<Visitor>(vis);
    };
};

template<class Subject, class Assignments, Integer Length, class Temp_Storage>
struct comp_initializer
{
    template<class Val, class Local_Storage>
    inline_initializer
    static void eval_comp(Local_Storage& ls, Temp_Storage* ev)
    {
        (void)ev;
        comp_initializer_impl<Subject,Assignments,Length>::eval<Val>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        comp_initializer_impl<Subject,Assignments,Length>::accept(vis);
    };
};

}}