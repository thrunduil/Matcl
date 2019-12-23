#pragma once

#include "mkgen/mkgen.h"
#include "matcl-fft/matcl_fft_config.h"
#include "matcl-fft/matcl_dct.h"
#include "cos_dct.h"

#include <iostream>

#define CALL_FUNC_dct1 0
#define CALL_FUNC_dct2 0
#define CALL_FUNC_dct4 0

template<int ... Vals>
struct list_int{};

namespace matcl { namespace fft
{

namespace mk = mkgen;

//---------------------------------------------------------------------
//                              TAGS
//---------------------------------------------------------------------
template<int Lev, int Num>
static void print_tmp(std::ostream& os, const std::string& name)
{
    os << name << "<" << Lev << "," << Num << ">";
};

template<class Tag, int Lev, int Num>
struct temp : mk::temp_tag_base
{
    static void print(std::ostream& os, int) { print_tmp<Lev,Num>(os, get_string<Tag>::eval()); }
};

template<class Tag, int Lev, int Num>
struct expr_tag{};

template<Integer M>
struct tag_virt {};

template<Integer M>
struct tag_tmp1 : mk::temp_tag_base
{     static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp1"); }
};
template<Integer M>
struct tag_tmp2 : mk::temp_tag_base
{ 
    static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp2"); }
};
template<Integer M>
struct tag_tmp3 : mk::temp_tag_base
{ 
    static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp3"); }
};

template<Integer M>
struct tag_tmp4 : mk::temp_tag_base
{ 
    static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp4"); }
};

template<Integer M>
struct tmp_dct4_1 : mk::temp_tag_base
{ 
    static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp_dct4_1"); }
};
template<Integer M>
struct tmp_dct4_2 : mk::temp_tag_base
{ 
    static void print(std::ostream& os,int) { print_tmp<M>(os, "tmp_dct4_2"); }
};

template<Integer Code, Integer Lev>
struct tag_test
{
    static std::string name()
    {
        std::ostringstream out;
        out << "tmp_" << Code << "_" << Lev;
        return out.str();
    };
};

template<class Seed, Integer Lev>
struct add_level
{};

template<int... Vals, Integer Lev>
struct add_level<list_int<Vals...>, Lev>
{
    using type = list_int<Vals...,Lev>;
};

template<int ... Vals>
struct print_list
{};
template<int Val, int ... Vals>
struct print_list<Val,Vals...>
{
    static void eval(std::ostream& os)
    {
        os << Val;
        if constexpr(sizeof...(Vals) > 0)
        {
            os << ",";
        };
        print_list<Vals...>::eval(os);
    };
};
template<>
struct print_list<>
{
    static void eval(std::ostream& os)
    {
        (void)os;
    };
};


template<class Seed>
struct get_string
{};
template<int ... Vals>
struct get_string<list_int<Vals...>>
{
    static std::string eval()
    {
        std::ostringstream os;
        os << "T[";
        print_list<Vals...>::eval(os);
        os << "]";
        return os.str();
    };
};

template<Integer M, class Scal>
struct tag_cos_dct1
{
    static const bool           is_continuous   = true; 
    using                       root_align_type = mk::align_full; 
    static const Integer        step            = 1; 

    static constexpr Integer get_offset(Integer Row, Integer Col)
    {
        return Row - 1 + (Col - 1) * M;
    };

    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "cos_dct1" << "<" << M << ">";
    };

    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        (void)dp;
        return cos_table_dct1<M,Scal>::get().m_table + Offset;
    }
};

template<Integer M, class Scal>
struct tag_cos_dct2
{
    static const bool           is_continuous   = true; 
    using                       root_align_type = mk::align_full; 
    static const Integer        step            = 1; 

    static constexpr Integer get_offset(Integer Row, Integer Col)
    {
        return Row - 1 + (Col - 1) * M;
    };

    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "cos_dct2" << "<" << M << ">";
    };

    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        (void)dp;
        return cos_table_dct2<M,Scal>::get().m_table + Offset;
    }
};

template<Integer M, class Scal>
struct tag_cos_dct3
{
    static const bool           is_continuous   = true; 
    using                       root_align_type = mk::align_full; 
    static const Integer        step            = 1; 

    static constexpr Integer get_offset(Integer Row, Integer Col)
    {
        return Row - 1 + (Col - 1) * M;
    };

    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "cos_dct3" << "<" << M << ">";
    };

    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        (void)dp;
        return cos_table_dct3<M,Scal>::get().m_table + Offset;
    }
};

template<Integer M, class Scal>
struct tag_cos_dct4
{
    static const bool           is_continuous   = true; 
    using                       root_align_type = mk::align_full; 
    static const Integer        step            = 1; 

    static constexpr Integer get_offset(Integer Row, Integer Col)
    {
        return Row - 1 + (Col - 1) * M;
    };

    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "cos_dct4" << "<" << M << ">";
    };

    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        (void)dp;
        return cos_table_dct4<M,Scal>::get().m_table + Offset;
    }
};

//---------------------------------------------------------------------
//                              HELPERS
//---------------------------------------------------------------------
template<class Mat, Integer Off, Integer Min_Length>
struct allow_recursive
{
    static const Integer length = Mat::rows + Off;
    static const bool value = (length % 2) == 0 && (length > Min_Length);
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, 
        bool Rec = (allow_recursive<Matrix_Type,-1,Config::min_rec_length_dct1>::value)> 
struct eval_dct1;

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, 
        bool Rec = (allow_recursive<Matrix_Type,0,Config::min_rec_length_dct2>::value) >
struct eval_dct2;

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, 
        bool Rec = (allow_recursive<Matrix_Type,0,Config::min_rec_length_dct3>::value) >
struct eval_dct3;

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, 
        bool Rec = (allow_recursive<Matrix_Type,0,Config::min_rec_length_dct4>::value)>
struct eval_dct4;

template<Integer M, class Scal> struct tag_dct1     { static const bool is_continuous = false; };
template<Integer M, class Scal> struct tag_dct2     { static const bool is_continuous = false; };
template<Integer M, class Scal> struct tag_dct3     { static const bool is_continuous = false; };
template<Integer M, class Scal> struct tag_dct4     { static const bool is_continuous = false; };
template<Integer M>             struct tag_dct4_2   { static const bool is_continuous = true; };

/*
//TODO: reimplement
template<Integer M, class Scal, class Tag, Integer Row, Integer Col>
struct mk::make_expr_ufunc<tag_dct1<M,Scal>, mkd::element<Tag, Row, Col>>
{
    //C_k     = sum_{n=1}^{N-2} x_n*cos(pi*n*k/(N-1))
    static const Integer ind        = (Row-1) * (Col-1);
    static const Integer step       = 1;
    using info                      = matcl::fft::cos_info<ind, M-1>;

    static const Integer base_ind   = info::base_index;
    static const bool negate        = info::negate;
    static const bool zero_v        = info::zero;
    static const bool one_v         = info::one;
    static const bool mone_v        = info::mone;
    static const bool half_v        = info::half;

    using type                      = typename matcl::fft::get_element_cos<M-1, base_ind, step,
                                        negate, zero_v, one_v, mone_v, half_v, Scal>::type;
};

//TODO: reimplement
template<Integer M, class Scal, class Tag, Integer Row, Integer Col>
struct mk::make_expr_ufunc<tag_dct2<M,Scal>,mkd::element<Tag, Row, Col>>
{
    //sum_{n=0}^{N-1} x_n*cos( pi*k*(2*n+1)/(2*N) )
    static const Integer ind        = (Row-1) * (2*(Col-1)+1);
    static const Integer step       = 1;
    using info                      = matcl::fft::cos_info<ind, 2*M>;

    static const Integer base_ind   = info::base_index;
    static const bool negate        = info::negate;
    static const bool zero_v        = info::zero;
    static const bool one_v         = info::one;
    static const bool mone_v        = info::mone;
    static const bool half_v        = info::half;

    using type                      = typename matcl::fft::get_element_cos<2*M, base_ind, step,
                                        negate, zero_v, one_v, mone_v, half_v, Scal>::type;
};

//TODO: reimplement
template<Integer M, class Scal, class Tag, Integer Row, Integer Col>
struct mk::make_expr_ufunc<tag_dct3<M,Scal>,mkd::element<Tag, Row, Col>>
{
    //sum_{n=0}^{N-1} x_n*cos( pi*n*(2*k+1)/(2*N) )
    static const Integer ind        = (Col-1) * (2*(Row-1)+1);
    static const Integer step       = 1;
    using info                      = matcl::fft::cos_info<ind, 2*M>;

    static const Integer base_ind   = info::base_index;
    static const bool negate        = info::negate;
    static const bool zero_v        = info::zero;
    static const bool one_v         = info::one;
    static const bool mone_v        = info::mone;
    static const bool half_v        = info::half;

    using type                      = typename matcl::fft::get_element_cos<2*M, base_ind, step,
                                        negate, zero_v, one_v, mone_v, half_v, Scal>::type;
};

template<Integer M, class Scal, class Tag, Integer Row, Integer Col>
struct mk::make_expr_ufunc<tag_dct4<M,Scal>,mkd::element<Tag, Row, Col>>
{
    //sum_{n=0}^{N-1} x_n*cos( pi*(2*k+1)*(2*n+1)/(4*N) )
    static const Integer ind        = (2*(Row-1)+1) * (2*(Col-1)+1);
    static const Integer step       = 2;
    using info                      = matcl::fft::cos_info<ind, 4*M>;

    static const Integer base_ind   = info::base_index;
    static const bool negate        = info::negate;
    static const bool zero_v        = info::zero;
    static const bool one_v         = info::one;
    static const bool mone_v        = info::mone;
    static const bool half_v        = info::half;

    using type                      = typename matcl::fft::get_element_cos<4*M, base_ind, step,
                                        negate, zero_v, one_v, mone_v, half_v, Scal>::type;
};

template<Integer M, class Tag, Integer Row, Integer Col>
struct mk::make_expr_ufunc<tag_dct4_2<M>,mkd::element<Tag, Row, Col>>
{
    //cos ( pi * (2*n+1) / 4N )
    static const Integer ind        = 2*(Row-1) + 1;
    static const Integer step       = 2;
    using info                      = matcl::fft::cos_info<ind, 4*M>;

    static const Integer base_ind   = info::base_index;
    static const bool negate        = info::negate;
    static const bool zero_v        = info::zero;
    static const bool one_v         = info::one;
    static const bool mone_v        = info::mone;
    static const bool half_v        = info::half;

    using type                      = typename matcl::fft::get_element_cos<4*M, base_ind, step,
                                            negate,zero_v, one_v, mone_v, half_v, two>::type;
};
*/

template<Integer Rows, class Scal>
using cos_dct1  = mk::func_unary<tag_dct1<Rows,Scal>>;

template<Integer Rows, class Scal>
using cos_dct2  = mk::func_unary<tag_dct2<Rows,Scal>>;

template<Integer Rows, class Scal>
using cos_dct3  = mk::func_unary<tag_dct3<Rows,Scal>>;

template<Integer Rows, class Scal>
using cos_dct4  = mk::func_unary<tag_dct4<Rows,Scal>>;

template<Integer Rows>
using cos_dct4_dct2  = mk::func_unary<tag_dct4_2<Rows>>;

//TODO
/*
template<class Scal, class Mat>
auto cos_dct1(Mat)      -> typename mk::func_unary<tag_dct1<Mat::rows,Scal>,Mat>::type;

template<class Scal, class Mat>
auto cos_dct2(Mat)      -> typename mk::func_unary<tag_dct2<Mat::rows,Scal>,Mat>::type;

template<class Scal, class Mat>
auto cos_dct3(Mat)      -> typename mk::func_unary<tag_dct3<Mat::rows,Scal>,Mat>::type;

template<class Scal, class Mat>
auto cos_dct4(Mat)      -> typename mk::func_unary<tag_dct4<Mat::rows,Scal>,Mat>::type;

template<class Mat>
auto cos_dct4_dct2(Mat) -> typename mk::func_unary<tag_dct4_2<Mat::rows>,Mat>::type;
*/

template<Integer M>
struct tag_t            
{ 
    //TODO: remove
    static constexpr Integer get_offset(Integer Row, Integer Col)
    {
        (void)Row;
        (void)Col;
        return 0;
    };

    //TODO: remove
    using                       root_align_type = mk::align_full; 
    static const Integer        step            = 1; 

    static void print(std::ostream& os,int) { os << "T"; }

    template<Integer Row, Integer Col>
    static mkd::element<tag_t,Row,Col> get_elem();
};

template<class Output, class Input>
struct dct1_evaler
{
    static void eval(double* output, const double* input);
};
template<class Output, class Input>
struct dct2_evaler
{
    static void eval(double* output, const double* input);
};
template<class Output, class Input>
struct dct3_evaler
{
    static void eval(double* output, const double* input);
};
template<class Output, class Input>
struct dct4_evaler
{
    static void eval(double* output, const double* input);
};

template<Integer Rows0, Integer Cols0>
struct rows_cols
{
    static const Integer rows   = Rows0;
    static const Integer cols   = Cols0;
};

template<class Scal>
struct func_dct1
{
    static std::string  name()  { return "dct1"; };

    template<Integer Rows, Integer Cols>
    using return_size           = rows_cols<Rows,Cols>;

    template<class Output, class Input>
    using evaler                = dct1_evaler<Output,Input>;
};

template<class Scal>
struct func_dct2
{
    static std::string  name()  { return "dct2"; };

    template<Integer Rows, Integer Cols>
    using return_size           = rows_cols<Rows,Cols>;

    template<class Output, class Input>
    using evaler                = dct2_evaler<Output,Input>;
};

template<class Scal>
struct func_dct3
{
    static std::string  name()  { return "dct3"; };

    template<Integer Rows, Integer Cols>
    using return_size           = rows_cols<Rows,Cols>;

    template<class Output, class Input>
    using evaler                = dct3_evaler<Output,Input>;
};

template<class Scal>
struct func_dct4
{
    static std::string  name()  { return "dct4"; };

    template<Integer Rows, Integer Cols>
    using return_size           = rows_cols<Rows,Cols>;

    template<class Output, class Input>
    using evaler                = dct4_evaler<Output,Input>;
};

//-----------------------------------------------------------------------
//                          dct2 static expression
//-----------------------------------------------------------------------
template<Integer M, class Config, bool Gen_Mat, class Scal>
struct make_cos_mat_dct2
{
    using T                     = mk::const_mat<M, M, tag_t<M>>;
    using type                  = decltype(cos_dct2<M,Scal>::eval(T()));
};

template<Integer M, class Config, class Scal>
struct make_cos_mat_dct2<M,Config,true, Scal>
{
    using type                  = mk::gen_mat<M,M,tag_cos_dct2<M,Scal>>;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, bool Recursive>
struct eval_dct2
{
    private:
        static const Integer M  = Matrix_Type::rows;
        static const bool gm    = (M <= Config::max_genmat_dct2);
        using cos_mat           = typename make_cos_mat_dct2<M, Config, gm, Scal>::type;
        using result            = decltype( cos_mat() * Matrix_Type() );

    public:
        using expression        = result;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal>
struct eval_dct2<Config, Matrix_Type, Tag_Seed, Scal, true>
{
    private:
        static const Integer M      = Matrix_Type::rows;
        static const bool force_1   = false;
        static const bool force_2   = false;

        using seed          = typename add_level<Tag_Seed, 2>::type;
        using tag_virt      = tag_virt<M>;
        using tmp1          = temp<seed, M, 1>;
        using tmp2          = temp<seed, M, 2>;
        using tmp3          = temp<seed, M, 3>;
        using tmp4          = temp<seed, M, 4>;

        using ys_1_expr     = decltype(Matrix_Type::sub(mk::colon2<1,M/2>()) 
                                       + Matrix_Type::sub(mk::colon3<M,-1,M/2+1>()) );
        using ys_2_expr     = decltype(Matrix_Type::sub(mk::colon2<1,M/2>()) 
                                       - Matrix_Type::sub(mk::colon3<M,-1,M/2+1>()) );

        using ys_1          = decltype( ys_1_expr().make_temp<tmp1>() );
        using ys_2          = decltype( ys_2_expr().make_temp<tmp2>() );

      #if CALL_FUNC_dct2 == 0
        using ev_dct2       = eval_dct2<Config, ys_1, seed, Scal>;
        using res_10        = typename ev_dct2::expression;
        using res_1         = decltype(res_10().make_temp<tmp3, force_1>() );
      #elif CALL_FUNC_dct2 == 1
        template<class Arg>
        using func_1        = eval_dct2<Config, Arg, seed, Scal>;
        using res_1         = decltype( call_inline<tmp3, func_1>(ys_1()) );
      #else
        using func_1        = func_dct2<Scal>;
        using res_1         = decltype( mk::call_external<tmp3, func_1>(ys_1()) );
      #endif
        
      #if CALL_FUNC_dct2 == 0
        using ev_dct4       = eval_dct4<Config, ys_2, seed, Scal>;
        using res_20        = typename ev_dct4::expression;
        using res_2         = decltype(res_20().make_temp<tmp4, force_2>() );
      #elif CALL_FUNC_dct2 == 1
        template<class Arg>
        using func_2        = eval_dct4<Config, Arg, seed, Scal>;
        using res_2         = decltype( call_inline<tmp4, func_2>(ys_2()) );
      #else
        using func_2        = func_dct4<Scal>;
        using res_2         = decltype( mk::call_external<tmp4, func_2>(ys_2()) );
      #endif

        using result        = mk::virtual_mat<M,1,tag_virt>;

        using part_res1     = decltype(   result::assign_1(mk::colon3<1,2,M-1>(), res_1() ) );
        using final_result  = decltype(part_res1::assign_1(mk::colon3<2,2,M>(),   res_2() ) );

    public:
        using expression    = final_result;
};

//-----------------------------------------------------------------------
//                          dct3 static expression
//-----------------------------------------------------------------------
template<Integer M, class Config, bool Gen_Mat, class Scal>
struct make_cos_mat_dct3
{
    using T                     = mk::const_mat<M,M,tag_t<M>>;
    using type                  = decltype(cos_dct3<M, Scal>::eval(T::sub(mk::colon_all(), mk::colon2<2,M>())) );
};
template<Integer M, class Config, class Scal>
struct make_cos_mat_dct3<M,Config,true, Scal>
{
    using type                  = mk::gen_mat<M,M-1,tag_cos_dct3<M,Scal>>;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, bool Recursive>
struct eval_dct3
{
    private:
        //C_k     = 0.5*x_0 + sum_{n=1}^{N-1} x_n*cos( pi*n*(k+0.5)/N)
        static const Integer M  = Matrix_Type::rows;
        static const bool gm    = (M <= Config::max_genmat_dct3);
        using cos_mat           = typename make_cos_mat_dct3<M, Config, gm, Scal>::type;

        using seed              = typename add_level<Tag_Seed, 3>::type;
        using tmp1              = temp<seed,M,1>;

        using first             = decltype(Matrix_Type::elem(mk::colon<1>()));
        using constant          = decltype((Scal() * mk::half() *first()).compute<tmp1>() );

        using result            = decltype( constant() + cos_mat() * Matrix_Type::sub( mk::colon2<2,M>() ) );

    public:
        using expression        = result;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal>
struct eval_dct3<Config, Matrix_Type, Tag_Seed, Scal, true>
{
    private:
        static const Integer M      = Matrix_Type::rows;

        using seed          = typename add_level<Tag_Seed, 2>::type;
        using tag_virt      = tag_virt<M>;
        using tmp1          = temp<seed, M, 1>;
        using tmp2          = temp<seed, M, 2>;
        using tmp3          = temp<seed, M, 3>;
        using tmp4          = temp<seed, M, 4>;

        using ys_1          = decltype(Matrix_Type::sub(mk::colon3<1,2,M>()) );
        using ys_2          = decltype(Matrix_Type::sub(mk::colon3<2,2,M>()) );

      #if CALL_FUNC_dct3 == 0
        using ev_dct3       = eval_dct3<Config, ys_1, seed, Scal>;
        using res_10        = typename ev_dct3::expression;
        using res_1         = decltype(res_10().make_temp<tmp1>() );
      #elif CALL_FUNC_dct3 == 1
        template<class Arg>
        using func_1        = eval_dct3<Config, Arg, seed, Scal>;
        using res_1         = decltype( call_inline<tmp3, func_1>(ys_1()) );
      #else
        using func_1        = func_dct3<Scal>;
        using res_1         = decltype( mk::call_external<tmp3, func_1>(ys_1()) );
      #endif
        
      #if CALL_FUNC_dct3 == 0
        using ev_dct4       = eval_dct4<Config, ys_2, seed, Scal>;
        using res_20        = typename ev_dct4::expression;
        using res_2         = decltype(res_20().make_temp<tmp2, true>() );
      #elif CALL_FUNC_dct3 == 1
        template<class Arg>
        using func_2        = eval_dct4<Config, Arg, seed, Scal>;
        using res_2         = decltype( call_inline<tmp4, func_2>(ys_2()) );
      #else
        using func_2        = func_dct4<Scal>;
        using res_2         = decltype( mk::call_external<tmp4, func_2>(ys_2()) );
      #endif

        using result        = mk::virtual_mat<M,1,tag_virt>;

        using part_res1     = decltype(   result::assign_1(mk::colon3<1,1,M/2>(), res_1() + res_2()) );
        using final_result  = decltype(part_res1::assign_1(mk::colon3<M,-1,M/2+1>(),res_1() - res_2()) );

    public:
        using expression    = final_result;
};

//-----------------------------------------------------------------------
//                          dct4 static expression
//-----------------------------------------------------------------------
template<Integer M, class Config, bool Gen_Mat, class Scal>
struct make_cos_mat_dct4
{
    using T                     = mk::const_mat<M,M,tag_t<M>>;
    using type                  = decltype(cos_dct4<M, Scal>::eval(T()));
};
template<Integer M, class Config, class Scal>
struct make_cos_mat_dct4<M,Config,true, Scal>
{
    using type                  = mk::gen_mat<M,M,tag_cos_dct4<M,Scal>>;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, bool Recursive>
struct eval_dct4
{
    private:
        static const Integer M  = Matrix_Type::rows;
        static const bool gm    = (M <= Config::max_genmat_dct4);
        using cos_mat           = typename make_cos_mat_dct4<M, Config, gm, Scal>::type;

        using result            = decltype( cos_mat() * Matrix_Type() );

    public:
        using expression        = result;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal>
struct eval_dct4<Config, Matrix_Type, Tag_Seed, Scal, true>
{
    private:
        static const Integer M  = Matrix_Type::rows;
        using T                 = mk::const_mat<M,1,tag_t<M>>;
        using cos               = decltype( cos_dct4_dct2(T()) ) ;

        using seed              = typename add_level<Tag_Seed, 4>::type;
        using tmp1              = temp<seed,M,1>;
        using tmp2              = temp<seed,M,2>;
        using tmp3              = temp<seed,M,3>;

        using y_cos             = decltype( mult_rows( Matrix_Type(), cos()).make_temp<tmp3>() );

      #if CALL_FUNC_dct4 == 0
        using ev_dct2           = eval_dct2<Config, y_cos,seed, Scal>;
        using res               = typename ev_dct2::expression;
        using res_temp          = decltype( res().make_temp<tmp1>() );
      #elif CALL_FUNC_dct4 == 1
        template<class Arg>
        using func              = eval_dct2<Config, Arg,seed, Scal>;
        using res_temp          = decltype( call_inline<tmp1, func>(y_cos()) );
      #else
        using func              = func_dct2<Scal>;
        using res_temp          = decltype( mk::call_external<tmp1, func>(y_cos()) );
      #endif

        using comp              = mk::computation<tmp2,res_temp> ;
        using comp_1            = decltype( comp().assign(mk::colon<1>(), 
                                            mk::operator/(comp().elem(mk::colon<1>()) , mk::two()) ));
        struct I{};

        template<class Subject, class Context>
        struct expr
        {
            static const Integer I = mk::get_from_context<Context,I>::value;
            using type          = decltype(Subject().assign(mk::colon<I>(),  Subject().elem(mk::colon<I>()) 
                                        - Subject().elem(mk::colon<I-1>()) ) ); 
        };

        using result_comp       = typename mk::for_expr<I,2,M,expr,comp_1>::type;
        using result            = typename mk::get_computation_result<result_comp>::type;

    public:
        using expression        = result;
};

//-----------------------------------------------------------------------
//                          dct1 static expression
//-----------------------------------------------------------------------
template<Integer M, class Config, bool Gen_Mat, class Scal>
struct make_cos_mat_dct1
{
    using T                     = mk::const_mat<M,M,tag_t<M>>;
    using type                  = decltype(cos_dct1<M, Scal>::eval(T::sub(mk::colon_all(), mk::colon2<2,M-1>())) );
};

template<Integer M, class Config, class Scal>
struct make_cos_mat_dct1<M,Config,true,Scal>
{
    using type                  = mk::gen_mat<M,M-2,tag_cos_dct1<M,Scal>>;
};

template<class Config, class Matrix_Type, class Tag_Seed, class Scal, bool Recursive>
struct eval_dct1
{
    private:
        //C_k     = 0.5*[x_0 + (-1)^k*x_{N-1}] + sum_{n=1}^{N-2} x_n*cos(pi*n*k/(N-1))
        static const Integer M      = Matrix_Type::rows;
        using seed                  = typename add_level<Tag_Seed, 1>::type;
        using tmp1                  = temp<seed,M,1>;
        using tmp2                  = temp<seed,M,2>;
        using tmp3                  = temp<seed,M,3>;
        using tmp4                  = temp<seed,M,4>;

        static const bool gm        = (M <= Config::max_genmat_dct1 || M == 4);
        using cos_mat               = typename make_cos_mat_dct1<M, Config, gm, Scal>::type;
        using last                  = decltype(Matrix_Type::elem(mk::colon<M>()));
        using first                 = decltype(Matrix_Type::elem(mk::colon<1>()));
        using constant              = decltype((Scal() * mk::half() * ( first() + last() ) ).compute<tmp3>() );
        using constant_last         = decltype((Scal() * last()).compute<tmp4>() );
        
        using res                   = decltype( constant() + cos_mat() * Matrix_Type::sub( mk::colon2<2,M-1>() ) );
        using res_temp              = decltype( res().make_temp<tmp2>() );

        using comp                  = mk::computation<tmp1,res_temp> ;
        using comp_1                = decltype( comp().assign(mk::colon3<2,2,M>(), 
                                            res_temp::sub(mk::colon3<2,2,M>()) - constant_last() ));

        using result                = typename mk::get_computation_result<comp_1>::type;

    public:
        using expression            = result;
};
template<class Config, class Matrix_Type, class Tag_Seed, class Scal>
struct eval_dct1<Config, Matrix_Type,Tag_Seed,Scal,true>
{
    private:
        static const Integer M      = Matrix_Type::rows;
        static const Integer M2     = (M-1)/2;

        static const bool force_2   = false;

        using seed                  = typename add_level<Tag_Seed, 11>::type;
        using tag_virt              = tag_virt<M>;
        using tmp1                  = temp<seed, M, 1>;
        using tmp2                  = temp<seed, M, 2>;
        using tmp3                  = temp<seed, M, 3>;
        using tmp4                  = temp<seed, M, 4>;

        using x_odd                 = decltype(Matrix_Type::sub(mk::colon3<1,2,M>()) );
        using x_even                = decltype(Matrix_Type::sub(mk::colon3<2,2,M>()) );

      #if CALL_FUNC_dct1 == 0
        using ev_dct1               = eval_dct1<Config, x_odd,seed,Scal>;
        using res_10                = typename ev_dct1::expression;
        using res_1                 = decltype( res_10().make_temp<tmp1>() );
      #elif CALL_FUNC_dct1 == 1
        template<class Arg>
        using func_1                = eval_dct1<Config, Arg, seed,Scal>;
        using res_1                 = decltype( call_inline<tmp1, func_1>(x_odd()) );
      #else
        using func_1                = func_dct1<Scal>;
        using res_1                 = decltype( mk::call_external<tmp1, func_1>(x_odd().make_temp<tmp3>() ) );
      #endif

      #if CALL_FUNC_dct1 == 0
        using ev_dct2               = eval_dct2<Config, x_even,seed,Scal>;
        using res_20                = typename ev_dct2::expression;
        using res_2                 = decltype( res_20().make_temp<tmp2,force_2>() );
      #elif CALL_FUNC_dct1 == 1
        template<class Arg>
        using func_2                = eval_dct2<Config, Arg, seed,Scal>;
        using res_2                 = decltype( call_inline<tmp2, func_2>(x_even()) );
      #else
        using func_2                = func_dct2<Scal>;
        using res_2                 = decltype( mk::call_external<tmp2, func_2>(x_even().make_temp<tmp4>()) );
      #endif

        using sum_1                 = decltype( res_1::sub(mk::colon2<1,M2>()) + res_2() );
        using sum_2                 = decltype( res_1::sub(mk::colon2<1,M2>()) - res_2() );

        using result                = mk::virtual_mat<M,1,tag_virt>;

        using part_res1             = decltype(   result::assign_1(mk::colon3<1,1,M2>(), sum_1() ) );
        using part_res2             = decltype(part_res1::assign_1(mk::colon<M2+1>(), res_1::sub(mk::colon<M2+1>()) ) );
        using part_res3             = decltype(part_res2::assign_1(mk::colon3<M,-1,M2+2>(),sum_2() ) );

        using final_result          = part_res3;

    public:
        using expression            = final_result;
};
template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
struct dct1_impl
{
    using seed                      = list_int<Config::seed,1>;
    struct tag{};

    struct data_provider
    {
        const double* restricted    input;
        double* restricted          output;
    };

    struct tag_y
    { 
        static const bool           is_continuous   = true; 
        using                       root_align_type = mk::align_full; 
        static const Integer        step            = 1;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return Row - 1;
        };

        static void                 print(std::ostream& os,int)             { os << "dct1_impl.y"; }

        template<class Val, Integer Offset>
        static const Val* restricted
                                    get_data_ptr(data_provider& dp)         { return dp.input + Offset; };
    };

    struct tag_ret  
    { 
        static const bool           is_continuous   = (Out_Step == 1); 
        using                       root_align_type = mk::align_none; 
        static const Integer        step            = Out_Step;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return (Row - 1) * Out_Step;
        };

        static void                 print(std::ostream& os,int)             { os << "dct1_impl.ret"; }

        template<class Val, Integer Offset>
        static Val* restricted      get_data_ptr(data_provider& dp)         { return dp.output + Offset * Out_Step; };
    };

    using y                         = mk::gen_mat<M, 1, tag_y>;
    using ret                       = mk::output_mat<M, 1, tag_ret>;

    using ev_dct1                   = eval_dct1<Config, y, seed, Scal>;
    using expression_type           = typename ev_dct1::expression;
    using code_gen                  = Code_Gen;
    using evaler_type               = mk::expr_evaler<tag,code_gen, ret, expression_type, double>;

    static void eval(double* output, const double* input);
    static void print();
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
struct dct2_impl
{
    using seed                      = list_int<Config::seed,2>;

    struct data_provider
    {
        const double* restricted    input;
        double* restricted          output;
    };

    struct tag_y
    { 
        static const bool           is_continuous   = true; 
        using                       root_align_type = mk::align_full; 
        static const Integer        step            = 1;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return Row - 1;
        };

        static void                 print(std::ostream& os,int)             { os << "dct2_impl.y"; }

        template<class Val, Integer Offset>
        static const Val* restricted
                                    get_data_ptr(data_provider& dp)         { return dp.input + Offset; };
    };

    struct tag_ret  
    { 
        static const bool           is_continuous   = (Out_Step == 1); 
        using                       root_align_type = mk::align_none; 
        static const Integer        step            = Out_Step;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return (Row - 1) * Out_Step;
        };

        static void                 print(std::ostream& os,int)             { os << "dct2_impl.ret"; }

        template<class Val, Integer Offset>
        static Val* restricted      get_data_ptr(data_provider& dp)         { return dp.output + Offset * Out_Step; };
    };

    using y                         = mk::gen_mat<M, 1, tag_y>;
    using ret                       = mk::output_mat<M, 1, tag_ret>;

    struct tag{};

    using ev_dct2                   = eval_dct2<Config, y, seed, Scal>;
    using expression_type           = typename ev_dct2::expression;
    using evaler_type               = mk::expr_evaler<tag, Code_Gen, ret, expression_type, double>;

    static void eval(double* output, const double* input);
    static void print();
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
struct dct3_impl
{
    using seed                      = list_int<Config::seed,3>;

    struct data_provider
    {
        const double* restricted    input;
        double* restricted          output;
    };

    struct tag_y
    { 
        static const bool           is_continuous   = true; 
        using                       root_align_type = mk::align_full; 
        static const Integer        step            = 1;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return Row - 1;
        };

        static void                 print(std::ostream& os,int)             { os << "dct3_impl.y"; }

        template<class Val, Integer Offset>
        static const Val* restricted
                                    get_data_ptr(data_provider& dp)         { return dp.input + Offset; };
    };

    struct tag_ret  
    { 
        static const bool           is_continuous   = (Out_Step == 1); 
        using                       root_align_type = mk::align_none; 
        static const Integer        step            = Out_Step;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return (Row - 1) * Out_Step;
        };

        static void                 print(std::ostream& os,int)             { os << "dct3_impl.ret"; }

        template<class Val, Integer Offset>
        static Val* restricted      get_data_ptr(data_provider& dp)         { return dp.output + Offset * Out_Step; };
    };

    using y                         = mk::gen_mat<M, 1, tag_y>;
    using ret                       = mk::output_mat<M, 1, tag_ret>;

    struct tag{};

    using ev_dct3                   = eval_dct3<Config, y, seed, Scal>;
    using expression_type           = typename ev_dct3::expression;
    using evaler_type               = mk::expr_evaler<tag, Code_Gen, ret, expression_type, double>;

    static void eval(double* output, const double* input);
    static void print();
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
struct dct4_impl
{
    using seed                      = list_int<Config::seed,4>;

    struct data_provider
    {
        const double* restricted    input;
        double* restricted          output;
    };

    struct tag_y
    { 
        static const bool           is_continuous   = true; 
        using                       root_align_type = mk::align_full; 
        static const Integer        step            = 1;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return Row - 1;
        };

        template<class Val, Integer Offset>
        static const Val* restricted
                                    get_data_ptr(data_provider& dp)         { return dp.input + Offset; };
        static void                 print(std::ostream& os,int)             { os << "dct4_impl.y"; }
    };

    struct tag_ret  
    { 
        static const bool           is_continuous   = (Out_Step == 1); 
        using                       root_align_type = mk::align_none; 
        static const Integer        step            = Out_Step;

        static constexpr Integer get_offset(Integer Row, Integer Col)
        {
            (void)Col;
            return (Row - 1) * Out_Step;
        };

        template<class Val, Integer Offset>
        static Val* restricted      get_data_ptr(data_provider& dp)         { return dp.output + Offset * Out_Step; };

        static void                 print(std::ostream& os,int)             { os << "dct4_impl.ret"; }
    };

    using y                         = mk::gen_mat<M, 1, tag_y>;
    using ret                       = mk::output_mat<M, 1, tag_ret>;

    struct tag{};

    using ev_dct4                   = eval_dct4<Config, y, seed, Scal>;
    using expression_type           = typename ev_dct4::expression;
    using evaler_type               = mk::expr_evaler<tag, Code_Gen, ret, expression_type, double>;

    static void eval(double* output, const double* input);
    static void print();
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
inline_lev_root
void dct1_impl<Config, M,Out_Step, Code_Gen, Scal>::eval(double* output, const double* input)
{
    data_provider dp{input,output};
    evaler_type ev;
    ev.eval(dp);
};
template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
void dct1_impl<Config, M,Out_Step, Code_Gen, Scal>::print()
{
    evaler_type::print(std::cout,0);
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
inline_lev_root
void dct2_impl<Config, M,Out_Step,Code_Gen,Scal>::eval(double* output, const double* input)
{
    data_provider dp{input,output};
    evaler_type ev;
    ev.eval(dp);
};
template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
void dct2_impl<Config, M,Out_Step, Code_Gen, Scal>::print()
{
    evaler_type::print(std::cout,0);
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
inline_lev_root
void dct3_impl<Config, M,Out_Step,Code_Gen,Scal>::eval(double* output, const double* input)
{
    data_provider dp{input,output};
    evaler_type ev;
    ev.eval(dp);
};
template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
void dct3_impl<Config, M,Out_Step, Code_Gen, Scal>::print()
{
    evaler_type::print(std::cout,0);
};

template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
inline_lev_root
void dct4_impl<Config, M,Out_Step,Code_Gen, Scal>::eval(double* output, const double* input)
{
    data_provider dp{input,output};
    evaler_type ev;
    ev.eval(dp);
};
template<class Config, Integer M, Integer Out_Step, class Code_Gen, class Scal>
void dct4_impl<Config, M, Out_Step, Code_Gen,Scal>::print()
{
    evaler_type::print(std::cout,0);
};

template<class Output, class Input>
void dct1_evaler<Output,Input>::eval(double* output, const double* input)
{
    static_assert(Input::rows == Output::rows, "invalid arguments");
    static_assert(Input::cols == Output::cols, "invalid arguments");
    static_assert(Input::cols == 1, "not implemented");

    static const Integer M          = Input::rows;
    static const Integer out_start  = Output::start;
    static const Integer out_step   = Output::step;

    return dct1_impl<M,out_step>::eval(output, input);
};

template<class Output, class Input>
void dct2_evaler<Output,Input>::eval(double* output, const double* input)
{
    static_assert(Input::rows == Output::rows, "invalid arguments");
    static_assert(Input::cols == Output::cols, "invalid arguments");
    static_assert(Input::cols == 1, "not implemented");

    static const Integer M          = Input::rows;
    static const Integer out_start  = Output::start;
    static const Integer out_step   = Output::step;

    return dct2_impl<M,out_step>::eval(output, input);
};

template<class Output, class Input>
void dct3_evaler<Output,Input>::eval(double* output, const double* input)
{
    static_assert(Input::rows == Output::rows, "invalid arguments");
    static_assert(Input::cols == Output::cols, "invalid arguments");
    static_assert(Input::cols == 1, "not implemented");

    static const Integer M          = Input::rows;
    static const Integer out_start  = Output::start;
    static const Integer out_step   = Output::step;

    return dct3_impl<M,out_step>::eval(output, input);
};

template<class Output, class Input>
void dct4_evaler<Output,Input>::eval(double* output, const double* input)
{
    static_assert(Input::rows == Output::rows, "invalid arguments");
    static_assert(Input::cols == Output::cols, "invalid arguments");
    static_assert(Input::cols == 1, "not implemented");

    static const Integer M          = Input::rows;
    static const Integer out_start  = Output::start;
    static const Integer out_step   = Output::step;

    return dct4_impl<M,out_step>::eval(output, input);
};

}};