#include "mkgen/mkgen.h"
#include "mkgen/expression/expressions.h"
#include "mkgen/details/matrix/matrix_arrays.h"

#include <iostream>

namespace matcl { namespace mkgen { namespace test
{

namespace mk = matcl::mkgen;

struct basic_codegen
{
    static const bool simd_enable                   = true;
    static const bool simd_allow_negative_step      = true;
    static const bool simd_allow_unaligned          = true;

    static const bool simd_half_allow               = true;
    static const bool simd_half_allow_negative_step = true;
    static const bool simd_half_allow_unaligned     = true;

    using simd_type                                 = simd::maximum_tag;
    static const int simd_bits                      = simd::maximum_bits;
};

struct empty_context
{
    using deps_temp     = mk::dps<>;
    using code_gen      = basic_codegen;
};

struct tag_ret
{ 
    template<Integer Row, Integer Col>
    using get_offset        = std::integral_constant<Integer, Row - 1>;

    /*
    static const bool           is_continuous   = true; 

    static void     print(std::ostream& os) { os << "B"; }
    
    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        return  dp.ptr_C + Offset;
    }
    */
};

void test_scalar()
{
    // Define expression:    
    using expr = decltype(integer_scalar<1>() + integer_scalar<2>());
    expr::print<empty_context>(std::cout, 0);

    /*
    using ret               = output_mat<1, 1, array<tag_ret>>;
    struct etag{};

    using evaler_type       = expr_evaler<etag, basic_codegen, ret, expr, double>;

    evaler_type::print(std::cout, 0);
    //evaler_type().eval(dp);
    */
};

struct tag_TM_1
{
    static constexpr Integer get_offset(Integer Row, Integer Col) { return Row + Col; };
    using root_align_type   = mk::align_full; 
    static constexpr Integer step   = 1;
    static const bool is_continuous = true;

    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "TM1";
    };
};

void test_matrix()
{
    struct tag_TM_2{};
    struct tag_TM_3{};

    using colon1    = mk::colon<2>;
    using colon2    = mk::colon_all;
    using colon3    = mk::colon2<2,4>;
    using colon4    = mk::colon3<2,2,4>;

    using mat   = mk::gen_mat<5, 5, tag_TM_1>;
    using e1    = decltype(mat::elem(colon<1>()));
    using e2    = decltype(mat::elem(colon<1>(), colon<2>()));

    using e3    = decltype(mat::sub(colon1()));
    using e4    = decltype(mat::sub(colon2()));
    using e5    = decltype(mat::sub(colon3()));
    using e6    = decltype(mat::sub(colon4()));

    using e11   = decltype(mat::sub(colon1(), colon1()));
    using e21   = decltype(mat::sub(colon2(), colon1()));
    using e31   = decltype(mat::sub(colon3(), colon1()));
    using e41   = decltype(mat::sub(colon4(), colon1()));

    using e12   = decltype(mat::sub(colon1(), colon2()));
    using e22   = decltype(mat::sub(colon2(), colon2()));
    using e32   = decltype(mat::sub(colon3(), colon2()));
    using e42   = decltype(mat::sub(colon4(), colon2()));

    using e13   = decltype(mat::sub(colon1(), colon3()));
    using e23   = decltype(mat::sub(colon2(), colon3()));
    using e33   = decltype(mat::sub(colon3(), colon3()));
    using e43   = decltype(mat::sub(colon4(), colon3()));

    using e14   = decltype(mat::sub(colon1(), colon4()));
    using e24   = decltype(mat::sub(colon2(), colon4()));
    using e34   = decltype(mat::sub(colon3(), colon4()));
    using e44   = decltype(mat::sub(colon4(), colon4()));

    using t1    = decltype(mat::make_temp<tag_TM_2, false>());
    using t2    = decltype(mat::make_temp<tag_TM_3, true>());

    e44::print<empty_context>(std::cout, 0);
};

struct tag_pi : mk::scal_data_const_value_tag<tag_pi>
              , mk::scal_data_value_tag<tag_pi>
              , mk::scal_data_gen_value_tag<tag_pi>
{
    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "pi";
    };

    template<class Val>
    static constexpr Val value()
    {
        return Val(3.1415);
    }

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage&)
    { 
        return value<Val>(); 
    };

};

struct tag_e : mk::scal_data_const_value_tag<tag_e>
{
    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "e";
    };

    template<class Val>
    static constexpr Val value()
    {
        return Val(2.7187);
    }
};

struct tag_one : mk::scal_data_const_value_tag<tag_one>
{
    static void print(std::ostream& os, int prior)
    {
        (void)prior;
        os << "1";
    };

    template<class Val>
    static constexpr Val value()
    {
        return Val(1.0);
    }
};

struct local_storage_dummy
{};

void test_matrix2()
{
    local_storage_dummy ls;

    {
        using val_pi    = mk::gen_scalar<tag_pi>;

        using res       = decltype(two() *( half()*(one() + val_pi()) - two()* val_pi()));
        res::print<empty_context>(std::cout, 0);
        std::cout << " " << res::eval<double>(ls) << "\n";
    };
    {
        using val_pi    = mk::gen_scalar<tag_pi>;

        using res       = decltype(two() *( half()*val_pi() - two()* val_pi()));
        res::print<empty_context>(std::cout, 0);
        std::cout << " " << res::eval<double>(ls) << "\n";
    };

    {
        using val_pi    = mk::const_value_scalar<tag_pi, double>;
        using val_1     = mk::const_value_scalar<tag_one, double>;
        using val_e     = mk::const_value_scalar<tag_e, double>;

        using res       = decltype((val_pi() * val_1() - val_e())/val_e());
        res::print<empty_context>(std::cout, 0);

        std::cout << res::data_type::value<double>() << "\n";
        static_assert(res::data_type::value<double>() != 1.0);
    };
};

//TODO
#if 0
template<class Tag>
struct elem;

struct context
{
    using deps_temp = empty_deps;
};

/*
void test_gemm()
{
    // form GEMM product: C = A*B + a*C
    //
    // matrix dimensions:
    
    static const Integer M = 4;
    static const Integer K = 4;
    static const Integer N = 4;

    // tags for matrices and scalars:
    struct tag_A{ static void print(std::ostream& os) { os << "A"; }};
    struct tag_B{ static void print(std::ostream& os) { os << "B"; }};
    struct tag_C{ static void print(std::ostream& os) { os << "C"; }};
    struct tag_a{ static void print(std::ostream& os) { os << "a"; }};

    // Define expression GEMM_expr:    
    using GEMM_expr = 
        mat_assign
        <
            ct_matrix<M, N, array<tag_C>, empty_deps>,
     
            mat_plus
            <
                mat_mult< ct_matrix<M,K, array<tag_A>, empty_deps>, 
                            ct_matrix<K,N, array<tag_B>, empty_deps> > :: type,
                mat_mult< ct_scalar<elem<tag_a>, empty_deps>, 
                            ct_matrix<M,N, array<tag_C>, empty_deps> > :: type
            > :: type
    
        >::type;
         
    GEMM_expr::print<context>(std::cout, 0);
};
*/

//data provider:
struct dp_GEMM
{
    double*         ptr_A;
    double*         ptr_B;
    double*         ptr_C;
    double          scal_a;
     
    template<class Tag, class Val, Integer row, Integer col, Integer pos>
    Val&    get_elem_write()    { return Tag::get_write<pos>(this); }
     
    template<class Tag, class Val, Integer row, Integer col, Integer pos>
    Val     get_elem_read()     { return Tag::get_read<pos>(this); }
     
    template<class Tag, class Val>
    Val     get_scalar()        { return Tag::get(this); }
};

    
struct tag_A
{ 
    static const bool           is_continuous   = true; 

    static void     print(std::ostream& os) { os << "A"; }        
        
    template<Integer pos>
    static double&  get_write(dp_GEMM* dp)  { return dp->ptr_A[pos]; }
        
    template<Integer pos>
    static double   get_read(dp_GEMM* dp)   { return dp->ptr_A[pos]; }
};

struct tag_B
{ 
    static const bool           is_continuous   = true; 

    static void     print(std::ostream& os) { os << "B"; }
    
    template<Integer pos>
    static double&  get_write(dp_GEMM* dp)  { return dp->ptr_B[pos]; }
        
    template<Integer pos>
    static double   get_read(dp_GEMM* dp)   { return dp->ptr_B[pos]; }
};

struct tag_C
{ 
    static const bool           is_continuous   = true; 

    static void     print(std::ostream& os) { os << "B"; }
    
    template<class T, Integer Offset, class DP>
    static const double* restricted get_data_ptr(const DP& dp)
    { 
        return  dp.ptr_C + Offset;
    }
};
     
struct tag_a
{ 
    static void     print(std::ostream& os) { os << "a"; }
    static double   get(dp_GEMM* dp)        { return dp->scal_a; };   
};

struct etag{};

void test_gemm_eval()
{
    // form GEMM product: C = A*B + a*C
    //
    // matrix dimensions:
    
    static const Integer M = 4;
    static const Integer K = 4;
    static const Integer N = 4;

    // Define expression GEMM_expr:    
    // output matrix
    using ret   = output_mat<M, N, array<tag_C>>;

    using GEMM_expr = 
        mat_assign
        <
            ret,
     
            mat_plus
            <
                mat_mult< ct_matrix<M,K, array<tag_A>, empty_deps>, 
                            ct_matrix<K,N, array<tag_B>, empty_deps> > :: type,
                mat_mult< ct_scalar<elem<tag_a>, empty_deps>, 
                            ct_matrix<M,N, array<tag_C>, empty_deps> > :: type
            > :: type
    
        >::type;

    //evaluation:
    //dp_GEMM dp;

    using evaler_type       = expr_evaler<etag, basic_codegen, ret, GEMM_expr, double>;

    evaler_type::print(std::cout, 0);
    //evaler_type().eval(dp);
};
#endif

}}};
