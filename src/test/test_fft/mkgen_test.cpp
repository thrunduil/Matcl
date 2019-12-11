#include "mkgen/mkgen.h"

#include <iostream>

namespace matcl { namespace mkgen { namespace test
{

struct tag_pi
{
    static constexpr double value = 3.14;
};

struct tag_2
{
    static constexpr double value = 2.0;
};

template<class Tag>
struct tag
{
    static constexpr double value = Tag::value;
};

template<class T, class S>
struct sum_tag
{
    static constexpr double value = T::value + S::value;
};

template<class Tag1, class Tag2>
auto operator+(const tag<Tag1>& a, const tag<Tag2>& b) -> tag<sum_tag<Tag1, Tag2>>;

void test()
{
    using type = decltype(tag<tag_pi>() + tag<tag_2>());
    std::cout << type::value << "\n";
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
