#include "matcl-fft/matcl_fft_config.h"
#include "matcl-fft/matcl_dct.h"

#include "mkgen/mkgen.h"
#include "static_dct.h"
#include "dct_kernels_test.h"

#define ENABLE_IACA 0

#if ENABLE_IACA == 1
    #include <iacaMarks.h>
#endif

#include <iostream>

namespace matcl { namespace fft
{

void MATCL_FFT_EXPORT test_dct_II();

struct dp_cheb;

struct default_dct_config
{
    static const Integer seed                       = -1;

    static const Integer min_rec_length_dct1        = 4;
    static const Integer min_rec_length_dct2        = 4;
    static const Integer min_rec_length_dct3        = 4;
    static const Integer min_rec_length_dct4        = 4;

    static const Integer max_genmat_dct1            = 0;
    static const Integer max_genmat_dct2            = 4;    
    static const Integer max_genmat_dct3            = 4;    
    static const Integer max_genmat_dct4            = 4;
};



struct tag_y
{ 
    static const bool           is_continuous   = true;
    using                       root_align_type = mk::align_full;
    static const Integer        step            = 1;

    template<Integer Row, Integer Col>
    using                       get_offset      = std::integral_constant<Integer,Row - 1>;

    static void                 print(std::ostream& os,int)             { os << "y"; }

    template<class T, Integer Offset>
    static const double* restricted
                                get_data_ptr(const dp_cheb& dp)         { return dp.ptr_y + Offset; }
};
struct tag_ret
{ 
    static const bool           is_continuous   = true;
    using                       root_align_type = mk::align_full;
    static const Integer        step            = 1;

    template<Integer Row, Integer Col>
    using                       get_offset      = std::integral_constant<Integer,Row - 1>;

    static void                 print(std::ostream& os,int) { os << "ret"; }
        
    template<class T, Integer Offset>
    static double* restricted   get_data_ptr(dp_cheb& dp)   { return dp.ptr_ret + Offset; }
};

struct dp_cheb
{    
    const Real* const restricted  ptr_y;
    Real* const restricted        ptr_ret;

    dp_cheb(const Real* y, Real* ret) : ptr_y(y), ptr_ret(ret){};
};

struct default_codegen;

template<Integer M>
struct dct_2
{
    using y                 = mk::gen_mat<M,1,tag_y>;
    using ret               = mk::output_mat<M,1,tag_ret>;
    using start_seed        = list_int<>;
    using ev_dct2           = eval_dct2<default_dct_config, y, start_seed, mk::one, false>;
    using expression_type   = typename ev_dct2::expression;

    using seed              = typename add_level<start_seed, 2>::type;
    using etag              = expr_tag<seed,M,1>;
    using evaler_type       = mk::expr_evaler<etag, default_codegen, ret, expression_type, double>;

    static void eval(const Real* restricted input, Real* restricted output,
                     Integer N, Integer in_ld, Integer out_ld)
    {
        for (Integer i = 0; i < N; ++i)
        {
            dp_cheb dp{input,output};
            evaler_type ev;
            ev.eval(dp);

            input   += in_ld;
            output  += out_ld;
        };
    };

    static void print(std::ostream& os)
    {
        os << "\n";
        evaler_type::print(os,0);
    };

    static mk::op_count get_op_count()
    {
        return evaler_type().get_op_count();
    };
};

template<Integer M>
struct dct_2_rec
{
    struct tag_temp { static void print(std::ostream& os,int)   { os << "temp"; }};

    using y                 = mk::gen_mat<M, 1, tag_y>;
    using ret               = mk::output_mat<M, 1, tag_ret>;
    using start_seed        = list_int<>;
    using ev_dct2           = eval_dct2<default_dct_config, y, start_seed, mk::one>;
    using expression_type   = typename ev_dct2::expression;

    using seed              = typename add_level<start_seed, 2>::type;
    using etag              = expr_tag<seed,M,2>;
    using evaler_type       = mk::expr_evaler<etag,default_codegen, ret, expression_type, double>;

    static void eval(const Real* restricted input, Real* restricted output,
                     Integer N, Integer in_ld, Integer out_ld)
    {
        #if ENABLE_IACA == 1
            IACA_START
        #endif

        for (Integer i = 0; i < N; ++i)
        {
            dp_cheb dp{input,output};
            evaler_type ev;
            ev.eval(dp);

            input   += in_ld;
            output  += out_ld;
        };

        #if ENABLE_IACA == 1
            IACA_END
        #endif
    };

    static void print(std::ostream& os)
    {
        os << "\n";
        evaler_type::print(os,0);
    };

    static mk::op_count get_op_count()
    {
        return evaler_type().get_op_count();
    };
};

template<Integer M, bool Rec>
struct dct_3
{
    using y                 = mk::gen_mat<M, 1, tag_y>;
    using ret               = mk::output_mat<M, 1, tag_ret>;
    using start_seed        = list_int<>;
    using ev_dct3           = eval_dct3<default_dct_config, y, start_seed, mk::one, Rec>;
    using expression_type   = typename ev_dct3::expression;

    using data_provider_type= dp_cheb;
    using seed              = typename add_level<start_seed, 3>::type;
    using etag              = expr_tag<seed,M,Rec>;

    using evaler_type       = mk::expr_evaler<etag, default_codegen, ret, expression_type, double>;

    static void eval(const Real* restricted input, Real* restricted output)
    {
        dp_cheb dp{input,output};
        evaler_type ev;
        ev.eval(dp);
    };

    static void eval(const Real* restricted input, Real* restricted output,
                     Integer N, Integer in_ld, Integer out_ld)
    {
        #if ENABLE_IACA == 1
            IACA_START
        #endif

        for (Integer i = 0; i < N; ++i)
        {
            dp_cheb dp{input,output};
            evaler_type ev;
            ev.eval(dp);

            input   += in_ld;
            output  += out_ld;
        };

        #if ENABLE_IACA == 1
            IACA_END
        #endif
    };

    static void print(std::ostream& os)
    {
        os << "\n";
        evaler_type::print(os,0);
    };

    static mk::op_count get_op_count()
    {
        return evaler_type().get_op_count();
    };
};

template<Integer M, bool Rec>
struct dct_1
{
    using y                 = mk::gen_mat<M, 1, tag_y>;
    using ret               = mk::output_mat<M, 1, tag_ret>;
    using start_seed        = list_int<>;
    using ev_dct1           = eval_dct1<default_dct_config, y, start_seed, mk::two, Rec>;
    using expression_type   = typename ev_dct1::expression;

    using data_provider_type= dp_cheb;
    using seed              = typename add_level<start_seed, 1>::type;
    using etag              = expr_tag<seed,M,Rec>;

    using evaler_type       = mk::expr_evaler<etag, default_codegen, ret, expression_type, double>;

    static void eval(const Real* restricted input, Real* restricted output)
    {
        dp_cheb dp{input,output};
        evaler_type ev;
        ev.eval(dp);
    };

    static void eval(const Real* restricted input, Real* restricted output,
                     Integer N, Integer in_ld, Integer out_ld)
    {
        #if ENABLE_IACA == 1
            IACA_START
        #endif

        evaler_type ev;
        for (Integer i = 0; i < N; ++i)
        {
            dp_cheb dp{input,output};            
            ev.eval(dp);

            input   += in_ld;
            output  += out_ld;
        };

        #if ENABLE_IACA == 1
            IACA_END
        #endif        
    };

    static void print(std::ostream& os)
    {
        os << "\n";
        evaler_type::print(os,0);
    };

    static mk::op_count get_op_count()
    {
        return evaler_type().get_op_count();
    };
};

template<Integer M, bool Rec>
struct dct_4
{
    using y                 = mk::gen_mat<M, 1, tag_y>;
    using ret               = mk::output_mat<M, 1, tag_ret>;
    using start_seed        = list_int<>;
    using ev_dct4           = eval_dct4<default_dct_config, y, start_seed, mk::one, Rec>;
    using expression_type   = typename ev_dct4::expression;

    using data_provider_type= dp_cheb;
    using seed              = typename add_level<start_seed, 4>::type;
    using etag              = expr_tag<seed,M,Rec>;

    using evaler_type       = mk::expr_evaler<etag, default_codegen, ret, expression_type, double>;

    static void eval(const Real* restricted input, Real* restricted output)
    {
        dp_cheb dp{input,output};
        evaler_type ev;
        ev.eval(dp);
    };
    static void eval(const Real* restricted input, Real* restricted output,
                     Integer N, Integer in_ld, Integer out_ld)
    {
        for (Integer i = 0; i < N; ++i)
        {
            dp_cheb dp{input,output};
            evaler_type ev;
            ev.eval(dp);

            input   += in_ld;
            output  += out_ld;
        };
    };

    static void print(std::ostream& os)
    {
        os << "\n";
        evaler_type::print(os,0);
    };

    static mk::op_count get_op_count()
    {
        return evaler_type().get_op_count();
    };
};

struct default_codegen
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

MATCL_FFT_EXPORT void fft::test_dct_II()
{
  #if 1
    using test_type = dct_kernels_test<DCT_SIZE>;

    test_type::print();
    test_type::check();
    test_type::eval_versions();
    test_type::test_perf();
  #endif

  #if 0
    static const Integer M      = 3;
    static const Integer M_OFF  = 1;
    #define DO_DEBUG 1

    using dct_rec       = dct_1<M, true && (M%2 == 1)>;
    using dct_nrec      = dct_1<M, false>;

    //using dct_rec     = dct_2_rec<M>;
    //using dct_nrec    = dct_2<M>;

    //using dct_rec     = dct_3<M,true>;
    //using dct_nrec    = dct_3<M,false>;

    //using dct_rec     = dct_4<M, true>;
    //using dct_nrec    = dct_4<M, false>;

  #if DO_DEBUG
    Integer N           = 1;
    Integer T           = 1;
  #else
    Integer N           = 100;
    Integer T           = 100000;
  #endif

    Matrix x            = randn(M+M_OFF,N);;
    const Real* X       = x.get_array<Real>();
    x                   = x(matcl::colon(1,M),matcl::colon());

  #if DO_DEBUG
    dct_rec::print(std::cout);
    dct_nrec::print(std::cout);
  #endif

    Real* restricted Y;
    Matrix y            = matcl::make_real_dense_noinit(M+M_OFF,N, Y);    
    y                   = y(matcl::colon(1,M),matcl::colon());

    Integer X_ld        = M+M_OFF;
    Integer Y_ld        = M+M_OFF;

  #if DO_DEBUG
    disp(x);    
  #endif

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Y;
        dct_nrec::eval(X2,Y2, N, X_ld, Y_ld);
    }
    Real t2 = matcl::toc();

  #if DO_DEBUG
    disp(y);
  #endif

    Real* restricted Z;
    Matrix z            = matcl::make_real_dense_noinit(M+M_OFF,N, Z);    
    z                   = z(matcl::colon(1,M),matcl::colon());

    matcl::tic();    
    for(int i = 0; i < T; ++i)
    {
        const Real* X2  = X;
        Real* Y2        = Z;

        dct_rec::eval(X2,Y2, N, X_ld, Y_ld);
    };    
    Real t1 = matcl::toc();

  #if DO_DEBUG
    disp(z);
    disp(div(y,z));
  #endif

  #if DO_DEBUG
    dct_nrec::get_op_count().print(std::cout);
    dct_rec::get_op_count().print(std::cout);
  #endif

  #if DO_DEBUG
    disp(odyn::fft::dct1(x,1));
    disp(div(odyn::fft::dct1(x,1),y));
  #endif

    matcl::tic();
    for(int i = 0; i < T; ++i)
    {
        Matrix Y        = odyn::fft::dct1(x,1);
    }
    Real t3 = matcl::toc();

    std::cout << M << ", " << N << ", " << T << "\n";
    std::cout   << "unrolled rec: " << t1 << "\n"
                << "unrolled: "     << t2 << ", rel: " << t2 / t1 << "\n"
                << ", dct2: "       << t3 << ", rel: " << t3 / t1 << "\n"
                << "\n";
  #endif
};

}};
