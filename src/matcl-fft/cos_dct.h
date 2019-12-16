#pragma once

#include "mkgen/mkgen.h"
#include "matcl-fft/matcl_fft_config.h"
#include "cos_dct_defs.h"

#include <iostream>

namespace matcl { namespace fft
{

namespace mk = mkgen;
namespace mkd = mk::details;

template<Integer M, Integer Step, class Scal>
struct cos_table
{
    public:        
        static alignas(VEC_ALIGN)
        Real                m_cos_table[M/Step + 1];
        static cos_table    instance;

        cos_table()         { initialize(); };        

        static cos_table    get() { return instance; };

        void initialize()
        {
            Real scal       = mk::get_scalar_value<Scal>::get<Real>();

            if constexpr(Step == 1)
            {
                m_cos_table[0]      = scal;

                for (Integer i = 1; i < M; ++i)
                {
                    m_cos_table[i]  = scal*cos(matcl::constants::pi() * double(i) / double(M));
                };
            }
            else
            {
                for (Integer i = 1, pos = 0; i < M; i += 2, ++pos)
                {
                    m_cos_table[pos]    = scal*cos(matcl::constants::pi() * double(i) / double(M));
                };
            };
        };
};
template<Integer M, Integer Step, class Scal>
Real cos_table<M,Step, Scal>::m_cos_table[] = {};

template<Integer M, Integer Step, class Scal>
cos_table<M,Step, Scal> cos_table<M,Step, Scal>::instance;

template<Integer M, class Scal>
struct cos_table_dct4
{
    public:        
        static alignas(VEC_ALIGN)
        Real                    m_table[M*M];
        static cos_table_dct4   instance;

        cos_table_dct4()        { initialize(); };        

        static cos_table_dct4   get() { return instance; };

        void initialize()
        {
            double scal         = mk::get_scalar_value<Scal>::get<double>();

            //sum_{n=0}^{N-1} x_n*cos( pi*(2*k+1)*(2*n+1)/(4*N) )
            for (Integer i = 0, k = 0; i < M; ++i)
            {
                for (Integer j = 0; j < M; ++j, ++k)
                {
                    double pos  = (2*j+1)*(2*i+1);
                    double val  = cos( matcl::constants::pi() * pos / double(4*M) );
                    m_table[k]  = val*scal;
                };
            };
        };
};
template<Integer M, class Scal>
Real cos_table_dct4<M,Scal>::m_table[];

template<Integer M,class Scal>
cos_table_dct4<M,Scal> cos_table_dct4<M,Scal>::instance;

template<Integer M, class Scal>
struct cos_table_dct1
{
    public:        
        static const Integer LD = ((M-1)/4 + 1) * 4;

        static alignas(VEC_ALIGN)
        Real                    m_table[LD*(M-2)];
        static cos_table_dct1   instance;

        cos_table_dct1()        { initialize(); };        

        static cos_table_dct1   get() { return instance; };

        void initialize()
        {
            double scal         = mk::get_scalar_value<Scal>::get<double>();

            //sum_{n=1}^{N-2} x_n*cos(pi*n*k/(N-1))

            Real* ptr           = m_table;

            for (Integer i = 1; i <= M - 2; ++i)
            {
                for (Integer j = 0; j < M; ++j)
                {
                    double pos  = j*i;
                    double val  = cos( matcl::constants::pi() * pos / double(M-1) );
                    ptr[j]      = scal*val;
                };

                ptr             += LD;
            };
        };
};
template<Integer M, class Scal>
Real cos_table_dct1<M,Scal>::m_table[];

template<Integer M,class Scal>
cos_table_dct1<M,Scal> cos_table_dct1<M,Scal>::instance;

template<Integer M, class Scal>
struct cos_table_dct2
{
    public:        
        static alignas(VEC_ALIGN)
        Real                    m_table[M*M];
        static cos_table_dct2   instance;

        cos_table_dct2()        { initialize(); };        

        static cos_table_dct2   get() { return instance; };

        void initialize()
        {
            double scal         = mk::get_scalar_value<Scal>::get<double>();

            //sum_{n=0}^{N-1} x_n*cos( pi*k*(2*n+1)/(2*N) )
            for (Integer i = 0, k = 0; i < M; ++i)
            {
                for (Integer j = 0; j < M; ++j, ++k)
                {
                    double pos  = j*(2*i+1);
                    double val  = cos( matcl::constants::pi() * pos / double(2*M) );
                    m_table[k]  = scal*val;
                };
            };
        };
};
template<Integer M,class Scal>
Real cos_table_dct2<M,Scal>::m_table[];

template<Integer M,class Scal>
cos_table_dct2<M,Scal> cos_table_dct2<M,Scal>::instance;

template<Integer M, class Scal>
struct cos_table_dct3
{
    public:        
        static alignas(VEC_ALIGN)
        Real                    m_table[M*(M-1)];
        static cos_table_dct3   instance;

        cos_table_dct3()        { initialize(); };        

        static cos_table_dct3   get() { return instance; };

        void initialize()
        {
            double scal         = mk::get_scalar_value<Scal>::get<double>();

            //sum_{n=1}^{N-1} x_n*cos( pi*n*(2*k+1)/(2*N) )
            for (Integer i = 1, k = 0; i < M; ++i)
            {
                for (Integer j = 0; j < M; ++j, ++k)
                {
                    double pos  = i*(2*j+1);
                    double val  = cos( matcl::constants::pi() * pos / double(2*M) );
                    m_table[k]  = scal*val;
                };
            };
        };
};

template<Integer M,class Scal>
Real cos_table_dct3<M,Scal>::m_table[];

template<Integer M,class Scal>
cos_table_dct3<M,Scal> cos_table_dct3<M,Scal>::instance;

template<Integer Base_Ind, Integer Step>
struct get_pos_cos_table
{
};
template<Integer Base_Ind>
struct get_pos_cos_table<Base_Ind,1>
{
    static const Integer value = Base_Ind;
};
template<Integer Base_Ind>
struct get_pos_cos_table<Base_Ind,2>
{
    //TODO: 
    //static_assert((Base_Ind-1)%2 == 0, "something is wrong");
    static const Integer value = (Base_Ind-1) / 2;
};

template<Integer M, Integer Base_Ind, Integer Step, class Scal>
struct cos_value_tag : mkd::scal_data_value_tag<cos_value_tag<M, Base_Ind, Step, Scal>>
{
    static const Integer pos    = get_pos_cos_table<Base_Ind,Step>::value;
    static const bool           is_continuous   = true; 
    using                       root_align_type = mk::align_full; 
    static const Integer        get_offset      = pos;

    using tag                   = cos_value_tag;

    template<class T, Integer Offset>
    inline_lev_1
    static double* get_data_ptr()
    {
        return cos_table<M,Step,Scal>::get().m_cos_table + Offset;
    };

    template<class Val>
    inline_lev_1
    static Val eval()
    {
        return Val(cos_table<M,Step,Scal>::get().m_cos_table[pos]);
    };

    static void print(std::ostream& os,int)
    {
        //os << "cos[" << pos << "]";
        os << cos_table<M, Step, Scal>::get().m_cos_table[pos];
    };

    //TODO
    /*
    // TODO: simplify this
    template<Integer Step2, class Arr_List>
    using get_arrays    = typename mk::push_back<Arr_List, mk::details::array_item<cos_value_tag,Step2,
                            mk::details::array_item_extern>> :: type;

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        vis.visit_load();
    };
    */
};

template<Integer Index, Integer M>
struct cos_info
{
    private:
        static const Integer ind                = Index % (2*M);
        static const Integer base_ind1          = ind % M;
        static const Integer base_ind2          = (2*base_ind1 > M)? M - base_ind1 : base_ind1;
        static const bool negate1               = ind > M;
        static const bool negate2               = 2*base_ind1 > M;

    public:
        static const Integer base_index         = base_ind2;
        static const bool negate                = negate1 != negate2;
        static const bool zero                  = (2*ind == M || 2*ind == 3*M);
        static const bool one                   = (ind == 0);
        static const bool mone                  = (one == false) && (base_ind1 == 0);
        static const bool half                  = (base_index * 3 == M);
};

template<Integer M, Integer Base_Ind, Integer Step, bool Negate, bool Zero, bool One, 
    bool Mone, bool Half, class Scal>
struct get_element_cos
{
    using tag       = cos_value_tag<M, Base_Ind, Step, Scal>;
    using scal_tag  = tag;
    using type      = mk::value_scalar<tag, double>;
};

template<Integer M, Integer Base_Ind, Integer Step, class Scal>
struct get_element_cos<M,Base_Ind,Step,true, false,false,false,false,Scal>
{
    using tag       = cos_value_tag<M, Base_Ind, Step, Scal>;
    using scal_tag  = tag;
    using type_cos  = mk::value_scalar<scal_tag, double>;
    using type      = typename mk::make_mult<mk::mone, type_cos>::type;
};

template<Integer M, Integer Base_Ind, Integer Step, bool Negate,class Scal>
struct get_element_cos<M,Base_Ind,Step,Negate,true,false,false,false,Scal>
{
    using type = mk::zero;
};

template<Integer M, Integer Base_Ind, Integer Step, bool Negate,class Scal>
struct get_element_cos<M,Base_Ind,Step,Negate,false,true,false,false,Scal>
{
    using type = typename mk::make_mult<mk::one,Scal>::type;
};

template<Integer M, Integer Base_Ind, Integer Step, bool Negate,class Scal>
struct get_element_cos<M,Base_Ind,Step,Negate,false,false,true,false,Scal>
{
    using type = typename mk::make_mult<mk::mone,Scal>::type;
};

template<Integer M, Integer Base_Ind, Integer Step, class Scal>
struct get_element_cos<M,Base_Ind,Step,false,false,false,false,true,Scal>
{
    using type = typename mk::make_mult<mk::half,Scal>::type;
};

template<Integer M, Integer Base_Ind, Integer Step, class Scal>
struct get_element_cos<M,Base_Ind,Step,true,false,false,false,true,Scal>
{
    using type = typename mk::make_mult<mk::rational_scalar<-1,2>,Scal>::type;
};

};};
