/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "matcl-matrep/func/raw/raw_vecfunc.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/func/raw/eval_vec_functor.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-matrep/func/raw/eval_vec_functor2.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-matrep/lib_functions/eval_functors.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-scalar/details/scalfunc_helpers.h"

namespace matcl { namespace raw { namespace details
{

namespace md    = matcl::details;
namespace mr    = matcl::raw;
namespace mrd   = matcl::raw::details;
namespace mdyf  = matcl::dynamic::functions;

//-----------------------------------------------------------
//              utils
//-----------------------------------------------------------
template<class T1, class T2>
inline bool less(const T1& A, const T2& B)        
{   
    if (mrd::isnan_helper<T1>::eval(A))
        return false;

    if (mrd::isnan_helper<T2>::eval(B))
        return true;

    return (bool)details::lt_helper<T1,T2>::eval(A,B); 
};

template<class T1, class T2>
inline bool greater(const T1& A, const T2& B)        
{
    if (mrd::isnan_helper<T1>::eval(A))
        return false;

    if (mrd::isnan_helper<T2>::eval(B))
        return true;

    return (bool)details::gt_helper<T1,T2>::eval(A,B); 
};

template<class V, class W>
static void mult_assign(V& ret, const W& in)
{ 
    ret *= in; 
};

template<class W>
static void mult_assign(Object& ret, const W& in) 
{ 
    ret = ret * in;
};

template<class T> struct min_value      {};
template<> struct min_value<Integer>    { static Integer eval() {return constants::min_int();}; };
template<> struct min_value<Real>       { static Real eval()    {return -constants::inf();};    };
template<> struct min_value<Float>      { static Float eval()   {return -constants::f_inf();};  };
template<> struct min_value<Complex>    { static Complex eval() {return -constants::inf();};    };
template<> struct min_value<Float_complex> 
                                        { static Float_complex eval() {return -constants::f_inf();};    };

template<class T> struct max_value      {};
template<> struct max_value<Integer>    { static Integer eval() {return constants::max_int();};    };
template<> struct max_value<Real>       { static Real eval()    {return constants::inf();};        };
template<> struct max_value<Float>      { static Float eval()   {return constants::f_inf();};        };
template<> struct max_value<Complex>    { static Complex eval() {return constants::inf();};        };
template<> struct max_value<Float_complex> 
                                        { static Float_complex eval() {return constants::f_inf();};};

//-----------------------------------------------------------
//              accumulators
//-----------------------------------------------------------

template<class in_type, class value_type>
class nnz_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  state_array;
        const value_type            Z;

        nnz_accumulator& operator=(const nnz_accumulator&) = delete;

    public:
        nnz_accumulator(ti::ti_type<value_type> ti)
            :state(md::default_value<value_type>(ti))
            ,Z(md::default_value<value_type>(ti))
            ,state_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer )         {};
        void reset()                            { state = Z;            };
        bool add_zero()                         { return false;         };
        void add_zero(Integer )                 {                       };

        const value_type& value()               { return state;         };        
        const value_type& value(Integer p)      { return state_array[p];};

        void reset_array(Integer s)
        { 
            state_array.resize(s,Z);
        };    

        bool add(const in_type& value)
        { 
            state += (mrd::is_zero(value)?0:1); 
            return false;
        };

        void add(Integer p,const in_type& val)
        { 
            state_array[p] += (mrd::is_zero(val)?0:1);
        };
};

template<class in_type, class value_type>
class sum_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  state_array;
        const value_type            Z;

        sum_accumulator& operator=(const sum_accumulator&) = delete;

    public:
        sum_accumulator(ti::ti_type<value_type> ti)
            :state(md::default_value<value_type>(ti))
            ,Z(md::default_value<value_type>(ti))
            ,state_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            //many operations, we require that ti does not chane
            return ti;
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer )         {};
        void reset()                            { state = Z;            };
        bool add_zero()                         { return false;         };
        void add_zero(Integer )                 {                       };

        const value_type& value()               { return state;         };
        const value_type& value(Integer p)      { return state_array[p];};        

        void reset_array(Integer s)                
        { 
            state_array.resize(s,Z);
        };    

        bool add(const in_type& value)
        { 
            state = mrd::plus_helper<value_type,in_type>::eval(state,value); 
            return false;
        };

        void add(Integer p,const in_type& val)
        { 
            state_array[p] = mrd::plus_helper<value_type,in_type>::eval(state_array[p],val);
        };
};

template<class in_type, class value_type>
struct prod_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  state_array;
        const value_type            Z;
        const value_type            one;

        prod_accumulator& operator=(const prod_accumulator&) = delete;

    public:
        prod_accumulator(ti::ti_type<value_type> ti)
            :state(md::default_value<value_type>(ti)) 
            , Z(md::default_value<value_type>(ti))
            , state_array(ti)
            , one(md::one_value<value_type>(ti))
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            //many operations, we require that ti does not chane
            return ti;
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer )         {};
        void reset()                            { state = one;              };
        bool add_zero()                         { state = Z; return true;   };
        void add_zero(Integer p)                { state_array[p] = Z;       };

        const value_type& value()               { return state;             };
        const value_type& value(Integer p)      { return state_array[p];    };        

        void reset_array(Integer s)
        { 
            state_array.resize(s,one);
        };    

        bool add(const in_type& value)
        { 
            mult_assign(state , value);
            return mrd::is_zero(state);
        };

        void add(Integer p,const in_type& val)
        { 
            mult_assign(state_array[p] , val);
        };
};

template<class in_type, class value_type>
class sumsq_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  state_array;
        const value_type            Z;

        sumsq_accumulator& operator=(const sumsq_accumulator&) = delete;

    public:
        sumsq_accumulator(ti::ti_type<value_type> ti)    
            : state(md::default_value<value_type>(ti)) 
            , Z(md::default_value<value_type>(ti))
            , state_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            //many operations, we require that ti does not chane
            return ti;
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer )         {};
        void reset()                            { state = Z;            };
        bool add_zero()                         { return false;         };
        void add_zero(Integer )                 {                       };

        const value_type& value()               { return state;         };
        const value_type& value(Integer p)      { return state_array[p];};        

        void reset_array(Integer s)                
        { 
            state_array.resize(s,Z);
        };

        bool add(const in_type& val)
        { 
            auto tmp    = val * val;
            state       = mrd::plus_helper<value_type,value_type>::eval(state,tmp); 
            return false;
        };

        void add(Integer p,const in_type& val)
        { 
            auto tmp        = val * val;
            state_array[p]  = mrd::plus_helper<value_type,value_type>::eval(state_array[p],tmp);
        };
};

template<class in_type, class value_type>
struct min_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        bool added;
        md::workspace2<value_type>  st_array;
        matcl::pod_workspace<bool>  st_array_flag;
        const value_type            Z;        

        min_accumulator& operator=(const min_accumulator&) = delete;

    public:
        min_accumulator(ti::ti_type<value_type> ti)
            :Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti))
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        void set_size(Integer)              { };
        void set_size(Integer,Integer )     {};
        void reset()                        { added = false;     };
        const value_type& value()           { return state;      };
        const value_type& value(Integer p)  { return st_array[p];};        

        bool add(const in_type& val)                
        { 
            if (added)  { state = min_nan(val,state);   }
            else        { state = val; added = true;    };

            return false;        
        };

        bool add_zero()                
        { 
            if (added)  { state = min_nan(Z,state); }
            else        { state = Z; added = true;  }
            
            return false;        
        };

        void reset_array(Integer s)            
        { 
            st_array.resize(s);    
            st_array_flag.resize(s,false);    
        };            

        void add(Integer p,const in_type& val)        
        { 
            if (st_array_flag[p])
            { 
                st_array[p] = min_nan(val,st_array[p]); 
            }
            else
            { 
                st_array[p] = val; 
                st_array_flag[p] = true; 
            };
        };        

        void add_zero(Integer p)            
        { 
            if (st_array_flag[p])
            { 
                st_array[p] = min_nan(Z,st_array[p]); 
            }
            else
            { 
                st_array[p] = Z; 
                st_array_flag[p] = true; 
            };        
        };

    private:
        value_type min_nan(const value_type& a, const value_type& b)
        {
            return mrd::min_helper<value_type,value_type>::eval(a,b);
        };
};

template<class in_type, class value_type>
struct max_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        bool                        added;
        md::workspace2<value_type>  st_array;
        matcl::pod_workspace<bool>  st_array_flag;
        const value_type            Z;

        max_accumulator& operator=(const max_accumulator&) = delete;

    public:
        max_accumulator(ti::ti_type<value_type> ti)    
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti))
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        void set_size(Integer )             {};
        void set_size(Integer,Integer )     {};
        void reset()                        { added = false;        };
        const value_type& value()           { return state;         };
        const value_type& value(Integer p)  { return st_array[p];   };

        bool add(const in_type& val)
        { 
            if (added)  { state = max_nan(val,state);   }
            else        { state = val; added = true;    };

            return false;        
        };

        bool add_zero()                        
        { 
            if (added)  { state = max_nan(Z,state); }
            else        { state = Z; added = true;  }
            
            return false;        
        };

        void reset_array(Integer s)            
        { 
            st_array.resize(s);    
            st_array_flag.resize(s,false);    
        };            

        void add(Integer p,const in_type& val)        
        { 
            if (st_array_flag[p])
            { 
                st_array[p] = max_nan(val,st_array[p]); 
            }
            else
            { 
                st_array[p]         = val;
                st_array_flag[p]    = true; 
            };
        };        

        void add_zero(Integer p)            
        { 
            if (st_array_flag[p])
            { 
                st_array[p] = max_nan(Z,st_array[p]); 
            }
            else
            { 
                st_array[p] = Z; 
                st_array_flag[p] = true; 
            };        
        };

    private:
        value_type max_nan(const value_type& a, const value_type& b)
        {
            return mrd::max_helper<value_type,value_type>::eval(a,b);
        };
};

template<class in_type, class value_type>
struct max_abs_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        bool                        added;
        md::workspace2<value_type>  st_array;
        matcl::pod_workspace<bool>  st_array_flag;
        const value_type            Z;

        max_abs_accumulator& operator=(const max_abs_accumulator&) = delete;

    public:
        max_abs_accumulator(ti::ti_type<value_type> ti)    
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti))
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {          
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        void set_size(Integer )             {};
        void set_size(Integer,Integer )     {};
        void reset()                        { added = false;        };
        const value_type& value()           { return state;         };
        const value_type& value(Integer p)  { return st_array[p];   };

        bool add(const in_type& val)                
        { 
            auto aval = mrd::abs_helper<in_type>::eval(val);

            if (added)  { state = max_nan(aval,state);   }
            else        { state = aval; added = true;    };

            return false;
        };

        bool add_zero()                        
        {
            if (added)  { state = max_nan(Z,state); }
            else        { state = Z; added = true;  }
            
            return false;
        };

        void reset_array(Integer s)            
        { 
            st_array.resize(s);    
            st_array_flag.resize(s,false);    
        };            

        void add(Integer p,const in_type& val)        
        { 
            auto aval = mrd::abs_helper<in_type>::eval(val);

            if (st_array_flag[p])
            { 
                st_array[p] = max_nan(aval,st_array[p]); 
            }
            else
            { 
                st_array[p] = aval; 
                st_array_flag[p] = true;
            };
        };        

        void add_zero(Integer p)
        {
            if (st_array_flag[p])
            { 
                st_array[p] = max_nan(Z,st_array[p]); 
            }
            else
            { 
                st_array[p] = Z; 
                st_array_flag[p] = true; 
            };        
        };

    private:
        value_type max_nan(const value_type& a, const value_type& b)
        {            
            return mrd::max_helper<value_type,value_type>::eval(a,b);
        };
};

template<class in_type, class value_type>
struct max_accumulator2
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  st_array;
        integer_dense               ind_array;
        Integer                     vec_pos;
        const value_type            Z;

        max_accumulator2& operator=(const max_accumulator2&) = delete;

    public:
        max_accumulator2(ti::ti_type<value_type> ti)    
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,ind_array(ti::ti_empty())
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        void set_size(Integer , Integer c, int dim)    
        { 
            Integer mr = 1, mc = 1;
            if (dim == 1)   mc = c;
            else            mr = c;

            ind_array.reset_unique(mr,mc);
            Integer* ptr_ind = ind_array.ptr();

            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = 0;
        };

        void reset()                                { };
        const value_type& value()                   { return state;         };
        const value_type& value(Integer p)          { return st_array[p];   };        
        const integer_dense& index_matrix() const   { return ind_array;     };
        void current_vector(Integer vec)            { vec_pos = vec;        };

        void reset_array(Integer s)
        { 
            st_array.resize(s);
        };    
        
        bool add(Integer v, const in_type& val)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || greater(val,state))
            {
                state = val;
                ind_array.ptr()[vec_pos] = v+1;                
            }
            return false;
        };

        void add(Integer p,Integer v, const in_type& val)        
        { 
            if (ind_array.ptr()[p] == 0 || greater(val,st_array[p]))
            {
                st_array[p] = val;        
                ind_array.ptr()[p] = v+1;
            }
        };

        bool add_zero(Integer v)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || greater(Z,state))
            {
                state = Z;
                ind_array.ptr()[vec_pos] = v+1;
            }
            else if(mrd::is_zero(state))
            {
                ind_array.ptr()[vec_pos] = std::min(v+1,ind_array.ptr()[vec_pos]);
            };
            return false;        
        };

        void add_zero(Integer p,Integer v)            
        { 
            if (ind_array.ptr()[p] == 0 || greater(Z,st_array[p]))
            {
                st_array[p] = Z;
                ind_array.ptr()[p] = v+1;
            }
            else if( mrd::is_zero(st_array[p]))
            {
                ind_array.ptr()[p] = std::min(v+1,ind_array.ptr()[p]);
            };
        };
};

template<class in_type, class value_type>
struct max_abs_accumulator2
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  st_array;
        integer_dense               ind_array;
        Integer                     vec_pos;
        const value_type            Z;

        max_abs_accumulator2&    operator=(const max_abs_accumulator2&) = delete;

    public:
        max_abs_accumulator2(ti::ti_type<value_type> ti)    
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,ind_array(ti::ti_empty())
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        void set_size(Integer , Integer c, int dim)    
        { 
            Integer mr = 1, mc = 1;
            if (dim == 1)   mc = c;
            else            mr = c;

            ind_array.reset_unique(mr,mc);
            Integer* ptr_ind = ind_array.ptr();

            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = 0;
        };

        void reset()                                { };        
        void current_vector(Integer vec)            { vec_pos = vec;        };
        const value_type& value()                   { return state;         };
        const value_type& value(Integer p)          { return st_array[p];   };        
        const integer_dense& index_matrix() const   { return ind_array;     };

        void reset_array(Integer s)
        { 
            st_array.resize(s);
        };    

        bool add(Integer v, const in_type& val)
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);
            if (ind_array.ptr()[vec_pos] == 0 || greater(aval,state))
            {
                state = aval;
                ind_array.ptr()[vec_pos] = v+1;                
            }
            return false;
        };

        void add(Integer p,Integer v, const in_type& val)        
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);
            if (ind_array.ptr()[p] == 0 || greater(aval,st_array[p]))
            {
                st_array[p] = aval;        
                ind_array.ptr()[p] = v+1;
            }
        };

        bool add_zero(Integer v)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || greater(Z,state))
            {
                state = Z;
                ind_array.ptr()[vec_pos] = v+1;
            }
            else if( mrd::is_zero(state))
            {
                ind_array.ptr()[vec_pos] = std::min(v+1,ind_array.ptr()[vec_pos]);
            };
            return false;        
        };

        void add_zero(Integer p,Integer v)            
        { 
            if (ind_array.ptr()[p] == 0 || greater(Z,st_array[p]))
            {
                st_array[p] = Z;
                ind_array.ptr()[p] = v+1;
            }
            else if( mrd::is_zero(st_array[p]))
            {
                ind_array.ptr()[p] = std::min(v+1,ind_array.ptr()[p]);
            };
        };
};

template<class in_type, class value_type>
class mean_accumulator 
{
    private:
        using in_ti         = ti::ti_type<in_type>;
        using ret_ti        = ti::ti_type<value_type>;
        using value_real    = typename md::real_type<value_type>::type;

    private:
        Integer                     rows;
        Integer                     cols;
        value_type                  state;
        value_type                  ret;
        md::workspace2<value_type>  st_array;
        const value_type            Z;

        mean_accumulator& operator=(const mean_accumulator&) = delete;

    public:
        mean_accumulator(ti::ti_type<value_type> ti)    
            :state(md::default_value<value_type>(ti)) 
            ,ret(md::default_value<value_type>(ti)) 
            ,Z(md::default_value<value_type>(ti))
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            //many sum operations, so we require that ti does not change
            //resulting ti is div(ti,Integer)
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::op_div::eval(), ti, 
                                                   ti::predefined::get_ti_int());
            return ret;
        };

        void set_size(Integer r)                { rows = r; cols = 1;     };
        void set_size(Integer r,Integer c)      { rows = r; cols = c;     };
        void reset()                            { state = Z;              };
        void reset_array(Integer s)             { st_array.resize(s,Z);   };            
        bool add_zero()                         { return false;           };
        void add_zero(Integer )                 {                         };

        bool add(const in_type& value)
        { 
            state = mrd::plus_helper<value_type,in_type>::eval(state,value);
            return false;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p]     = mrd::plus_helper<value_type,in_type>::eval(st_array[p],val);            
        };

        const value_type& value()
        { 
            value_real size = value_real(rows) * value_real(cols);
            ret             = mrd::div_helper<value_type,value_real>::eval(state,size);
            return ret;
        };

        const value_type& value(Integer p)
        { 
            value_real size = value_real(rows) * value_real(cols);

            ret             = mrd::div_helper<value_type,value_real>::eval(st_array[p],size);
            return ret;
        };
};

template<class in_type, class value_type>
class std_accumulator 
{
    private:
        using in_ti         = ti::ti_type<in_type>;
        using ret_ti        = ti::ti_type<value_type>;
        using value_real    = typename md::real_type_int_real<value_type>::type;

    private:
        struct state
        {
            value_type      sum;
            value_real      sum_sq;

            state(ret_ti ti) 
                : sum(md::default_value<value_type>(ti))
                , sum_sq(md::default_value<value_real>(ti)) 
            {};
        };
        
        using workspace     = md::workspace2<state,ret_ti>;

        Integer             rows;
        Integer             cols;
        state               m_state;
        value_real          m_ret;
        workspace           st_array;
        const value_type    Z;
        const state         ZS;
        bool                m_unbiased;

        std_accumulator& operator=(const std_accumulator&) = delete;

    public:
        std_accumulator(ret_ti ti, bool unbiased)
            : Z(md::default_value<value_type>(ti)) 
            , m_ret(md::default_value<value_real>(ti)) 
            , ZS(ti), m_state(ti), st_array(ti)
            ,m_unbiased(unbiased)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            //many different operations, so we require that ti does not change
            //resulting ti is div(ti,Integer)
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::op_div::eval(), ti, 
                                                   ti::predefined::get_ti_int());
            return ret;
        };

        void set_size(Integer r)                { rows = r; cols = 1;   };
        void set_size(Integer r,Integer c)      { rows = r; cols = c;   };
        void reset()                            { m_state = ZS;         };
        bool add_zero()                         { return false;         };
        void add_zero(Integer )                 {};

        void reset_array(Integer s)
        { 
            st_array.clear(); 
            st_array.resize(s,ZS);
        };    

        bool add(const in_type& val)
        { 
            m_state.sum     = mrd::plus_helper<value_type,in_type>::eval(m_state.sum,val);

            auto tmp        = mrd::abs_helper<value_type>::eval(val);
            m_state.sum_sq  = mrd::plus_helper<value_real, value_real>::eval(m_state.sum_sq,tmp);
            return false;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p].sum     = mrd::plus_helper<value_type,in_type>::eval(st_array[p].sum,val);

            auto tmp            = mrd::abs_helper<value_type>::eval(val);
            st_array[p].sum_sq  = mrd::plus_helper<value_real, value_real>::eval(st_array[p].sum_sq, tmp);
        };

        value_type value()
        { 
            if (rows <= 1 && cols <= 1)
                return Z;

            value_real tmp  = mrd::abs_helper<value_type>::eval(m_state.sum);
            value_real size = value_real(rows) * value_real(cols);

            tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,size);
            tmp     = mrd::minus_helper<value_real,value_real>::eval(m_state.sum_sq,tmp);

            if (m_unbiased)
            {
                value_real size_m1  = mrd::minus_helper<value_real,value_real>::eval(size,value_real(1));
                tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,size_m1);
            }
            else
                tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,size);

            m_ret   = mrd::sqrt_helper<value_real>::eval(tmp);

            return m_ret;
        };

        value_type value(Integer p)
        { 
            if (rows <= 1 && cols << 1)
                return Z;

            value_real size = value_real(rows) * value_real(cols);

            value_real tmp = mrd::abs_helper<value_type>::eval(st_array[p].sum);

            tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,size);
            tmp     = mrd::minus_helper<value_real,value_real>::eval(st_array[p].sum_sq,tmp);

            if (m_unbiased)
            {
                value_real size_m1  = mrd::minus_helper<value_real,value_real>::eval(size,value_real(1));
                tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,size_m1);
            }
            else
                tmp     = mrd::div_helper<value_real,value_real>::eval(tmp,value_real(size));

            m_ret   = mrd::sqrt_helper<value_real>::eval(tmp);

            return m_ret;
        };
};

template<class in_type, class value_type>
struct any_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        Integer                         state;
        matcl::pod_workspace<Integer>   st_array;

    public:
        any_accumulator(ti::ti_type<value_type> )    {};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer)          {};
        void reset()                            { state = false;        };
        bool add_zero()                         { return false;         };
        void add_zero(Integer )                 {                       };
        const value_type& value()               { return state;         };
        const value_type& value(Integer p)      { return st_array[p];   };

        void reset_array(Integer s)
        { 
            st_array.resize(s,false);
        };    

        bool add(const in_type& val)
        { 
            state = or_helper(state?true:false,val);
            return state?true:false;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p] = or_helper(st_array[p]?true:false,val);
        };

        template<class in_type_>
        bool or_helper(bool state,const in_type_& val)
        {
            return state || mrd::cast_bool_helper<in_type_>::eval(val);
        };
};

template<class in_type, class value_type>
struct any_accumulator_t 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        Integer                         state;
        matcl::pod_workspace<Integer>   st_array;
        const test_function&            test;
        const in_type                   Z;        

    public:
        any_accumulator_t(ti::ti_type<in_type> ti,const test_function& t)
            : test(t) 
            , Z(md::default_value<in_type>(ti)) 
        {};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer)          {};
        void reset()                            { state = false;    };
        bool add_zero()                         { return add(Z);    };
        void add_zero(Integer p)                { return add(p,Z);  };

        const value_type& value()               { return state;      };
        const value_type& value(Integer p)      { return st_array[p];};

        void reset_array(Integer s)
        { 
            st_array.resize(s,false);
        };    

        bool add(const in_type& val)
        { 
            state = (state||test.eval(val));
            return state?true:false;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p] = (st_array[p] || test.eval(val));
        };

    private:
        any_accumulator_t(const any_accumulator_t&) = delete;
        any_accumulator_t& operator=(const any_accumulator_t&) = delete;
};

template<class in_type, class value_type>
struct all_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        Integer                         state;
        matcl::pod_workspace<Integer>   st_array;

    public:
        all_accumulator(ti::ti_type<value_type> ){};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer)          {};
        void reset()                            { state = true;                 };
        bool add_zero()                         { state = false; return true;   };
        void add_zero(Integer p)                { st_array[p] = false;          };

        const value_type& value()               { return state;                 };
        const value_type& value(Integer p)      { return st_array[p];           };

        void reset_array(Integer s)
        { 
            st_array.resize(s,true);
        };    

        bool add(const in_type& val)
        { 
            state = and_helper(state?true:false,val);
            return !state;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p] = and_helper(st_array[p]?true:false,val);
        };

        template<class in_type_>
        bool and_helper(bool state,const in_type_& val)
        {
            return state && mrd::cast_bool_helper<in_type_>::eval(val);
        };
};

template<class in_type, class value_type>
struct all_accumulator_t
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        Integer                         state;
        matcl::pod_workspace<Integer>   st_array;
        const test_function&            test;
        const in_type                   Z;        

    public:
        all_accumulator_t(ti::ti_type<in_type> ti,const test_function& t)
            : test(t)
            , Z(md::default_value<in_type>(ti)) 
        {};

        static ret_ti get_ret_ti(in_ti)
        {
            return ret_ti();
        };

        void set_size(Integer )                 {};
        void set_size(Integer,Integer)          {};
        void reset()                            { state = true;         };
        bool add_zero()                         { return add(Z);        };
        void add_zero(Integer p)                { return add(p,Z);      };

        const value_type& value()               { return state;         };
        const value_type& value(Integer p)      { return st_array[p];   };
        
        void reset_array(Integer s)
        { 
            st_array.resize(s,true);
        };    

        bool add(const in_type& val)
        { 
            state = (state&&test.eval(val));
            return !state;
        };

        void add(Integer p,const in_type& val)
        { 
            st_array[p] = (st_array[p] && test.eval(val));
        };

    private:
        all_accumulator_t(const all_accumulator_t&) = delete;
        all_accumulator_t& operator=(const all_accumulator_t&) = delete;
};

template<class in_type, class value_type>
struct min_abs_accumulator 
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        bool                        added;
        md::workspace2<value_type>  st_array;
        matcl::pod_workspace<bool>  st_array_flag;
        const value_type            Z;        

        min_abs_accumulator& operator=(const min_abs_accumulator&) = delete;

    public:
        min_abs_accumulator(ti::ti_type<value_type> ti)
            :Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti))
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        void set_size(Integer)              { };
        void set_size(Integer,Integer)      {};
        void reset()                        { added = false;    };
        const value_type& value()           { return state;     };
        const value_type& value(Integer p)  { return st_array[p];};        

        bool add(const in_type& val)                
        { 
            auto aval = mrd::abs_helper<in_type>::eval(val);

            if (added)  { state = min_nan(aval,state);   }
            else        { state = aval; added = true;    };

            return false;
        };

        bool add_zero()                        
        { 
            if (added)  { state = min_nan(Z,state); }
            else        { state = Z; added = true;  }
            
            return true;
        };

        void reset_array(Integer s)
        { 
            st_array.resize(s);    
            st_array_flag.resize(s,false);
        };            

        void add(Integer p,const in_type& val)    
        { 
            auto aval = mrd::abs_helper<in_type>::eval(val);

            if (st_array_flag[p])
            { 
                st_array[p] = min_nan(aval,st_array[p]); 
            }
            else
            { 
                st_array[p] = aval;
                st_array_flag[p] = true;
            };
        };        

        void add_zero(Integer p)
        { 
            if (st_array_flag[p])
            { 
                st_array[p] = min_nan(Z,st_array[p]); 
            }
            else
            { 
                st_array[p] = Z; 
                st_array_flag[p] = true; 
            };
        };

    private:
        value_type min_nan(const value_type& a, const value_type& b)
        {
            return mrd::min_helper<value_type,value_type>::eval(a,b);
        };
};

template<class in_type, class value_type>
struct min_accumulator2
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  st_array;
        integer_dense               ind_array;
        Integer                     vec_pos;
        const value_type            Z;

        min_accumulator2& operator=(const min_accumulator2&) = delete;

    public:
        min_accumulator2(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,ind_array(ti::ti_empty())
            ,st_array(ti)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        void set_size(Integer , Integer c, Integer dim)    
        { 
            Integer mr = 1, mc = 1;
            if (dim == 1)   mc = c;
            else            mr = c;

            ind_array.reset_unique(mr,mc);
            Integer* ptr_ind = ind_array.ptr();

            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = 0;
        };

        void reset()                                { };
        void current_vector(Integer vec)            { vec_pos = vec;        };
        const value_type& value()                   { return state;         };
        const value_type& value(Integer p)          { return st_array[p];   };        
        const integer_dense& index_matrix() const   { return ind_array;     };

        void reset_array(Integer s)
        { 
            st_array.resize(s);
        };    
        
        bool add(Integer v, const in_type& val)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || less(val,state))
            {
                state = val;
                ind_array.ptr()[vec_pos] = v+1;
            }
            return false;
        };

        bool add_zero(Integer v)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || less(Z,state))
            {
                state = Z;
                ind_array.ptr()[vec_pos] = v+1;
            }
            else if (mrd::is_zero(state))
            {
                ind_array.ptr()[vec_pos] = std::min(v+1,ind_array.ptr()[vec_pos]);
            };
            return false;        
        };

        void add(Integer p,Integer v, const in_type& val)        
        { 
            if (ind_array.ptr()[p] == 0 || less(val,st_array[p]))
            {
                st_array[p] = val;        
                ind_array.ptr()[p] = v+1;
            };
        };

        void add_zero(Integer p,Integer v)            
        { 
            if (ind_array.ptr()[p] == 0 || less(Z,st_array[p]))
            {
                st_array[p] = Z;
                ind_array.ptr()[p] = v+1;
            }
            else if(mrd::is_zero(st_array[p]))
            {
                ind_array.ptr()[p] = std::min(v+1,ind_array.ptr()[p]);
            };
        }        
};

template<class in_type, class value_type>
struct min_accumulator2_vec
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type              state;
        Integer                 col_pos;
        Integer                 row_pos;
        Integer                 col_pos_best;
        const value_type        Z;

        min_accumulator2_vec& operator=(const min_accumulator2_vec&) = delete;

    public:
        min_accumulator2_vec(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,row_pos(0), col_pos(0), col_pos_best(0)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        Integer index_row() const           { return row_pos;       };
        Integer index_col() const           { return col_pos_best;  };

        void reset()                        { };
        void set_size(Integer, Integer)     { };
        void current_vector(Integer vec)    { col_pos = vec+1;      };
        const value_type& value()           { return state;         };

        bool add(Integer r, const in_type& val)
        { 
            if (row_pos == 0 || less(val,state))
            {
                state           = val;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;
        };

        bool add_zero(Integer r)
        { 
            if (row_pos == 0 || less(Z,state))
            {
                state           = Z;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;        
        };
};

template<class in_type, class value_type>
struct max_accumulator2_vec
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type              state;
        Integer                 col_pos;
        Integer                 col_pos_best;
        Integer                 row_pos;
        const value_type        Z;

        max_accumulator2_vec& operator=(const max_accumulator2_vec&) = delete;

    public:
        max_accumulator2_vec(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,row_pos(0), col_pos(0), col_pos_best(0)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            return ti;
        };

        Integer index_row() const           { return row_pos;       };
        Integer index_col() const           { return col_pos_best;  };

        void reset()                        { };
        void set_size(Integer, Integer)     { };
        void current_vector(Integer vec)    { col_pos = vec+1;      };
        const value_type& value()           { return state;         };

        bool add(Integer r, const in_type& val)
        { 
            if (row_pos == 0 || greater(val,state))
            {
                state           = val;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;
        };

        bool add_zero(Integer r)
        { 
            if (row_pos == 0 || greater(Z,state))
            {
                state           = Z;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;        
        };
};

template<class in_type, class value_type>
struct min_abs_accumulator2
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        md::workspace2<value_type>  st_array;
        integer_dense               ind_array;
        Integer                     vec_pos;
        const value_type            Z;

        min_abs_accumulator2& operator=(const min_abs_accumulator2&) = delete;

    public:
        min_abs_accumulator2(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,ind_array(ti::ti_empty())
            ,st_array(ti)
        {};
        
        static ret_ti get_ret_ti(in_ti ti)
        {
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        void set_size(Integer , Integer c, Integer dim)    
        { 
            Integer mr = 1, mc = 1;
            if (dim == 1)    mc = c;
            else            mr = c;

            ind_array.reset_unique(mr,mc);            
            Integer* ptr_ind = ind_array.ptr();

            for (Integer i = 0; i < c; ++i)
                ptr_ind[i] = 0;
        };

        void reset()                                { };        
        void current_vector(Integer vec)            { vec_pos = vec;        };
        const value_type& value()                   { return state;         };
        const value_type& value(Integer p)          { return st_array[p];   };        
        const integer_dense& index_matrix() const   { return ind_array;     };

        void reset_array(Integer s)
        { 
            st_array.resize(s);    
        };    

        bool add(Integer v, const in_type& val)
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);
            if (ind_array.ptr()[vec_pos] == 0 || less(aval,state))
            {
                state = aval;
                ind_array.ptr()[vec_pos] = v+1;
            }
            return false;
        };

        bool add_zero(Integer v)
        { 
            if (ind_array.ptr()[vec_pos] == 0 || less(Z,state))
            {
                state = Z;
                ind_array.ptr()[vec_pos] = v+1;
            }
            else if (mrd::is_zero(state))
            {
                ind_array.ptr()[vec_pos] = std::min(v+1,ind_array.ptr()[vec_pos]);
            };
            return true;    
        };

        void add(Integer p,Integer v, const in_type& val)        
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);

            if (ind_array.ptr()[p] == 0 || less(aval,st_array[p]))
            {
                st_array[p] = aval;        
                ind_array.ptr()[p] = v+1;
            };
        };

        void add_zero(Integer p,Integer v)            
        { 
            if (ind_array.ptr()[p] == 0 || less(Z,st_array[p]))
            {
                st_array[p] = Z;
                ind_array.ptr()[p] = v+1;
            }
            else if(mrd::is_zero(st_array[p]))
            {
                ind_array.ptr()[p] = std::min(v+1,ind_array.ptr()[p]);
            };
        }        
};

template<class in_type, class value_type>
struct min_abs_accumulator2_vec
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        Integer                     col_pos;
        Integer                     row_pos;
        Integer                     col_pos_best;
        const value_type            Z;

        min_abs_accumulator2_vec& operator=(const min_abs_accumulator2_vec&) = delete;

    public:
        min_abs_accumulator2_vec(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,row_pos(0), col_pos(0), col_pos_best(0)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        Integer index_row() const           { return row_pos;       };
        Integer index_col() const           { return col_pos_best;  };

        void reset()                        { };
        void set_size(Integer, Integer)     { };
        void current_vector(Integer vec)    { col_pos = vec+1;      };
        const value_type& value()           { return state;         };

        bool add(Integer r, const in_type& val)
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);

            if (row_pos == 0 || less(aval,state))
            {
                state           = aval;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;
        };

        bool add_zero(Integer r)
        { 
            if (row_pos == 0 || less(Z,state))
            {
                state           = Z;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;        
        };
};

template<class in_type, class value_type>
struct max_abs_accumulator2_vec
{
    private:
        using in_ti     = ti::ti_type<in_type>;
        using ret_ti    = ti::ti_type<value_type>;

    private:
        value_type                  state;
        Integer                     col_pos;
        Integer                     row_pos;
        Integer                     col_pos_best;
        const value_type            Z;

        max_abs_accumulator2_vec&   operator=(const max_abs_accumulator2_vec&) = delete;

    public:
        max_abs_accumulator2_vec(ti::ti_type<value_type> ti)        
            : Z(md::default_value<value_type>(ti)) 
            ,state(md::default_value<value_type>(ti)) 
            ,row_pos(0), col_pos(0), col_pos_best(0)
        {};

        static ret_ti get_ret_ti(in_ti ti)
        {
            ret_ti ret = ti::get_return_ti<ret_ti>(mdyf::abs::eval(), ti);
            return ret;
        };

        Integer index_row() const           { return row_pos;       };
        Integer index_col() const           { return col_pos_best;  };

        void reset()                        { };
        void set_size(Integer, Integer)     { };
        void current_vector(Integer vec)    { col_pos = vec+1;      };
        const value_type& value()           { return state;         };

        bool add(Integer r, const in_type& val)
        { 
            value_type aval = mrd::abs_helper<in_type>::eval(val);

            if (row_pos == 0 || greater(aval,state))
            {
                state           = aval;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;
        };

        bool add_zero(Integer r)
        { 
            if (row_pos == 0 || greater(Z,state))
            {
                state           = Z;
                row_pos         = r+1;
                col_pos_best    = col_pos;
            }
            return false;        
        };
};

//-----------------------------------------------------------
//              vec_manip_helper
//-----------------------------------------------------------
template<class M>
void vec_manip_helper<M>::eval_sum(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_sum::value_type;
    eval_vec_functor<ret_type_sum,M, sum_accumulator<in_type,out_type>>::eval(ret, m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_nnz(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_nnz::value_type;

    eval_vec_functor<ret_type_nnz,M, nnz_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_prod(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_prod::value_type;
    eval_vec_functor<ret_type_prod, M, prod_accumulator<in_type,out_type>>::eval(ret, m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_cumsum(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_cumsum::value_type;
    eval_vec_functor<ret_type_cumsum,M, sum_accumulator<in_type,out_type>>::eval_cum(ret, m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_cumprod(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_cumprod::value_type;
    eval_vec_functor<ret_type_cumprod,M, prod_accumulator<in_type,out_type>>::eval_cum(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_sumsq(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_sumsq::value_type;
    eval_vec_functor<ret_type_sumsq,M, sumsq_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_min(matcl::Matrix& ret, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);

        if (dim == 1)
            ret = matcl::Matrix(ret_type_min(ret_ti,std::min(r,1),c), false);
        else
            ret = matcl::Matrix(ret_type_min(ret_ti,r,std::min(c,1)),false);

        return;
    };

    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_min::value_type;
    eval_vec_functor<ret_type_min,M, min_accumulator<in_type,out_type>>::eval(ret, m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_min_abs(matcl::Matrix& ret, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);
        if (dim == 1)
            ret = matcl::Matrix(ret_type_min_abs(ret_ti,std::min(r,1),c),false);
        else
            ret = matcl::Matrix(ret_type_min_abs(ret_ti,r,std::min(c,1)), false);

        return;
    };

    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_min_abs::value_type;
    eval_vec_functor<ret_type_min_abs,M, min_abs_accumulator<in_type,out_type>>
                ::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_max(matcl::Matrix& ret, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);
        if (dim == 1)
            ret = matcl::Matrix(ret_type_max(ret_ti,std::min(r,1),c),false);
        else
            ret = matcl::Matrix(ret_type_max(ret_ti,r,std::min(c,1)), false);

        return;
    };

    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_max::value_type;
    eval_vec_functor<ret_type_max,M, max_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_max_abs(matcl::Matrix& ret, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);
        if (dim == 1)
            ret = matcl::Matrix(ret_type_max_abs(ret_ti,std::min(r,1),c),false);
        else
            ret = matcl::Matrix(ret_type_max_abs(ret_ti,r,std::min(c,1)), false);

        return;
    };

    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_max_abs::value_type;
    eval_vec_functor<ret_type_max_abs,M, max_abs_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_min2(matcl::Matrix& rx, matcl::Matrix& ri, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);

        if (dim == 1)
        {
            integer_dense ind(ti::ti_empty(),std::min(r,1),c);
            ret_type_min x(ret_ti,std::min(r,1),c);

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            return;
        }
        else
        {
            integer_dense ind(ti::ti_empty(),r,std::min(c,1));
            ret_type_min x(ret_ti,r,std::min(c,1));

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            return;
        };
    };

    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_min::value_type;
    using accum_type    = min_accumulator2<in_type,out_type>;

    accum_type accum(ret_ti);
    mr::eval_vec_functor2<ret_type_min,M, accum_type>::eval(rx,m,dim,accum);
    const integer_dense& ind = accum.index_matrix();

    ri  = matcl::Matrix(ind,true);
    return;
};

template<class M>
void vec_manip_helper<M>::eval_min_abs2(matcl::Matrix& rx, matcl::Matrix& ri, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    using ret_ti_type   = typename ti::get_ti_type<M>::type;
    ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::abs2::eval(), m.get_type());

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);

        if (dim == 1)
        {
            integer_dense ind(ti::ti_empty(),std::min(r,1),c);
            ret_type_min_abs x(ret_ti,std::min(r,1),c);

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        }
        else
        {
            integer_dense ind(ti::ti_empty(),r,std::min(c,1));
            ret_type_min_abs x(ret_ti,r,std::min(c,1));

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        };
    };

    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_min_abs::value_type;
    using accum_type    = min_abs_accumulator2<in_type,out_type>;

    accum_type accum(ret_ti);
    mr::eval_vec_functor2<ret_type_min_abs,M, accum_type>::eval(rx, m,dim,accum);

    const integer_dense& ind = accum.index_matrix();
    ri = matcl::Matrix(ind,false);
    return;
};

template<class M>
void vec_manip_helper<M>::eval_max2(matcl::Matrix& rx, matcl::Matrix& ri, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);
        
        if (dim == 1)
        {
            integer_dense ind(ti::ti_empty(),std::min(r,1),c);
            ret_type_max x(ret_ti,std::min(r,1),c);

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        }
        else
        {
            integer_dense ind(ti::ti_empty(),r,std::min(c,1));
            ret_type_max x(ret_ti,r,std::min(c,1));

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        };
    };

    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_max::value_type;
    using accum_type    = max_accumulator2<in_type,out_type>;

    accum_type accum(ret_ti);

    mr::eval_vec_functor2<ret_type_max,M, accum_type>::eval(rx, m,dim,accum);
    const integer_dense& ind = accum.index_matrix();

    ri  = matcl::Matrix(ind,true);
    return;
};

template<class M>
void vec_manip_helper<M>::eval_max_abs2(matcl::Matrix& rx, matcl::Matrix& ri, const M& m, int dim)
{
    Integer r = m.rows();
    Integer c = m.cols();

    using ret_ti_type   = typename ti::get_ti_type<M>::type;
    ret_ti_type ret_ti  = ti::get_return_ti<ret_ti_type>(mdyf::abs2::eval(), m.get_type());

    if (r == 0 || c == 0)
    {
        error::check_dim(dim);
        if (dim == 1)
        {
            integer_dense ind(ti::ti_empty(),std::min(r,1),c);
            ret_type_max_abs x(ret_ti,std::min(r,1),c);

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        }
        else
        {
            integer_dense ind(ti::ti_empty(),r,std::min(c,1));
            ret_type_max_abs x(ret_ti,r,std::min(c,1));

            rx  = matcl::Matrix(x,false);
            ri  = matcl::Matrix(ind,false);
            
            return;
        };
    };

    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_max_abs::value_type;
    using accum_type    = max_abs_accumulator2<in_type,out_type>;

    accum_type accum(ret_ti);

    mr::eval_vec_functor2<ret_type_max_abs,M, accum_type>::eval(rx, m,dim,accum);
    const integer_dense& ind = accum.index_matrix();

    ri = matcl::Matrix(ind,true);
    return;
};

template<class M>
void vec_manip_helper<M>::eval_mean(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_mean::value_type;
    eval_vec_functor<ret_type_mean,M, mean_accumulator<in_type,out_type>>::eval(ret, m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_std(matcl::Matrix& ret, const M& m, int dim, bool unbiased)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_std::value_type;

    using ret_ti_type   = typename ti::get_ti_type<ret_type_std>::type;
    using accum_type    = std_accumulator<in_type,out_type>;

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));

    accum_type accum(ret_ti,unbiased);
    eval_vec_functor<ret_type_std,M, accum_type>::eval(ret, m,dim,accum);
};

template<class M>
void vec_manip_helper<M>::eval_any(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_any::value_type;
    eval_vec_functor<ret_type_any,M, any_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_any(matcl::Matrix& ret, const M& m, const test_function& t, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_any::value_type;
    using accum     = any_accumulator_t<in_type,out_type>;

    accum a(m.get_type(),t);
    eval_vec_functor<ret_type_any,M,accum>::eval(ret,m,dim,a);
};

template<class M>
void vec_manip_helper<M>::eval_all(matcl::Matrix& ret, const M& m, int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_all::value_type;
    eval_vec_functor<ret_type_all,M, all_accumulator<in_type,out_type>>::eval(ret,m,dim);
};

template<class M>
void vec_manip_helper<M>::eval_all(matcl::Matrix& ret, const M& m, const test_function& t,int dim)
{
    using in_type   = typename M::value_type;
    using out_type  = typename ret_type_all::value_type;
    using accum     = all_accumulator_t<in_type,out_type>;

    accum a(m.get_type(),t);
    eval_vec_functor<ret_type_all,M, accum>::eval(ret,m,dim,a);
};

template<class M>
void vec_manip_helper<M>::eval_nnz_vec(Integer& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = Integer;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret = 0;
        return;
    };

    eval_vec_functor_vec<Integer, M, nnz_accumulator<in_type,out_type>>::eval(ret,m);
};

template<class M>
void vec_manip_helper<M>::eval_sum_vec(matcl::Matrix& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = in_type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();
    out_type tmp_ret = md::default_value<out_type>(ret_ti);
    eval_vec_functor_vec<out_type, M, sum_accumulator<in_type,out_type>>::eval(tmp_ret, m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_prod_vec(matcl::Matrix& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = in_type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        ret = md::one_value<out_type>(ret_ti);
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();
    out_type tmp_ret = md::default_value<out_type>(ret_ti);
    eval_vec_functor_vec<out_type, M, prod_accumulator<in_type,out_type>>::eval(tmp_ret, m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_sumsq_vec(matcl::Matrix& ret, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = in_type;
    using accum_type    = sumsq_accumulator<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        ret = md::default_value<out_type>(ret_ti);
        return;
    };
    
    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, accum_type>::eval(tmp_ret,m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_min_vec(matcl::Matrix& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = in_type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();
    out_type tmp_ret = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, min_accumulator<in_type,out_type>>::eval(tmp_ret, m);
    ret = tmp_ret;
}

template<class M>
void vec_manip_helper<M>::eval_max_vec(matcl::Matrix& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = in_type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();
    out_type tmp_ret = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, max_accumulator<in_type,out_type>>::eval(tmp_ret,m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_min_abs_vec(matcl::Matrix& ret, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_min_abs::value_type;
    using accum_type    = min_abs_accumulator<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, accum_type>::eval(tmp_ret,m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_max_abs_vec(matcl::Matrix& ret, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_max_abs::value_type;
    using accum_type    = max_abs_accumulator<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, accum_type>::eval(tmp_ret,m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_mean_vec(matcl::Matrix& ret, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_mean::value_type;
    using accum_type    = mean_accumulator<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    eval_vec_functor_vec<out_type,M, accum_type>::eval(tmp_ret,m);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_std_vec(matcl::Matrix& ret, const M& m, bool unbiased)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_std::value_type;
    using ret_ti_type   = typename ti::get_ti_type<ret_type_std>::type;
    using accum_type    = std_accumulator<in_type,out_type>;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        ret = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    accum_type accum(ret_ti,unbiased);    
    eval_vec_functor_vec<out_type,M, accum_type>::eval(tmp_ret, m, accum);
    ret = tmp_ret;
};

template<class M>
void vec_manip_helper<M>::eval_any_vec(bool& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = Integer;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret = false;
        return;
    };

    out_type tmp_ret;
    eval_vec_functor_vec<out_type,M, any_accumulator<in_type,out_type>>::eval(tmp_ret,m);
    ret = (tmp_ret == 0) ? false : true;
}

template<class M>
void vec_manip_helper<M>::eval_all_vec(bool& ret, const M& m)
{
    using in_type   = typename M::value_type;
    using out_type  = Integer;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret = true;
        return;
    };

    out_type tmp_ret;
    eval_vec_functor_vec<out_type,M, all_accumulator<in_type,out_type>>::eval(tmp_ret,m);
    ret = (tmp_ret == 0) ? false : true;
};

template<class M>
void vec_manip_helper<M>::eval_any_vec(bool& ret, const M& m, const test_function& t)
{
    using in_type   = typename M::value_type;
    using out_type  = Integer;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret = false;
        return;
    };

    using accum     = any_accumulator_t<in_type,out_type>;

    accum a(m.get_type(),t);

    out_type tmp_ret;
    eval_vec_functor_vec<out_type,M,accum>::eval(tmp_ret,m,a);
    ret = (tmp_ret == 0) ? false : true;
};

template<class M>
void vec_manip_helper<M>::eval_all_vec(bool& ret, const M& m, const test_function& t)
{
    using in_type   = typename M::value_type;
    using out_type  = Integer;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret = true;
        return;
    };

    using accum     = all_accumulator_t<in_type,out_type>;

    accum a(m.get_type(),t);

    out_type tmp_ret;
    eval_vec_functor_vec<out_type, M, accum>::eval(tmp_ret,m,a);
    ret = (tmp_ret == 0) ? false : true;
};

template<class M>
void vec_manip_helper<M>::eval_min2_vec(matcl::Matrix& x, Integer& i, Integer& j, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = in_type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        x = md::default_value<out_type>(ret_ti);
        i   = 0;
        j   = 0;
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();
    using accum_type    = min_accumulator2_vec<in_type,out_type>;

    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    accum_type accum(ret_ti);
    mr::eval_vec_functor_vec2<out_type,M, accum_type>::eval(tmp_ret,m,accum);
    x   = tmp_ret;
    i   = accum.index_row();
    j   = accum.index_col();
    return;
};

template<class M>
void vec_manip_helper<M>::eval_max2_vec(matcl::Matrix& x, Integer& i, Integer& j, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = in_type;
    using accum_type    = max_accumulator2_vec<in_type,out_type>;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        typename ti::get_ti_type<M>::type ret_ti = m.get_type();
        x   = md::default_value<out_type>(ret_ti);
        i   = 0;
        j   = 0;
        return;
    };

    typename ti::get_ti_type<M>::type ret_ti = m.get_type();

    out_type tmp_ret    = md::default_value<out_type>(ret_ti);

    accum_type accum(ret_ti);
    mr::eval_vec_functor_vec2<out_type,M, accum_type>::eval(tmp_ret,m,accum);
    x   = tmp_ret;
    i   = accum.index_row();
    j   = accum.index_col();
    return;
};

template<class M>
void vec_manip_helper<M>::eval_min_abs2_vec(matcl::Matrix& x, Integer& i, Integer& j, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_min_abs::value_type;
    using accum_type    = min_abs_accumulator2_vec<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        x = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);
    accum_type accum(ret_ti);

    mr::eval_vec_functor_vec2<out_type,M, accum_type>::eval(tmp_ret,m,accum);
    x   = tmp_ret;
    i   = accum.index_row();
    j   = accum.index_col();
    return;    
};

template<class M>
void vec_manip_helper<M>::eval_max_abs2_vec(matcl::Matrix& x, Integer& i, Integer& j, const M& m)
{
    using in_type       = typename M::value_type;
    using out_type      = typename ret_type_max_abs::value_type;
    using accum_type    = max_abs_accumulator2_vec<in_type,out_type>;
    using ret_ti_type   = typename ti::get_ti_type<out_type>::type;

    Integer r = m.rows();
    Integer c = m.cols();

    if (r == 0 || c == 0)
    {
        ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
        x = md::default_value<out_type>(ret_ti);
        return;
    };

    ret_ti_type ret_ti  = accum_type::get_ret_ti(ti::get_ti(m));
    out_type tmp_ret    = md::default_value<out_type>(ret_ti);
    accum_type accum(ret_ti);

    mr::eval_vec_functor_vec2<out_type,M, accum_type>::eval(tmp_ret,m,accum);
    x   = tmp_ret;
    i   = accum.index_row();
    j   = accum.index_col();
    return;       
};

}}};

MACRO_INSTANTIATE_G_1(matcl::raw::details::vec_manip_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::vec_manip_helper)
