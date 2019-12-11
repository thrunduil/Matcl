#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

struct align_none{};
struct align_half{};
struct align_full{};

template<class Align_1, class Align_2>
struct link_alignment;

template<> struct link_alignment<align_none, align_none> { using type = align_none; };
template<> struct link_alignment<align_none, align_half> { using type = align_none; };
template<> struct link_alignment<align_none, align_full> { using type = align_none; };
template<> struct link_alignment<align_half, align_none> { using type = align_none; };
template<> struct link_alignment<align_half, align_half> { using type = align_half; };
template<> struct link_alignment<align_half, align_full> { using type = align_half; };
template<> struct link_alignment<align_full, align_none> { using type = align_none; };
template<> struct link_alignment<align_full, align_half> { using type = align_half; };
template<> struct link_alignment<align_full, align_full> { using type = align_full; };

template<bool Full, bool Half>
struct make_offset_align_type;

template<bool Half>
struct make_offset_align_type<true,Half>
{
    using type = align_full;
};
template<>
struct make_offset_align_type<false,true>
{
    using type = align_half;
};
template<>
struct make_offset_align_type<false,false>
{
    using type = align_none;
};

template<class Val, Integer Offset, Integer Step, Integer Alignment>
struct check_offset_aligned
{
    static const Integer size   = sizeof(Val);
    static const Integer off    = Alignment / size - 1;    
    static const bool value = (Step > 0 &&  (size * Offset ) % Alignment == 0)
                            || (Step < 0 && (size * (Offset-off) ) % Alignment == 0)
                            || (Step == 0);
};
template<class Val, Integer Offset, Integer Step>
struct get_offset_alignment
{
    static const bool value1    = check_offset_aligned<Val,Offset,Step,VEC_ALIGN>::value;
    static const bool value2    = check_offset_aligned<Val,Offset,Step,VEC_ALIGN/2>::value;
    using type                  = typename make_offset_align_type<value1,value2>::type;
};

template<class Align>   struct is_align_full                { static const bool value = false; };
template<>              struct is_align_full<align_full>    { static const bool value = true; };

template<class Align>   struct is_align_half                { static const bool value = false; };
template<>              struct is_align_half<align_full>    { static const bool value = true; };
template<>              struct is_align_half<align_half>    { static const bool value = true; };

template<class Align_Type>
struct align_type_size;

template<>              struct align_type_size<align_none>  { static const Integer value = 1; };
template<>              struct align_type_size<align_half>  { static const Integer value = 2; };
template<>              struct align_type_size<align_full>  { static const Integer value = 4; };


template<class Align_Type, Integer Vec_Size>
struct is_aligned
{ 
    static const bool value = (Vec_Size <= align_type_size<Align_Type>::value); 
};

}}