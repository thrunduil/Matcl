#pragma once

#include "mkgen/TODO/utils/utils.h"

namespace matcl { namespace mkgen
{

template<Integer Pos, class Colon>
struct get_pos_colon
{
    static_assert(details::dependent_false<Colon>::value,
                  "this type should not be instantiated");
};
template<Integer Pos>
struct get_pos_colon<Pos,colon_all>
{
    static const Integer value = Pos;
};
template<Integer Pos, Integer Start, Integer Step, Integer End>
struct get_pos_colon<Pos, colon3<Start,Step,End>>
{
    static const Integer value = Start + Step * (Pos - 1);
};
template<Integer Pos, Integer Start, Integer End>
struct get_pos_colon<Pos, colon2<Start, End>>
{
    static const Integer value = Start + (Pos - 1);
};
template<Integer Pos, Integer Start>
struct get_pos_colon<Pos, colon<Start>>
{
    static_assert(Pos == 1, "invalid element");
    static const Integer value = Start;
};

template<class Colon, Integer M>
struct get_size_colon
{
    static_assert(details::dependent_false<Colon>::value, "unknown colon type");
};
template<class Colon, Integer M>
struct get_first_colon
{
    static_assert(details::dependent_false<Colon>::value, "unknown colon type");
};
template<class Colon, Integer M>
struct get_last_colon
{
    static_assert(details::dependent_false<Colon>::value, "unknown colon type");
};
template<class Colon>
struct get_offset_colon
{
    static_assert(details::dependent_false<Colon>::value, "unknown colon type");
};
template<class Colon>
struct get_step_colon
{
    static_assert(details::dependent_false<Colon>::value, "unknown colon type");
};

template<Integer M>
struct get_size_colon<colon_all, M>
{
    static const Integer value  = M;
};
template<Integer M>
struct get_first_colon<colon_all, M>
{
    static const Integer value  = 1;
};
template<Integer M>
struct get_last_colon<colon_all, M>
{
    static const Integer value  = M;
};
template<>
struct get_offset_colon<colon_all>
{
    static const Integer value  = 0;
};
template<>
struct get_step_colon<colon_all>
{
    static const Integer value  = 1;
};

template<Integer M, Integer M1>
struct get_size_colon<colon<M1>, M>
{
    static const Integer value  = 1;
};
template<Integer M, Integer M1>
struct get_first_colon<colon<M1>, M>
{
    static const Integer value  = M1;
};
template<Integer M, Integer M1>
struct get_last_colon<colon<M1>, M>
{
    static const Integer value  = M1;
};
template<Integer M1>
struct get_offset_colon<colon<M1>>
{
    static const Integer value  = M1 - 1;
};
template<Integer M1>
struct get_step_colon<colon<M1>>
{
    static const Integer value  = 1;
};

template<Integer M, Integer Start, Integer End>
struct get_size_colon<colon2<Start,End>, M>
{
    static const Integer value0 = End - Start + 1;
    static const Integer value  = (value0 < 0 ? 0 : value0);
};
template<Integer M, Integer Start, Integer End>
struct get_first_colon<colon2<Start,End>, M>
{
    static const Integer value  = Start;
};
template<Integer M, Integer Start, Integer End>
struct get_last_colon<colon2<Start,End>, M>
{
    static const Integer value  = End;
};
template<Integer Start, Integer End>
struct get_offset_colon<colon2<Start,End>>
{
    static const Integer value  = Start - 1;
};
template<Integer Start, Integer End>
struct get_step_colon<colon2<Start,End>>
{
    static const Integer value  = 1;
};

template<Integer M, Integer Start, Integer Step, Integer End>
struct get_size_colon<colon3<Start, Step, End>, M>
{
    static const Integer d      = (End - Start) / Step;
    static const Integer value  = d + 1;
};
template<Integer M, Integer Start, Integer Step, Integer End>
struct get_first_colon<colon3<Start, Step, End>, M>
{
    static const Integer value  = Start;
};
template<Integer M, Integer Start, Integer Step, Integer End>
struct get_last_colon<colon3<Start, Step, End>, M>
{
    static const Integer size   = get_size_colon<colon3<Start, Step, End>,M>::value;
    static const Integer value  = Start + (size-1) * Step;
};
template<Integer Start, Integer Step, Integer End>
struct get_offset_colon<colon3<Start, Step, End>>
{
    static const Integer value  = Start - 1;
};
template<Integer Start, Integer Step, Integer End>
struct get_step_colon<colon3<Start, Step, End>>
{
    static const Integer value  = Step;
};

}}
