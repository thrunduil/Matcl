/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-dynamic/register_function.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-dynamic/predefined_functions.h"
#include "special_functions.h"
#include "matcl-dynamic/details/object.inl"

namespace matcl { namespace dynamic
{

// we consider as standard conversions:
//      - floating point promotions (float -> real, Float_complex -> Complex)
//      - integer -> real promotion
//      - equivalent conversions I (float -> Float_complex, Real->Complex)
//      - equivalent conversions II (Integer->Complex, float->Complex)
//      - derived-to-base conversions. 
//      - decay conversions I (a conversion with precision lost, 
//          Real->Float, Complex->Float_complex)
//      - decay conversions II (Real->Float_complex)
//      - int-float conversion (Integer -> Float, Integer -> Float_complex)
//      - an explicit cast (floating point values to integer, complex values
//          to floating point values and integers)

struct numeric_standard_promotion_1   
{ 
    static function_name eval() {return special_functions::convert_numeric_promotion_1();}; 
};
struct numeric_standard_promotion_2
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_promotion_2();
    }; 
};
struct numeric_convert_equivalent_1
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_equivalent_1();
    }; 
};
struct numeric_convert_equivalent_2
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_equivalent_2();
    }; 
};

struct numeric_convert_decay_1
{ 
    static function_name eval() 
    {
        return special_functions::convert_numeric_decay_1();
    }; 
};
struct numeric_convert_decay_2
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_decay_2();
    }; 
};
struct numeric_convert_int_float
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_int_float();
    }; 
};

struct numeric_convert_explicit
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_explicit();
    }; 
};

struct numeric_convert_cast
{ 
    static function_name eval()
    {
        return special_functions::convert_numeric_cast();
    }; 
};
struct any_standard_convert
{ 
    static function_name eval()
    {
        return special_functions::convert_any();
    }; 
};

struct unit_convert
{ 
    static function_name eval()
    {
        return special_functions::convert_unit();
    }; 
};

struct convert_id
{ 
    static function_name eval()
    {
        return special_functions::convert_id();
    }; 
};

//----------------------------------------------------------------------------
//                  from Integer
//----------------------------------------------------------------------------
struct conv_int_float : register_convert<conv_int_float, numeric_convert_int_float>
{
	static OFloat eval(const OInteger& val) { return OFloat(Float(val.get())); };
};

struct conv_int_real : register_convert<conv_int_real, numeric_standard_promotion_2>
{
	static OReal eval(const OInteger& val)  { return OReal(Real(val.get())); };
};

struct conv_int_compl : register_convert<conv_int_compl, numeric_convert_equivalent_2>
{
	static OComplex eval(const OInteger& v) { return OComplex(Complex(v.get())); };
};

struct conv_int_fcompl : register_convert<conv_int_fcompl, numeric_convert_int_float>
{
	static OFloat_complex eval(const OInteger& val)
    { 
        return OFloat_complex(Float_complex(Float(val.get()))); 
    };
};

//----------------------------------------------------------------------------
//                  from Float
//----------------------------------------------------------------------------
struct conv_float_int : register_convert<conv_float_int, numeric_convert_cast>
{
	static OInteger eval(const OFloat& val)
    { 
        return OInteger(Integer(val.get())); 
    };
};

struct conv_float_real : register_convert<conv_float_real, numeric_standard_promotion_1>
{
	static OReal eval(const OFloat& val)
    { 
        return OReal(Real(val.get())); 
    };
};

struct conv_float_fcompl : register_convert<conv_float_fcompl, numeric_convert_equivalent_1>
{
	static OFloat_complex eval(const OFloat& val)
    { 
        return OFloat_complex(Float_complex(val.get())); 
    };
};

struct conv_float_compl : register_convert<conv_float_compl, numeric_convert_equivalent_2>
{
	static OComplex eval(const OFloat& val)
    { 
        return OComplex(Complex(val.get())); 
    };
};

//----------------------------------------------------------------------------
//                  from Real
//----------------------------------------------------------------------------
struct conv_real_int : register_convert<conv_real_int, numeric_convert_cast>
{
	static OInteger eval(const OReal& val)
    { 
        return OInteger(Integer(val.get())); 
    };
};

struct conv_real_float : register_convert<conv_real_float, numeric_convert_decay_1>
{
	static OFloat eval(const OReal& val) 
    { 
        return OFloat(Float(val.get())); 
    };
};

struct conv_real_compl : register_convert<conv_real_compl, numeric_convert_equivalent_1>
{
	static OComplex eval(const OReal& val)
    { 
        return OComplex(Complex(val.get())); 
    };
};

struct conv_real_fcompl : register_convert<conv_real_fcompl, numeric_convert_decay_2>
{
	static OFloat_complex eval(const OReal& val) 
    { 
        return OFloat_complex(Float_complex(Float(val.get()))); 
    };
};

//----------------------------------------------------------------------------
//                  from Complex
//----------------------------------------------------------------------------
struct conv_compl_int : register_convert<conv_compl_int, numeric_convert_cast>
{
	static OInteger eval(const OComplex& val)    
    { 
        return OInteger(Integer(real(val.get()))); 
    };
};

struct conv_compl_float : register_convert<conv_compl_float, numeric_convert_cast>
{
	static OFloat eval(const OComplex& val)
    { 
        return OFloat(Float(real(val.get()))); 
    };
};

struct conv_compl_real : register_convert<conv_compl_real, numeric_convert_cast>
{
	static OReal eval(const OComplex& val)
    { 
        return OReal(Real(real(val.get())));
    };
};

struct conv_compl_fcompl : register_convert<conv_compl_fcompl, numeric_convert_decay_1>
{
	static OFloat_complex eval(const OComplex& val)
    { 
        return OFloat_complex(Float_complex(val.get())); 
    };
};

//----------------------------------------------------------------------------
//                  from Float_complex
//----------------------------------------------------------------------------
struct conv_fcompl_int : register_convert<conv_fcompl_int, numeric_convert_cast>
{
	static OInteger eval(const OFloat_complex& val)
    { 
        return OInteger(Integer(real(val.get()))); 
    };
};

struct conv_fcompl_float : register_convert<conv_fcompl_float, numeric_convert_cast>
{
	static OFloat eval(const OFloat_complex& val)
    { 
        return OFloat(Float(real(val.get()))); 
    };
};

struct conv_fcompl_real : register_convert<conv_fcompl_real, numeric_convert_cast>
{
	static OReal eval(const OFloat_complex& val)
    { 
        return OReal(Real(real(val.get()))); 
    };
};

struct conv_fcompl_compl : register_convert<conv_fcompl_compl, numeric_standard_promotion_1>
{
	static OComplex eval(const OFloat_complex& val)
    { 
        return OComplex(Complex(val.get())); 
    };
};

//----------------------------------------------------------------------------
//                  to Any
//----------------------------------------------------------------------------
struct conv_to_any : register_convert<conv_to_any, any_standard_convert>
{
	static Any eval(const object& val)  { return Any(any_type(val)); };
};

//----------------------------------------------------------------------------
//                  to Unit
//----------------------------------------------------------------------------
struct conv_to_unit : register_convert<conv_to_unit, unit_convert>
{
	static Unit eval(const object& val)  { return Unit(unit_type(val)); };
};

//----------------------------------------------------------------------------
//                  from Bool
//----------------------------------------------------------------------------

struct conv_bool_int : register_convert<conv_bool_int, convert_explicit>
{
	static OInteger eval(const OBool& val) { return OInteger(Integer(val.get())); };
};

//----------------------------------------------------------------------------
//                  identity
//----------------------------------------------------------------------------
struct conv_id : register_convert<conv_id, convert_id>
{
	static object eval(const object& val)  { return val; };
};

//----------------------------------------------------------------------------
//                  unifiers
//----------------------------------------------------------------------------
struct unif_real : register_unifier<unif_real>
{
    static OReal eval() { return OReal();};
};

struct unif_compl : register_unifier<unif_compl>
{
    static OComplex eval() {return OComplex();}
};

};};