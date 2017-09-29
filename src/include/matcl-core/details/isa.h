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

#pragma once

#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/details/mpl.h"
#include "matcl-core/matrix/enums.h"

namespace matcl { namespace details
{

template<class T>	
            struct is_scalar					: is_object<T>{};
template<>  struct is_scalar<bool>				{	static const bool value = false; };
template<>  struct is_scalar<unsigned char>		{	static const bool value = true;	 };
template<>  struct is_scalar<unsigned short>	{	static const bool value = true;	 };
template<>  struct is_scalar<unsigned int>		{	static const bool value = true;	 };
template<>  struct is_scalar<unsigned long>		{	static const bool value = true;	 };
template<>  struct is_scalar<signed char>		{	static const bool value = true;	 };
template<>  struct is_scalar<signed short>		{	static const bool value = true;	 };
template<>  struct is_scalar<signed int>		{	static const bool value = true;	 };
template<>  struct is_scalar<signed long>		{	static const bool value = true;	 };
template<>  struct is_scalar<float>				{	static const bool value = true;	 };
template<>  struct is_scalar<double>			{	static const bool value = true;	 };
template<>  struct is_scalar<long double>		{	static const bool value = true;	 };
template<>  struct is_scalar<Complex>           {	static const bool value = true;	 };
template<>  struct is_scalar<Float_complex>     {	static const bool value = true;	 };
template<>  struct is_scalar<Object>            {	static const bool value = true;	 };

template<class T>	
            struct is_matcl_scalar              {	static const bool value = false; };
template<>	struct is_matcl_scalar<Integer>     {   static const bool value = true;	 };
template<>	struct is_matcl_scalar<Float>       {   static const bool value = true;	 };
template<>	struct is_matcl_scalar<Real>        {   static const bool value = true;	 };
template<>	struct is_matcl_scalar<Complex>     {   static const bool value = true;	 };
template<>	struct is_matcl_scalar<Float_complex>{  static const bool value = true;	 };
template<>	struct is_matcl_scalar<Object>      {   static const bool value = true;	 };

template<class T> 
struct is_any_scalar
{	
    static const bool value = (is_external_scalar<T>::value || is_scalar<T>::value)
                            && is_object<T>::value == false; 
};

template<class T>
            struct is_float_scalar              {	static const bool value = false; };
template<>	struct is_float_scalar<float>		{	static const bool value = true;	 };
template<>	struct is_float_scalar<double>		{	static const bool value = true;	 };
template<>	struct is_float_scalar<long double>	{	static const bool value = true;	 };
template<>	struct is_float_scalar<Complex>     {	static const bool value = true;	 };
template<>	struct is_float_scalar<Float_complex>{	static const bool value = true;	 };

template<class T>
            struct is_float_real_scalar                 {	static const bool value = false; };
template<>	struct is_float_real_scalar<float>          {	static const bool value = true;	 };
template<>	struct is_float_real_scalar<double>         {	static const bool value = true;	 };
template<>	struct is_float_real_scalar<long double>    {	static const bool value = true;	 };

template<class T>
            struct is_real_scalar                       {	static const bool value = false; };
template<>	struct is_real_scalar<Integer>              {	static const bool value = true;	 };
template<>	struct is_real_scalar<float>                {	static const bool value = true;	 };
template<>	struct is_real_scalar<double>               {	static const bool value = true;	 };
template<>	struct is_real_scalar<long double>          {	static const bool value = true;	 };

template<class T>
            struct is_single_precision				    {	static const bool value = false; };
template<>	struct is_single_precision<Integer>         {   static const bool value = false; };
template<>	struct is_single_precision<Float>           {   static const bool value = true;	 };
template<>	struct is_single_precision<Real>            {   static const bool value = false; };
template<>	struct is_single_precision<Complex>         {   static const bool value = false; };
template<>	struct is_single_precision<Float_complex>   {   static const bool value = true;  };
template<>	struct is_single_precision<Object>          {   static const bool value = false; };

template<class T>
            struct is_double_precision				    {	static const bool value = false; };
template<>	struct is_double_precision<Integer>         {   static const bool value = false; };
template<>	struct is_double_precision<Float>           {   static const bool value = false; };
template<>	struct is_double_precision<Real>            {   static const bool value = true;  };
template<>	struct is_double_precision<Complex>         {   static const bool value = true;  };
template<>	struct is_double_precision<Float_complex>   {   static const bool value = false; };
template<>	struct is_double_precision<Object>          {   static const bool value = false; };

};};