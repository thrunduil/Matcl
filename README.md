# MATCL - Matrix Computation Library

Matrix computation and linear algebra library written in C++.
This library is a part of the Morfa project. See morfalang.org.

## matcl-core

Library contains utility functions used by other libraries including
memory allocators, input-output functions, pretty printing, option
handling, error handling, management of global objects, definition of complex type,
debugging utilities.

## matcl-dynamic

Library defines the object type, which can store values of any type, and
implements powerful multiple dispatching (or multimethods). Arbitrary function 
can be dynamically dispatched based on type info of supplied arguments (of object type).


Object of given type and functions defined on this object must be registered
first, but can be registered in other libraries. Thanks to implicit conversion rules
only a small set of overloads needs to be defined in order to call a function 
with arbitrary combination of object types.

Primary purpose of this library is ilustrated by the following example:
```cpp
#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-scalar/matcl_scalar.h"

using namespace matcl;
using namespace matcl::dynamic;

// define function exp(x - y) for real values
static Real expm(Real x, Real y)
{
    return exp(x - y);
};

// define function exp(x - y) for complex values
static Complex expm(const Complex& x, const Complex& y)
{
    return exp(x - y);
};

// define function exp(x - y) for objects
object expm(const object& x, const object& y)
{
    static function_name f("expm");

    object ret;
    eval_function::eval(f, ret, x, y);
    return ret;
};

// register functions expm defined for real and complex values
MATCL_DEFINE_FUNCTION_NAME(expm)

MATCL_REGISTER_BIN_FUNC(tag1, expm, Real, Real, 
                MATCL_FUNCTION_NAME(expm))

MATCL_REGISTER_BIN_FUNC(tag2, expm, Complex, Complex,
                MATCL_FUNCTION_NAME(expm))

void example_object()
{
    object x = object(2.0f);
    object y = object(Integer(3));
    object z = object(Float_complex(1.0f, 2.0f));

    object ret1 = expm(x, y);
    object ret2 = expm(y, z);

    disp(ret1);
    disp(ret2);
};

Output:
 0.36788
 -3.0749 - 6.7188i
```            

## matcl-mp

Library defines multiprecision integer, rational, floating point
and complex types as well as functions operating on these types.

For integer, rational, and floating point types this is just a wraper
over MPFR and MPIR library with unified interface.
Matcl-mp defines however functions for multiprecision complex type
accurate up to 1 ulp.

## matcl-mp-obj

Library defines objects storing multiprecision types from matcl-mp
library. This library does not export new functions. 

## matcl-scalar

Library defines elementary functions operating on predefined scalar
types (i.e. Integer, Float, Real, Float_complex, and Complex) and 
registers these functions for objects.

Library also defines new IO functions, timing functions, and
random number generators.

## matcl-simd

Higher level abstraction of single instruction, multiple data (SIMD) instructions
for single and double precision floating point numbers and complex numbers.

## Licence

Copyright (C) 2017  Paweł Kowal

This library is published under GPL licence.