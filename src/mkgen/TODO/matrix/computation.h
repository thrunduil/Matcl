#pragma once

namespace matcl { namespace mkgen
{

namespace mkd = matcl::mkgen::details;

template<class Tag, class Matrix, class Assignments_List>
struct computation
{
    static_assert(is_temporary_mat<Matrix,true>::value,"matrix must be temporary");

    public:
        using subject       = Matrix;
        using deps          = typename subject::dps_type;
        using assignments   = Assignments_List;

    public:
        template<Integer Pos>
        static auto elem(colon<Pos>)        -> ct_scalar<mkd::scalar_mat_elem_1<Matrix, Pos>, deps>;

        template<Integer Row, Integer Col>
        static auto elem(colon<Row>, colon<Col>) 
                                            -> ct_scalar<mkd::scalar_mat_elem_2<Matrix, Row, Col>, deps>;

        template<class Colon_1, class Mat>
        static auto assign(Colon_1, Mat)    -> typename comp_assign_1<computation, Mat, Colon_1>::type;
};

template<class Computation>
struct get_computation_result
{
    using type = typename make_comp_result<Computation>::type;
};

}}