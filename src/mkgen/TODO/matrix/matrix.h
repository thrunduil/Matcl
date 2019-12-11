#pragma once

namespace matcl { namespace mkgen
{

namespace mkd = matcl::mkgen::details;

//stores statically known data, accessible through tag argument
template<Integer M, Integer N, class Tag>
using const_mat = ct_matrix<M,N, const_array<Tag>,empty_deps>;

//stores generic data unknown statically, usually supplied by some data_provider
template<Integer M, Integer N, class Tag>
using gen_mat = ct_matrix<M,N,array<Tag>,extern_deps<Tag>>;

//stores results
template<Integer M, Integer N, class Tag>
using output_mat = ct_matrix<M,N,output_array<Tag>,extern_deps<Tag>>;

template<Integer M, Integer N, class Tag>
using temp_output_mat = ct_matrix<M,N,temp_output_array<Tag,M,N>, return_deps<Tag,M,N>>;

//matrix for storing local computation results, without creating
//runtime buffers.
template<Integer M, Integer N, class Tag>
using virtual_mat = ct_matrix<M,N,virtual_array<Tag>,empty_deps>;

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
        static auto elem(colon<Pos>)        -> ct_scalar<mkd::scalar_data<mkd::scalar_mat_elem_1<Matrix, Pos>>,deps>;

        template<Integer Row, Integer Col>
        static auto elem(colon<Row>, colon<Col>) 
                                            -> ct_scalar<mkd::scalar_data<mkd::scalar_mat_elem_2<Matrix, Row, Col>>,
                                                    deps>;

        template<class Colon_1, class Mat>
        static auto assign(Colon_1, Mat)    -> typename comp_assign_1<computation, Mat, Colon_1>::type;
};

template<class Computation>
struct get_computation_result
{
    using type = typename make_comp_result<Computation>::type;
};

}}