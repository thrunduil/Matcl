/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/general/config.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-core/matrix/enums.h"

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl 
{

namespace details
{
    struct MATCL_MATREP_EXPORT struct_flag_register_impl
    {
        using ptr_type = std::shared_ptr<user_flag_config>;
        static const size_t* get_code(const ptr_type&, const type_info& cl_name);
    };

    template<class User_struct>
    struct struct_flag_register
    {
        static const size_t* code;
    };
};

// predefined structures
enum class predefined_struct_type
{
    id,         // identity matrix
    diag,       // diagonal matrix
    tril,       // lower triangular matrix
    triu,       // upper triangular matrix
    qtril,      // lower quasi triangular matrix
    qtriu,      // upper quasi triangular matrix
    hessl,      // lower hessenberg matrix
    hessu,      // upper hessenberg matrix
    sym,        // symmetric matrix, it is preferred to mark real
                // symmetric matrix as symmetric but not hermitian
    her,        // hermitian matrix
};

// value classification used by struct_flag
enum class value_struct_class
{
    vc_zero,    // value represents zero scalar
    vc_one,     // value represents one scalar
    vc_pos_real,// value represents positive finite real scalar
    vc_neg_real,// value represents negative finite real scalar
    vc_real,    // value represents other real scalar
    vc_general  // unclassified values
};

// base class of all user-defined struct flags; user-defined struct flags should also
// be derived from user_flag_config
class MATCL_MATREP_EXPORT user_flag
{
    private:
        const size_t*   m_code;

    private:
        explicit user_flag(const size_t* code);

        friend details::registered_user_flags;

    public:
        // create empty flag, i.e. flag not representing any structure
        user_flag();

        // true if given flag is empty, i.e. not representing any structure
        bool                is_empty() const;

        // constructor of user-defined struct flag should call this function
        // to get valid code
        template<class Derived>
        static user_flag get()
        {            
            const size_t* code = details::struct_flag_register<Derived>::code;
            return user_flag(code);
        };

    //internal use
    public:
        size_t              get_code() const;
};

// set of functions that must be implemented for user-defined struct flags
class MATCL_MATREP_EXPORT user_flag_config
{    
    public:
        //-----------------------------------------------------------------
        //                  config functions
        //-----------------------------------------------------------------

        // check if given matrix has given structure
        virtual bool        test(const matcl::Matrix& mat) const = 0;
        
        // unique string representing given structure, tag cannot contain
        // spaces and underscores ("_"), tag is case insensitive
        virtual std::string tag() const = 0;

        //-----------------------------------------------------------------
        //            operations on user struct flag
        //-----------------------------------------------------------------
        // on default empty flag is returned in al functions in this group

        // functons conj() was called on nonreal matrix
        virtual user_flag   conj(struct_flag sf_mat) const;

        // functons abs() was called on X; is is_square = true, then X is a 
        // square matrix
        virtual user_flag   abs(struct_flag sf_mat, bool is_square) const;

        // functons real() was called on nonreal matrix 
        virtual user_flag   real(struct_flag sf_mat) const;

        // functons -X was called
        virtual user_flag   uminus(struct_flag sf_mat) const;

        // functons trans() was called
        virtual user_flag   trans(struct_flag sf_mat) const;
        
        // functons ctrans() was called
        virtual user_flag   ctrans(struct_flag sf_mat) const;
        
        // precision is lost during conversion
        virtual user_flag   precision_lost(struct_flag sf_mat) const;

        // precision is increased during conversion from single to double
        virtual user_flag   precision_increased(struct_flag sf_mat) const;

        // multiplication by scalar value of given class
        virtual user_flag   scal(struct_flag sf_mat, value_struct_class vc) const;

        // evaluation of trans(X,t1) * trans(X,t2) on a matrix with value code vc; additionally
        // function mult_both was called, thus this function allows for setting additional
        // structures resulting for multiplication the same matrix; if is_square = true, then
        // resulting matrix is square
        virtual user_flag   seft_mult(struct_flag sf_mat, trans_type t1, trans_type t2, 
                                value_code vc, bool is_square) const;

        // function kron(X,Y); both operands have given structure
        virtual user_flag   kron_both(struct_flag sf_left, struct_flag sf_right) const;

        // function kron(X,Y); first operand has given structure, second does not
        virtual user_flag   kron_1(struct_flag sf_left, struct_flag sf_right) const;

        // function kron(X,Y); second operand has given structure, first does not
        virtual user_flag   kron_2(struct_flag sf_left, struct_flag sf_right) const;

        // operator*; both operands have given structure; if resulting matrix is square, then
        // is_square = true
        virtual user_flag   mult_both(struct_flag sf_left, struct_flag sf_right, bool is_square) const;
        
        // operator*; first operand has given structure, second does not; if resulting matrix is square,
        // then is_square = true
        virtual user_flag   mult_1(struct_flag sf_left, struct_flag sf_right, bool is_square) const;
        
        // operator*; second operand has given structure, first does not; if resulting matrix is square,
        // then is_square = true
        virtual user_flag   mult_2(struct_flag sf_left, struct_flag sf_right, bool is_square) const;

        // operator+; both operands have given structure
        virtual user_flag   plus_both(struct_flag sf_left, struct_flag sf_right) const;
        
        // operator+; first operand has given structure, second does not
        virtual user_flag   plus_1(struct_flag sf_left, struct_flag sf_right) const;
        
        // operator+; second operand has given structure, first does not
        virtual user_flag   plus_2(struct_flag sf_left, struct_flag sf_right) const;

        // operator-; both operands have given structure
        virtual user_flag   minus_both(struct_flag sf_left, struct_flag sf_right) const;
        
        // operator-; first operand has given structure, second does not
        virtual user_flag   minus_1(struct_flag sf_left, struct_flag sf_right) const;
        
        // operator-; second operand has given structure, first does not
        virtual user_flag   minus_2(struct_flag sf_left, struct_flag sf_right) const;

        //-----------------------------------------------------------------
        //                  actions on other flags
        //-----------------------------------------------------------------
        // functions in the group 'operations' are called only if a matrix has given user flag;
        // functions in this group allows to set a new flag in case of some matrix operations;
        // for example positive_diagonal flag can be set if one multiplies identity matrix by
        // positive scalar; given user flag must register itself to check operations of given
        // type

        // return true if given user flag can be set if conj(X) is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_conj() const;

        // conj(X) is evaluated; return true if given user flag can be set in this case;
        // on default return false
        virtual bool        visit_conj(struct_flag mat_flags) const;

        // return true if given user flag can be set if abs(X) is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_abs() const;

        // abs(X) is evaluated; if is_square = true, then X is a square matrix; return true if
        // given user flag can be set in this case; on default return false
        virtual bool        visit_abs(struct_flag mat_flags, bool is_square) const;

        // return true if given user flag can be set if real(X) is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_real() const;

        // real(X) is evaluated; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_real(struct_flag mat_flags) const;

        // return true if given user flag can be set if -X is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_uminus() const;

        // -X is evaluared; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_uminus(struct_flag mat_flags) const;

        // return true if given user flag can be set if trans(X) is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_trans() const;

        // trans(X) is evaluated; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_trans(struct_flag mat_flags) const;

        // return true if given user flag can be set if ctrans(X) is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_ctrans() const;

        // ctrans(X) is evaluared; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_ctrans(struct_flag mat_flags) const;

        // return true if given user flag can be set if alpha * X is evaluated, but X does not 
        // have given flag; on default return false
        virtual bool        register_visit_scal() const;

        // multiplication by scalar value of given class; return true if given user flag can be set
        // in this case; on default empty flag is returned
        virtual bool        visit_scal(struct_flag mat_flags, value_struct_class vc) const;

        // return true if given user flag can be set if function trans(X,t1) * trans(X,t2) is
        // evaluated, but X does not have given flag; on default return false
        virtual bool        register_visit_self_mult() const;

        // evaluation of trans(X,t1) * trans(X,t2) on a matrix with value code vc; 
        // if is_square = true, then resulting matrix is square; return true if given user flag
        // can be set in this case; on default return false
        virtual bool        visit_self_mult(struct_flag mat_flags, trans_type t1, trans_type t2, 
                                    value_code vc, bool is_square) const;

        // return true if given user flag can be set if function kron(X,Y) is evaluated, 
        // but X and Y do not have given flag; on default return false
        virtual bool        register_visit_kron() const;

        // kron(X,Y) is evaluated; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_kron(struct_flag sf_X, struct_flag sf_Y) const;

        // return true if given user flag can be set if function mmul(X,Y) (directly or 
        // indirectly) is evaluated, but X and Y do not have given flag; on default return false
        virtual bool        register_visit_mmul() const;

        // return true if given user flag can be set if function mmul(X,Y); when function 
        // mmul(X,Y,t1, t2) was called, then transpositions of struct_flags are already evaluated;
        // on default return false
        virtual bool        visit_mmul(struct_flag sf_X, struct_flag sf_Y) const;

        // return true if given user flag can be set if function X + Y is evaluated, but X and
        // Y do not have given flag; on default return false
        virtual bool        register_visit_plus() const;

        // X + Y is evaluated; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_plus(struct_flag sf_X, struct_flag sf_Y) const;

        // return true if given user flag can be set if function X - Y is evaluated, but X and
        // Y do not have given flag; on default return false
        virtual bool        register_visit_minus() const;

        // X - Y is evaluated; return true if given user flag can be set in this case; 
        // on default return false
        virtual bool        visit_minus(struct_flag sf_X, struct_flag sf_Y) const;
};

// class representing additional structures assigned to a matrix
class MATCL_MATREP_EXPORT struct_flag
{
    public:
        // structure of upper or lower triangular part of a matrix
        enum diag_type
        {
            general = 0,        //unknwon number of sub- or superdiagonals
            one     = 1,        //one sub- or superdiagonals 
            qtriang = 2,        //quasi-triangular            
            zero    = 3,        //zero sub- or superdiagonals            
        };

    private:
        struct impl_type
        {
            struct id_tag{};

            size_t          m_ldiags    : 2;
            size_t          m_udiags    : 2;
            size_t          m_sym       : 1;
            size_t          m_her       : 1;
            size_t          m_id        : 1;
            size_t          m_user      : sizeof(size_t)*8 - 7;

            impl_type();
            impl_type(size_t code);
            impl_type(id_tag);
            impl_type(diag_type ld, diag_type ud, bool sym, bool her);

            impl_type       add(impl_type other) const;
            impl_type       add_user(impl_type other) const;
        };

    private:
        mutable impl_type   m_flag;

    public:
        // create empty struct flag (i.e. general matrix flag)
        struct_flag();

        // create one of predefined struct type
        struct_flag(predefined_struct_type type);

        // create user struct type
        struct_flag(user_flag uf);

        // change struct type 
        void                set(const struct_flag& t) const { m_flag = t.m_flag; };
        
        // add struct type, resulting struct_flag has all properties from this and
        // new struct_flag
        void                add(const struct_flag& t) const;
        
        // add user struct type, resulting struct_flag has all properties from this and
        // user propertied from new struct_flag
        void                add_user(const struct_flag& t) const;
        
        // reset struct type to general matrix
        void                reset();

        // clear all value dependent flags (i.e sym/her, id) and user flags
        void                reset_value();

        // remove all user defined flags
        void                reset_user()                    { m_flag.m_user = 0; };

        // set id flag
        void                set_id(bool is_id)              { m_flag.m_id = is_id; };
        
        // set symmetric flag; it is preferred to mark real symmetric matrices as
        // symmetric but not hermitian
        void                set_sym(bool is_sym)            { m_flag.m_sym = is_sym; };
        
        // set hermitian flag; it is preferred to mark only complex matrices as 
        // hermitian matrix matrix
        void                set_her(bool is_her)            { m_flag.m_her = is_her; };
        
        // set structure of lower triangular part
        void                set_ldiags(diag_type ld)        { m_flag.m_ldiags = ld; };
        
        // set structure of upper triangular part
        void                set_udiags(diag_type ud)        { m_flag.m_udiags = ud; };
        
        // set user defined structure 
        void                set_user(user_flag uf);

        // add symmetric flag; it is preferred to mark real symmetric matrices as
        // symmetric but not hermitian
        void                add_sym(bool is_sym)            { m_flag.m_sym |= is_sym; };
        
        // add hermitian flag; it is preferred to mark only complex matrices as 
        // hermitian matrix matrix
        void                add_her(bool is_her)            { m_flag.m_her |= is_her; };
        
        // add structure of lower triangular part
        void                add_ldiags(diag_type ld);
        
        // add structure of upper triangular part
        void                add_udiags(diag_type ud);
        
        // test for equality
        bool                operator==(const struct_flag& other) const;
        
        // test for not equal
        bool                operator!=(const struct_flag& other) const;

        // get all properties as integer; resulting value cannot be modified since
        // representation of struct_flag is for internal use only
        size_t              to_int() const;
        
        // create struct_flag from rep; rep must be created by to_int() function
        static struct_flag  from_int(size_t rep);

        // get all user properties as integer; resulting value cannot be modified since
        // representation of struct_flag is for internal use only
        size_t              to_int_user() const;
        
        // create struct_flag from representation of user flags; rep must be created by 
        // to_int_user() function
        static struct_flag  from_int_user(size_t rep);

        // easy to read string representation of struct_flag
        std::string         to_string() const;

        // serialize function
        void                save(oarchive_impl & ar, const unsigned int ver) const;
        
        // deserialize function
        void                load(iarchive_impl & ar, const unsigned int ver);
        
        // string encoding of struct type
        std::string         save_as_string() const;
        
        // load struct type from string encoding returned from save_as_string
        void                load_from_string(const std::string& name);

        // structure of lower triangular part of a matrix
        diag_type           get_ldiags() const      { return (diag_type)m_flag.m_ldiags; };
        
        // structure of upper triangular part of a matrix
        diag_type           get_udiags() const      { return (diag_type)m_flag.m_udiags; };
        
        // true if given user defined structure is set
        bool                get_user(user_flag uf) const;

        // true if matrix is diagonal
        bool                is_diag() const         { return is_tril() && is_triu();};
        
        // true if matrix is lower triangular
        bool                is_tril() const         { return m_flag.m_udiags == 3;};
        
        // true if matrix is upper triangular
        bool                is_triu() const         { return m_flag.m_ldiags == 3;};
        
        // true if matrix is lower quasi-triangular
        bool                is_qtril() const        { return m_flag.m_udiags >= 2;};
        
        // true if matrix is upper quasi-triangular
        bool                is_qtriu() const        { return m_flag.m_ldiags >= 2;};
        
        // true if matrix is lower hessenberg
        bool                is_hessl() const        { return m_flag.m_udiags >= 1;};
        
        // true if matrix is upper hessenberg
        bool                is_hessu() const        { return m_flag.m_ldiags >= 1;};
        
        // true if matrix has symmetric flag
        bool                has_sym_flag() const    { return m_flag.m_sym == 1; };
        
        // true if matrix has hermitian flag
        bool                has_her_flag() const    { return m_flag.m_her == 1; };
        
        // true if matrix has symmetric or hermitian flag
        bool                has_symher_flag() const { return has_sym_flag() || has_her_flag(); };

        // true if trans(A) == A for a matrix with A this flag
        bool                is_symmetric(bool is_square, bool is_real) const;

        // true if conj_trans(A) == A for a matrix A with this flag
        bool                is_hermitian(bool is_square, bool is_real) const;
        
        // true if identity matrix
        bool                is_id() const           { return m_flag.m_id == 1; };
        
        // true if no struct type is set
        bool                is_general() const      { return to_int() == 0; };
};

// check structures assigned to the matrix
void MATCL_MATREP_EXPORT check_struct(const matcl::Matrix& mat);

namespace details
{
    template<class User_struct>
    const size_t* struct_flag_register<User_struct>::code
        = struct_flag_register_impl
            ::get_code(struct_flag_register_impl::ptr_type(new User_struct()),
                        typeid(User_struct));
};

};

#pragma warning(pop)
