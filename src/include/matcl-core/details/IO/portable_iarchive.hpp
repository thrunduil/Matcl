/*****************************************************************************/
/*****************************************************************************/
 
// this version of EOS library is taken from Bayeux library (https://supernemo.org/Bayeux/)
// TODO: use default implementation from https://epa.codeplex.com/ when possible

 #pragma once
 
 #include <istream>
 
#ifndef BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #define BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #define BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
#endif

 // basic headers
 #include <boost/version.hpp>
 #include <boost/utility/enable_if.hpp>
 #include <boost/archive/basic_binary_iprimitive.hpp>
 #include <boost/archive/basic_binary_iprimitive.hpp>
 #include <boost/archive/basic_binary_iarchive.hpp>
 
 // funny polymorphics
 #include <boost/archive/detail/polymorphic_iarchive_route.hpp>
 
 // endian and fpclassify
 #include <boost/endian/conversion.hpp>
 
 #if BOOST_VERSION < 106900
 #include <boost/spirit/home/support/detail/math/fpclassify.hpp>
 // namespace alias for fp utilities
 namespace fp = boost::spirit::math;
 #else
 #if !defined(BOOST_MATH_DISABLE_STD_FPCLASSIFY)
 #error "You must build using -DBOOST_MATH_DISABLE_STD_FPCLASSIFY!"
 #endif
 #include <boost/math/special_functions/fpclassify.hpp>
 // namespace alias for fp utilities
 namespace fp = boost::math;
 # endif
 
 #ifndef BOOST_NO_STD_WSTRING
 // used for wstring to utf8 conversion
 #include <boost/program_options/config.hpp>
 #include <boost/program_options/detail/convert.hpp>
 #endif
 
 // generic type traits for numeric types
 #include <boost/type_traits/is_integral.hpp>
 #include <boost/type_traits/is_unsigned.hpp>
 #include <boost/type_traits/is_arithmetic.hpp>
 #include <boost/type_traits/is_floating_point.hpp>
 
 #include "portable_archive_exception.hpp"
 
 namespace eos {
 
 using namespace boost::archive;

         // forward declaration
         class portable_iarchive;
 
         typedef basic_binary_iprimitive <
                 portable_iarchive
                 , std::istream::char_type
                 , std::istream::traits_type
         > portable_iprimitive;
 
         class portable_iarchive : public portable_iprimitive
 
                 // Robert's example derives from common_oarchive but that lacks the
                 // load_override functions so we chose to stay one level higher
                 , public basic_binary_iarchive<portable_iarchive>
         {
                 // only needed for Robert's hack in basic_binary_iarchive::init
                 friend class basic_binary_iarchive<portable_iarchive>;
 
                 // workaround for gcc: use a dummy struct
                 // as additional argument type for overloading
                 template <int> struct dummy { dummy(int) {} };
 
                 // loads directly from stream
                 inline signed char load_signed_char()
                 {
                         signed char c;
                         portable_iprimitive::load(c);
                         return c;
                 }
 
                 // archive initialization
                 void init(unsigned flags)
                 {
                         archive_version_type input_library_version(3);
 
                         // it is vital to have version information!
                         // if we don't have any we assume boost 1.33
                         if (flags & no_header)
                                 set_library_version(input_library_version);
 
                         // extract and check the magic byte header
                         else if (load_signed_char() != magic_byte)
                                 throw archive_exception(archive_exception::invalid_signature);
 
                         else
                         {
                                 // extract version information
                                 operator>>(input_library_version);
 
                                 // throw if file version is newer than we are
                                 if (input_library_version > archive_version)
                                         throw archive_exception(archive_exception::unsupported_version);
 
                                 // else set the library version accordingly
                                 else set_library_version(input_library_version);
                         }
                 }
 
         public:
                 portable_iarchive(std::istream& is, unsigned flags = 0)
                         : portable_iprimitive(*is.rdbuf(), flags & no_codecvt)
                         , basic_binary_iarchive<portable_iarchive>(flags)
                 {
                         init(flags);
                 }
 
                 portable_iarchive(std::streambuf& sb, unsigned flags = 0)
                         : portable_iprimitive(sb, flags & no_codecvt)
                         , basic_binary_iarchive<portable_iarchive>(flags)
                 {
                         init(flags);
                 }
 
                 void load(std::string& s)
                 {
                         portable_iprimitive::load(s);
                 }
 
 #ifndef BOOST_NO_STD_WSTRING
 
                 void load(std::wstring& s)
                 {
                         std::string utf8;
                         load(utf8);
                         s = boost::from_utf8(utf8);
                 }
 #endif
 
                 void load(bool& b)
                 {
                         switch (signed char c = load_signed_char())
                         {
                         case 0: b = false; break;
                         case 1: b = load_signed_char(); break;
                         default: throw portable_archive_exception(c);
                         }
                 }
 
                 template <typename T>
                 typename boost::enable_if<boost::is_integral<T> >::type
                         load(T & t, dummy<2> = 0)
                 {
                         // get the number of bytes in the stream
                         if (signed char size = load_signed_char())
                         {
                                 // check for negative value in unsigned type
                                 if (size < 0 && boost::is_unsigned<T>::value)
                                         throw portable_archive_exception();
 
                                 // check that our type T is large enough
                                 else if ((unsigned)abs(size) > sizeof(T))
                                         throw portable_archive_exception(size);
 
                                 // reconstruct the value
                                 T temp = size < 0 ? -1 : 0;
                                 load_binary(&temp, abs(size));
                                 // load the value from little endian - it is then converted
                                 // to the target type T and fits it because size <= sizeof(T)
                                 t = boost::endian::little_to_native(temp);
                         }
 
                         else t = 0; // zero optimization
                 }
 
                 template <typename T>
                 typename boost::enable_if<boost::is_floating_point<T> >::type
                         load(T & t, dummy<3> = 0)
                 {
                         typedef typename fp::detail::fp_traits<T>::type traits;
 
                         // if you end here there are three possibilities:
                         // 1. you're serializing a long double which is not portable
                         // 2. you're serializing a double but have no 64 bit integer
                         // 3. your machine is using an unknown floating point format
                         // after reading the note above you still might decide to 
                         // deactivate this static assert and try if it works out.
                         typename traits::bits bits;
                         BOOST_STATIC_ASSERT(sizeof(bits) == sizeof(T));
                         BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_iec559);
 
                         load(bits);
                         traits::set_bits(t, bits);
 
                         // if the no_infnan flag is set we must throw here
                         if (get_flags() & no_infnan && !fp::isfinite(t))
                                 throw portable_archive_exception(t);
 
                         // if you end here your floating point type does not support 
                         // denormalized numbers. this might be the case even though 
                         // your type conforms to IEC 559 (and thus to IEEE 754)
                         if (std::numeric_limits<T>::has_denorm == std::denorm_absent
                                 && fp::fpclassify(t) == (int)FP_SUBNORMAL) // GCC4
                                 throw portable_archive_exception(t);
                 }
 
                 // in boost 1.44 version_type was splitted into library_version_type and
                 // item_version_type, plus a whole bunch of additional strong typedefs.
                 template <typename T>
                 typename boost::disable_if<boost::is_arithmetic<T> >::type
                         load(T& t, dummy<4> = 0)
                 {
                         // we provide a generic load routine for all types that feature
                         // conversion operators into an unsigned integer value like those
                         // created through BOOST_STRONG_TYPEDEF(X, some unsigned int) ie.
                         // library_version_type, collection_size_type, item_version_type,
                         // class_id_type, object_id_type, version_type and tracking_type
                         load((typename boost::uint_t<sizeof(T)*CHAR_BIT>::least&)(t));
                 }
         };
 
         // polymorphic portable binary iarchive typedef
         typedef detail::polymorphic_iarchive_route<portable_iarchive> polymorphic_portable_iarchive;
 
 }
 
 // this is required by export which registers all of your
 // classes with all the inbuilt archives plus our archive.
 BOOST_SERIALIZATION_REGISTER_ARCHIVE(eos::portable_iarchive)
 BOOST_SERIALIZATION_REGISTER_ARCHIVE(eos::polymorphic_portable_iarchive)

#ifdef BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
    #undef BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #undef BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
#endif