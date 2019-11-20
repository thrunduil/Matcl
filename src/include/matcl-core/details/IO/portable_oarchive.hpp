 /*****************************************************************************/
 /*****************************************************************************/
 
// this version of EOS library is taken from Bayeux library (https://supernemo.org/Bayeux/)
// TODO: use default implementation from https://epa.codeplex.com/ when possible

 #pragma once
 
 #include <ostream>
 
#ifndef BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #define BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #define BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
#endif

 // basic headers
 #include <boost/version.hpp>
 #include <boost/utility/enable_if.hpp>
 #include <boost/archive/basic_binary_oprimitive.hpp>
 #include <boost/archive/basic_binary_oarchive.hpp>
 #include <boost/archive/detail/polymorphic_oarchive_route.hpp>
 
 // endian and fpclassify
 #include <boost/endian/conversion.hpp>
 
 // Boost Spirit does not provide FP tools from version 1.69
 #if BOOST_VERSION < 106900
 #include <boost/spirit/home/support/detail/math/fpclassify.hpp>
 // namespace alias for fp utilities
 namespace fp = boost::spirit::math;
 #else
 // So we use here the Boost Math FP stuff but this implies
 // to build with BOOST_MATH_DISABLE_STD_FPCLASSIFY otherwise
 // the std fpclassify misses some required typedefs
 // required from boost::math::detail::fp_traits (bits...)
 // Also the boost::math::detail::fp_traits_non_native declare
 // the static constant 'significand' in place of the spirir version 'mantissa';
 // this implies a workaround (see below).
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
 #include <boost/type_traits/is_signed.hpp>
 #include <boost/type_traits/is_arithmetic.hpp>
 #include <boost/type_traits/is_floating_point.hpp>
 
 #include "portable_archive_exception.hpp"
 
 namespace eos {
 
 using namespace boost::archive;

         // forward declaration
         class portable_oarchive;
 
         typedef basic_binary_oprimitive <
                 portable_oarchive
                 , std::ostream::char_type
                 , std::ostream::traits_type
         > portable_oprimitive;
 
         class portable_oarchive : public portable_oprimitive
 
                 // Robert's example derives from common_oarchive but that lacks the
                 // save_override functions so we chose to stay one level higher
                 , public basic_binary_oarchive<portable_oarchive>
         {
                 // workaround for gcc: use a dummy struct
                 // as additional argument type for overloading
                 template<int> struct dummy { dummy(int) {} };
 
                 // stores a signed char directly to stream
                 inline void save_signed_char(const signed char& c)
                 {
                         portable_oprimitive::save(c);
                 }
 
                 // archive initialization
                 void init(unsigned flags)
                 {
                         // it is vital to have version information if the archive is
                         // to be parsed with a newer version of boost::serialization
                         // therefor we create a header, no header means boost 1.33
                         // for backwards compatibility
                         if (flags & no_header)
                                 BOOST_ASSERT(archive_version == 3);
                         else
                         {
                                 // write our minimalistic header (magic byte plus version)
                                 // the boost archives write a string instead - by calling
                                 // boost::archive::basic_binary_oarchive<derived_t>::init()
                                 save_signed_char(magic_byte);
 
                                 // write current version
 //                              save<unsigned>(archive_version);
                                 operator<<(archive_version);
                         }
                 }
 
         public:
                 portable_oarchive(std::ostream& os, unsigned flags = 0)
                         : portable_oprimitive(*os.rdbuf(), flags & no_codecvt)
                         , basic_binary_oarchive<portable_oarchive>(flags)
                 {
                         init(flags);
                 }
 
                 portable_oarchive(std::streambuf& sb, unsigned flags = 0)
                         : portable_oprimitive(sb, flags & no_codecvt)
                         , basic_binary_oarchive<portable_oarchive>(flags)
                 {
                         init(flags);
                 }
 
                 void save(const std::string& s)
                 {
                         portable_oprimitive::save(s);
                 }
 
 #ifndef BOOST_NO_STD_WSTRING
 
                 void save(const std::wstring& s)
                 {
                         save(boost::to_utf8(s));
                 }
 #endif
 
                 void save(const bool& b)
                 {
                         save_signed_char(b);
                         if (b) save_signed_char('T');
                 }
 
                 template <typename T>
                 typename boost::enable_if<boost::is_integral<T> >::type
                         save(const T & t, dummy<2> = 0)
                 {
                         if (T temp = t)
                         {
                                 // examine the number of bytes
                                 // needed to represent the number
                                 signed char size = 0;
                                 do { temp >>= CHAR_BIT; ++size; } while (temp != 0 && temp != (T)-1);
 
                                 // encode the sign bit into the size
                                 save_signed_char(t > 0 ? size : -size);
                                 BOOST_ASSERT(t > 0 || boost::is_signed<T>::value);
 
                                 // we choose to use little endian because this way we just
                                 // save the first size bytes to the stream and skip the rest
                                 temp = boost::endian::native_to_little(t);
 
                                 save_binary(&temp, size);
                         }
                         // zero optimization
                         else save_signed_char(0);
                 }
 
                 template <typename T>
                 typename boost::enable_if<boost::is_floating_point<T> >::type
                         save(const T & t, dummy<3> = 0)
                 {
                         typedef typename fp::detail::fp_traits<T>::type traits;
 
                         // if the no_infnan flag is set we must throw here
                         if (get_flags() & no_infnan && !fp::isfinite(t))
                                 throw portable_archive_exception(t);
 
                         // if you end here there are three possibilities:
                         // 1. you're serializing a long double which is not portable
                         // 2. you're serializing a double but have no 64 bit integer
                         // 3. your machine is using an unknown floating point format
                         // after reading the note above you still might decide to 
                         // deactivate this static assert and try if it works out.
                         typename traits::bits bits;
                         BOOST_STATIC_ASSERT(sizeof(bits) == sizeof(T));
                         BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_iec559);
                         
                         // examine value closely
                         switch (fp::fpclassify(t))
                         {
                                 //case FP_ZERO: bits = 0; break; 
 #if defined(BOOST_MATH_FP_TRAITS_HPP)
                         // Using the traits::significand constant as the mantissa mask from Boost/Math
                         case FP_NAN: bits = traits::exponent | traits::significand; break;
 #else
 #if defined(BOOST_SPIRIT_MATH_FP_TRAITS_HPP) 
                         // Using the traits::mantissa constant as the mantissa mask from Boost/Spirit
                         case FP_NAN: bits = traits::exponent | traits::mantissa; break;
 #endif
 #endif
                         case FP_INFINITE: bits = traits::exponent | (t < 0) * traits::sign; break;
                         case FP_SUBNORMAL: assert(std::numeric_limits<T>::has_denorm); // pass
                         case FP_ZERO: // note that floats can be ?0.0
                         case FP_NORMAL: traits::get_bits(t, bits); break;
                         default: throw portable_archive_exception(t);
                         }
 
                         save(bits);
                 }
 
                 // in boost 1.44 version_type was splitted into library_version_type and
                 // item_version_type, plus a whole bunch of additional strong typedefs.
                 template <typename T>
                 typename boost::disable_if<boost::is_arithmetic<T> >::type
                         save(const T& t, dummy<4> = 0)
                 {
                         // we provide a generic save routine for all types that feature
                         // conversion operators into an unsigned integer value like those
                         // created through BOOST_STRONG_TYPEDEF(X, some unsigned int) like
                         // library_version_type, collection_size_type, item_version_type,
                         // class_id_type, object_id_type, version_type and tracking_type
                         save((typename boost::uint_t<sizeof(T)*CHAR_BIT>::least)(t));
                 }
         };
 
         // polymorphic portable binary oarchive typedef
         typedef detail::polymorphic_oarchive_route<portable_oarchive> polymorphic_portable_oarchive;
 
 }
 
 // required by export
 BOOST_SERIALIZATION_REGISTER_ARCHIVE(eos::portable_oarchive)
 BOOST_SERIALIZATION_REGISTER_ARCHIVE(eos::polymorphic_portable_oarchive)

#ifdef BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
    #undef BOOST_MATH_DISABLE_STD_FPCLASSIFY
    #undef BOOST_MATH_DISABLE_STD_FPCLASSIFY_ADDED
#endif