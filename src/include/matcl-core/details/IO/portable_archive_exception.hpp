/*****************************************************************************/
/****************************************************************************/

// this version of EOS library is taken from Bayeux library (https://supernemo.org/Bayeux/)
// TODO: use default implementation from https://epa.codeplex.com/ when possible
 
 #pragma once
 
 #include <boost/lexical_cast.hpp>
 #include <boost/archive/basic_archive.hpp>
 #include <boost/archive/archive_exception.hpp>
 
 // hint from Johan Rade: on VMS there is still support for
 // the VAX floating point format and this macro detects it
 #if defined(__vms) && defined(__DECCXX) && !__IEEE_FLOAT
 #error "VAX floating point format is not supported!"
 #endif
 
 namespace boost {
   namespace archive {
 
     // this value is written to the top of the stream
     const signed char magic_byte = 127;
 
     // flag for fp serialization
     const unsigned no_infnan = 64;
 
     // integral type for the archive version
     typedef library_version_type archive_version_type;
 
     // version of the linked boost archive library
     const archive_version_type archive_version(BOOST_ARCHIVE_VERSION());
 
     class portable_archive_exception : public archive_exception
     {
       std::string msg;
 
     public:
       portable_archive_exception(signed char invalid_size)
         : archive_exception(other_exception)
         , msg("requested integer size exceeds type size: ")
       {
         msg += lexical_cast<std::string, int>(invalid_size);
       }
 
       portable_archive_exception()
         : archive_exception(other_exception)
         , msg("cannot read a negative number into an unsigned type")
       {
       }
 
       template <typename T>
       portable_archive_exception(const T& abnormal)
         : archive_exception(other_exception)
         , msg("serialization of illegal floating point value: ")
       {
         msg += lexical_cast<std::string>(abnormal);
       }
 
       const char* what() const throw() { return msg.c_str(); }
       ~portable_archive_exception() throw() {}
     };
 
   } // namespace archive
 } // namespace boost
