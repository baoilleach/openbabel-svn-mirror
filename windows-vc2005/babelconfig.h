#ifndef OB_BCONFIG_H
#define OB_BCONFIG_H

//For use with Visual C++ (version 8)
#define HAVE_CONIO_H 1
#define HAVE_CLOCK_T 1
#define HAVE_IOSTREAM	1
#define HAVE_FSTREAM 1
#define HAVE_SSTREAM 1
#define HAVE_SNPRINTF 1
#define HAVE_STRNCASECMP 1
#define HAVE_LIBZ 1 //Compression enabled...
#define DISABLE_WRITE_COMPRESSION 1 //...but not for output
#define BABEL_VERSION  "2.2"

#define BABEL_DATADIR "."
#define snprintf _snprintf
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#define rindex(a,b) strrchr((a),(b))

#define TESTDATADIR "../../test/files/"

// Supress warning on deprecated functions
#pragma warning(disable : 4996)
// Supress warning that compiler is ignoring C++ exception specification
#pragma warning( disable : 4290 )
// Supress warning on signed/unsigned comparison with < or > (harmless, but maybe should be fixed)
#pragma warning( disable : 4018 )
// Supress warning on unreferenced formal parameter
#pragma warning( disable : 4100 )

#ifndef __FUNCTION__
#define __FUNCTION__ __FILE__
#endif

//Makes unix input files compatible with VC++6
#define ALL_READS_BINARY

//Do not use newlinebuf. (All varieties of line endings seem to be handled ok anyway)
#define NO_NEWLINEBUF

// Visual Studio Express does not have the classes in namespace std::tr1, particularly
// shared_ptr. Use the Boost library instead. It is necessary to have download Boost
// and put its top level folder in the Include section of Visual Studio's VC++ Directories.
#define USE_BOOST

/* Export of functions and global variables from DLLs
*/

#if defined (SINGLE_DLL) 
  //builds with a single DLL OpenBabel.dll
  // Supress warning on wrong dll-interface
  /*Exporting templated classes and classes with templated member variable
  requires some care. See http://support.microsoft.com/kb/168958*/
  #pragma warning(disable : 4251)
  #define OBAPI __declspec(dllexport)
  #define OBCONV __declspec(dllexport)
  #define OBERROR __declspec(dllexport)
  #define OBCOMMON __declspec(dllexport)
  #define OBFPRT __declspec(dllexport)
  #define EXTERN __declspec(dllexport) extern
#elif defined(USING_DYNAMIC_LIBS)
  //builds with multiple DLLs like OBDLL.dll, OBConv.dll , *.obf
  #pragma warning (disable : 4251) //no dll interface for some templated classes

  #if defined(OBDLL_EXPORTS)
    //OBDLL being built
    #define OBAPI __declspec(dllexport)
    #define EXTERN __declspec(dllexport) extern
  #else
    #define OBAPI __declspec(dllimport)
    #define EXTERN __declspec(dllimport) extern
  #endif

  #ifdef OBCONV_EXPORTS
    #define OBCONV __declspec(dllexport)
  #else
    #define OBCONV __declspec(dllimport)
  #endif

  #if defined(OBFPRT_EXPORTS)
    //OBFPRT being built
    #define OBFPRT __declspec(dllexport)
  #else
    #define OBFPRT __declspec(dllimport)
  #endif

  #if defined(OBDESC_EXPORTS)
    //OBDESC being built
    #define OBDESC __declspec(dllexport)
  #else
    #define OBDESC __declspec(dllimport)
  #endif

#if defined(OBCOMMON_EXPORTS)
   //OBCommonFormats being built
    #define OBCOMMON __declspec(dllexport)
  #else
    #define OBCOMMON __declspec(dllimport)
  #endif

  #if defined(OBERROR_EXPORTS)
    //OBError being built
    #define OBERROR __declspec(dllexport)
  #else
    #define OBERROR __declspec(dllimport)
  #endif
#else
  //builds without DLLs
  #define EXTERN extern
  #define OBAPI
  #define OBCONV
  #define OBERROR
  #define OBCOMMON
  #define OBFPRT

#endif
#endif //OB_BCONFIG_H