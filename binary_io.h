#ifndef OE_IO_BINARY_INCLUDED
#define OE_IO_BINARY_INCLUDED

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <vector>
#include <string>

#ifndef __sgi
#include <iostream>
#else
#include <iostream.h>
#endif

using namespace std;

namespace OpenEye {

/* generic binary readers for basic data types */
unsigned int OE_io_read_binary(char *ccc, char *x, unsigned int size, unsigned int count);
unsigned int OE_io_read_binary(istream &, char *x, unsigned int size, unsigned int count);

/* generic binary writers for basic data types */
unsigned int OE_io_write_binary(char *ccc, const char *x, unsigned int size, unsigned int count);
unsigned int OE_io_write_binary(ostream &, const char *x, unsigned int size, unsigned int count);

/* binary readers for STL strings */
unsigned int OE_io_read_binary(char* ccc, string& str);
unsigned int OE_io_read_binary(istream& ifs, string& str);

/* binary writers for STL strings */
unsigned int OE_io_write_binary(char* ccc, const string& str);
unsigned int OE_io_write_binary(ostream& ofs, const string& str);

/* binary writers for compressed data */
unsigned int OE_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Writes 1+(bits*Nvalues)/8 bytes
unsigned int OE_io_write_binary_compressed(ostream&, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Writes 1+(bits*Nvalues)/8 bytes
unsigned int OE_io_write_binary_compressed(char*ccc, float *x, unsigned int bits, unsigned int Nvalues); //Writes 9+(bits*Nvalues)/8 bytes
unsigned int OE_io_write_binary_compressed(ostream&, float *x, unsigned int bits, unsigned int Nvalues); //Writes 9+(bits*Nvalues)/8 bytes

/* binary readers for compressed data */
unsigned int OE_io_read_binary_compressed(char*ccc, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Reads 1+(bits*Nvalues)/8 bytes
unsigned int OE_io_read_binary_compressed(istream&, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Reads 1+(bits*Nvalues)/8 bytes
unsigned int OE_io_read_binary_compressed(char*ccc, float *x, unsigned int bits, unsigned int Nvalues); //Reads 9+(bits*Nvalues)/8 bytes
unsigned int OE_io_read_binary_compressed(istream&, float *x, unsigned int bits, unsigned int Nvalues); //Reads 9+(bits*Nvalues)/8 bytes

/* utilities for reading/writing compressed data */
unsigned int OE_io_util_calc_NumBits(unsigned int *x, unsigned int N);
unsigned int OE_io_util_calc_NumBits(float *x, unsigned int N, float res);

/*!
**\brief An abstract base class with binary read and write function.
**Members of this class understand how to read and write themselves
**to and from platform independent binary streams and arrays using
**the functions of this class.
*/
class OEBinaryIO {
    public:
        virtual ~OEBinaryIO() {}
        virtual unsigned int WriteBinary(char *ccc) const = 0;
        virtual unsigned int ReadBinary(char *ccc) = 0;
        virtual unsigned int WriteBinary(ostream& ostr) const = 0;
        virtual unsigned int ReadBinary(istream& istr) = 0;
        virtual unsigned int BinarySize() const = 0;
  };

/*!
**\brief Writes a OEBinaryIO object to a file stream, preappending
**the record with the size of the OEBinaryIO objects binary record.
*/
ostream& operator<<(ostream& ostr, const OEBinaryIO& obj);

/*!
**\brief Reads in a OEBinaryIO object that has the size of the 
**binary record preappended to the record.
*/
istream& operator>>(istream& istr, OEBinaryIO& obj);

}//End namespace OpenEye
#endif
