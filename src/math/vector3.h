/**********************************************************************
vector3.h - Handle 3D coordinates.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_VECTOR_H
#define OB_VECTOR_H

#include <ostream>
#include <math.h>

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

namespace OpenBabel
{

  class matrix3x3; // declared in math/matrix3x3.h
  class OBRandom; // declared in obutil.h

  // class introduction in vector3.cpp
  class	OBAPI vector3
  {
  private :
    double _vx, _vy, _vz ;

  public :
    //! Constructor
    vector3 (const double inX=0.0, const double inY=0.0, const double inZ=0.0):
      _vx(inX), _vy(inY), _vz(inZ)
      {}
    //! Copy Constructor
    vector3 (const vector3& v):
      _vx(v._vx), _vy(v._vy), _vz(v._vz)
        { }

    //! set x,y and z-component of a vector
    void Set(const double inX, const double inY, const double inZ)
    {
      _vx = inX;
      _vy = inY;
      _vz = inZ;
    }
    //! set x,y and z-component of a vector from c[0]..c[2]
    void Set(const double *c)
    {
      _vx = c[0];
      _vy = c[1];
      _vz = c[2];
    }
    //! access function to get the x-coordinate of the vector
    void SetX(const double inX)
    {
      _vx = inX;
    }
    //! access function to get the y-coordinate of the vector
    void SetY(const double inY)
    {
      _vy = inY;
    }
    //! access function to get the z-coordinate of the vector
    void SetZ(const double inZ)
    {
      _vz = inZ;
    }
    //! set c[0]..c[2] to the components of the vector
    void Get(double *c)
    {
      c[0]=_vx;
      c[1]=_vy;
      c[2]=_vz;
    }
    //! access function to x: [0], y: [1], and z[2]
    double operator[] ( unsigned int i);

    //! assignment
    vector3& operator= ( const vector3& v)
      {
        _vx = v._vx;
        _vy = v._vy;
        _vz = v._vz;
        return *this;
      }

    //! returns the vector as a const double *.
    const double *AsArray()
    {
      return &_vx;
    }

    //! Vector addition (returns *this + v)
    vector3& operator+= ( const vector3& v)
      {
        _vx += v._vx;
        _vy += v._vy;
        _vz += v._vz;
        return *this;
      };
    //! Vector subtraction (returns *this - v)
    vector3& operator-= ( const vector3& v)
      {
        _vx -= v._vx;
        _vy -= v._vy;
        _vz -= v._vz;
        return *this;
      };
    //! Scalar addition
    vector3& operator+= ( const double* f)
      {
        _vx += f[0];
        _vy += f[1];
        _vz += f[2];
        return *this;
      };
    //! Scalar subtraction
    vector3& operator-= ( const double* f)
      {
        _vx -= f[0];
        _vy -= f[1];
        _vz -= f[2];
        return *this;
      };
    //! Scalar multiplication
    vector3& operator*= ( const double& c)
      {
        _vx *= c;
        _vy *= c;
        _vz *= c;
        return *this;
      };

    //! Scalar division
    vector3& operator/= ( const double& c)
      {
        double inv = 1.0 / c;
        return( (*this) *= inv );
      };
    //! multiplication of matrix and vector
    vector3& operator*= ( const matrix3x3 &);

    //! create a random unit vector
    void randomUnitVector(OBRandom *oeRand= 0L);

    //  Member Functions

    //! scales a vector to give it length one.
    vector3& normalize () ;

    //! tests whether a vector can be normalized
    bool CanBeNormalized () const;

    //! vector length squared
    inline double length_2 () const
    {
      return _vx*_vx + _vy*_vy + _vz*_vz;
    };
    //! vector length
    double length () const
    {
      return sqrt( length_2() );
    };
    //! access function to get the x-coordinate of the vector
    const double & x () const
    {
      return _vx ;
    } ;
    //! access function to get the y-coordinate of the vector
    const double & y () const
    {
      return _vy ;
    } ;
    //! access function to get the z-coordinate of the vector
    const double & z () const
    {
      return _vz ;
    } ;
    //! access function to set the x-coordinate of the vector
    double & x ()
    {
      return _vx ;
    } ;
    //! access function to set the y-coordinate of the vector
    double & y ()
    {
      return _vy ;
    } ;
    //! access function to set the z-coordinate of the vector
    double & z ()
    {
      return _vz ;
    } ;

    //! Comparison Methods
    // @{
    //! Equivalence of vectors
    //! \deprecated This method uses unreliable floating point == comparisons
    //!    Use vector3::IsApprox() instead.
    //! \return true if every component is equal
    int operator== ( const vector3& ) const;
    //! \deprecated This method uses unreliable floating point == comparisons
    //!    Use vector3::IsApprox() instead.
    //! \return true if at least one component of the two vectors is !=
    int operator!= ( const vector3& other ) const
    {
      return ! ( (*this) == other );
    }
    //! Safe comparison for floating-point vector3
    //! \return true if the vector *this is approximately equal to the vector
    //!         @p other, to the precision @p precision. More specifically,
    //!         this method works exactly like the OpenBabel::IsApprox()
    //!         function, replacing the absolute value for doubles by the norm
    //!         for vectors.
    //! \param other The vector for comparison
    //! \param precision This parameter plays the same role as in
    //!        OpenBabel::IsApprox().
    bool IsApprox( const vector3 & other, const double & precision ) const;
    //! }@

    //! square to the distance between *this and vv
    /*! equivalent to length_2(*this-vv)
     */
    double distSq(const vector3 &vv) const
    {
      double dx = x() - vv.x();
      double dy = y() - vv.y();
      double dz = z() - vv.z();
      return( dx*dx + dy*dy + dz*dz );
    }

    //! creates a vector of length one, orthogonal to *this.
    bool createOrthoVector(vector3 &v) const;

  };

  //! prints a representation of the vector as a row vector of the form "<0.1,1,2>"
  OBAPI std::ostream& operator<< ( std::ostream&, const vector3& );

  //  Sum, Difference, Scalar Product
  //! vector addition
  inline OBAPI vector3 operator+ ( const vector3& v1, const vector3& v2)
  {
    return vector3(v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z());
  };
  //! vector subtraction
  inline OBAPI vector3 operator- ( const vector3& v1, const vector3& v2)
  {
    return vector3(v1.x()-v2.x(), v1.y()-v2.y(), v1.z()-v2.z());
  };
  //! unary minus
  inline OBAPI vector3 operator- ( const vector3& v)
  {
    return vector3(-v.x(), -v.y(), -v.z());
  };
  //! multiplication with a scalar
  inline OBAPI vector3 operator* ( const double& c, const vector3& v)
    {
      return vector3( c*v.x(), c*v.y(), c*v.z());
    };
  //! multiplication with a scalar
  inline OBAPI vector3 operator* ( const vector3& v, const double& c)
    {
      return vector3( c*v.x(), c*v.y(), c*v.z());
    };
  //! division by a scalar
  inline OBAPI vector3 operator/ ( const vector3& v, const double& c)
  {
    return vector3( v.x()/c, v.y()/c, v.z()/c);
  };
  // @removed@ misleading operation
  // friend vector3 operator* ( const vector3 &,const vector3 &);

  //vector and matrix ops
  // @removed@ misleading operation; matrix multiplication is not commutitative
  //     friend vector3 operator *(const vector3 &v,const matrix3x3 &m);

  //! multiplication of matrix and vector
  OBAPI vector3 operator *(const matrix3x3 &m, const vector3 &v);

  //! dot product of two vectors
  inline OBAPI double dot ( const vector3& v1, const vector3& v2 )
  {
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z() ;
  }
  //! cross product of two vectors
  OBAPI vector3 cross ( const vector3&, const vector3& );

  //! Calculate angle between vectors
  OBAPI double vectorAngle ( const vector3& v1, const vector3& v2 );

  //! Calculate the torsion angle between vectors
  OBAPI double CalcTorsionAngle(const vector3 &a, const vector3 &b,
                                        const vector3 &c, const vector3 &d);

  //! Calculate the distance of point a to the plane determined by b,c,d
  OBAPI double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d);

  //  The global constant vector3 objects
  //! The zero vector: <0.0, 0.0, 0.0>
  extern OBAPI const vector3 VZero;
  //! The x unit vector: <1.0, 0.0, 0.0>
  extern OBAPI const vector3 VX;
  //! The y unit vector: <0.0, 1.0, 0.0>
  extern OBAPI const vector3 VY;
  //! The z unit vector: <0.0, 0.0, 1.0>
  extern OBAPI const vector3 VZ;

}

#endif // OB_VECTOR_H

//! \file
//! \brief Handle 3D coordinates.
