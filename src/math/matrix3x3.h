/**********************************************************************
matrix3x3.cpp - Handle 3D Rotation matrix.
 
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

#ifndef OB_MATRIX3x3_H
#define OB_MATRIX3x3_H

#include <ostream>
#include <math.h>

#include "math/vector3.h"
#include "oberror.h"

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

namespace OpenBabel
{
  class OBRandom; // class introduction in obutil.h

  // class introduction in matrix3x3.cpp
  class OBAPI matrix3x3
    {
      //! Elements of the matrix
      /*! This array holds the matrix. The first index refers to the
        row, the second the column. */
      double ele[3][3];

    public:
      //! constructs the zero-matrix
      matrix3x3(void)
        {
          ele[0][0] = 0.0;
          ele[0][1] = 0.0;
          ele[0][2] = 0.0;
          ele[1][0] = 0.0;
          ele[1][1] = 0.0;
          ele[1][2] = 0.0;
          ele[2][0] = 0.0;
          ele[2][1] = 0.0;
          ele[2][2] = 0.0;
        }

      //! constructs s times the unit matrix
      matrix3x3(double s)
        {
          ele[0][0] = s;
          ele[0][1] = 0.0;
          ele[0][2] = 0.0;
          ele[1][0] = 0.0;
          ele[1][1] = s;
          ele[1][2] = 0.0;
          ele[2][0] = 0.0;
          ele[2][1] = 0.0;
          ele[2][2] = s;
        }

      //! constructs a matrix from row vectors
      matrix3x3(vector3 row1,vector3 row2,vector3 row3)
        {
          ele[0][0] = row1.x();
          ele[0][1] = row1.y();
          ele[0][2] = row1.z();
          ele[1][0] = row2.x();
          ele[1][1] = row2.y();
          ele[1][2] = row2.z();
          ele[2][0] = row3.x();
          ele[2][1] = row3.y();
          ele[2][2] = row3.z();
        }

      //! constructs a matrix from a 3x3-array of doubles
      /*! constructs a matrix from a 3x3-array of doubles. The first
        index represents the row, the second index the column */
      matrix3x3(double d[3][3])
        {
          ele[0][0] = d[0][0];
          ele[0][1] = d[0][1];
          ele[0][2] = d[0][2];
          ele[1][0] = d[1][0];
          ele[1][1] = d[1][1];
          ele[1][2] = d[1][2];
          ele[2][0] = d[2][0];
          ele[2][1] = d[2][1];
          ele[2][2] = d[2][2];
        }

      //! access function
      /*! writes the matrix into the 1-dimensional array m, row by
        row. The array must be able to hold 9 doubles, otherwise your
        program will segfault. */
      void GetArray(double *m)
        {
          m[0] = ele[0][0];
          m[1] = ele[0][1];
          m[2] = ele[0][2];
          m[3] = ele[1][0];
          m[4] = ele[1][1];
          m[5] = ele[1][2];
          m[6] = ele[2][0];
          m[7] = ele[2][1];
          m[8] = ele[2][2];
        }

      /*! Returns a constant reference to an element of the matrix.
          row and column must be between 0 and 2. No check is done. */
      const double & operator() (int row, int column ) const
      {
        return ele[row][column];
      }

      /*! Returns a non-constant reference to an element of the matrix.
          row and column must be between 0 and 2. No check is done. */
      double & operator() (int row, int column )
      {
        return ele[row][column];
      }

      //! Calculates the inverse of a matrix.
      matrix3x3 inverse(void) const
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! Calculates the transpose of a matrix.
      matrix3x3 transpose(void) const;

      //! generates a matrix for a random rotation
      void randomRotation(OBRandom &rnd);

      //! returns the determinant of the matrix
      double determinant() const;

      //! Checks if a matrix is symmetric
      bool isSymmetric(void) const;

      //! Checks if a matrix is orthogonal
      /*! This method checks if a matrix is orthogonal, i.e.
        if all column vectors are normalized and
        are mutually orthogonal. A matrix is orthogonal if,
        and only if the transformation it describes is orthonormal.
        An orthonormal transformation is a
        transformation that preserves length and angle.

        The check is performed using the method isUnitMatrix() to
        check if
        \code
        *this * transpose()
        \endcode
        is a unit matrix. The criterion is therefore numerically quite
        tight. */
      bool isOrthogonal(void) const
        {
          return (*this * transpose()).isUnitMatrix();
        };

      //! Checks if a matrix is diagonal
      bool isDiagonal(void) const;

      //! Checks if a matrix is the unit matrix
      bool isUnitMatrix(void) const;

      //! access function
      /*! \warning row or column are not in the range 0..2, zero is returned
       *! \deprecated use the constant operator() instead
       */
      double Get(int row,int column) const
        {
#ifdef OB_OLD_MATH_CHECKS
          if (row >= 0 && row <= 2 && column >= 0 && column <= 2)
            return(ele[row][column]);
          else
            return 0.0f;
#else
          return(ele[row][column]);
#endif
        }

      //! access function
      /*! \warning if row or column are not in the range 0..2, nothing will happen
       *! \deprecated use the non-constant operator() instead
       */
      void Set(int row,int column, double v)
        {
#ifdef OB_OLD_MATH_CHECKS
          if (row >= 0 && row <= 2 && column >= 0 && column <= 2)
            ele[row][column]= v;
#else
          ele[row][column]= v;
#endif
        }

      //! access function
      /*! \warning If column is not in the range 0..2, the vector
        remains unchanged and an exception is thrown. */
      void SetColumn(int column, const vector3 &v)
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! access function
      /*! \warning If column is not in the range 0..2, the vector
        remains unchanged and an exception is thrown. */
      void SetRow(int row, const vector3 &v)
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! access function
      /*! \warning If col is not in the range 0..2, an exception is
        thrown. */
      vector3 GetColumn(unsigned int col) const
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! access function
      /*! \warning If row is not in the range 0..2, an exception is
        thrown. */
      vector3 GetRow(unsigned int row) const
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! multiplies all entries of the matrix by a scalar c
      matrix3x3 &operator*=(const double &c)
      {
        for( int i = 0; i < 3; i++ )
          for( int j = 0; j < 3; j++ )
            ele[i][j] *= c;
        return *this;
      }

      //! divides all entries of the matrix by a scalar c
      matrix3x3 &operator/=(const double &c)
      {
        return( (*this) *= ( 1.0 / c ) );
      }

      //! Calculate a rotation matrix for rotation about the x, y, and z
      //! axes by the angles specified (in degrees)
      void SetupRotMat(double x, double y, double z);

      //! calculates a matrix that represents reflection on a plane
      void PlaneReflection(const vector3 &norm);

      //! Calculates a rotation matrix, rotating around the specified axis by
      //! the specified angle (in degrees)
      void RotAboutAxisByAngle(const vector3 &axis, const double angle);

      //! Calculate an orthogonalisation matrix for a unit cell
      //! specified by the parameters alpha, beta, gamma, a, b, c
      //! where alpha, beta, and gamma are the cell angles (in degrees)
      //! and a, b, and c are the cell vector lengths
      //! Used by OBUnitCell
      void FillOrth(double alpha, double beta, double gamma,
                    double a, double b, double c);

      //! find the eigenvalues and -vectors of a symmetric matrix
      matrix3x3 findEigenvectorsIfSymmetric(vector3 &eigenvals) const
#ifdef OB_OLD_MATH_CHECKS
  throw(OBError)
#endif
      ;

      //! matrix-vector multiplication
      friend OBAPI vector3 operator *(const matrix3x3 &,const vector3 &);

      //! matrix-matrix multiplication
      friend OBAPI matrix3x3 operator *(const matrix3x3 &,const matrix3x3 &);

      friend OBAPI std::ostream& operator<< ( std::ostream&, const matrix3x3 & ) ;

      //! eigenvalue calculation
      static void jacobi(unsigned int n, double *a, double *d, double *v);
    };

#ifndef SWIG
  OBAPI vector3 center_coords(double*,int);
#endif
}

#endif // OB_MATRIX3x3_H

//! \file matrix3x3.h
//! \brief Handle 3D Rotation matrix.
