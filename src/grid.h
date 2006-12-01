/**********************************************************************
grid.h - Handle grids of values.
 
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

#ifndef OB_GRID_H
#define OB_GRID_H

#include "babelconfig.h"

#include <iosfwd>
#include <algorithm>
#include <vector>
#include <string>

#ifndef OBPolarGrid
#define OBPolarGrid 0x01 /* polar interactions? */
#endif //OBPolarGrid

#ifndef OBLipoGrid
#define OBLipoGrid 0x02 /* lipophilicity? */
#endif //OBLipoGrid

namespace OpenBabel
{

  //! A grid for determining the proximity of a given point to atoms in an OBMol
  class OBAPI OBProxGrid
    {
      int _gridtype;
      int _nxinc,_nyinc,_nzinc,_maxinc;
      double _xmin,_xmax,_ymin,_ymax,_zmin,_zmax,_inc;
      std::vector<std::vector<int> > cell;

    public:

      OBProxGrid(int gridtype=0)
        {
          _gridtype=gridtype;
        }
      ~OBProxGrid()
        {}
      void Setup(OBMol &mol,OBMol &box, double cutoff,double resolution = 0.5);
      void Setup(OBMol &mol,OBMol &box, double cutoff,
                 std::vector<bool> &use,double resolution = 0.5);
      std::vector<int> *GetProxVector(double,double,double);
      std::vector<int> *GetProxVector(double*);

      bool LipoGrid()
        {
          return((_gridtype&OBLipoGrid) ? true : false);
        }
      bool PolarGrid()
        {
          return(_gridtype&OBPolarGrid);
        }
      void SetGridType(int gridtype)
        {
          _gridtype = gridtype;
        }
      bool PointIsInBox(double x,double y,double z)
        {
          return (x>=_xmin) && (x<=_xmax) &&
            (y>=_ymin) && (y<=_ymax) &&
            (z>=_zmin) && (z<=_zmax);
        }
    };

  //! Handle floating-point 3D grids (i.e., charge density around an OBMol)
  class OBAPI OBFloatGrid
    {
    protected:
      double *_val;             //!< doubleing point values
      int   *_ival;            //!< for integer values (deprecated)
      double _midz,_midx,_midy;   //!< center of grid in world coordinates
      int _ydim,_xdim,_zdim;     //!< grid dimensions
      double _spacing,_inv_spa;  //!< spacing between grid points and its invers
      double _xmin,_xmax,_ymin,_ymax,_zmin,_zmax; //!< the min/max values in XYZ axes (i.e., the box)
      double _halfSpace;         //!< half of the grid spacing */

    public:

      OBFloatGrid() : _halfSpace(0.0), _val(NULL), _ival(NULL) {}
      ~OBFloatGrid()
        {
          if (_ival) delete [] _ival;
          if (_val)  delete [] _val;
        }
      //! Initialize the grid using this molecule as a box (plus a padding)
      //!  with the supplied spacing between points
      void Init(OBMol &box,double spacing, double pad= 0.0);
      //! \return whether the supplied XYZ coordinates fall within the box
      bool PointIsInBox(double x,double y,double z)
        {
          return (x>=_xmin) && (x<=_xmax) &&
            (y>=_ymin) && (y<=_ymax) &&
            (z>=_zmin) && (z<=_zmax);
        }
      //! \return true if the point falls within the box 
      bool PointIsInBox(double *c)
        {
          return (c[0]>=_xmin) && (c[0]<=_xmax) &&
            (c[1]>=_ymin) && (c[1]<=_ymax) &&
            (c[2]>=_zmin) && (c[2]<=_zmax);
        }
      double GetXmin() const    { return(_xmin);    }
      double GetYmin() const    { return(_ymin);    }
      double GetZmin() const    { return(_zmin);    }
      double GetXmax() const    { return(_xmax);    }
      double GetYmax() const    { return(_ymax);    }
      double GetZmax() const    { return(_zmax);    }
      void GetMin(double *a)
        {
          a[0]=_xmin;
          a[1]=_ymin;
          a[2]=_zmin;
        }
      void GetMax(double *a)
        {
          a[0]=_xmax;
          a[1]=_ymax;
          a[2]=_zmax;
        }

      double GetSpacing() const { return(_spacing); }
      void GetSpacing(double &s)
        {
          s=_spacing;
        }
      double GetScale() const   { return(_inv_spa); }
      double GetHalfSpace() const {return(_halfSpace);}
      int GetXdim() const       { return(_xdim);    }
      int GetYdim() const       { return(_ydim);    }
      int GetZdim() const       { return(_zdim);    }
      void GetDim(int *a)
        {
          a[0]=_xdim;
          a[1]=_ydim;
          a[2]=_zdim;
        }
      vector3 GetMidpointVector()
        {
          vector3 v;
          v.Set(_midx,_midy,_midz);
          return(v);
        }
      double *GetVals()    {        return(_val);    }
      void SetVals(double *ptr)    {  _val = ptr;    }
      vector3 Center()
        {
          return vector3(_midx,_midy,_midz);
        }

      friend std::ostream& operator<< ( std::ostream&, const OBFloatGrid& ) ;
      friend std::istream& operator>> ( std::istream&,OBFloatGrid& ) ;

      //! \return the value at the given point (rounding as needed)
      double Inject(double x,double y,double z);

      void IndexToCoords(int idx, double &x, double &y, double &z);
      void CoordsToIndex(int*,double*);
      int CoordsToIndex(double &x, double &y, double &z);
      //! \return the interpolated value for the given point
      double Interpolate(double,double,double);
      //! \return the interpolated value for the given point and set the dx, dy, dz derivatives
      double InterpolateDerivatives(double,double,double,double *derivatives);
    };

  // scoring function used: PLP = Piecewise Linear Potential or ChemScore algorithm
  typedef enum { Undefined = -1, PLP, ChemScore } score_t;

  //! A base class for scoring docking interactions between multiple molecules
  class OBAPI OBScoreGrid
    {
    protected:

      score_t gridtype;
      bool verbose;

    public:

      double score;

      OBScoreGrid(void)
        {
          verbose = false;
        }
      virtual ~OBScoreGrid(void)
        {}

      void    SetVerbose(bool v)
        {
          verbose = v;
        }
      void    SetType(score_t type)
        {
          gridtype = type;
        }
      score_t GetType(void)
        {
          return gridtype;
        }

      virtual void   Clear(void)
        { }
      virtual double  Eval(double *)
        {
          return -1;
        }
      virtual double  Eval(OBMol &mol)
        {
          return Eval(mol.GetCoordinates());
        }
      virtual void   Init(OBMol &, OBMol &, std::string &, double)
        {}
      virtual void   Setup(OBMol &)
        {}
      virtual void   Setup(OBMol &, std::vector<int> &)
        {}
      virtual void   Setup(std::vector<int> &)
        {}
      virtual void   Config(std::string)
        {}
      virtual bool   Read(std::string)
        {
          return false;
        }
      virtual bool   Write(std::string)
        {
          return false;
        }
      virtual vector3 Center()
        {
          return VZero;
        }
      virtual vector3 CenterMol(OBMol &)
        {
          return VZero;
        }
    };

} // end namespace OpenBabel

#endif // OB_GRID_H

//! \file grid.h
//! \brief Handle grids of values.
