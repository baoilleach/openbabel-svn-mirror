/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef _GRID_H
#define _GRID_H

#ifdef __sgi
#include <iostream.h>
#else
#include <iostream>
#endif

#include <algorithm>
#include <vector>
#include <string>
using namespace std;

#ifndef OEPolarGrid
#define OEPolarGrid 0x01
#endif //OEPolarGrid

#ifndef OELipoGrid
#define OELipoGrid 0x02
#endif //OELipoGrid

namespace OpenBabel {

class OEProxGrid
{
  int _gridtype;
  int _nxinc,_nyinc,_nzinc,_maxinc;
  float _xmin,_xmax,_ymin,_ymax,_zmin,_zmax,_inc;
  vector<vector<int> > cell;

public:

  OEProxGrid(int gridtype=0){_gridtype=0;}
  ~OEProxGrid(){}
  void Setup(OEMol &,OEMol &,float,float res=0.5);
  void Setup(OEMol &,OEMol &,float,vector<bool>&,float res=0.5);
  vector<int> *GetProxVector(float,float,float);
  vector<int> *GetProxVector(float*);
  // For HasFlag force return type to bool so VC6.0 doesn't complain
  bool LipoGrid() {return((_gridtype&OELipoGrid) ? true : false);}
  bool PolarGrid() {return(_gridtype&OEPolarGrid);}
  void SetGridType(int gridtype) {_gridtype = gridtype;}
  bool PointIsInBox(float x,float y,float z)
    {
      if (x>_xmin && x<_xmax && y>_ymin && y<_ymax &&
	  z>_zmin && z<_zmax)
	return(true);
      return(false);
    }
};

class OEFloatGrid
{
  float *_val;             /* floating point values */
  int   *_ival;            /* for integer values */
  float _midz,_midx,_midy;   /* center of grid in world coordinates */
  int _ydim,_xdim,_zdim;     /* grid dimensions */
  float _spacing,_inv_spa;  /* spacing between grid points and its inverse*/
  float _xmin,_xmax,_ymin,_ymax,_zmin,_zmax;
  float _halfSpace;         /* half of the grid spacing */        

public:

    OEFloatGrid() : _halfSpace(0.0f) {_val=NULL;_ival=NULL;}
  ~OEFloatGrid() {if (_ival) delete [] _ival;if (_val) delete [] _val;}
  void Init(OEMol &,float, float pad= 0.0f); //initialized using boxfile
  bool PointIsInBox(float x,float y,float z)
    {
      if (x < _xmin || x > _xmax) return(false);
      if (y < _ymin || y > _ymax) return(false);
      if (z < _zmin || z > _zmax) return(false);
      return(true);
    }
  bool PointIsInBox(float *c)
    {
      if (c[0] < _xmin || c[0] > _xmax) return(false);
      if (c[1] < _ymin || c[1] > _ymax) return(false);
      if (c[2] < _zmin || c[2] > _zmax) return(false);
      return(true);
    }
  float GetXmin() const {return(_xmin);}
  float GetYmin() const {return(_ymin);}
  float GetZmin() const {return(_zmin);}
  float GetXmax() const {return(_xmax);}
  float GetYmax() const {return(_ymax);}
  float GetZmax() const {return(_zmax);}
  float GetSpacing() const {return(_spacing);}
  float GetScale() const {return(_inv_spa);}
  float GetHalfSpace() const {return(_halfSpace);}
  int GetXdim() const {return(_xdim);}
  int GetYdim() const {return(_ydim);}
  int GetZdim() const {return(_zdim);}
  void GetMin(float *a) {a[0]=_xmin;a[1]=_ymin;a[2]=_zmin;}
  void GetMax(float *a) {a[0]=_xmax;a[1]=_ymax;a[2]=_zmax;}
  void GetDim(int *a)   {a[0]=_xdim;a[1]=_ydim;a[2]=_zdim;}
  void GetSpacing(float &s) {s=_spacing;}
  Vector GetMidpointVector() 
    {Vector v; v.Set(_midx,_midy,_midz); return(v);}
  float *GetVals() {return(_val);}
  void SetVals(float *ptr) {_val = ptr;}
  Vector Center() { return Vector(_midx,_midy,_midz); } //added by jjc
  friend ostream& operator<< ( ostream&, const OEFloatGrid& ) ;
  friend istream& operator>> ( istream&,OEFloatGrid& ) ;

  float Inject(float x,float y,float z)
    {
      if((x<=_xmin)||(x>=_xmax)) return(0.0);
      if((y<=_ymin)||(y>=_ymax)) return(0.0);
      if((z<=_zmin)||(z>=_zmax)) return(0.0);

      int gx=(int)((x-_xmin)*_inv_spa);
      int gy=(int)((y-_ymin)*_inv_spa);
      int gz=(int)((z-_zmin)*_inv_spa);

      return(_val[(gz*_ydim*_xdim)+(gy*_xdim)+gx]);
    }


  void IndexToCoords(int idx, float &x, float &y, float &z);
  void CoordsToIndex(int*,float*);
  int CoordsToIndex(float &x, float &y, float &z);
  float Interpolate(float,float,float);
	float InterpolateDerivatives(float,float,float,float *derivatives);
};

typedef enum { Undefined = -1, PLP, ChemScore } score_t;

class OEScoreGrid
{
protected:

  score_t gridtype;
  bool verbose;

public:

  float score;

  OEScoreGrid(void) { verbose = false; }
  virtual ~OEScoreGrid(void) {}

  void    SetVerbose(bool v)    { verbose = v; }
  void    SetType(score_t type) { gridtype = type; }
  score_t GetType(void)         { return gridtype; }

  virtual void   Clear(void) { }
  virtual float  Eval(float *)    { return -1; }
  virtual float  Eval(OEMol &mol) { return Eval(mol.GetCoordinates()); }
  virtual void   Init(OEMol &, OEMol &, string &, float) {}
  virtual void   Setup(OEMol &) {}
  virtual void   Setup(OEMol &, vector<int> &) {}
  virtual void   Setup(vector<int> &) {}
  virtual void   Config(string) {}
  virtual bool   Read(string)       { return false; }
  virtual bool   Write(string)      { return false; }
  virtual Vector Center()           { return VZero; }
  virtual Vector CenterMol(OEMol &) { return VZero; }
};

}

#endif //_GRID_H
