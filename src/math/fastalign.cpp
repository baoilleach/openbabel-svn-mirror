/**********************************************************************
alignfast.cpp - Align two molecules or vectors of vector3
 
Copyright (C) 2010 by Noel M. O'Boyle
 
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

#include <openbabel/babelconfig.h>

#include <vector>
#include <climits> // UINT_MAX

#include <openbabel/math/fastalign.h>
#include <openbabel/graphsym.h>
#include <openbabel/math/vector3.h>

#include <Eigen/Dense>

using namespace std;

vector<double> CalcQuarticCoeffs(const double* coords1, const double* coords2, const int len)
{
  vector<double> coeff(4);

  double          Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
  double          Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                  SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                  SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                  SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
  double          x1, x2, y1, y2, z1, z2;
  int             i;

  Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
  for (i = 0; i < len*3; i+=3)
  {
      x1 = coords1[i];
      y1 = coords1[i+1];
      z1 = coords1[i+2];
      x2 = coords2[i];
      y2 = coords2[i+1];
      z2 = coords2[i+2];
 
      Sxx += (x1 * x2);
      Sxy += (x1 * y2);
      Sxz += (x1 * z2);

      Syx += (y1 * x2);
      Syy += (y1 * y2);
      Syz += (y1 * z2);

      Szx += (z1 * x2);
      Szy += (z1 * y2);
      Szz += (z1 * z2);  
  }

  Sxx2 = Sxx * Sxx;
  Syy2 = Syy * Syy;
  Szz2 = Szz * Szz;

  Sxy2 = Sxy * Sxy;
  Syz2 = Syz * Syz;
  Sxz2 = Sxz * Sxz;

  Syx2 = Syx * Syx;
  Szy2 = Szy * Szy;
  Szx2 = Szx * Szx;

  SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
  Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  /* coeff[4] = 1.0; */
  /* coeff[3] = 0.0; */
  coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
  coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

  SxzpSzx = Sxz+Szx;
  SyzpSzy = Syz+Szy;
  SxypSyx = Sxy+Syx;
  SyzmSzy = Syz-Szy;
  SxzmSzx = Sxz-Szx;
  SxymSyx = Sxy-Syx;
  SxxpSyy = Sxx+Syy;
  SxxmSyy = Sxx-Syy;
  Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
           + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
           + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
           + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
           + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
           + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
  
  return coeff;
}

namespace OpenBabel
{
  OBFastAlign::OBFastAlign(bool includeH, bool symmetry) {
    _ready = false;
    _symmetry = symmetry;
    _includeH = includeH;
    _prefmol = 0;
    _mref = 0;
    _mtarget = 0;
  }
/*
  OBFastAlign::OBFastAlign(const vector<vector3> &ref, const vector<vector3> &target)
  {
    SetRef(ref);
    SetTarget(target);
    _symmetry = false;
    _prefmol = 0;
  }

  OBFastAlign::OBFastAlign(const OBMol &refmol, const OBMol &targetmol, bool includeH, bool symmetry) {
    _symmetry = symmetry;
    _includeH = includeH;
    SetRefMol(refmol);
    SetTargetMol(targetmol);
  }

  void OBFastAlign::VectorsToMatrix(const vector<vector3> *pcoords, Eigen::MatrixXd &coords) {
    
    vector<vector3>::size_type N = pcoords->size();
    coords.resize(3, N);

    // Create a 3xN matrix of the coords
    vector<vector3>::const_iterator it;
    vector<vector3>::size_type colm;
    for (colm=0,it=pcoords->begin();colm<N;++colm,++it)
      coords.col(colm) = Eigen::Vector3d( it->AsArray() );
  }
*/
  vector<double> OBFastAlign::MoveToOrigin(const double** origcoords, double* &newcoords) {

    int             i;
    double          xsum, ysum, zsum;
    //double         *x = coords[0], *y = coords[1], *z = coords[2];

    const double* orig = (*origcoords);
    xsum = ysum = zsum = 0.0;
    for (i = 0; i < _len*3; i+=3)
    {
        xsum += orig[i];
        ysum += orig[i+1];
        zsum += orig[i+2];
    }

    xsum /= _len;
    ysum /= _len;
    zsum /= _len;

    for (i = 0; i < _len*3; i+=3)
    {
        newcoords[i] = orig[i] - xsum;
        newcoords[i+1] = orig[i+1] - ysum;
        newcoords[i+2] = orig[i+2] - zsum; 
    }

    vector<double> centroid(3);
    centroid[0] = xsum;
    centroid[1] = ysum;
    centroid[2] = zsum;

    return centroid;
  }

  void OBFastAlign::SetRef(const double* ref, int len) {
    _pref = &ref;
    _len = len;
    
    if (_mref) {
      delete[] _mref;
      _mref = 0;
    }
    _mref = new double [_len * 3];
    if (_mtarget) {
      delete[] _mtarget;
      _mtarget = 0;
    }
    _mtarget = new double [_len * 3];

    _ref_centr = MoveToOrigin(_pref, _mref);

    _ready = false;
  }

  void OBFastAlign::SetTarget(const double* target) {
    _ptarget = &target;
    _target_centr = MoveToOrigin(_ptarget, _mtarget);

    _ready = false;
  }
/*
  void OBFastAlign::SetRefMol(const OBMol &refmol) {
    _prefmol = &refmol;

    // Set up the BitVec for the hydrogens and store the refmol coords
    _frag_atoms.Clear();
    _frag_atoms.Resize(refmol.NumAtoms() + 1);
    _refmol_coords.resize(0);
    OBAtom* atom;
    int delta = 1;
    _newidx.resize(0);

    for (int i=1; i<=refmol.NumAtoms(); ++i) {
      atom = refmol.GetAtom(i);
      if (_includeH || !atom->IsHydrogen()) {
        _frag_atoms.SetBitOn(i);
        _newidx.push_back(i - delta);
        _refmol_coords.push_back(atom->GetVector());
      }
      else {
        delta++;
        _newidx.push_back(UINT_MAX);
      }
    }
    SetRef(_refmol_coords);

    if (_symmetry)
      _aut = FindAutomorphisms((OBMol*)&refmol, _frag_atoms);
  }

  void OBFastAlign::SetTargetMol(const OBMol &targetmol) {
    _ptargetmol = &targetmol;
    _targetmol_coords.resize(0);
    OBAtom const *atom;
    for (int i=1; i<=targetmol.NumAtoms(); ++i) {
      atom = targetmol.GetAtom(i);
      if (_includeH || !atom->IsHydrogen())
        _targetmol_coords.push_back(atom->GetVector());
    }
    SetTarget(_targetmol_coords);
  }

  void OBFastAlign::SetMethod(OBFastAlign::AlignMethod method) {
    _method = method;
  }
*/
/* Evaluates the Newton-Raphson correction for the Horn quartic.
   only 11 FLOPs */
  inline static double eval_horn_NR_corrxn(const vector<double> &c, const double x)
  {
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
  }

  /* Newton-Raphson root finding */
  inline static double QCProot(const vector<double> &coeff, double guess, const double delta)
  {
    int             i;
    double          oldg;
    double initialg = guess;

    for (i = 0; i < 50; ++i)
    {
        oldg = guess;
        /* guess -= (eval_horn_quart(coeff, guess) / eval_horn_quart_deriv(coeff, guess)); */
        guess -= eval_horn_NR_corrxn(coeff, guess);
    
        if (fabs(guess - oldg) < fabs(delta*guess))
            return(guess);
    }
    
    return initialg + 1.0; // Failed to converge!
  }

/* Calculate the inner product of some coordinates.
   This is the same as the squared radius of gyration without normalization
   for the number of atoms. */
static double
CoordsInnerProd(const double *coords, const int len)
{
    int             i;
    double          sum, tmp;

    sum = 0.0;
    for (i = 0; i < len * 3; ++i)
    {
      tmp = coords[i];
      sum += tmp * tmp;
    }

    return(sum);
}


  void OBFastAlign::SimpleAlign(const double* mtarget)
  {
    double innerprod = CoordsInnerProd(_mref, _len) + CoordsInnerProd(_mtarget, _len);

    vector<double> coeffs = CalcQuarticCoeffs(_mref, _mtarget, _len);
    double lambdamax = QCProot(coeffs, 0.5 * innerprod, 1e-6);
    if (lambdamax > (0.5 * innerprod))
      _fail = true;
    else {
      double sqrdev = innerprod - (2.0 * lambdamax);
      _rmsd = sqrt(sqrdev / _len);
    }
  }

  /*
  void OBFastAlign::SimpleAlign(const Eigen::MatrixXd &mtarget)
  {
    // Covariance matrix C = X times Y(t)
    Eigen::Matrix3d C = _mref * mtarget.transpose();
    
    // Singular Value Decomposition of C into USV(t)
    Eigen::SVD<Eigen::Matrix3d> svd(C);

    // Prepare matrix T
    double sign = (C.determinant() > 0) ? 1. : -1.; // Sign of determinant
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    T(2,2) = sign;

    // Optimal rotation matrix, U, is V T U(t)
    _rotMatrix = svd.matrixV() * T * svd.matrixU().transpose();
    
    // Rotate target using rotMatrix
    _result = _rotMatrix.transpose() * mtarget;

    Eigen::MatrixXd deviation = _result - _mref;
    Eigen::MatrixXd sqr = deviation.cwise().square();
    double sum = sqr.sum();
    _rmsd = sqrt( sum / sqr.cols() );

  }
*/
  bool OBFastAlign::Align()
  {
    //vector<vector3>::size_type N = _ptarget->size();
/*
    if (_pref->size() != N) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot align the reference and target as they are of different size" , obError);
      return false;
    }
*/
    if (!_symmetry || _aut.size() == 1) {
        SimpleAlign(_mtarget);
    }/*
    else {  // Iterate over the automorphisms
   
      // ...for storing the results from the lowest rmsd to date
      double min_rmsd = DBL_MAX;
      Eigen::MatrixXd result, rotMatrix;

      // Try all of the symmetry-allowed permutations
      OBIsomorphismMapper::Mappings::const_iterator cit;
      Eigen::MatrixXd mtarget(_mtarget.rows(), _mtarget.cols());
      
      for (int k = 0; k < _aut.size(); ++k) {
        // Rearrange columns of _mtarget for this permutation
        int i=0;
        for (int j=1; j<=_prefmol->NumAtoms(); ++j) {
          if (_frag_atoms.BitIsSet(j)) {
            mtarget.col(i) = _mtarget.col(_newidx[_aut[k][j - 1]]);
            i++;
          }
        }
        if (_method == OBFastAlign::Kabsch)
          SimpleAlign(mtarget);
        else
          //TheobaldAlign(mtarget);
          ;
        if (_rmsd < min_rmsd) {
          min_rmsd = _rmsd;
          result = _result;
          rotMatrix = _rotMatrix;
        }
      }

      // Restore the best answer from memory
      _rmsd = min_rmsd;
      _result = result;
      _rotMatrix = rotMatrix;
    }*/

    _ready = true;
    return true;
  }
/*
  vector<vector3> OBFastAlign::GetAlignment() {
    vector<vector3> aligned_coords;
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "Alignment not available until you call Align()" , obError);
      return aligned_coords;
    }

    if (!_prefmol || _includeH) {
      // Add back the centroid of the reference and convert to vv3
      Eigen::Vector3d tmp;
      aligned_coords.reserve(_result.cols());
      for (int i=0; i<_result.cols(); ++i) {
        tmp = _result.col(i) + _ref_centr;
        aligned_coords.push_back(vector3(tmp(0), tmp(1), tmp(2)));
      }
    }
    else { // Need to deal with the case where hydrogens were excluded
      vector<vector3> target_coords;
      for (int i=1; i<=_ptargetmol->NumAtoms(); ++i)
        target_coords.push_back(_ptargetmol->GetAtom(i)->GetVector());
      Eigen::MatrixXd mtarget;
      VectorsToMatrix(&target_coords, mtarget);

      // Subtract the centroid of the non-H atoms
      for (vector<vector3>::size_type i=0; i<mtarget.cols(); ++i)
        mtarget.col(i) -= _target_centr;

      // Rotate
      Eigen::MatrixXd result = mtarget.transpose() * _rotMatrix;
      result.transposeInPlace();

      // Add back the centroid of the reference and convert to vv3
      Eigen::Vector3d tmp;
      aligned_coords.reserve(_result.cols());
      for (int i=0; i<result.cols(); ++i) {
        tmp = result.col(i) + _ref_centr;
        aligned_coords.push_back(vector3(tmp(0), tmp(1), tmp(2)));
      }
    }

    return aligned_coords;
  }

  bool OBFastAlign::UpdateCoords(OBMol* target) {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "Alignment not available until you call Align()" , obError);
      return false;
    }

    vector<vector3> newcoords = GetAlignment();
    if (newcoords.size() != target->NumAtoms()) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot update the target molecule with the alignment coordinates as they are of different size" , obError);
      return false;
    }

    int i = 0;
    FOR_ATOMS_OF_MOL(a, *target) {
      a->SetVector(newcoords.at(i));
      i++;
    }

    return true;
  }
*/
  double OBFastAlign::GetRMSD() {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "RMSD not available until you call Align()" , obError);
      return (double) NULL;
    }
    
    return _rmsd;
  }

} // namespace OpenBabel

//! \file align.cpp
//! \brief Handle 3D coordinates.
