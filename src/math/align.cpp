/**********************************************************************
align.cpp - Align two molecules or vectors of vector3
 
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

#include <openbabel/math/align.h>
#include <openbabel/graphsym.h>
#include <openbabel/permutation.h>
#include <openbabel/math/vector3.h>
#include <Eigen/Dense>

//#define VERBOSE_DEBUG

using namespace std;

namespace OpenBabel
{
  OBAlign::OBAlign() {
    _ready = false;
    _symmetry = false;
    _includeH = false;
  }

  OBAlign::OBAlign(const vector<vector3> &ref, const vector<vector3> &target)
  {
    SetRef(ref);
    SetTarget(target);
    _symmetry = false;
  }

  OBAlign::OBAlign(const OBMol &refmol, const OBMol &targetmol, bool includeH, bool symmetry) {
    _symmetry = symmetry;
    _includeH = includeH;
    SetRefMol(refmol);
    SetTargetMol(targetmol);
  }

  void OBAlign::VectorsToMatrix(const vector<vector3> *pcoords, Eigen::MatrixXd &coords) {
    
    vector<vector3>::size_type N = pcoords->size();
    coords.resize(3, N);

    // Create a 3xN matrix of the coords
    vector<vector3>::const_iterator it;
    vector<vector3>::size_type colm;
    for (colm=0,it=pcoords->begin();colm<N;++colm,++it)
      coords.col(colm) = Eigen::Vector3d( it->AsArray() );
  }

  Eigen::Vector3d OBAlign::MoveToOrigin(Eigen::MatrixXd &coords) {

    vector<vector3>::size_type N = coords.cols();

    // Find the centroid
    Eigen::Vector3d centroid;
    centroid = coords.rowwise().sum() / N;

    // Subtract the centroids
    for (vector<vector3>::size_type i=0; i<N; ++i)
      coords.col(i) -= centroid;
    return centroid;
  }

  void OBAlign::SetRef(const vector<vector3> &ref) {
    _pref = &ref;
    VectorsToMatrix(_pref, _mref);
    _ref_centr = MoveToOrigin(_mref);

    _ready = false;
  }

  void OBAlign::SetTarget(const vector<vector3> &target) {
    _ptarget = &target;
    VectorsToMatrix(_ptarget, _mtarget);
    _target_centr = MoveToOrigin(_mtarget);

    _ready = false;
  }

  void OBAlign::SetRefMol(const OBMol &refmol) {
    _prefmol = &refmol;
    _refmol_coords.resize(0);
    OBAtom const *atom;
    for (int i=1; i<=refmol.NumAtoms(); ++i) {
      atom = refmol.GetAtom(i);
      if (_includeH || !atom->IsHydrogen())
        _refmol_coords.push_back(atom->GetVector());
    }
    SetRef(_refmol_coords);
  }

  void OBAlign::SetTargetMol(const OBMol &targetmol) {
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

  void OBAlign::SimpleAlign(Eigen::MatrixXd &mtarget)
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
    _result = mtarget.transpose() * _rotMatrix;
    _result.transposeInPlace(); // To give 3xN matrix

    Eigen::MatrixXd deviation = _result - _mref;
    Eigen::MatrixXd sqr = deviation.cwise().square();
    double sum = sqr.sum();
    _rmsd = sqrt( sum / sqr.size() );

  }

  bool OBAlign::Align()
  {
    vector<vector3>::size_type N = _ptarget->size();

    if (_pref->size() != N) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot align the reference and target as they are of different size" , obError);
      return false;
    }

    if (!_symmetry) {
      SimpleAlign(_mtarget);
    }
    else {  
      // Find the automorphisms of the Reference Molecule
      OBMol workmol = *_prefmol; // OBGraphSym requires non-const OBMol

      OBBitVec* frag_atoms = new OBBitVec(workmol.NumAtoms());
      FOR_ATOMS_OF_MOL(a, workmol)
        if (_includeH || !a->IsHydrogen())
          frag_atoms->SetBitOn(a->GetIdx());
      
      OBGraphSym gs(&workmol, frag_atoms); 
      vector<unsigned int> sym_classes;
      gs.GetSymmetry(sym_classes, true);
      // delete frag_atoms; // Necessary?

      // Does any symmetry class occur twice?
      vector<unsigned int> x = sym_classes;
      sort(x.begin(), x.end());
      if (adjacent_find( x.begin(), x.end() ) == x.end()) { // No duplicate symmetry classes
        SimpleAlign(_mtarget);
      }
      else { // Get the isomorphisms and iterate over them
    
        // ...for storing the results from the lowest rmsd to date
        double min_rmsd = DBL_MAX;
        Eigen::MatrixXd result, rotMatrix;

        // Try all of the symmetry-allowed permutations
        PermutationGroup pg = OpenBabel::findAutomorphisms(&workmol, sym_classes);
        std::vector<Permutation>::const_iterator cit;
        Eigen::MatrixXd mtarget;
        
        for (cit = pg.permutations.begin(); cit != pg.permutations.end(); ++cit) {
          mtarget = _mtarget*cit->matrix().cast<double>(); // Permute the columns
          SimpleAlign(mtarget);
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
      }
    }

    _ready = true;
    return true;
  }

  vector<vector3> OBAlign::GetAlignment() {
    vector<vector3> aligned_coords;
    if (!_ready) {
      // Warn!
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

  double OBAlign::GetRMSD() {
    if (!_ready) {
      obErrorLog.ThrowError(__FUNCTION__, "RMSD not available until you call Align()" , obError);
      return (double) NULL;
    }
    
    return _rmsd;
  }

} // namespace OpenBabel

//! \file align.cpp
//! \brief Handle 3D coordinates.
