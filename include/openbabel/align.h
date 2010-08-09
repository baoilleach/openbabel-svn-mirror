/**********************************************************************
align.h - Align two molecules or vectors of vector3
 
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

#ifndef OB_ALIGN_H
#define OB_ALIGN_H

#include <openbabel/mol.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/permutation.h>
#include <Eigen/Core>

using namespace std;

namespace OpenBabel
{
  class OBAPI OBAlign {
  public: 
    OBAlign(bool includeH=false, bool symmetry=true);
    OBAlign(const OBMol &refmol, const OBMol &targetmol, bool includeH=false, bool symmetry=true);
    //OBAlign(const OBMol &refmol, const OBMol &targetmol, const vector<double> wts);
    OBAlign(const vector<vector3> &ref, const vector<vector3> &target);
    //OBAlign(const vector<vector3> &ref, const vector<vector3> &target, const vector<double> wts);

    // Partial Setup
    void SetRef(const vector<vector3> &ref);
    void SetTarget(const vector<vector3> &target);
    void SetRefMol(const OBMol &refmol);
    void SetTargetMol(const OBMol &targetmol);

    // Run the algorithm
    bool Align();

    // Accessor methods
    double GetRMSD();
    vector<vector3> GetAlignment();
    bool UpdateCoords(OBMol* target);
    inline size_t NumAutos() {return _pg.permutations.size();}

  private:
    bool _ready;
    bool _symmetry;
    bool _includeH;
    double _rmsd;
    OBBitVec _frag_atoms;
    PermutationGroup _pg;
    const OBMol* _prefmol;
    const OBMol* _ptargetmol;
    Eigen::MatrixXd _rotMatrix;
    Eigen::Vector3d _ref_centr, _target_centr;
    const vector<vector3> *_pref;
    const vector<vector3> *_ptarget;
    vector<vector3> _refmol_coords;
    vector<vector3> _targetmol_coords;
    Eigen::MatrixXd _result;
    Eigen::MatrixXd _mref, _mtarget;
    void VectorsToMatrix(const vector<vector3> *pcoords, Eigen::MatrixXd &coords);
    Eigen::Vector3d MoveToOrigin(Eigen::MatrixXd &coords);
    void SimpleAlign(const Eigen::MatrixXd &mtarget);
    void GetAutomorphisms();
    // Generate a mapping from the permutation map to the index of
    // correct column in _mtarget. Need to handle the fact that the
    // permutation group contains non-fragment atoms.
    // For example, map(213465) will be converted to newidx(102354).
    // If the atom with Idx=3 is not in the fragment, it will be
    // converted to newidx(10X243) instead.
    vector<unsigned int> _newidx;
  };
}

#endif // OB_ALIGN_H
