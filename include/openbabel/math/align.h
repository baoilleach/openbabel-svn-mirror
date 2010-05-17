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
#include <Eigen/Core>

using namespace std;

namespace OpenBabel
{
  class OBAPI OBAlign {
  public: 
    OBAlign();
    OBAlign(const OBMol &refmol, const OBMol &targetmol);
    //OBAlign(const OBMol &refmol, const OBMol &targetmol, const vector<double> wts);
    OBAlign(const vector<vector3> &ref, const vector<vector3> &target);
    //OBAlign(const vector<vector3> &ref, const vector<vector3> &target, const vector<double> wts);

    // Partial Setup
    void OBAlign::SetRef(const vector<vector3> &ref);
    void OBAlign::SetTarget(const vector<vector3> &target);
    void OBAlign::SetRefMol(const OBMol &refmol);
    void OBAlign::SetTargetMol(const OBMol &targetmol);

    // Run the algorithm
    bool Align();

    // Accessor methods
    //void GetRotMatrix(matrix3x3 *rotMatrix);
    double GetRMSD();
    vector<vector3> GetAlignment();

  private:
    bool _ready;
    double _rmsd;
    Eigen::MatrixXd _rotMatrix;
    Eigen::Vector3d _ref_centr;
    const vector<vector3> *_pref;
    const vector<vector3> *_ptarget;
    Eigen::MatrixXd _result;
    Eigen::MatrixXd _mref, _mtarget;
    void VectorsToMatrix(const vector<vector3> *pcoords, Eigen::MatrixXd &coords);
    Eigen::Vector3d MoveToOrigin(Eigen::MatrixXd &coords);
  };
}

#endif // OB_ALIGN_H
