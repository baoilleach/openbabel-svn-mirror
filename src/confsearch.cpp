/**********************************************************************
confsearch.cpp - Conformer Searching Routines (see also forcefield.cpp)

Copyright (C) 2010 Noel O'Boyle <baoilleach@gmail.com>

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

#include <openbabel/forcefield.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/confsearch.h>
#include <openbabel/align.h>
#include <openbabel/tree/tree.hh>
#include <openbabel/tree/tree_util.hh>
#include <openbabel/math/vector3.h>

#include <float.h> // For DBL_MAX
#include <algorithm> // For min

namespace OpenBabel
{


  ostream& operator<<( ostream &out, OBMol &mol_b) {
    if (mol_b.NumAtoms() == 0)
      out << "Root";
    else
      out << mol_b.GetTorsion(mol_b.GetAtom(1), mol_b.GetAtom(2), mol_b.GetAtom(3), mol_b.GetAtom(4));
    return out;
  }

  OBDiversePoses::OBDiversePoses(double RMSD) {
    n_rmsd = 0;
    cutoff = RMSD;

    static const double arr[] = {3.0, 2.0, 1.5, 1.0, 0.5, 0.25};
    std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    vec.erase(std::remove_if(vec.begin(), vec.end(), std::bind2nd(std::less<double>(), (cutoff + 0.1) )), vec.end());
    vec.push_back(cutoff);

    levels = vec;

    poses.insert(poses.begin(), OBMol()); // Add a dummy top node
  }

  bool OBDiversePoses::AddPose(const OBMol &mol) {
    OBMol cmol = mol; // Store a copy of the molecule in the tree

    Tree_it node = poses.begin();
    int level = 0;
    bool first_time = true;
    align.SetRefMol(mol);

    while(true) {

      // Find whether the molecule is similar to any of the children of this node.
      // - min_node will hold the result of this search
      
      Tree_it min_node = poses.end(); // Point it at the end iterator (will test its value later)
      double rmsd;

      Tree_sit sib = poses.begin(node);
      // Skip the first child after the first time through this loop
      // - it will already have been tested against at the end of the previous loop
      if (!first_time)
        ++sib;
      for (; sib != poses.end(node); ++sib) { // Iterate over children of node
        align.SetTargetMol(*sib);
        align.Align();
        rmsd = align.GetRMSD();
        n_rmsd++;
        if (rmsd < levels.at(level)) {
          if (rmsd < cutoff)
            return false;
          min_node = sib;
          break; // Exit as soon as one is found
        }
      } // end of for loop

      if (min_node == poses.end()) {
        // No similar molecule found, so append it the children
        // and add it as the first child all the way down through the levels
        
        for(int k = level; k < levels.size(); ++k) {
          node = poses.append_child(node, cmol);
        }
        //kptree::print_tree_bracketed(poses);
        return true;
      }

      // If we reach here, then a similar molecule was found
      node = min_node;
      level++;
      while (rmsd < levels.at(level)) {
        node = poses.child(node, 0); // Get the first child
        level++;
      }

      first_time = false;

    } // end of while loop

    return true;
  }

  size_t OBDiversePoses::GetSize() {
    return poses.size() - 1; // Remove the dummy
  }

////////////////////////////////////////////////////////////////////////////////////////////////////

  OBDiversePosesB::OBDiversePosesB(const OBMol &ref, double RMSD) {
    natoms = ref.NumAtoms();
    align.SetRefMol(ref);
    n_rmsd = 0;
    cutoff = RMSD;

    static const double arr[] = {3.0, 2.0, 1.5, 1.0, 0.5, 0.25};
    std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    vec.erase(std::remove_if(vec.begin(), vec.end(), std::bind2nd(std::less<double>(), (cutoff + 0.1) )), vec.end());
    vec.push_back(cutoff);

    levels = vec;
    vector<vector3> pdummy;
    poses.insert(poses.begin(), PosePair(pdummy, 0.0)); // Add a dummy top node
  }

  bool OBDiversePosesB::AddPose(double* coords, double energy) {
    Tree_it node = poses.begin();
    int level = 0;
    bool first_time = true;

    // Convert coords to vector<vector3>
    vector<vector3> vcoords;
    for (unsigned int a = 0; a < natoms; ++a)
      vcoords.push_back(vector3(coords[a*3], coords[a*3+1], coords[a*3+2]));
    align.SetRef(vcoords);

    while(true) {

      // Find whether the molecule is similar to any of the children of this node.
      // - min_node will hold the result of this search
      
      Tree_it min_node = poses.end(); // Point it at the end iterator (will test its value later)
      double rmsd;

      Tree_sit sib = poses.begin(node);
      // Skip the first child after the first time through this loop
      // - it will already have been tested against at the end of the previous loop
      if (!first_time)
        ++sib;
      for (; sib != poses.end(node); ++sib) { // Iterate over children of node
        align.SetTarget((*sib).first);
        align.Align();
        rmsd = align.GetRMSD();
        n_rmsd++;
        if (rmsd < levels.at(level)) {
          if (rmsd < cutoff)
            return false;
          min_node = sib;
          break; // Exit as soon as one is found
        }
      } // end of for loop

      if (min_node == poses.end()) {
        // No similar molecule found, so append it the children
        // and add it as the first child all the way down through the levels
        
        for(int k = level; k < levels.size(); ++k) {
          node = poses.append_child(node, PosePair(vcoords, energy));
        }
        //kptree::print_tree_bracketed(poses);
        return true;
      }

      // If we reach here, then a similar molecule was found
      node = min_node;
      level++;
      while (rmsd < levels.at(level)) {
        node = poses.child(node, 0); // Get the first child
        level++;
      }

      first_time = false;

    } // end of while loop

    return true;
  }

  size_t OBDiversePosesB::GetSize() {
    return poses.size() - 1; // Remove the dummy
  }

  bool sortpred(const OBDiversePosesB::PosePair *a, const OBDiversePosesB::PosePair *b) {
    return (a->second < b->second);
  }

  vector<double> OBDiversePosesB::UpdateConformers(OBMol &mol) {

    // Initialise a vector of pointers to the nodes of the tree
    vector <PosePair*> confs;
    for (Tree_it node = poses.begin(); node != poses.end(); ++node)
      if (node->first.size() > 0) // Don't include the dummy head node
        confs.push_back(&*(node));

    // Sort the confs by energy (lowest first)
    sort(confs.begin(), confs.end(), sortpred);

    // For the moment, store the results in a vector
    typedef vector<PosePair*> vpp;
    vpp chosen_confs;
    
    // Loop through the confs and chose a diverse bunch
    for (vpp::iterator conf = confs.begin(); conf!=confs.end(); ++conf) {
      align.SetRef((*conf)->first);

      bool keepconf = true;
      for (vpp::iterator chosen = chosen_confs.begin(); chosen!=chosen_confs.end(); ++chosen) {
        align.SetTarget((*chosen)->first);
        align.Align();
        double rmsd = align.GetRMSD();
        if (rmsd < cutoff) {
          keepconf = false;
          break;
        }
      }

      if (keepconf)
        chosen_confs.push_back(*conf);
      
    }

    cout << "Tree size " << GetSize() << "Len chosen confs = " << chosen_confs.size();

    // Add confs to the molecule's conformer data and return the energies (these will be added by the calling function)
    vector<double> energies;
    for (vpp::iterator chosen = chosen_confs.begin(); chosen!=chosen_confs.end(); ++chosen) {
      energies.push_back((*chosen)->second);
      
      // To avoid making copies of vectors or vector3s, I am using pointers throughout
      vector<vector3> *tmp = &((*chosen)->first);
      double *confCoord = new double [natoms * 3];
      for(unsigned int a = 0; a<natoms; ++a) {
        vector3* pv3 = &(*tmp)[a];
        confCoord[a*3] = pv3->x();
        confCoord[a*3 + 1] = pv3->y();
        confCoord[a*3 + 2] = pv3->z();
      }
      mol.AddConformer(confCoord);
    }

    return energies;
  }

  int OBForceField::FastRotorSearch(bool permute)
  {
    if (_mol.NumRotors() == 0)
      return 0;
  
    int origLogLevel = _loglvl;

    // Remove all conformers (e.g. from previous conformer generators) except for current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    _energies.clear(); // Wipe any energies from previous conformer generators

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    // Start with all of the rotors in their 0 position
    // (perhaps instead I should set them randomly?)
    std::vector<int> init_rotorKey(rl.Size() + 1, 0);
    std::vector<int> rotorKey(init_rotorKey);

    unsigned int j, minj;
    double currentE, minE, best_minE;

    double *verybestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date
    double *bestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date in the current permutation
    double *minconf = new double [_mol.NumAtoms() * 3];  // store the best conformer for the current rotor
    memcpy((char*)bestconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());

    double energy_offset;
    // Can take shortcut later, as 4 components of the energy will be constant
    rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
    SetupPointers();
    energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);

    // This function relies on the fact that Rotors are ordered from the most
    // central to the most peripheral (due to CompareRotors in rotor.cpp)
    std::vector<OBRotor *> vrotors;
    OBRotor *rotor;
    for (rotor = rl.BeginRotor(ri); rotor; rotor = rl.NextRotor(ri))
      vrotors.push_back(rotor);

    // The permutations are ordered so that the first 2 permutations cover the
    // combinations of 2, and the first 6 permutations cover the combinations of 3
    const char permutations[24*4] = {0,1,2,3, 1,0,2,3, 0,2,1,3, 1,2,0,3, 2,0,1,3, 2,1,0,3,
                                     0,1,3,2, 0,2,3,1, 0,3,1,2, 0,3,2,1, 1,0,3,2, 1,2,3,0,
                                     1,3,0,2, 1,3,2,0, 2,0,3,1, 2,1,3,0, 2,3,0,1, 2,3,1,0,
                                     3,0,1,2, 3,0,2,1, 3,1,0,2, 3,1,2,0, 3,2,0,1, 3,2,1,0};
    const char factorial[5] = {0, 1, 2, 6, 24};

    char num_rotors_to_permute, num_permutations;
    if (permute)
      num_rotors_to_permute = std::min<size_t> (4, vrotors.size());
    else
      num_rotors_to_permute = 1; // i.e. just use the original order
    num_permutations = factorial[num_rotors_to_permute];

    // Initialize reordered_rotors - the order in which to test rotors
    std::vector<unsigned int> reordered_rotors(vrotors.size());
    for (int i=0; i<vrotors.size(); ++i)
      reordered_rotors[i] = i;

    std::set<unsigned int> seen;
    best_minE = DBL_MAX;  
    for (int N=0; N<num_permutations; ++N) {
      for (int i=0; i<num_rotors_to_permute; ++i)
        reordered_rotors.at(i) = *(permutations + N*4 + i);

      rotorKey = init_rotorKey;
      _mol.SetCoordinates(store_initial);
      bool quit = false;

      for (int i=0; i<reordered_rotors.size(); ++i) {
        unsigned int idx = reordered_rotors[i];
        rotor = vrotors.at(idx);
        
        minE = DBL_MAX;  
        
        for (j = 0; j < rotor->GetResolution().size(); j++) { // For each rotor position
          // Note: we could do slightly better by skipping the rotor position we already
          //       tested in the last loop (position 0 at the moment). Note that this
          //       isn't as simple as just changing the loop starting point to j = 1.
          _mol.SetCoordinates(bestconf);
          rotorKey[idx + 1] = j;
          rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
          SetupPointers();
          
          currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);
          
          if (currentE < minE) {
            minE = currentE;
            minj = j;
            memcpy((char*)minconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
          }
        } // Finished testing all positions of this rotor
        rotorKey[idx + 1] = minj;

        if (i==4) { // Check whether this rotorKey has already been chosen
          // Create a hash of the rotorKeys (given that the max value of any rotorKey is 11 from torlib.txt)
          unsigned int hash = rotorKey[1] + rotorKey[2]*12 + rotorKey[3]*12*12 + rotorKey[4]*12*12*12;
          
          if (seen.find(hash) == seen.end()) // Not seen before
            seen.insert(hash);
          else { // Already seen - no point continuing
            quit = true;
            break;
          }
        }

        memcpy((char*)bestconf,(char*)minconf,sizeof(double)*3*_mol.NumAtoms());
//#ifdef _DEBUG
//        std::cout << "Energy now is " << energy_offset + minE << std::endl;
//#endif
      } // end of this permutation
      if (!quit) {
          if (minE < best_minE) {
            best_minE = minE;
            memcpy((char*)verybestconf,(char*)bestconf,sizeof(double)*3*_mol.NumAtoms());
          }
      }

//#ifdef _DEBUG
//      if (!quit) {
//        std::cout << "Final energy is " << energy_offset + minE << " with rotor keys: ";
//        for (int k = 1; k < rotorKey.size(); ++k)
//          std::cout << rotorKey.at(k) << " ";
//        std::cout << std::endl;
//      }
//#endif
    } // end of final permutation

    _mol.SetCoordinates(verybestconf);
    SetupPointers();

    delete [] store_initial;
    delete [] bestconf;
    delete [] verybestconf;
    delete [] minconf;

//#ifdef _DEBUG
//    std::cout << "Very final energy is " << Energy() << std::endl;
//#endif
    return true;
  }

template < class V >
std::ostream& operator<< (std::ostream &out, std::vector<V> &v) {
  copy(v.begin(), v.end(), ostream_iterator<V>(out, ", ") );
  return out;
}


 int OBForceField::DiverseConfGen()
  {
    if (_mol.NumRotors() == 0)
      return 0;
  
    int origLogLevel = _loglvl;

    // Remove all conformers (e.g. from previous conformer generators) except for current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    _energies.clear(); // Wipe any energies from previous conformer generators

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    double *bestconf = new double [_mol.NumAtoms() * 3]; // store the best conformer to date
    memcpy((char*)bestconf,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());

    double energy_offset;
    // Can take shortcut later, as 4 components of the energy will be constant
    SetupPointers();
    energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);

    OBRotorKeys rotorKeys;
    OBRotor* rotor = rl.BeginRotor(ri);
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) // foreach rotor
      rotorKeys.AddRotor(rotor->GetResolution().size());
    
    // Main loop over rotamers
    OBDiversePosesB divposes(_mol, 0.25);
    do {
      _mol.SetCoordinates(store_initial);
      cout << rotorKeys.GetKey() << ": ";
      rotamerlist.SetCurrentCoordinates(_mol, rotorKeys.GetKey());
      SetupPointers();
      double currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);
      bool added = divposes.AddPose(_mol.GetCoordinates(), currentE);
      cout << currentE << " " << added << endl;
    } while (rotorKeys.Next());

    // Get results from tree
    _energies = divposes.UpdateConformers(_mol);
    // Add back the energy offset
    transform(_energies.begin(), _energies.end(), _energies.begin(), bind2nd(plus<double>(), energy_offset));

    // Clean up
    delete [] store_initial;

    return 0;
 }

} // end of namespace OpenBabel

//! \file confsearch.cpp
//! \brief Conformer searching routines
