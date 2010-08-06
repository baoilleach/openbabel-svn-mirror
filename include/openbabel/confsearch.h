#ifndef OB_CONFSEARCH_H
#define OB_CONFSEARCH_H

#include <openbabel/tree/tree.hh>
#include <openbabel/align.h>
//#include <vld.h>

namespace OpenBabel
{
  class OBMol;

  class OBAPI OBDiversePoses {
    public:
      OBDiversePoses(double RMSD);
      bool AddPose(const OBMol &mol);
      typedef tree<OBMol> Tree;
      typedef tree<OBMol>::iterator Tree_it;
      typedef tree<OBMol>::sibling_iterator Tree_sit;
      size_t GetSize();
      inline int GetNRMSD() {
        return n_rmsd;
      }
      
    private:
      Tree poses;
      std::vector<double> levels;
      OBAlign align;
      double cutoff;
      int n_rmsd;
    };

  class OBAPI OBDiversePosesB {
    public:
      OBDiversePosesB(const OBMol &ref, double RMSD, bool percise=false);
      bool AddPose(double* coords, double energy);
      bool AddPose(vector<vector3> coords, double energy);
      vector<double> UpdateConformers(OBMol &mol);
      typedef std::pair<vector<vector3>, double> PosePair;
      typedef tree<PosePair> Tree;
      Tree* GetTree() { return &poses; }
      typedef tree<PosePair>::iterator Tree_it;
      typedef tree<PosePair>::sibling_iterator Tree_sit;
      size_t GetSize();
      inline int GetNRMSD() {
        return n_rmsd;
      }
      
    private:
      bool _percise;
      vector<vector3> GetHeavyAtomCoords(const vector<vector3> &all_coords);
      int natoms;
      Tree poses;
      std::vector<double> levels;
      OBAlign align;
      double cutoff;
      int n_rmsd;
      OBBitVec hydrogens;
    };

}

#endif // OB_CONFSEARCH_H