#ifndef OB_CONFSEARCH_H
#define OB_CONFSEARCH_H

#include <openbabel/tree/tree.hh>
#include <openbabel/math/align.h>
//#include <vld.h>

namespace OpenBabel
{
  class OBMol;

  class OBAPI LFSR
  {
    /**
     * Usage:
     *     LFSR lfsr(N); // where N < 2^31
     *     unsigned int d;
     *     do {
     *       d = lfsr.GetNext();
     *     } while (d != 1);
     **/
  public:
    LFSR::LFSR(unsigned int range);
    unsigned int GetNext(); // Return 0 when finished
  private:
    unsigned int _range, _lfsr, _poly;
    static const unsigned int _polynomials[31];
  };

  class OBAPI OBDiversePosesB {
    public:
      OBDiversePosesB(const OBMol &ref, double RMSD, bool percise=false);
      ~OBDiversePosesB() {
        delete palign; // Allocated with 'new'
      }
      bool AddPose(double* coords, double energy);
      bool AddPose(vector<vector3> coords, double energy);
      //vector<double> UpdateConformers(OBMol &mol);
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
      OBAlign* palign;
      const double cutoff;
      int n_rmsd;
      OBBitVec hydrogens;
    };

}

#endif // OB_CONFSEARCH_H