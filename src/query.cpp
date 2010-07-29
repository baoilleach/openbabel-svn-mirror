#include <openbabel/query.h>
#include <openbabel/obconversion.h>

using namespace std;

namespace OpenBabel {
 
  OBQuery* CompileMoleculeQuery(OBMol *mol, const OBBitVec &mask)
  {
    // set all atoms to 1 if the mask is empty
    OBBitVec mask2 = mask;
    if (!mask2.CountBits())
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        mask2.SetBitOn(i);

    OBQuery *query = new OBQuery;
    FOR_ATOMS_OF_MOL (obatom, mol) {
      if (!mask2.BitIsSet(obatom->GetIndex()))
        continue;
      query->AddAtom(new OBQueryAtom(obatom->GetAtomicNum()));
    }
    FOR_BONDS_OF_MOL (obbond, mol) {
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      if (!mask2.BitIsSet(beginIndex) || !mask2.BitIsSet(endIndex))
        continue;
 
      query->AddBond(new OBQueryBond(query->GetAtoms()[beginIndex], query->GetAtoms()[endIndex],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;  
  }

  OBQuery* CompileSmilesQuery(const std::string &smiles, const OBBitVec &mask)
  {
    OBConversion conv;
    conv.SetInFormat("smi");
    OBMol mol;
    conv.ReadString(&mol, smiles);
    return CompileMoleculeQuery(&mol, mask);  
  }












}
