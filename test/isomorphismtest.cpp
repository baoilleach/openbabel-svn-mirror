#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

void testIsomorphism1()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "CC1CCC(C)CC1");

  OBQuery *query = CompileMoleculeQuery(&mol);
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps = mapper->MapAll(&mol);

  OB_ASSERT( maps.size() == 4 );

  delete query;
  delete mapper;

  query = CompileSmilesQuery("C1(C)CCC(C)CC1");
  mapper = OBIsomorphismMapper::GetInstance(query);
  
  OB_ASSERT( mapper->MapFirst(&mol).size() == 8 );
  OB_ASSERT( mapper->MapUnique(&mol).size() == 1 );
  OB_ASSERT( mapper->MapAll(&mol).size() == 4 );

  delete query;
  delete mapper;
}

void testIsomorphism2()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "Cc1ccc(C)cc1");

  OBQuery *query = CompileSmilesQuery("C1=CC=CC=C1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps = mapper->MapUnique(&mol);

  cout << maps.size() << endl;

  OB_ASSERT( maps.size() == 1 );

  delete query;
  delete mapper;
}

void testIsomorphismMask()
{
  // read file: 3 6-rings
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/ \/ \/
  //
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("cml");
  std::ifstream ifs(GetFilename("isomorphism1.cml").c_str());
  OB_REQUIRE( ifs );
  conv.Read(&mol, &ifs);

  OBQuery *query = CompileSmilesQuery("C1CCCCC1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  
  // no mask
  OBIsomorphismMapper::Mappings maps = mapper->MapUnique(&mol);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 3 );

  // mask first ring
  OBBitVec mask;
  for (int i = 0; i < 6; ++i)
    mask.SetBitOn(i+1);
  maps = mapper->MapUnique(&mol, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 1 );

  // mask second ring also
  for (int i = 6; i < 10; ++i)
    mask.SetBitOn(i+1);
  maps = mapper->MapUnique(&mol, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 2 );

  // just mask last ring (atomIds 7-8, 10-13)
  mask.Clear();
  for (int i = 10; i < 14; ++i)
    mask.SetBitOn(i+1);
  mask.SetBitOn(7 + 1); mask.SetBitOn(8 + 1);
  maps = mapper->MapUnique(&mol, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 1 ); // Should be same result as masking just the first ring

  delete query;
  delete mapper;
}

void testAutomorphismMask() {
  // read file: 3 6-rings
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/ \/ \/
  //
  cout <<  "testAutomorphismMask" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("cml");
  std::ifstream ifs(GetFilename("isomorphism1.cml").c_str());
  OB_REQUIRE( ifs );
  conv.Read(&mol, &ifs);

  OBIsomorphismMapper::Mappings maps; 

  // First of all, how many automorphisms are there without any mask?
  // This takes about 20 seconds, so you may want to comment this out while debugging
  maps = FindAutomorphisms(&mol);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 4 );

  // Now, let's remove the bridge (atomId 6) of the central ring.
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/    \/
  // both rings can be flipped around exocyclic bond, the whole molecule can be mirrored
  // horizontally, this results in 2 x 2 x 2 = 8 automorphisms
  OBBitVec mask;
  mask.SetRangeOn(1, mol.NumAtoms());
  mask.SetBitOff(6+1);
  maps = FindAutomorphisms(&mol, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 8 );

  // Verify that atom Id 6 does not occur anywhere in the mappings
  OBIsomorphismMapper::Mappings::const_iterator a;
  OBIsomorphismMapper::Mapping::const_iterator b;
  for (a = maps.begin(); a != maps.end(); ++a)
    for (b = a->begin(); b!= a->end(); ++b) {
      OB_ASSERT( b->first != 6 );
      OB_ASSERT( b->second != 6 );
    }
}

void testAutomorphismMask2()
{
  // The test molecule is progesterone, a steroid (four fused non-planar rings)
  cout <<  "testAutomorphismMask2" << endl;
  OBMol mol;
  OBConversion conv;

  conv.SetInFormat("sdf");
  std::ifstream ifs(GetFilename("progesterone.sdf").c_str());
  OB_REQUIRE( ifs );
  OB_REQUIRE( conv.Read(&mol, &ifs) );
  
  vector<map<unsigned int, unsigned int> > _aut;
  OBBitVec _frag_atoms;
  FOR_ATOMS_OF_MOL(a, mol) {
    if(!(a->IsHydrogen()))
      _frag_atoms.SetBitOn(a->GetIdx());
  }
  _aut = FindAutomorphisms((OBMol*)&mol, _frag_atoms);
  OB_ASSERT( _aut.size() == 1 );

}

void testAutomorphismPreMapping()
{
  cout <<  "testAutomorphismPreMapping" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "c1(C)c(C)c(C)c(C)c(C)c1");

  vector<map<unsigned int, unsigned int> > aut;
  aut = FindAutomorphisms((OBMol*)&mol);
  cout << aut.size() << endl;
  OB_ASSERT( aut.size() == 2 );
}

int main() 
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  testIsomorphism1();
  testIsomorphism2();
  testIsomorphismMask();
  testAutomorphismMask();
  testAutomorphismMask2(); 
  testAutomorphismPreMapping();

  return 0;
}

                
