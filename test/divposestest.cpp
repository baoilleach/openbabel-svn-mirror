#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/align.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/builder.h>
#include <openbabel/confsearch.h>
#include <openbabel/mol.h>

using namespace std;
using namespace OpenBabel;

typedef vector<vector3> vv3;

std::ostream& OpenBabel::operator<<( std::ostream &out, OBMol &mol) {
    out << "Hey!";
    return out;
  }

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;

  return path;
}

void test_basic()
{
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );

  OBMol mol;
  OB_REQUIRE( conv.ReadString(&mol, "ClCCC(=O)Cl") );

  OBBuilder builder;
  OB_REQUIRE( builder.Build(mol) );

  OBDiversePoses poses(0.25);
  OB_ASSERT( poses.GetSize() == 0 );

  poses.AddPose(mol);
  size_t s = poses.GetSize();
  OB_ASSERT( s == 6 );

  poses.AddPose(mol);
  s = poses.GetSize();
  OB_ASSERT( s == 6 );

  OBMol mol_b = mol;
  mol_b.SetTorsion(mol_b.GetAtom(1), mol_b.GetAtom(2), mol_b.GetAtom(3), mol_b.GetAtom(4), 2.055);
  poses.AddPose(mol_b);
  s = poses.GetSize();
  OB_ASSERT( s == 7 );

  OBMol mol_c = mol;
  mol_c.SetTorsion(mol_c.GetAtom(1), mol_c.GetAtom(2), mol_c.GetAtom(3), mol_c.GetAtom(4), 0);
  poses.AddPose(mol_c);
  s = poses.GetSize();
  OB_ASSERT( s == 9 );
}

int main()
{
//  test_tree();
  test_basic();

  return 0;
}
