#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/math/align.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/builder.h>

using namespace std;
using namespace OpenBabel;

typedef vector<vector3> vv3;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;

  return path;
}

void test_simpleAlign()
{
  //          |
  //  c       |        b                   e
  //          |
  //          |
  // ---------+------------------------------
  //          |
  //          |
  //  a       |        d
  //          |
  //

  vector3 a(-1, -1, 0), b(1, 1, 0), c(-1, 1, 0), d(1, -1, 0), e(sqrt(8.) + 1, 1, 0);

  vv3 ref(2), target(2), result;
  ref[0] = d; ref[1] = c;

  // Align dc to dc
  target[0] = d; target[1] = c;
  OBAlign align(ref, target);

  align.Align();
  result = align.GetAlignment();
  OB_ASSERT( result[0].IsApprox(ref[0], 1.0E-08) );
  OB_ASSERT( result[1].IsApprox(ref[1], 1.0E-08) );
  OB_ASSERT( fabs(align.GetRMSD()) < 1.0E-08 );

  // Align ab to dc
  target[0] = a; target[1] = b;
  align.SetTarget(target);
  align.Align();
  result = align.GetAlignment();
  OB_ASSERT( result[0].IsApprox(ref[0], 1.0E-08) );
  OB_ASSERT( result[1].IsApprox(ref[1], 1.0E-08) );
  OB_ASSERT( fabs(align.GetRMSD()) < 1.0E-08 );

  // Align be to dc
  target[0] = b; target[1] = e;
  align.SetTarget(target);
  align.Align();
  result = align.GetAlignment();
  OB_ASSERT( result[0].IsApprox(ref[0], 1.0E-08) );
  OB_ASSERT( result[1].IsApprox(ref[1], 1.0E-08) );
  OB_ASSERT( fabs(align.GetRMSD()) < 1.0E-08 );

  // Align bd to ac
  ref[0] = a; ref[1] = c;
  target[0] = b; target[1] = d;
  align.SetRef(ref);
  align.SetTarget(target);
  align.Align();
  result = align.GetAlignment();
  OB_ASSERT( result[0].IsApprox(ref[0], 1.0E-08) );
  OB_ASSERT( result[1].IsApprox(ref[1], 1.0E-08) );
  OB_ASSERT( fabs(align.GetRMSD()) < 1.0E-08 );
}

void test_RMSD()
{
  //  c       |
  //          |        b                   e
  //          |
  //          |
  // ---------+------------------------------
  //          |
  //          |
  //  a       |        d
  //          |
  //

  vector3 a(-1, -1, 0), b(1, 1, 0), c(-1, 1.1, 0), d(1, -1, 0), e(sqrt(8.) + 1, 1, 0);

  double rmsd;
  vv3 ref(2), target(2);
  ref[0] = d; ref[1] = b;

  // Align ac to db
  target[0] = a; target[1] = c;
  OBAlign align(ref, target);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd - 0.05) < 1.0E-06 );
}

void test_alignMol(){
  OBConversion conv;
  bool success = conv.SetInFormat("xyz");
  OB_REQUIRE( success );

  OBMol mol;
  success = conv.ReadFile(&mol, TESTDATADIR + string("test3d.xyz"));
  OB_REQUIRE( success );

  // Align molecule to itself (not using symmetry)
  OBAlign align = OBAlign(mol, mol, true, false);
  align.Align();
  double rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );

  // Rotate molecule and align it to itself
  OBMol mol_b = mol;
  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(1.0, -0.3, 0.23), 67);
  double rot_array[9];
  rot.GetArray(rot_array);
  mol_b.Rotate(rot_array);

  // Assert that rotation has occured
  OB_ASSERT( !mol_b.GetAtom(1)->GetVector().IsApprox(mol.GetAtom(1)->GetVector(), 1.0E-8) );

  align.SetTargetMol(mol_b);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );
}

void test_alignMolWithSym(){
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );

  OBMol mol;
  OB_REQUIRE( conv.ReadString(&mol, "ClC(=O)Cl") );

  OBBuilder builder;
  OB_REQUIRE( builder.Build(mol) );

  // Offset Atom#1
  OBAtom *patom = mol.GetAtom(1);
  patom->SetVector( patom->GetVector() + vector3(.1, .1, .1) );

  OBMol mol_b = mol;

  // Align mol to mol_b
  OBAlign align = OBAlign(mol, mol_b, true, true);
  align.Align();
  double rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );

  // Swap atom #1 and #4 in mol_b, and align again (also with symmetry)
  vector<int> a(4);
  a[0] = 4; a[1] = 2; a[2] = 3; a[3] = 1;
  mol_b.RenumberAtoms(a);
  align.SetTargetMol(mol_b);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );
  
  // Now align without symmetry
  align = OBAlign(mol, mol_b, true, false);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) > 1.0E-2 );

}

void test_alignWithoutHydrogens() {
  OBConversion conv;
  bool success = conv.SetInFormat("xyz");
  OB_REQUIRE( success );

  OBMol mol;
  success = conv.ReadFile(&mol, TESTDATADIR + string("test3d.xyz"));
  OB_REQUIRE( success );

  // Align molecule to itself without hydrogens
  OBAlign align = OBAlign(mol, mol, false, false);
  align.Align();
  double rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );

  // Move one of the hydrogens and rotate molecule
  OBMol clone = mol;
  OBAtom *atom = clone.GetAtom(8);
  OB_REQUIRE( atom->IsHydrogen() );
  atom->SetVector(atom->GetVector() + vector3(0.1, 0.1, 0.1));

  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(1.0, -0.3, 0.23), 67);
  double rot_array[9];
  rot.GetArray(rot_array);
  clone.Rotate(rot_array);

  // Assert that rotation has occured
  OB_ASSERT( !clone.GetAtom(1)->GetVector().IsApprox(mol.GetAtom(1)->GetVector(), 1.0E-8) );

  // Align molecule to clone, with hydrogens
  align = OBAlign(mol, clone, true, false);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) > 1.0E-3 );
  vector<vector3> result = align.GetAlignment();
  OB_ASSERT( result.size() == mol.NumAtoms() );

  // Align molecule to clone, without hydrogens
  align = OBAlign(mol, clone, false, false);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );
  result = align.GetAlignment();
  OB_ASSERT( result.size() == mol.NumAtoms() );
  OB_ASSERT( result.at(0).IsApprox( mol.GetAtom(1)->GetVector(), 1.0E-8 ) );

  // Align molecule to clone, without hydrogens but with sym
  align = OBAlign(mol, clone, false, true);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 1.0E-6 );
  result = align.GetAlignment();
  OB_ASSERT( result.size() == mol.NumAtoms() );
  OB_ASSERT( result.at(0).IsApprox( mol.GetAtom(1)->GetVector(), 1.0E-8 ) );
}

void test_alignWithSymWithoutHydrogens() {
  OBConversion conv;
  bool success = conv.SetInFormat("smi");
  OB_REQUIRE( success );

  OBMol mol;
  success = conv.ReadString(&mol, "BrCC(Cl)(Cl)Cl");
  OB_REQUIRE( success );

  OBBuilder builder;
  OB_REQUIRE( builder.Build(mol) );
  mol.AddHydrogens();

  // Rotate the CCl3
  OBMol clone = mol;
  double ang = mol.GetTorsion(1, 2, 3, 4);
  clone.SetTorsion( clone.GetAtom(1), clone.GetAtom(2), clone.GetAtom(3), clone.GetAtom(4), (ang + 120) * DEG_TO_RAD );

  // Align molecule to clone with hydrogens and without sym
  OBAlign align = OBAlign(mol, clone, true, false);
  align.Align();
  double rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) > 1.2 );

  // Align molecule to clone with hydrogens and with sym
  align = OBAlign(mol, clone, true, true);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 0.017 );

  // Align molecule to clone without hydrogens and with sym
  align = OBAlign(mol, clone, false, true);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 0.019 );
}

int main()
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  test_simpleAlign();

  test_RMSD();

  test_alignMol();

  test_alignMolWithSym();

  test_alignWithoutHydrogens();

  test_alignWithSymWithoutHydrogens();

  return 0;
}
