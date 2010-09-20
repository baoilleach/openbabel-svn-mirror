#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/math/fastalign.h>
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
/*
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
  OB_ASSERT( fabs(rmsd - 0.0288675) < 1.0E-06 );
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
  OB_ASSERT( fabs(rmsd) > 0.5 );

  // Align molecule to clone with hydrogens and with sym
  align = OBAlign(mol, clone, true, true);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 0.010 ); // It's actually around 0.0096 but this is as good as it gets it seems

  // Align molecule to clone without hydrogens and with sym
  align = OBAlign(mol, clone, false, true);
  align.Align();
  rmsd = align.GetRMSD();
  OB_ASSERT( fabs(rmsd) < 0.011 ); // It's actually around 0.0105
}

void testQCP()
{
    vector<vector<double> > frag_a;
    vector<vector<double> > frag_b;
    for (int i=0; i<3; ++i) {
      frag_a.push_back(vector<double>(7));
      frag_b.push_back(vector<double>(7));
      
    }
    
    frag_a[0][0] =  -2.803;
    frag_a[1][0] = -15.373;
    frag_a[2][0] =  24.556;
    frag_a[0][1] =   0.893;
    frag_a[1][1] = -16.062;
    frag_a[2][1] =  25.147;
    frag_a[0][2] =   1.368;
    frag_a[1][2] = -12.371;
    frag_a[2][2] =  25.885;
    frag_a[0][3] =  -1.651;
    frag_a[1][3] = -12.153;
    frag_a[2][3] =  28.177;
    frag_a[0][4] =  -0.440;
    frag_a[1][4] = -15.218;
    frag_a[2][4] =  30.068;
    frag_a[0][5] =   2.551;
    frag_a[1][5] = -13.273;
    frag_a[2][5] =  31.372;
    frag_a[0][6] =   0.105;
    frag_a[1][6] = -11.330;
    frag_a[2][6] =  33.567;

    frag_b[0][0] = -14.739;
    frag_b[1][0] = -18.673;
    frag_b[2][0] =  15.040;
    frag_b[0][1] = -12.473;
    frag_b[1][1] = -15.810;
    frag_b[2][1] =  16.074;
    frag_b[0][2] = -14.802;
    frag_b[1][2] = -13.307;
    frag_b[2][2] =  14.408;
    frag_b[0][3] = -17.782;
    frag_b[1][3] = -14.852;
    frag_b[2][3] =  16.171;
    frag_b[0][4] = -16.124;
    frag_b[1][4] = -14.617;
    frag_b[2][4] =  19.584;
    frag_b[0][5] = -15.029;
    frag_b[1][5] = -11.037;
    frag_b[2][5] =  18.902;
    frag_b[0][6] = -18.577;
    frag_b[1][6] = -10.001;
    frag_b[2][6] =  17.996;

    vector<vector3> ref, target;
    for(int i=0; i<7; ++i) {
      ref.push_back(vector3(frag_a[0][i], frag_a[1][i], frag_a[2][i]));
      target.push_back(vector3(frag_b[0][i], frag_b[1][i], frag_b[2][i]));
    }

    OBAlign align(ref, target);
    align.Align();
    double rmsd = align.GetRMSD();
    cout << "RMSD = " << rmsd << endl;


}
*/
double randno() {
  return (rand() % 400) / 100.;
}

void speedtest() {
  const int N = 15;
  double* target = new double(N*3);
  double* ref = new double(N*3);

  OBFastAlign align(true, false);

  srand(0);
  for (int i = 0; i < N*3; ++i) {
    ref[i] = randno();
  }
  align.SetRef(ref, N);

  for (int j = 0; j < 100000; ++j) { // 10000000
    srand(j);
    for (int i = 0; i < N*3; ++i) {
      target[i] = randno();
    }
    align.SetTarget(target);
    align.Align();   
  }
  cout << align.GetRMSD() << "\n";
  delete[] target;
  delete[] ref;

}

int main()
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  
/*
  test_simpleAlign();

  test_RMSD();

  test_alignMol();

  test_alignMolWithSym();

  test_alignWithoutHydrogens();

  test_alignWithSymWithoutHydrogens();
*/
  //testQCP();
  speedtest();
  return 0;
}

