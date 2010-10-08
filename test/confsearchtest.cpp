#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/builder.h>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;

  return path;
}

void test_simple() 
{
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("sdf") );
  
  OBMol mol;
  OB_REQUIRE( conv.ReadFile(&mol, GetFilename("tmp.sdf")) );

  cout << mol.NumAtoms() << endl;

  OBBuilder builder;
  builder.SetKeepRings();
  builder.Build(mol);
  

  mol.AddHydrogens();

  OBForceField* pff = OBForceField::FindType("mmff94");
  pff->Setup(mol);
  
  pff->FastRotorSearch(false);
  //pff->SteepestDescent(50);
  //pff->FastRotorSearch(0, 0);
  //pff->FastRotorSearch(0, 0);

  
  cout << mol.NumAtoms() << endl;
}

int main()
{

  test_simple();

  return 0;
}
