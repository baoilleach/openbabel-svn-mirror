#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>

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
  OB_REQUIRE( conv.ReadFile(&mol, GetFilename("8_min.sdf")) );

  OBForceField* pff = OBForceField::FindType("mmff94");
  pff->Setup(mol);
  pff->FastRotorSearch();
}

int main()
{

  test_simple();

  return 0;
}
