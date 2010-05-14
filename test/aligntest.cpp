#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/math/align.h>

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
  OB_ASSERT( fabs(rmsd - 0.0288675) < 1.0E-06 );
}

int main()
{

  test_simpleAlign();

  test_RMSD();

  return 0;
}
