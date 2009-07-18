#include "obtest.h"
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

const int From2D = 2;
const int From3D = 3;

std::string GetFilename(const std::string &filename)
{
  #ifdef TESTDATADIR
    string testdatadir = TESTDATADIR;
    string path = testdatadir + filename;
  #else
    string path = "files/" + filename;
  #endif

  return path;
}

std::string test_singleTetrahedral(const std::string &file, 
    const OBTetrahedralStereo::Config &correct, int fromDims = From3D)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("sdf");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return std::string();
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  // perceive stereochemistry from 2D/3D
  if (fromDims == From3D)
    StereoFrom3D(&mol);
  else
    StereoFrom2D(&mol);

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 1 );

  // compare the stereochemistry
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
      OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
      OB_ASSERT( ts->GetConfig() == correct );
      if ( ts->GetConfig() != correct ) {
        cout << ts->GetConfig() << endl;
        cout << correct << endl;
      }

      OBTetrahedralStereo::Config cfg = ts->GetConfig();
      // change refs 
      OBStereo::Permutate(cfg.refs, 1, 2);
      OB_ASSERT( cfg != correct );
      // change winding
      cfg.winding = (cfg.winding == OBStereo::Clockwise) ? OBStereo::AntiClockwise : OBStereo::Clockwise;
      OB_ASSERT( cfg == correct );
      // change view
      cfg.view = (cfg.view == OBStereo::ViewFrom) ? OBStereo::ViewTowards : OBStereo::ViewFrom;
      OB_ASSERT( cfg != correct );
      cfg.view = (cfg.view == OBStereo::ViewFrom) ? OBStereo::ViewTowards : OBStereo::ViewFrom;
      OB_ASSERT( cfg == correct );
      // change center
      cfg.center = 3994;
      OB_ASSERT( cfg != correct );


    }
  }

  return conv.WriteString(&mol);
}

std::string test_singleCisTrans(const std::string &file, 
    const OBCisTransStereo::Config &correct, int fromDims = From3D)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("sdf");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return std::string();
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  // perceive stereochemistry from 3D
  if (fromDims == From3D)
    StereoFrom3D(&mol);
  else
    StereoFrom2D(&mol);


  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 1 );

  // compare the stereochemistry
  for (std::vector<OBGenericData*>::iterator data = stereoData.begin(); data != stereoData.end(); ++data) {
    if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
      OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
      OB_ASSERT( ct->GetConfig() == correct );
      if ( ct->GetConfig() != correct ) {
        cout << ct->GetConfig() << endl;
        cout << correct << endl;
      }
    }
  }

  return conv.WriteString(&mol);
}

void test_noStereo(const std::string &file)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("sdf");
  conv.SetOutFormat("can");

  const std::string filename = GetFilename(file);
  ifstream ifs(filename.c_str());
  if (!ifs) {
    cerr << "Could not open " << filename << endl;
    return;
  }

  conv.Read(&mol, &ifs);
  ifs.close();

  // perceive stereochemistry from 3D
  StereoFrom3D(&mol);

  std::vector<OBGenericData *> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
  OB_ASSERT( stereoData.size() == 0 );
} 


int main()
{
  // There are several cases to test:

  //////////////////////////////////////////////////////////////////////////////
  // 1      StereoFrom3D for tetrahedral atoms
  //

  // 1.1    Input molecule with 4 refs (explicit H)
 
  // 1.1.1  234 == 234  @
  string smiles3D_1 = test_singleTetrahedral("tetrahedral1.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));
  // 1.1.2  234 == 234  @@
  string smiles3D_2 = test_singleTetrahedral("tetrahedral2.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::Clockwise));
  
  // 1.1.3 compare using ImplicitId
  string smiles3D_5 = test_singleTetrahedral("tetrahedral1.sdf",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::AntiClockwise));
   
  // 1.2    Input molecule with 3 refs (implicit H)
  
  // 1.2.1  23H == 23H  @@
  string smiles3D_3 = test_singleTetrahedral("tetrahedral3.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::Clockwise));
  // 1.2.2  23H == 23H  @
  string smiles3D_4 = test_singleTetrahedral("tetrahedral4.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), OBStereo::AntiClockwise));
  
  // 1.2.3 compare using explicit id 
  string smiles3D_6 = test_singleTetrahedral("tetrahedral4.sdf",
      OBTetrahedralStereo::Config(0, 1, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise));

  OB_ASSERT( smiles3D_1 == smiles3D_2 );
  OB_ASSERT( smiles3D_1 == smiles3D_3 );
  OB_ASSERT( smiles3D_1 == smiles3D_4 );
  OB_ASSERT( smiles3D_1 == smiles3D_5 );
  OB_ASSERT( smiles3D_1 == smiles3D_6 );

  //////////////////////////////////////////////////////////////////////////////
  // 2      StereoFrom3D for cis/trans bonds
  //

  // 2.1    Input molecule with 4 refs (explicit H)
 
  // 2.1.1  F      H   2      5
  //         \    /
  //          C==C       1  3
  //         /    \     
  //        H      F   0      4
  string smiles7 = test_singleCisTrans("cistrans1.sdf",
      OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5), OBStereo::ShapeZ));
  // implicit refs
  test_singleCisTrans("cistrans1.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitId, 2, 4, 5), OBStereo::ShapeZ));
  test_singleCisTrans("cistrans1.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitId), OBStereo::ShapeZ));
  test_singleCisTrans("cistrans1.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitId, 2, 4, OBStereo::ImplicitId), OBStereo::ShapeZ));
 
  // 2.1.2  F      F   2      4
  //         \    /
  //          C==C       1  3
  //         /    \     
  //        H      H   0      5
  string smiles8 = test_singleCisTrans("cistrans4.sdf",
      OBCisTransStereo::Config(1, 3, OBStereo::MakeRefs(0, 2, 4, 5), OBStereo::ShapeU));
  // implicit refs
  test_singleCisTrans("cistrans4.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitId, 2, 4, 5), OBStereo::ShapeU));
  test_singleCisTrans("cistrans4.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitId), OBStereo::ShapeU));
  test_singleCisTrans("cistrans4.sdf", OBCisTransStereo::Config(1, 3, 
        OBStereo::MakeRefs(OBStereo::ImplicitId, 2, 4, OBStereo::ImplicitId), OBStereo::ShapeU));
 
  // 2.2    Input molecule with 3 refs (implicit H)

  // 2.2.1  F          1      
  //         \    
  //          C==C       0  2
  //              \     
  //               F          3
  string smiles9 = test_singleCisTrans("cistrans2.sdf", OBCisTransStereo::Config(0, 2, 
      OBStereo::MakeRefs(1, OBStereo::ImplicitId, 3, OBStereo::ImplicitId), OBStereo::ShapeU));

  // 2.2.1  F      F   2      4
  //         \    /
  //          C==C       1  3
  //         /
  //        H          0
  string smiles10 = test_singleCisTrans("cistrans5.sdf", OBCisTransStereo::Config(1, 3, 
      OBStereo::MakeRefs(0, 2, 4, OBStereo::ImplicitId), OBStereo::ShapeU));

  OB_ASSERT( smiles7 == smiles9 );
  OB_ASSERT( smiles8 == smiles10 );

  // 3      No stereochemistry

  // 3.1    H      H   3      2
  //         \    /
  //          C==C       0  1
  //              \     
  //               F          4
  test_noStereo("cistrans3.sdf");

  //////////////////////////////////////////////////////////////////////////////
  // 3      StereoFrom2D for tetrahedral atoms
  //
  
  // 3.1    Input molecule with 4 refs (2x in plane, one behind plane, one in front of plane)

  // 3.1.1  Input molecule with 4 refs (2x in plane bond, hash & wedge bond)
  //
  //         F   I         C-F : hash (from C to F)
  //          \ /          C-I : wedge (from C to F)
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_1 = test_singleTetrahedral("tetrahedral2D1.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise), From2D);
 
  // 3.1.2  Input molecule with 4 refs (2x in plane bond, 'real' hash & 'inverted' hash bond)
  //
  //         F   I         C-F : 'real' hash (from C to F)
  //          \ /          C-I : 'inverted' hash (from I to C) = 'real' wedge
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_2 = test_singleTetrahedral("tetrahedral2D4.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise), From2D);
 
  // 3.1.3  Input molecule with 4 refs (2x in plane bond, 'real' wedge & 'inverted' wedge bond)
  //
  //         F   I         C-F : 'inverted' wedge (from F to C) = 'real' hash
  //          \ /          C-I : 'real' wedge (from C to I)
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_3 = test_singleTetrahedral("tetrahedral2D5.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, 4), OBStereo::AntiClockwise), From2D);

  OB_ASSERT( smiles2D_1 == smiles2D_2 );
  OB_ASSERT( smiles2D_1 == smiles2D_3 );
  
  // 3.2    Input molecule with 3 refs (2x in plane, one behind plane or in front of plane)
 
  // 3.2.1  Input molecule with 3 refs (2x in plane bond, real wedge bond)
  //
  //           F           C-F : 'real' wedge (from C to F)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_4 = test_singleTetrahedral("tetrahedral2D2.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), 
      OBStereo::AntiClockwise), From2D);
 
  // 3.2.2  Input molecule with 3 refs (2x in plane bond, inverted wedge bond)
  //
  //           F           C-F : 'inverted' wedge (from F to C)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_5 = test_singleTetrahedral("tetrahedral2D6.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), 
      OBStereo::Clockwise), From2D);
 
  // 3.2.3  Input molecule with 3 refs (2x in plane bond, real hash bond)
  //
  //           F           C-F : 'real' hash (from C to F)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_6 = test_singleTetrahedral("tetrahedral2D3.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), 
      OBStereo::Clockwise), From2D);
 
  // 3.2.3  Input molecule with 3 refs (2x in plane bond, real hash bond)
  //
  //           F           C-F : 'inverted' hash (from F to C)
  //           |           
  //           C
  //          / \
  //        Br   Cl
  string smiles2D_7 = test_singleTetrahedral("tetrahedral2D7.mol",
      OBTetrahedralStereo::Config(1, 0, OBStereo::MakeRefs(2, 3, OBStereo::ImplicitId), 
      OBStereo::AntiClockwise), From2D);

 /*
  OB_ASSERT( smiles2D_4 == smiles2D_7 );
  OB_ASSERT( smiles2D_5 == smiles2D_6 );
*/

  return 0;
}
