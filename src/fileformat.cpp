/**********************************************************************
Copyright (C) 2000 by Geoffrey Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "fileformat.h"
#include "mol.h"

using namespace std;

namespace OpenBabel {

bool OBFileFormat::ReadMolecule(istream &ifs, OBMol &mol, char *title)
{
  switch(mol.GetInputType())
    {
    case ALCHEMY:   ReadAlchemy(ifs,mol,title);	   break;
    case BALLSTICK: ReadBallAndStick(ifs,mol,title);break;
    case BIOSYM:    ReadBiosymCAR(ifs,mol,title);  break;
    case BOX:       ReadBox(ifs,mol,title);        break;
    case CACAO:	    ReadCaccrt(ifs,mol,title);	   break;
    case CCC:       ReadCCC(ifs,mol,title);        break;
    case CHEM3D1:   ReadChem3d1(ifs,mol,title);    break;
    case CHEM3D2:   ReadChem3d2(ifs,mol,title);    break;
    case CML:       ReadCML(ifs,mol,title);        break;
    case DMOL:      ReadDMol(ifs,mol,title);       break;
    case FEATURE:   ReadFeat(ifs,mol,title);	   break;
    case GAMESSOUT: ReadGAMESS(ifs,mol,title);	   break;
    case GHEMICAL:  ReadGhemical(ifs,mol,title);   break; 
    case HIN:	    ReadHIN(ifs,mol,title);	   break;
    case NWCHEMOUT: ReadNWChem(ifs, mol, title);   break;
    case MMD:       ReadMacroModel(ifs,mol,title); break;
    case MMADS:     ReadMmads(ifs,mol,title);      break;
    case MOL2:      ReadMol2(ifs,mol,title);       break;
    case MOPACOUT:  ReadMOPAC(ifs,mol,title);	   break;
    case MOPACCART: ReadMOPACCartesian(ifs,mol,title);break;
    case MPQC:      ReadMPQC(ifs,mol,title);	   break;
    case OEBINARY:  ReadBinary(ifs,mol); 	   break;
    case PDB:       ReadPDB(ifs,mol,title);        break;
    case PREP:	    ReadAmberPrep(ifs,mol,title);  break;
    case JAGUAROUT: ReadJaguar(ifs,mol,title);     break;
    case QCHEMOUT:  ReadQChem(ifs,mol,title);	   break;
    case SDF:       ReadSDFile(ifs,mol,title);     break;
    case SMI:       ReadSmiles(ifs,mol,title);     break;
    case UNICHEM:   ReadUnichem(ifs,mol,title);	   break;
    case VIEWMOL:   ReadViewMol(ifs,mol,title);	   break;
    case XYZ:	    ReadXYZ(ifs,mol,title);	   break;

    default:
      ThrowError("Input type not defined");
    }
  
  return((ifs) ? true : false);
}

bool OBFileFormat::WriteMolecule(ostream &ofs,OBMol &mol, 
				 char *dimension, char *options)
{
  switch(mol.GetOutputType())
    {
    case ALCHEMY:   WriteAlchemy(ofs,mol);		break;
    case BALLSTICK: WriteBallAndStick(ofs,mol);		break;
    case BGF:	    WriteBGF(ofs,mol);			break;
    case CACAO:     WriteCaccrt(ofs,mol);		break;
    case CACAOINT:  WriteCacaoInternal(ofs,mol);	break;
    case CACHE:     WriteCache(ofs,mol);		break;
    case CHEMDRAW:  WriteChemDraw(ofs,mol);		break;
    case CHEM3D1:   WriteChem3d1(ofs,mol);		break;
    case CHEM3D2:   WriteChem3d2(ofs,mol);		break;
    case CML:       WriteCML(ofs,mol,dimension, options);break;
    case CSR:       WriteCSR(ofs,mol);			break;
    case CSSR:      WriteCSSR(ofs,mol);			break;
    case DMOL:      WriteDMol(ofs,mol);			break;
    case DELPDB:    WriteDelphiPDB(ofs,mol);  		break;
    case FEATURE:   WriteFeat(ofs,mol);			break;
    case FH:	    WriteFenskeZmat(ofs,mol);		break;
    case FIX:       WriteFixFile(ofs,mol);    		break;
    case GAMESSIN:  WriteGAMESS(ofs,mol);		break;
    case GHEMICAL:  WriteGhemical(ofs,mol);   		break;
    case GROMOS96A: WriteGromos96A(ofs,mol);		break;
    case GROMOS96N: WriteGromos96N(ofs,mol);		break;
    case GAUSSIANCART:WriteGaussianCart(ofs,mol);	break;
    case HIN:	    WriteHIN(ofs,mol);			break;
    case JAGUARIN:  WriteJaguar(ofs,mol);		break;
    case OEBINARY:  WriteBinary(ofs,mol);     		break;
    case NWCHEMIN:  WriteNWChem(ofs, mol);		break;
    case MMD:       WriteMacroModel(ofs,mol); 		break;
    case MMADS:     WriteMmads(ofs,mol); 		break;
    case MOL2:      WriteMol2(ofs,mol,dimension);  	break;
    case MOPACCART: WriteMOPACCartesian(ofs,mol);	break;
    case PDB:       WritePDB(ofs,mol);			break;
    case QCHEMIN:   WriteQChem(ofs,mol);		break;
    case REPORT:    WriteReport(ofs,mol);		break;
    case SDF:       WriteSDFile(ofs,mol,dimension);     break;
    case SMI:       WriteSmiles(ofs,mol);     		break;
    case TINKER:    WriteTinker(ofs,mol);		break;
    case TITLE:	    WriteTitles(ofs,mol); 	        break;
    case UNICHEM:   WriteUnichem(ofs,mol);		break;
    case VIEWMOL:   WriteViewMol(ofs,mol);		break;
    case XED:	    WriteXED(ofs,mol);			break;
    case XYZ:	    WriteXYZ(ofs,mol);			break;

    default:
      ThrowError("Output type not defined");
    }

  return((ofs) ? true : false);
}

}
