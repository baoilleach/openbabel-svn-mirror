/**********************************************************************
reportformat.cpp - Report information about the molecule: charge, distance
             matrix angle, chiral info, etc.

Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/obmolecformat.h>

using namespace std;
namespace OpenBabel
{

  class ReportFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ReportFormat()
    {
      OBConversion::RegisterFormat("report",this);
    }

    virtual const char* Description() //required
    {
      return
        "Open Babel report format\n"
        "A detailed report on the geometry of a molecule\n"
"The report format presents a report of various molecular information,\n"
"including:\n\n"

    "* Filename / molecule title\n"
    "* Molecular formula\n"
    "* Mass\n"
    "* Exact mass (i.e., for high-resolution mass spectrometry, the mass of the most abundant elements)\n"
    "* Total charge (if not electrically neutral)\n"
    "* Total spin (if not singlet)\n"
    "* Interatomic distances\n"
    "* Atomic charges\n"
    "* Bond angles\n"
    "* Dihedral angles\n"
    "* Chirality information (including which atoms are chiral)\n"
    "* Additional comments in the input file \n\n"
        "Example for benzene::\n\n"

" FILENAME: benzene.report\n"
" FORMULA: C6H6\n"
" MASS: 78.1118\n"
" EXACT MASS: 78.0469502\n"
" INTERATOMIC DISTANCES\n\n"

"               C   1      C   2      C   3      C   4      C   5      C   6\n"
"               ------------------------------------------------------------------\n"
"    C   1    0.0000 \n"
"    C   2    1.3958     0.0000 \n"
"    C   3    2.4176     1.3958     0.0000 \n"
"    C   4    2.7916     2.4176     1.3958     0.0000 \n"
"    C   5    2.4176     2.7916     2.4176     1.3958     0.0000 \n"
"    C   6    1.3958     2.4176     2.7916     2.4176     1.3958     0.0000 \n"
"    H   7    1.0846     2.1537     3.4003     3.8761     3.4003     2.1537 \n"
"    H   8    2.1537     1.0846     2.1537     3.4003     3.8761     3.4003 \n"
"    H   9    3.4003     2.1537     1.0846     2.1537     3.4003     3.8761 \n"
"    H  10    3.8761     3.4003     2.1537     1.0846     2.1537     3.4003 \n"
"    H  11    3.4003     3.8761     3.4003     2.1537     1.0846     2.1537\n"
"    H  12    2.1537     3.4003     3.8761     3.4003     2.1537     1.0846 \n\n"

"               H   7      H   8      H   9      H  10      H  11      H  12\n"
"               ------------------------------------------------------------------\n"
"    H   7    0.0000 \n"
"    H   8    2.4803     0.0000 \n"
"    H   9    4.2961     2.4804     0.0000 \n"
"    H  10    4.9607     4.2961     2.4803     0.0000 \n"
"    H  11    4.2961     4.9607     4.2961     2.4803     0.0000 \n"
"    H  12    2.4803     4.2961     4.9607     4.2961     2.4804     0.0000\n\n"

" ATOMIC CHARGES\n"
"    C   1   -0.1000000000\n"
"    C   2   -0.1000000000\n"
"    C   3   -0.1000000000\n"
"    C   4   -0.1000000000\n"
"    C   5   -0.1000000000\n"
"    C   6   -0.1000000000\n"
"    H   7    0.1000000000\n"
"    H   8    0.1000000000\n"
"    H   9    0.1000000000\n"
"    H  10    0.1000000000\n"
"    H  11    0.1000000000\n"
"    H  12    0.1000000000\n\n"

" BOND ANGLES\n"
"    7    1    2   HC  Car  Car    120.000\n"
"    1    2    3  Car  Car  Car    120.000\n"
"    1    2    8  Car  Car   HC    120.000\n"
"    8    2    3   HC  Car  Car    120.000\n"
"    2    3    4  Car  Car  Car    120.000\n"
"    2    3    9  Car  Car   HC    120.000\n"
"    9    3    4   HC  Car  Car    120.000\n"
"    3    4    5  Car  Car  Car    120.000\n"
"    3    4   10  Car  Car   HC    120.000\n"
"   10    4    5   HC  Car  Car    120.000\n"
"    4    5    6  Car  Car  Car    120.000\n"
"    4    5   11  Car  Car   HC    120.000\n"
"   11    5    6   HC  Car  Car    120.000\n"
"    5    6    1  Car  Car  Car    120.000\n"
"    5    6   12  Car  Car   HC    120.000\n"
"   12    6    1   HC  Car  Car    120.000\n"
"    6    1    2  Car  Car  Car    120.000\n"
"    6    1    7  Car  Car   HC    120.000\n"
"    2    1    7  Car  Car   HC    120.000\n"
"    3    2    8  Car  Car   HC    120.000\n"
"    4    3    9  Car  Car   HC    120.000\n"
"    5    4   10  Car  Car   HC    120.000\n"
"    6    5   11  Car  Car   HC    120.000\n"
"    1    6   12  Car  Car   HC    120.000\n\n"

" TORSION ANGLES\n"
"    6    1    2    3      0.026\n"
"    6    1    2    8   -179.974\n"
"    7    1    2    3    179.974\n"
"    7    1    2    8     -0.026\n"
"    1    2    3    4     -0.026\n"
"    1    2    3    9   -179.974\n"
"    8    2    3    4    179.974\n"
"    8    2    3    9      0.026\n"
"    2    3    4    5      0.026\n"
"    2    3    4   10    179.974\n"
"    9    3    4    5    179.974\n"
"    9    3    4   10     -0.026\n"
"    3    4    5    6     -0.026\n"
"    3    4    5   11    179.974\n"
"   10    4    5    6   -179.974\n"
"   10    4    5   11      0.026\n"
"    4    5    6    1      0.026\n"
"    4    5    6   12    179.974\n"
"   11    5    6    1   -179.974\n"
"   11    5    6   12     -0.026\n"
"    5    6    1    2     -0.026\n"
"    5    6    1    7   -179.974\n"
"   12    6    1    2    179.974\n"
"   12    6    1    7      0.026\n\n"

".. seealso::\n\n"

"  :ref:`Open_Babel_molecule_report`\n\n"
        ;
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    static void WriteCharges(ostream &ofs,OBMol &mol);
    static void WriteDistanceMatrix(ostream &ofs,OBMol &mol);
    static void WriteTorsions(ostream &ofs,OBMol &mol);
    static void WriteAngles(ostream &ofs,OBMol &mol);
    static void WriteChiral(ostream &ofs,OBMol &mol);

  };

  //Make an instance of the format class
  ReportFormat theReportFormat;

  ////////////////////////////////////////////////////////////////

  bool ReportFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    ofs << "FILENAME: " << mol.GetTitle() << "\n";
    ofs << "FORMULA: " << mol.GetFormula() << "\n";
    ofs << "MASS: ";
    snprintf(buffer, BUFF_SIZE, "%5.4f\n", mol.GetMolWt());
    ofs << buffer;
    ofs << "EXACT MASS: ";
    snprintf(buffer, BUFF_SIZE, "%5.7f", mol.GetExactMass());
    ofs << buffer << "\n";
    if (mol.GetTotalCharge() != 0)
      {
        ofs << "TOTAL CHARGE: ";
        snprintf(buffer, BUFF_SIZE, "%d", mol.GetTotalCharge());
        ofs << buffer << "\n";
      }
    if (mol.GetTotalSpinMultiplicity() != 1)
      {
        ofs << "TOTAL SPIN: ";
        snprintf(buffer, BUFF_SIZE, "%d", mol.GetTotalSpinMultiplicity());
        ofs << buffer << "\n";
      }
    ofs << "INTERATOMIC DISTANCES" << "\n";
    WriteDistanceMatrix(ofs, mol);
    ofs << "\n" << "\n" << "ATOMIC CHARGES" << "\n";
    WriteCharges(ofs, mol);
    ofs << "\n" << "\n" << "BOND ANGLES" << "\n";
    WriteAngles(ofs, mol);
    ofs << "\n" << "\n" << "TORSION ANGLES" << "\n";
    WriteTorsions(ofs, mol);
    if (mol.IsChiral())
      {
        ofs << "\n" << "\n" << "CHIRAL ATOMS" << "\n";
        WriteChiral(ofs, mol);
      }
    if (mol.HasData(OBGenericDataType::CommentData))
      {
        ofs << "\n" << "\n" << "COMMENTS" << "\n";
        OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
        ofs << cd->GetData() << "\n";
      }
    ofs << "\n" << "\n";
    return(true);
  }

  /////////////////////////////////////////////////////////////
  void ReportFormat::WriteCharges(ostream &ofs,OBMol &mol)
  {
    unsigned int i;
    OBAtom *atom;
    char buffer[BUFF_SIZE];

    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "%4s%4d   % 2.10f",
                etab.GetSymbol(atom->GetAtomicNum()),
                i,
                atom->GetPartialCharge());

        ofs << buffer << "\n";
      }
  }

  void ReportFormat::WriteDistanceMatrix(ostream &ofs,OBMol &mol)
  {
    int columns = 7;
    unsigned int max, min = 1;
    unsigned int i,j;
    string type;
    OBAtom *atom, *atom2;
    char buffer[BUFF_SIZE];
    double dst;

    max = columns;
    while (max <= mol.NumAtoms() + columns)
      {
        ofs << "\n";
        if (min > mol.NumAtoms())
          break;
        atom = mol.GetAtom(min);

        snprintf(buffer,BUFF_SIZE,"%15s%4d",
                etab.GetSymbol(atom->GetAtomicNum()),
                min);
        ofs << buffer;

        for (i = min + 1; ((i < max) && (i <= mol.NumAtoms())); i++)
          if (i <= mol.NumAtoms())
            {
              atom = mol.GetAtom(i);
              snprintf(buffer,BUFF_SIZE, "%7s%4d",
                      etab.GetSymbol(atom->GetAtomicNum()),
                      i);
              ofs << buffer;
            }
        ofs << "\n";

        snprintf(buffer, BUFF_SIZE, "%14s","");
        ofs << buffer;
        for (i = min; i < max; i++)
          if (i <= mol.NumAtoms())
            {
              ofs << "-----------";
            }

        ofs << "\n";
        for (i = min; i <= mol.NumAtoms(); i++)
          {
            atom = mol.GetAtom(i);
            snprintf(buffer, BUFF_SIZE, "%4s%4d",
                    etab.GetSymbol(atom->GetAtomicNum()),
                    i);
            ofs << buffer;
            for (j = min; j < max; j++)
              if (j <= i)
                {
                  atom2 = mol.GetAtom(j);
                  dst = SQUARE(atom->GetX() - atom2->GetX());
                  dst += SQUARE(atom->GetY() - atom2->GetY());
                  dst += SQUARE(atom->GetZ() - atom2->GetZ());
                  dst = sqrt(dst);
                  snprintf(buffer, BUFF_SIZE, "%10.4f ",dst);
                  ofs << buffer;
                }
            ofs << "\n";
          }
        max += columns - 1;
        min += columns - 1;
      }
    ofs << "\n";
  }

  void ReportFormat::WriteTorsions(ostream &ofs,OBMol &mol)
  {
    vector<OBBond*>::iterator bi1,bi2,bi3;
    OBBond* bond;
    OBAtom *a,*b,*c,*d;
    char buffer[BUFF_SIZE];

    //loop through all bonds generating torsions
    for(bond = mol.BeginBond(bi1); bond; bond = mol.NextBond(bi1))
      {
        b = bond->GetBeginAtom();
        c = bond->GetEndAtom();

        for(a = b->BeginNbrAtom(bi2);a;a = b->NextNbrAtom(bi2))
          {
            if(a == c)
              continue;

            for(d = c->BeginNbrAtom(bi3);d;d = c->NextNbrAtom(bi3))
              {
                if(d == b)
                  continue;

                snprintf(buffer, BUFF_SIZE, "%4d %4d %4d %4d %10.3f",
                        a->GetIdx(), b->GetIdx(),c->GetIdx(),d->GetIdx(),
                        CalcTorsionAngle(a->GetVector(), b->GetVector(),
                                         c->GetVector(), d->GetVector()));
                ofs << buffer << "\n";
              }
          }
      }
  }

  void ReportFormat::WriteAngles(ostream &ofs,OBMol &mol)
  {
    // Alas, we still need to sort these to only list unique entries...
    vector3 v1, v2;
    OBAtom *a, *b, *c, *d;
    OBBond *bond1, *bond2, *bond3;
    vector<OBBond*>::iterator i, j, k;
    char buffer[BUFF_SIZE];

    for (bond1 = mol.BeginBond(i); bond1; bond1 = mol.NextBond(i))
      {
        b = bond1->GetBeginAtom();
        c = bond1->GetEndAtom();

        for (bond2 = b->BeginBond(j); bond2; bond2 = b->NextBond(j))
          {
            if (bond2->GetEndAtomIdx() != c->GetIdx()
                && bond2->GetEndAtomIdx() != b->GetIdx())
              {
                a = bond2->GetEndAtom();

                v1 = a->GetVector() - b->GetVector();
                v2 = c->GetVector() - b->GetVector();

                snprintf(buffer, BUFF_SIZE, "%4d %4d %4d %4s %4s %4s %10.3f",
                        a->GetIdx(),b->GetIdx(),c->GetIdx(),
                        a->GetType(),b->GetType(),c->GetType(),
                        vectorAngle(v1, v2));
                ofs << buffer << "\n";

                for (bond3 = c->BeginBond(k); bond3; bond3 = c->NextBond(k))
                  if (bond3->GetEndAtomIdx() != b->GetIdx()
                      && bond3->GetEndAtomIdx() != c->GetIdx())
                    {
                      d = bond3->GetEndAtom();

                      v1 = b->GetVector() - c->GetVector();
                      v2 = d->GetVector() - c->GetVector();

                      snprintf(buffer, BUFF_SIZE, "%4d %4d %4d %4s %4s %4s %10.3f",
                              b->GetIdx(),c->GetIdx(),d->GetIdx(),
                              b->GetType(),c->GetType(),d->GetType(),
                              vectorAngle(v1, v2));
                      ofs << buffer << "\n";
                    }
              }
          }
      }
  }

  void ReportFormat::WriteChiral(ostream &ofs,OBMol &mol)
  {
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    char buffer[BUFF_SIZE];

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if (atom->IsChiral())
          {
            /* @todo
            snprintf(buffer, BUFF_SIZE, "%4s %5d is chiral: %s",
                    etab.GetSymbol(atom->GetAtomicNum()),
                    atom->GetIdx(),
                    (atom->IsClockwise() ? "clockwise" : "counterclockwise"));
                    */

            ofs << buffer << "\n";
          }
      }
  }

} //namespace OpenBabel
