/**********************************************************************
Copyright  (C) 2010 by Jens Thomas
 
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

#include <iomanip>

namespace OpenBabel
{

  class DlpolyInputReader
  {
    /*
     *  Base class for the CONFIG and HISTORY file parsers
     */

  public:

    // Parse the header; levcfg, imcon & unitcell
    bool ParseHeader( std::istream &ifs, OBMol &mol );

    // Parse the unit cell info
    bool ParseUnitCell( std::istream &ifs, OBMol &mol );

    // Read the data associated with a single atom
    bool ReadAtom( std::istream &ifs, OBMol &mol );
    
    /**
     * Converts a string to a numerical type
     * This purloined from: http://www.codeguru.com/forum/showthread.php?t=231054
     */
    template <class T>
    bool from_string(T& t, const std::string& s,
                     std::ios_base& (*f)(std::ios_base&))
    {
      std::istringstream iss(s);
      return !(iss >> f >> t).fail();
    }
    
    int LabelToAtomicNumber(std::string label);
    
    std::stringstream errorMsg;

    char buffer[BUFF_SIZE];
    std::string line; // For convenience so we can refer to lines from the iterator as 'line'
    std::vector<std::string> tokens; // list of lines and list of tokens on a line

    // These are data that are accessed by several subroutines, so we keep them as class variables
    int levcfg,imcon;
    std::vector< vector3 > forces;

    
  }; // End DlpolyInputReader
  
  
  int DlpolyInputReader::LabelToAtomicNumber(std::string label)
  {
    /*
     * Given a string with the label for an atom return the atomic number
     * As we are using the GetAtomicNum function case is not important
     */
    
    // See if the first 2 characters give us a valid atomic number
    int Z=etab.GetAtomicNum(label.substr(0,2).c_str());
    
    // If not try the first one
    if (Z==0) Z=etab.GetAtomicNum(label.substr(0,1).c_str());
    
    if (Z==0){
      // Houston...
      errorMsg << "LabelToAtomicNumber got bad Label: " << label << std::endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    return Z;
    
  } //End LabelToAtomicNumber
  
  bool DlpolyInputReader::ParseHeader( std::istream &ifs, OBMol &mol )
  {

    // Title line
    if ( ! ifs.getline(buffer,BUFF_SIZE) )
      {
        obErrorLog.ThrowError(__FUNCTION__, "Problem reading title line", obWarning);
        return false;
      }  
    mol.BeginModify();
    mol.SetTitle(buffer);
    mol.EndModify();
  
    // levcfg, imcon & poss natms
    if ( ! ifs.getline(buffer,BUFF_SIZE) )
      {
        line=buffer;
        line="Problem reading levcfg line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }
    
    tokenize(tokens, buffer, " \t\n");
    if ( ! tokens.size() >= 2 || ! ( from_string<int>(levcfg, tokens.at(0), std::dec) 
                                     && from_string<int>(imcon, tokens.at(1), std::dec) ) )
      {
        line=buffer;
        line="Problem reading keytrj line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }

    return true;

  } // End ParseHeader

  bool DlpolyInputReader::ParseUnitCell( std::istream &ifs, OBMol &mol )
  {

    bool ok;
    double x,y,z;
    ifs.getline(buffer,BUFF_SIZE);
    tokenize(tokens, buffer, " \t\n");
    ok = from_string<double>(x, tokens.at(0), std::dec);
    ok = from_string<double>(y, tokens.at(1), std::dec);
    ok = from_string<double>(z, tokens.at(2), std::dec);
    vector3 vx = vector3( x, y, z );

    ifs.getline(buffer,BUFF_SIZE);
    tokenize(tokens, buffer, " \t\n");
    ok = from_string<double>(x, tokens.at(0), std::dec);
    ok = from_string<double>(y, tokens.at(1), std::dec);
    ok = from_string<double>(z, tokens.at(2), std::dec);
    vector3 vy = vector3( x, y, z );

    ifs.getline(buffer,BUFF_SIZE);
    tokenize(tokens, buffer, " \t\n");
    ok = from_string<double>(x, tokens.at(0), std::dec);
    ok = from_string<double>(y, tokens.at(1), std::dec);
    ok = from_string<double>(z, tokens.at(2), std::dec);
    vector3 vz = vector3( x, y, z );

    // Add the Unit Cell to the molecule
    OBUnitCell * unitcell = new OBUnitCell();
    unitcell->SetData( vx, vy, vz );
    //std::cout << "Set unit cell " << vx << vy << vz << std::endl;
    mol.BeginModify();
    mol.SetData( unitcell );
    mol.EndModify();

    return true;

  } // End ParseUnitCell


  bool DlpolyInputReader::ReadAtom( std::istream &ifs, OBMol &mol )
  {

    std::string AtomLabel;
    int AtomIndex,AtomicNumber=-1;
    double x,y,z;
    OBAtom *atom;
    bool ok;

    // Line with AtomLabel, AtomIndex & AtomicNumber - only AtomLabel required
    if ( ! ifs.getline(buffer,BUFF_SIZE) ) return false;
    //std::cout << "Got Atom line " << buffer << std::endl;

    tokenize(tokens, buffer, " \t\n");
    AtomLabel = tokens.at(0);
        
    // Currently we ignore atom index as it is optional - assume atoms are in order
    if ( tokens.size() >= 2 ) ok = from_string<int>(AtomIndex, tokens.at(1), std::dec);

    if ( tokens.size() == 3 )
      {
        ok = from_string<int>(AtomicNumber, tokens.at(2), std::dec);
        if ( ! ok ) AtomicNumber=-1;
      }
        
    // Got data so read in  coordinates on next line
    if ( !ifs.getline(buffer,BUFF_SIZE) ) return false;
    tokenize(tokens, buffer, " \t\n");
    ok = from_string<double>(x, tokens.at(0), std::dec);
    ok = from_string<double>(y, tokens.at(1), std::dec);
    ok = from_string<double>(z, tokens.at(2), std::dec);
               
    // Check if we've read in the atomic charge or need to try and extract it from the label
    if ( AtomicNumber == -1 ) {
      AtomicNumber = LabelToAtomicNumber( AtomLabel );
    }

    atom = mol.NewAtom();
    atom->SetAtomicNum(AtomicNumber);
    atom->SetVector(x, y, z); //set coordinates
        
    // Reset Atomic Number
    AtomicNumber=-1;
        
    // levcfg > 0, next line will be velocities - skip for now
    if (levcfg > 0 )
      {
        if ( !ifs.getline(buffer,BUFF_SIZE) ) return false;
      }
    
    // levcfg > 1, next line will be forces
    if ( levcfg > 1 )
      {
        if ( !ifs.getline(buffer,BUFF_SIZE) ) return false;
        tokenize(tokens, buffer, " \t\n");
        ok =  from_string<double>(x, tokens.at(0), std::dec);
        ok =  from_string<double>(y, tokens.at(1), std::dec);
        ok =  from_string<double>(z, tokens.at(2), std::dec);
        forces.push_back( vector3( x,y,z ) );
        //std::cout << "ADDING FORCE " << x << ":" << y << ":" << z << std::endl;
      }

    return true;

  } // End ReadAtom
  
  
  class DlpolyConfigFormat : public OBMoleculeFormat, public DlpolyInputReader
  {
  public:
    //Register this format type ID
    DlpolyConfigFormat()
    {
      OBConversion::RegisterFormat("CONFIG",this);
    }
    
    virtual const char* Description() //required
    {
      return "DL-POLY CONFIG\n";
    };
    
    virtual const char* SpecificationURL()
    { 
      return "http://www.cse.scitech.ac.uk/ccg/software/DL_POLY";
    };
    
    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return WRITEONEONLY;
    };
    
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    
  };
  
  //Make an instance of the format class
  DlpolyConfigFormat theDlpolyConfigFormat;
  
  bool DlpolyConfigFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    
    bool ok;

    // Reset data
    levcfg=0;
    imcon=0;
    forces.clear();

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;
    
    //Define some references so we can use the old parameter names
    std::istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    
    if ( ! ParseHeader( ifs, mol ) ) return false;

    // If imcon > 0 then there are 3 lines with the cell vectors
    if ( imcon > 0 ) ParseUnitCell( ifs, mol );

    mol.BeginModify();
    ok = true;
    while ( ok )
      {
        ok = ReadAtom( ifs, mol );
      } 

    // Add forces as conformer data
    if ( levcfg > 1 && forces.size() )
      {
        OBConformerData * conformer = new OBConformerData();
        std::vector< std::vector< vector3 > > conflist;
        conflist.push_back( forces );
        conformer->SetForces( conflist );
        mol.SetData( conformer );
      }

    mol.EndModify();

    if ( mol.NumAtoms() == 0 )
      return(false);
    else
      return(true);

  } // End ReadMolecule

  
  bool DlpolyConfigFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {

    /**
     * Write a DL-POLY CONFIG file. Ints are 10 chars wide, floats are formatted as 20.15f
     */

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;
    
    //Define some references so we can use the old parameter names
    std::ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    
    // We only print the coordinates and labels
    levcfg=0;
    imcon=0;
    
    int idx=0; // For counting molecule index
    std::string title = mol.GetTitle();

    // 80 char title
    ofs << title.substr(0,80) << std::endl;
    
    // Set levcfg & imcon
    ofs << std::setw(10) << levcfg << std::setw(10) << imcon << std::endl;
    
    // Now loop through molecule
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        
        ofs << std::setw(8) << etab.GetSymbol(atom->GetAtomicNum()) << std::setw(10) << ++idx << std::setw(10) << atom->GetAtomicNum() << std::endl;
        snprintf(buffer, BUFF_SIZE, "%20.15f %20.15f %20.15f\n",
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ()
                 );
        ofs << buffer;
      }

    return(true);
    
  } //End WriteMolecule
  

  class DlpolyHISTORYFormat : public OBMoleculeFormat, public DlpolyInputReader
{
public:
  //Register this format type ID
  DlpolyHISTORYFormat()
  {
    OBConversion::RegisterFormat("HISTORY",this);
  }
  
  virtual const char* Description() //required
  {
    return "DL-POLY HISTORY\n";
  };
  
  virtual const char* SpecificationURL()
  { 
    return "http://www.cse.scitech.ac.uk/ccg/software/DL_POLY";
  };
  
  //Flags() can return be any the following combined by | or be omitted if none apply
  // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
  virtual unsigned int Flags()
  {
    return NOTWRITABLE;
  };
  
  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

};
    
  //Make an instance of the format class
  DlpolyHISTORYFormat theDlpolyHISTORYFormat;

  bool DlpolyHISTORYFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    // gotHeader used so we don't parse the header when called again
    static bool gotHeader=false;
    static std::string origTitle;
    std::string title;
    bool ok;
    int nstep,natms=0;

    levcfg=0;
    imcon=0;
    forces.clear();
  
    OBAtom *atom;
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;
    
    //Define some references so we can use the old parameter names
    std::istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

  
    if ( ! gotHeader )
      {
        if ( ! ParseHeader( ifs, mol ) ) return false;
        gotHeader=true;
        // We get the title from the first line of the file, so we need to remember
        // it for subsequent calls
        origTitle=mol.GetTitle();
      } // End of header parsing


    /*
     * Read the trajectory line - this tells us how many atoms we read in and
     * also the timestep which we use to set the title
     */
    if ( ! ifs.getline(buffer,BUFF_SIZE) ) return false;
    tokenize(tokens, buffer, " \t\n");
    if ( ! tokens.size() >= 6  )
      {
        line=buffer;
        line="Problem reading trajectory line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }

    ok = from_string<int>(nstep, tokens.at(1), std::dec);
    ok = from_string<int>(natms, tokens.at(2), std::dec);
    // It ain't gonna work if we don't have natms
    if ( ! ok  )
      {
        line=buffer;
        line="Problem reading natoms on trajectory line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }
    // Get imcon as it could change?
    ok = from_string<int>(levcfg, tokens.at(3), std::dec);
    ok = from_string<int>(imcon, tokens.at(4), std::dec);

    // Set the title
    title=origTitle + ": timestep=" + tokens.at(1);
    mol.SetTitle( title );

    // If imcon > 0 then there are 3 lines with the cell vectors 
    if ( imcon > 0 ) ParseUnitCell( ifs, mol );

    // Start of coordinates - just loop through reading in data
    int atomsRead=0;
    mol.BeginModify();
    while ( true )
      {
        
        if ( ! ReadAtom( ifs, mol ) ) break;
        atomsRead++;
        if ( atomsRead >= natms) break;

      } // End while reading loop
    
    // Add forces as conformer data
    if ( levcfg > 1 && forces.size() )
      {
        OBConformerData * conformer = new OBConformerData();
        std::vector< std::vector< vector3 > > conflist;
        conflist.push_back( forces );
        conformer->SetForces( conflist );
        mol.SetData( conformer );
      }

    mol.EndModify();
    
    if ( mol.NumAtoms() == 0 )
      return(false);
    else
      return(true);

  } // End ReadMolecule

} //namespace OpenBabel
