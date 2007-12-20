/**********************************************************************
forcefield.h - Handle OBForceField class.
 
Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#ifndef OB_FORCEFIELD_H
#define OB_FORCEFIELD_H

#include <vector>
#include <string>
#include <map>

#include <list>
#include <set>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/plugin.h>

#ifndef OBFPRT
#define OBFPRT
#endif

namespace OpenBabel
{
  // log levels
#define OBFF_LOGLVL_NONE	0   //!< no output
#define OBFF_LOGLVL_LOW		1   //!< SteepestDescent progress... (no output from Energy())
#define OBFF_LOGLVL_MEDIUM	2   //!< individual energy terms
#define OBFF_LOGLVL_HIGH	3   //!< individual calculations and parameters

  // terms
#define OBFF_ENERGY		(1 << 0)   //!< all terms
#define OBFF_EBOND		(1 << 1)   //!< bond term
#define OBFF_EANGLE		(1 << 2)   //!< angle term
#define OBFF_ESTRBND		(1 << 3)   //!< strbnd term
#define OBFF_ETORSION		(1 << 4)   //!< torsion term
#define OBFF_EOOP		(1 << 5)   //!< oop term
#define OBFF_EVDW		(1 << 6)   //!< vdw term
#define OBFF_EELECTROSTATIC	(1 << 7)   //!< electrostatic term

  // constraint types
#define OBFF_CONST_IGNORE	(1 << 0)   //!< ignore the atom while setting up calculations
#define OBFF_CONST_ATOM		(1 << 1)   //!< fix the atom position
#define OBFF_CONST_ATOM_X	(1 << 2)   //!< fix the x coordinate of the atom position
#define OBFF_CONST_ATOM_Y	(1 << 3)   //!< fix the y coordinate of the atom position
#define OBFF_CONST_ATOM_Z	(1 << 4)   //!< fix the z coordinate of the atom position
#define OBFF_CONST_DISTANCE	(1 << 5)   //!< constrain distance length
#define OBFF_CONST_ANGLE	(1 << 6)   //!< constrain angle
#define OBFF_CONST_TORSION	(1 << 7)   //!< constrain torsion

  // mode arguments for SteepestDescent, ConjugateGradients, ...
#define OBFF_NUMERICAL_GRADIENT   (1 << 0)  //!< use numerical gradients
#define OBFF_ANALYTICAL_GRADIENT	(1 << 1)  //!< use analytical gradients

#define KCAL_TO_KJ	4.1868

  // inline if statements for logging.
#define IF_OBFF_LOGLVL_LOW    if(_loglvl >= OBFF_LOGLVL_LOW)
#define IF_OBFF_LOGLVL_MEDIUM if(_loglvl >= OBFF_LOGLVL_MEDIUM)
#define IF_OBFF_LOGLVL_HIGH   if(_loglvl >= OBFF_LOGLVL_HIGH)
  
  //! \class OBFFParameter forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold forcefield parameters
  class OBFPRT OBFFParameter {
  public:
    //! Used to store integer atom types
    int         a, b, c, d;
    //! used to store string atom types
    std::string _a, _b, _c, _d; 

    //! Used to store integer type parameters (bondtypes, multiplicity, ...)
    std::vector<int>    _ipar;
    //! Used to store double type parameters (force constants, bond lengths, angles, ...)
    std::vector<double> _dpar;

    //! Assignment 
    OBFFParameter& operator=(const OBFFParameter &ai) 
    {
      if (this != &ai) {
        a = ai.a;
        b = ai.b;
        c = ai.c;
        d = ai.d;
        _a = ai._a;
        _b = ai._b;
        _c = ai._c;
        _d = ai._d;
        _ipar = ai._ipar;
        _dpar = ai._dpar;
      }
        
      return *this;
    }

    //! Reset the atom types and set all parameters to zero
    void clear () 
    {
      a = b = c = d = 0;
      _ipar.clear();
      _dpar.clear();
    }
  }; // class OBFFParameter
  
  // specific class introductions in forcefieldYYYY.cpp (for YYYY calculations)
  //! \class OBFFCalculation forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold energy and gradient calculations on specific force fields
  class OBFPRT OBFFCalculation
  {
    public:
      //! Used to store the energy for this OBFFCalculation
      double energy;
      //! Used to store the gradients for this OBFFCalculation
      vector3 grada, gradb, gradc, gradd;
      //! Used to store the atoms for this OBFFCalculation
      OBAtom *a, *b, *c, *d;

      //! Constructor
      OBFFCalculation() 
        {
	  a = NULL;
	  b = NULL;
	  c = NULL;
	  d = NULL;
	  energy = 0.0;
          grada = VZero;
          gradb = VZero;
          gradc = VZero;
          gradd = VZero;
        }
      //! Destructor
      virtual ~OBFFCalculation()
        {
        }
      
      //! Compute the energy and gradients for this OBFFCalculation
      virtual void Compute(bool gradients = true) 
        {
        }
      //! \return Energy for this OBFFCalculation (call Compute() first)
      virtual double GetEnergy() 
      {
        if (!energy)
	  Compute(false);

        return energy; 
      }
      //! \return Gradient for this OBFFCalculation with respect to coordinates of atom (call Compute() first)
      virtual vector3 GetGradient(OBAtom *atom) 
      {
        if (atom == a)
          return grada;
        else if (atom == b)
          return gradb;
        else if (atom == c)
          return gradc;
        else if (atom == d)
          return gradd;
        else 
          return  VZero;
      }
  };
  
  //! \class OBFFConstraint forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to hold constraints
  class OBFPRT OBFFConstraint
  {
    public:
      //! Used to store the contraint energy for this OBFFConstraint
      double factor, constraint_value;
      double rab0, rbc0;
      //! Used to store the contraint type for this OBFFConstraint
      int type, ia, ib, ic, id;
      //! Used to store the atoms for this OBFFCostraint
      OBAtom *a, *b, *c, *d;
      //! Used to store the gradients for this OBFFCalculation
      vector3 grada, gradb, gradc, gradd;

      //! Constructor
      OBFFConstraint() 
      {
	a = b = c = d = NULL;
	ia = ib = ic = id = 0;
        constraint_value = 0.0;
	factor = 0.0;
      }
      //! Destructor
      ~OBFFConstraint()
      {
      }
      
      vector3 GetGradient(int a) 
      {
        if (a == ia)
          return grada;
        else if (a == ib)
          return gradb;
        else if (a == ic)
          return gradc;
        else if (a == id)
          return gradd;
        else 
          return  VZero;
      }
  };

  //! \class OBFFConstraints forcefield.h <openbabel/forcefield.h>
  //! \brief Internal class for OBForceField to handle constraints
  class OBFPRT OBFFConstraints
  {
    public:
      //! Constructor
      OBFFConstraints();
      //! Destructor
      ~OBFFConstraints()
      {
        _constraints.clear();
      }
      //! Clear all constraints
      void Clear();
      //! Get the constraint energy
      double GetConstraintEnergy();
      //! Get the constraint gradient for atom with index a
      vector3 GetGradient(int a);
      //! Get the constrain gradient for the atom
      OBFFConstraints& operator=(const OBFFConstraints &ai) 
      {
        if (this != &ai) {
          _constraints = ai._constraints;
        }
        return *this;
      }

      /*! Translate indices to OBAtom* objects, this function is called from OBForceField::Setup,
       *  this function doesn't have to be called from anywhere else.
       */
      void Setup(OBMol &mol);

      /////////////////////////////////////////////////////////////////////////
      // Set Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to set constraints
      //@{
      //! Set Constraint factor
      void SetFactor(double factor);
      //! Ignore the atom while setting up calculations
      void AddIgnore(int a);
      //! Fix the position of an atom
      void AddAtomConstraint(int a);
      //! Fix the x coordinate of the atom position
      void AddAtomXConstraint(int a);
      //! Fix the y coordinate of the atom position
      void AddAtomYConstraint(int a);
      //! Fix the z coordinate of the atom position
      void AddAtomZConstraint(int a);
      //! Constrain the bond length a-b
      void AddBondConstraint(int a, int b, double length);
      //! Constrain the angle a-b-c
      void AddAngleConstraint(int a, int b, int c, double angle);
      //! Constrain the torsion angle a-b-c-d
      void AddTorsionConstraint(int a, int b, int c, int d, double torsion);
      //! Delete a constraint
      //! \par index constraint index
      void DeleteConstraint(int index);
      //@}
      /////////////////////////////////////////////////////////////////////////
      // Get Constraints                                                     //
      /////////////////////////////////////////////////////////////////////////
      //! \name Methods to get information about set constraints
      //@{
      //! Get Constraint factor
      double GetFactor();
      //! \returns the number of set constraints
      int Size() const;
      /*! The following constraint types are known: OBFF_CONST_IGNORE (ignore 
       *  the atom while setting up calculations, forcefield implementations 
       *  need to check this value in their setup function), OBFF_CONST_ATOM
       *  (fix atom position), OBFF_CONST_ATOM_X (fix x coordinate), 
       *  OBFF_CONST_ATOM_Y (fix y coordinate), OBFF_CONST_ATOM_Z (fix z 
       *  coordinate), OBFF_CONST_BOND (constrain bond length), OBFF_CONST_ANGLE
       *  (constrain angle), OBFF_CONST_TORSION (constrain torsion angle)
       *  \return the constraint type
       */
      int GetConstraintType(int index) const;
      /*! \return The constraint value, this can be a bond length, angle or 
       *   torsion angle depending on the constraint type.
       */
      double GetConstraintValue(int index) const;
      //! \return The constraint atom a (or fixed atom)
      //! \par index constraint index
      int GetConstraintAtomA(int index) const;
      //! \return The constraint atom b
      //! \par index constraint index
      int GetConstraintAtomB(int index) const;
      //! \return The constraint atom c
      //! \par index constraint index
      int GetConstraintAtomC(int index) const;
      //! \return The constraint atom d
      //! \par index constraint index
      int GetConstraintAtomD(int index) const;
      //! \return true if this atom is ignored
      //! \par a atom index
      bool IsIgnored(int a);
      //! \return true if this atom is fixed
      //! \par a atom index
      bool IsFixed(int a);
      //! \return true if the x coordinate for this atom is fixed
      //! \par a atom index
      bool IsXFixed(int a);
      //! \return true if the y coordinate for this atom is fixed
      //! \par a atom index
      bool IsYFixed(int a);
      //! \return true if the z coordinate for this atom is fixed
      //! \par a atom index
      bool IsZFixed(int a);
      //@}
 
    private:
      std::vector<OBFFConstraint> _constraints;
      double _factor;
  };
 
  // Class OBForceField
  // class introduction in forcefield.cpp
  class OBFPRT  OBForceField : public OBPlugin
  {
  
    MAKE_PLUGIN(OBForceField)
  
    public:
    //!Clone the current instance. May be desirable in multithreaded environments,
    //!Should be deleted after use
    virtual OBForceField* MakeNewInstance()=0;

    protected:

    /*! 
      Get the correct OBFFParameter from a OBFFParameter vector.
       
      \code vector<OBFFParameter> parameters; \endcode
      
      this vector is filled with entries (as OBFFParameter) from 
      a parameter file. This happens in the Setup() function.
      
      \code GetParameter(a, 0, 0, 0, parameters); \endcode
        
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a (pa = parameter.a)
      
      use: vdw parameters, ...
      
      \code GetParameter(a, b, 0, 0, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b      (ab)
      or: pa = b & pb = a      (ba)
            
      use: bond parameters, vdw parameters (pairs), ...
      
      \code GetParameter(a, b, c, 0, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b & pc = c     (abc)
      or: pa = c & pb = b & pc = a     (cba)
      
      use: angle parameters, ...
      
      \code GetParameter(a, b, c, d, parameters); \endcode
      
      returns the first OBFFParameter from vector<OBFFParameter> 
      parameters where: pa = a & pb = b & pc = c & pd = d    (abcd)
      or: pa = d & pb = b & pc = c & pd = a    (dbca)
      or: pa = a & pb = c & pc = b & pd = d    (acbd)
      or: pa = d & pb = c & pc = b & pd = a    (dcba)
      
      use: torsion parameters, ...
    */
    OBFFParameter* GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
    //! see GetParameter(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
    OBFFParameter* GetParameter(const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);
    //! Get index for vector<OBFFParameter> ...
    int GetParameterIdx(int a, int b, int c, int d, std::vector<OBFFParameter> &parameter);
           
    //! Calculate the potential energy function derivative numerically with repect to the coordinates of atom with index a (this vector is the gradient)
    /*!
     * \param a  provides coordinates
     * \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     * \return the negative gradient of atom a
     */
    vector3 NumericalDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    //! OB 3.0
    vector3 NumericalSecondDerivative(OBAtom *a, int terms = OBFF_ENERGY);
    /*! Calculate the potential energy function derivative analyticaly with repect to the coordinates of atom with index a (this vector is the gradient)
     *
     *  If the currently selected forcefield doesn't have analytical gradients, 
     *  we can still call this function which will return the result of 
     *  NumericalDerivative()
     *  \param a  provides coordinates
     *  \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
     *  \return the negative gradient of atom a
     */
    virtual vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY) 
    { 
      return -NumericalDerivative(a, terms); 
    }
    /*! Check if two atoms are in the same ring. [NOTE: this function uses SSSR, 
     *  this means that not all rings are found for bridged rings. This causes 
     *  some problems with the MMFF94 validation.]
     *  \param a atom a
     *  \param b atom b
     *  \return true if atom a and b are in the same ring
     */
    bool IsInSameRing(OBAtom* a, OBAtom* b);
 
    OBMol _mol; //!< Molecule to be evaluated or minimized
    bool _init; //!< Used to make sure we only parse the parameter file once, when needed

    OBFFConstraints _constraints; //!< Constraints

    std::ostream* _logos; //! Output for logfile
    char _logbuf[BUFF_SIZE]; //!< Temporary buffer for logfile output
    int _loglvl; //!< Log level for output
    int _origLogLevel;
    
    
    int _current_conformer; //! used to hold i for current conformer (needed by UpdateConformers)
    std::vector<double> _energies; //! used to hold the energies for all conformers

    double _econv, _e_n1; //! Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    int _method, _cstep, _nsteps; //! Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    double *_grad1, *_dir1; //! Used for conjugate gradients and steepest descent(Initialize and TakeNSteps)
    int _ncoords; //!< Number of coordinates for conjugate gradients

    double _timestep; //! Molecular dynamics time step in picoseconds
    double _T; //! Molecular dynamics temperature in Kelvin
    double *_velocityPtr; //! pointer to the velocities
  
  public:
    //! Destructor
    virtual ~OBForceField()
      {
        if (_grad1 != NULL) {
          delete [] _grad1;
          _grad1 = NULL;
        }
        if (_dir1 != NULL) {
          delete [] _dir1;
          _dir1 = NULL;
        }
      }
    const char* TypeID()
      {
        return "forcefields";
      }

    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const std::string& ID)
    { 
      return FindType(ID.c_str());
    } 
    /*! \param ID forcefield id (Ghemical, ...)
     *  \return A pointer to a forcefield (the default if ID is empty), or NULL if not available
     */
    static OBForceField* FindForceField(const char *ID)
    {
      return FindType(ID);
    }
    //! \return The unit (kcal/mol, kJ/mol, ...) in which the energy is expressed as std::string
    virtual std::string GetUnit() { return std::string("au"); }
    /* Does this force field have analytical gradients defined for all
     * calculation components (bonds, angles, non-bonded, etc.)
     * If this is true, code should default to using OBFF_ANALYTICAL_GRADIENT
     * for SteepestDescent() or ConjugateGradients()
     * \return true if all analytical gradients are implemented.
     */
    virtual bool HasAnalyticalGradients() { return false; }
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.). Keep current constraints.
     *  \param mol the OBMol object that contains the atoms and bonds
     *  \return True if succesfull
     */
    bool Setup(OBMol &mol); 
    /*! Setup the forcefield for mol (assigns atom types, charges, etc.). Use new constraints. 
     *  \param mol the OBMol object that contains the atoms and bonds
     *  \param constraints the OBFFConstraints object that contains the constraints
     *  \return True if succesfull
     */
    bool Setup(OBMol &mol, OBFFConstraints &constraints);
    /*! Load the parameters (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup())
     */
    virtual bool ParseParamFile() { return false; }
    /*! Set the atom types (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup())
     */
    virtual bool SetTypes() { return false; }
    /*! Set the formal charges (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup())
     */
    virtual bool SetFormalCharges() { return false; }
    /*! Set the partial charges (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup())
     */
    virtual bool SetPartialCharges() { return false; }
    /*! Setup the calculations (this function is overloaded by the individual forcefields,
     *  and is called autoamically from OBForceField::Setup())
     */
    virtual bool SetupCalculations() { return false; }
    /*! Compare the internal forcefield OBMol object to mol. If the two have the
     *  same number of atoms and bonds, and all atomic numbers are the same, 
     *  this function returns false, and no call to Setup is needed.
     *  \return true if Setup needs to be called
     */
    bool IsSetupNeeded(OBMol &mol);
    /*! Get coordinates for current conformer
     *  \param mol the OBMol object to copy the coordinates to (from OBForceField::_mol)
     *  \return true if succesfull
     */
    bool GetCoordinates(OBMol &mol);
    bool UpdateCoordinates(OBMol &mol) {return GetCoordinates(mol); } // = GetCoordinates, depricated
    /*! Get coordinates for all conformers
     *  \param mol the OBMol object to copy the coordinates to (from OBForceField::_mol)
     *  \return true if succesfull
     */
    bool GetConformers(OBMol &mol);
    bool UpdateConformers(OBMol &mol) { return GetConformers(mol); } // = GetConformers, depricated
    /*! Set coordinates for current conformer
     *  \param mol the OBMol object to copy the coordinates from (to OBForceField::_mol)
     *  \return true if succesfull
     */
    bool SetCoordinates(OBMol &mol);
    /*! Set coordinates for all conformers
     *  \param mol the OBMol object to copy the coordinates from (to OBForceField::_mol)
     *  \return true if succesfull
     */
    bool SetConformers(OBMol &mol);
    
    /////////////////////////////////////////////////////////////////////////
    // Energy Evaluation                                                   //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for energy evaluation
    //@{
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Total energy
     *   \par Output to log:
     *    OBFF_LOGLVL_NONE:   none \n
     *    OBFF_LOGLVL_LOW:    none \n
     *    OBFF_LOGLVL_MEDIUM: energy for indivudual energy terms \n
     *    OBFF_LOGLVL_HIGH:   energy for individual energy interactions \n
     */
    virtual double Energy(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Bond stretching energy
     *   \par Output to log:
     *    see Energy()
     */
    virtual double E_Bond(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Angle bending energy
     *  \par Output to log:
     *   see Energy()
     */
    virtual double E_Angle(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Stretch bending energy
     *   \par Output to log:
     *    see Energy()
     */ 
    virtual double E_StrBnd(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Torsional energy
     *    \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Torsion(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Out-Of-Plane bending energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_OOP(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Van der Waals energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_VDW(bool gradients = true) { return 0.0f; }
    /*! \param gradients Set to true when the gradients need to be calculated (needs to be done before calling GetGradient())
     *  \return Electrostatic energy
     *   \par Output to log:
     *	  see Energy()
     */ 
    virtual double E_Electrostatic(bool gradients = true) { return 0.0f; }
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for logging
    //@{
    //! Print the atom types 
    void PrintTypes();
    //! Print the formal charges (atom.GetPartialCharge(), MMFF94 FC's are not always int) 
    void PrintFormalCharges();
    //! Print the partial charges
    void PrintPartialCharges();
    //! Print the velocities
    void PrintVelocities();
    /*! Set the stream for logging (can also be &cout for logging to screen)
     *  \param pos stream
     *  \return True if succesfull
     */
    bool SetLogFile(std::ostream *pos);
    /*!
      Set the log level (OBFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH)

      Inline if statements for logging are available: 
	
      \code
      #define IF_OBFF_LOGLVL_LOW    if(_loglvl >= OBFF_LOGLVL_LOW)
      #define IF_OBFF_LOGLVL_MEDIUM if(_loglvl >= OBFF_LOGLVL_MEDIUM)
      #define IF_OBFF_LOGLVL_HIGH   if(_loglvl >= OBFF_LOGLVL_HIGH)
      \endcode

      example:
      \code
      SetLogLevel(OBFF_LOGLVL_MEDIUM);
  
      IF_OBFF_LOGLVL_HIGH {
        OBFFLog("this text will NOT be logged...\n");
      }
   
      IF_OBFF_LOGLVL_LOW {
        OBFFLog"this text will be logged...\n");
      }
  
      IF_OBFF_LOGLVL_MEDIUM {
        OBFFLog("this text will also be logged...\n");
      }
      \endcode
    */
    bool SetLogLevel(int level);
    //! \return log level
    int GetLogLevel() { return _loglvl; }
    /*! Print msg to the logfile
     *  \param msg the message
     */
    void OBFFLog(std::string msg)
    {
      if (!_logos)
        return;
      
      *_logos << msg;
    }
    /*! Print msg to the logfile
     *  \param msg the message
     */
    void OBFFLog(const char *msg)
    {
      if (!_logos)
        return;
      
      *_logos << msg;
    }
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Structure Generation                                                //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for structure generation
    //@{
    //! Generate coordinates for the molecule (distance geometry). (OB 3.0)
    void DistanceGeometry();
    /*! Generate conformers for the molecule (systematicaly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. SystematicRotorSearch works by rotating around
     *  the rotatable bond in a molecule (see OBRotamerList class). This rotating generates 
     *  multiple conformers. The energy for all these conformers is then evaluated and the 
     *  lowest energy conformer is selected.
     *
     * \param geomSteps the number of steps to take during geometry optimization
     *        
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n 
     */
    void SystematicRotorSearch(unsigned int geomSteps = 2500);
    /*! Generate conformers for the molecule by systematicaly rotating torsions. To be used in combination with 
     *  SystematicRotorSearchNexConformer().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class 
     *  pFF->SystematicRotorSearchInitialize(300);
     *  while (pFF->SystematicRotorSearchNextConformer(300)) {
     *    // do some updating in your program (show last generated conformer, ...)
     *  }
     *  \endcode
     * 
     *  If you don't need any updating in your program, SystematicRotorSearch() is recommended.
     *
     *  \param geomSteps the number of steps to take during geometry optimization
     *  \return the number of conformers
     */
    int SystematicRotorSearchInitialize(unsigned int geomSteps = 2500);
    //! \return true if there are more conformers
    bool SystematicRotorSearchNextConformer(unsigned int geomSteps = 2500);
    /*! Generate conformers for the molecule (randomly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. RandomRotorSearch works by randomly rotating around
     *  the rotatable bonds in a molecule (see OBRotamerList class). This rotating generates 
     *  multiple conformers. The energy for all these conformers is then evaluated and the 
     *  lowest energy conformer is selected.
     *
     * \param conformers the number of random conformers to consider during the search
     * \param geomSteps the number of steps to take during geometry optimization for each conformer
     *        
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n 
     */
    void RandomRotorSearch(unsigned int conformers, unsigned int geomSteps = 2500);
    /*! Generate conformers for the molecule by randomly rotating torsions. To be used in combination with 
     *  RandomRotorSearchNexConformer().
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class 
     *  pFF->RandomRotorSearchInitialize(300);
     *  while (pFF->RandomRotorSearchNextConformer(300)) {
     *    // do some updating in your program (show last generated conformer, ...)
     *  }
     *  \endcode
     * 
     *  If you don't need any updating in your program, RandomRotorSearch() is recommended.
     *
     *  \param conformers the number of random conformers to consider during the search
     *  \param geomSteps the number of steps to take during geometry optimization
     *  \return the number of conformers
     */
    void RandomRotorSearchInitialize(unsigned int conformers, unsigned int geomSteps = 2500);
    //! \return true if there are more conformers
    bool RandomRotorSearchNextConformer(unsigned int geomSteps = 2500);
    /*! Generate conformers for the molecule (randomly rotating torsions).
     *  
     *  The initial starting structure here is important, this structure should be
     *  minimized for the best results. WeightedRotorSearch works by randomly rotating around
     *  the rotatable bonds in a molecule (see OBRotamerList class). Unlike RandomRotorSearch()
     *  the random choice of torsions is reweighted based on the energy of the generated conformer.
     *  Over time, the generated conformers for each step should become increasingly better.
     *  The lowest energy conformer is selected.
     *
     * \param conformers the number of random conformers to consider during the search
     * \param geomSteps the number of steps to take during geometry optimization for each conformer
     *        
     *	\par Output to log:
     *  This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
     *	too much information about the energy calculations needed for this function will interfere with the output for
     *	this function. \n\n
     *  OBFF_LOGLVL_NONE:   none \n
     *  OBFF_LOGLVL_LOW:    number of rotatable bonds, energies for the conformers, which one is the lowest, ... \n
     *  OBFF_LOGLVL_MEDIUM: see note above \n
     *  OBFF_LOGLVL_HIGH:   see note above \n 
     */
    void WeightedRotorSearch(unsigned int conformers, unsigned int geomSteps);

    /////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for energy minimization
    //@{
    /*! Perform a linesearch starting at atom in direction direction
      \deprecated Current code should use LineSearch(double *, double*) instead
      \param atom start coordinates
      \param direction the search direction

      \return vector which starts at atom and stops at the minimum (the length is the ideal stepsize)

      \par Output to log:
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    none \n
        OBFF_LOGLVL_MEDIUM: none \n
        OBFF_LOGLVL_HIGH:   none \n
    */
    vector3 LineSearch(OBAtom *atom, vector3 &direction);

    /*! Perform a linesearch for the entire molecule in direction direction
      \param currentCoords start coordinates
      \param direction the search direction

      \return alpha, the scale of the step we moved along the direction vector

      \par Output to log:
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    none \n
        OBFF_LOGLVL_MEDIUM: none \n
        OBFF_LOGLVL_HIGH:   none \n
    */
    double LineSearch(double *currentCoords, double *direction);

    /*! Perform steepest descent optimalization for steps steps or until convergence criteria is reached.
      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
        OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps and first step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n 
    */
    void SteepestDescent(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize steepest descent optimalization, to be used in combination with SteepestDescentTakeNSteps().
      
      example:
      \code
      // pFF is a pointer to a OBForceField class 
      pFF->SteepestDescentInitialize(100, 1e-5f);
      while (pFF->SteepestDescentTakeNSteps(5)) {
        // do some updating in your program (redraw structure, ...)
      }
      \endcode
      
      If you don't need any updating in your program, SteepestDescent() is recommended.

      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    void SteepestDescentInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a steepestdescent optimalization that was previously initialized with SteepestDescentInitialize().
      \param n the number of steps to take
      
      \return false if convergence or the number of steps given by SteepestDescentInitialize() has been reached 
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    bool SteepestDescentTakeNSteps(int n);
    /*! Perform conjugate gradient optimalization for steps steps or until convergence criteria is reached.
      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS

      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    information about the progress of the minimization \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
 
    */
    void ConjugateGradients(int steps, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Initialize conjugate gradient optimalization and take the first step, to be used in combination with ConjugateGradientsTakeNSteps().
      
      example:
      \code
      // pFF is a pointer to a OBForceField class 
      pFF->ConjugateGradientsInitialize(100, 1e-5f);
      while (pFF->ConjugateGradientsTakeNSteps(5)) {
        // do some updating in your program (redraw structure, ...)
      }
      \endcode
      
      If you don't need any updating in your program, ConjugateGradients() is recommended.

      \param steps the number of steps 
      \param econv energy convergence criteria (defualt is 1e-6)
      \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    header including number of steps and first step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
    */
    void ConjugateGradientsInitialize(int steps = 1000, double econv = 1e-6f, int method = OBFF_ANALYTICAL_GRADIENT);
    /*! Take n steps in a conjugate gradient optimalization that was previously initialized with ConjugateGradientsInitialize().
      \param n the number of steps to take
      
      \return false if convergence or the number of steps given by ConjugateGradientsInitialize() has been reached 
      
      \par Output to log:
        This function should only be called with the log level set to OBFF_LOGLVL_NONE or OBFF_LOGLVL_LOW. Otherwise
	too much information about the energy calculations needed for the minimization will interfere with the list 
	of energies for succesive steps. \n\n
	OBFF_LOGLVL_NONE:   none \n
        OBFF_LOGLVL_LOW:    step number, energy and energy for the previous step \n
        OBFF_LOGLVL_MEDIUM: see note above \n
        OBFF_LOGLVL_HIGH:   see note above \n
    */
    bool ConjugateGradientsTakeNSteps(int n);
    //@}
    
    /////////////////////////////////////////////////////////////////////////
    // Molecular Dynamics                                                  //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for molecular dynamics
    //@{
    /*! Generate starting velocities with a Maxwellian distribution.
     */
    void GenerateVelocities();
    /*! Correct the velocities so that the following is true:
     *  
     *  \code
     *        3N
     *       ----
     *  0.5  \    m_i * v_i^2 = 0.5 * Ndf * kB * T = E_kin
     *       /
     *       ----
     *       i=1
     *  
     *  E_kin : kinetic energy
     *  m_i : mass of atom i
     *  v_i : velocity of atom i
     *  Ndf : number of degrees of freedom (3 * number of atoms)
     *  kB : Boltzmann's constant
     *  T : temperature
     *  \endcode
     *  
     */
    void CorrectVelocities();
    /*! Take n steps at tempetature T. If no velocities are set, they will be generated.
     *
     *  example:
     *  \code
     *  // pFF is a pointer to a OBForceField class 
     *  while (pFF->MolecularDynamicsTakeNSteps(5, 300)) {
     *    // do some updating in your program (redraw structure, ...)
     *  }
     * \endcode
     *
     *  \param T absolute temperature in Kelvin
     *  \param timestep the time step in picoseconds (10e-12)
        \param method OBFF_ANALYTICAL_GRADIENTS (default) or OBFF_NUMERICAL_GRADIENTS
     */
    void MolecularDynamicsTakeNSteps(int n, double T, double timestep = 0.001, int method = OBFF_ANALYTICAL_GRADIENT);
    //@}

    /////////////////////////////////////////////////////////////////////////
    // Constraints                                                         //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for constraints
    //@{
    //! Get the constraints 
    OBFFConstraints& GetConstraints() { return _constraints; }
    //! Set the constraints 
    void SetConstraints(OBFFConstraints& constraints) 
    { 
      _constraints = constraints; 
      if (_mol.NumAtoms())
        _constraints.Setup(_mol); 
    }
    //@}

 
    /////////////////////////////////////////////////////////////////////////
    // Validation                                                          //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for forcefield validation
    //@{
    //! (debugging)
    bool DetectExplosion();
    //! (debugging)
    vector3 ValidateLineSearch(OBAtom *atom, vector3 &direction);
    //! (debugging)
    void ValidateSteepestDescent(int steps);
    //! (debugging)
    void ValidateConjugateGradients(int steps);
    //! Validate the force field implementation (debugging)
    virtual bool Validate() { return false; }
    /*! 
      Validate the analytical gradients by comparing them to numerical ones. This function has to
      be implemented force field specific. (debugging)
    */
    virtual bool ValidateGradients() { return false; }
    /*! 
      Calculate the error of the analytical gradient (debugging)
      \return  error = fabs(numgrad - anagrad) / anagrad * 100% 
   */
    vector3 ValidateGradientError(vector3 &numgrad, vector3 &anagrad);
    //@}
     
    /////////////////////////////////////////////////////////////////////////
    // Vector Analysis                                                     //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for vector analysis (used by OBFFXXXXCalculationYYYY)
    //@{
    /*! Calculate the derivative of a vector length. The vector is given by a - b, 
     * the length of this vector rab = sqrt(ab.x^2 + ab.y^2 + ab.z^2).
     * \param a atom a (coordinates), will be changed to -drab/da
     * \param b atom b (coordinates), will be changed to -drab/db
     * \return The distance between a and b (bondlength for bond stretching, separation for vdw, electrostatic)
     */
    static double VectorLengthDerivative(vector3 &a, vector3 &b);
    /*! Calculate the derivative of a angle a-b-c. The angle is given by dot(ab,cb)/rab*rcb.
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \return The angle between a-b-c
     */
    static double VectorAngleDerivative(vector3 &a, vector3 &b, vector3 &c);
    /*! Calculate the derivative of a OOP angle a-b-c-d. b is the central atom, and a-b-c is the plane. 
     * The OOP angle is given by 90° - arccos(dot(corss(ab,cb),db)/rabbc*rdb).
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \param d atom d (coordinates), will be changed to -dtheta/dd
     * \return The OOP angle for a-b-c-d
     */
    static double VectorOOPDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
 
    /*! Calculate the derivative of a torsion angle a-b-c-d. The torsion angle is given by arccos(dot(corss(ab,bc),cross(bc,cd))/rabbc*rbccd).
     * \param a atom a (coordinates), will be changed to -dtheta/da
     * \param b atom b (coordinates), will be changed to -dtheta/db
     * \param c atom c (coordinates), will be changed to -dtheta/dc
     * \param d atom d (coordinates), will be changed to -dtheta/dd
     * \return The tosion angle for a-b-c-d
     */
    static double VectorTorsionDerivative(vector3 &a, vector3 &b, vector3 &c, vector3 &d);
    //@}

  }; // class OBForceField

}// namespace OpenBabel

#endif   // OB_FORCEFIELD_H

//! \file forcefield.h
//! \brief Handle forcefields
