/**********************************************************************
  stereo.h - OBStereo & OBStereoBase

  Copyright (C) 2009 by Tim Vandermeersch
 
  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#ifndef OB_STEREO_H
#define OB_STEREO_H

#include <openbabel/base.h> // OBGenericData
#include <vector>
#include <map>
#include <climits> // UINT_MAX

namespace OpenBabel {
  
  ///@addtogroup stereo Stereochemistry 
  ///@{

  /**
   * @brief Placeholder for enums & Ref/Refs related functions.
   *
   * The OBStereo struct contains a number of enums with predefined values. 
   * These are OBStereo::BondDirection, OBStereo::Type, OBStereo::Shape, 
   * OBStereo::View, OBStereo::Winding. There are enums 
   * which only apply to certain types of stereochemistry but having them 
   * in 1 place makes it easier to remember. 
   *
   * The OBStereo struct also contains typedefs and functions which 
   * are crucial to fully understand how to use the OBStereoBase derived 
   * classes (i.e. OBTetrahedralStereo, OBCisTransStereo, ...). Ref variables
   * and Refs lists are a way to uniquely reference atoms in the molecule. In
   * most cases these Ref variables are the same as the unique atom ids 
   * (OBAtom::GetId). However, 2 special cases are provided:
   *
   * - OBStereo::NoRef: An initial value for Ref variables. The constructors of 
   *   the various Config structs set all refs to NoRef (Refs lists will remain 
   *   empty). This value is considered invalid when comparing stereochemistry.
   * - OBStereo::ImplicitRef: Can be used to replace implicit hydrogen ids. 
   *   Even with explicit hydrogens, it is still valid to replace the Ref of 
   *   the hydrogen with ImplicitRef. This flexibility also applies to 
   *   comparing stereochemistry. It is possible to compare stereochemistry
   *   (e.g. OBTetrahedral::Config::operator==) between a Config struct with 
   *   an ImplicitRef and a Config struct where the Ref is the atom id of the 
   *   explicit hydrogen. See actual documentation for details about operator==.
   *
   * There are utility functions which make it easier to handle Refs. The most
   * frequently used one is OBStereo::MakeRefs to create lists containing 3 or
   * 4 Ref values. Formats and library use normally doesn't need the other 
   * functions.
   *
   * @sa OBStereoBase OBStereoFacade
   */
  struct OBAPI OBStereo 
  {
    /**
     * The various types of stereochemistry
     */
    enum Type {
      CisTrans            = (1<<0), //!< cis/trans double bond
      ExtendedCisTrans    = (1<<1), //!< allene, biphenyl, ...
      SquarePlanar        = (1<<2), //!< Square-planar stereochemistry
      Tetrahedral         = (1<<3), //!< tetrahedral
      ExtendedTetrahedral = (1<<4), //!< extended tetrahedral
      TrigonalBipyramidal = (1<<5), //!< Trigonal-bipyramidal stereochemistry
      Octahedral          = (1<<6)  //!< Octahedral stereochemistry
    };

    /**
     * Bond directions used by StereoFrom0D() to translate to
     * internal CisTransStereo representation.
     */
    enum BondDirection { // Values taken from MDL format
      NotStereo =   0,
      UpBond =      1,
      DownBond =    6,
      UnknownDir =  4
    };

    /**
     * Shapes used by OBTetraPlanarStereo subclasses for 
     * setting/getting reference ids.
     *
     * @image html shape.png
     * @sa OBTetraPlanarStereo
     */
    enum Shape {
      ShapeU = 1,
      ShapeZ = 2,
      Shape4 = 3
    };

    /**
     * Views used by OBTetraNonPlanarStereo subclasses for
     * setting/getting reference ids.
     * @sa OBTetraNonPlanarStereo
     */
    enum View
    {
      ViewFrom = 1, //!< view from the atom (id parameter) towards the center atom
      ViewTowards = 2, //!< view from center atom towards the atom (id paramater)
    };

    /**
     * Windings used by OBTetraNonPlanar subclasses for 
     * setting/getting reference ids.
     * @sa OBTetraNonPlanar
     */
    enum Winding {
      Clockwise = 1,     //!< Clockwise winding
      AntiClockwise = 2  //!< AntiClockiwe winding (or CounterClockwise
    };

    ///@name Ref & Refs types 
    //@{
    /**
     * All stereo classes work with variables of the type Ref to uniquely
     * identify atoms. In most cases these Ref variables are the same as the
     * unique atom ids. 
     *
     * @sa OBAtom::GetId() OBStereo
     */
    typedef unsigned long Ref;
    /**
     * Special case Ref values.
     */
    enum {
      NoRef = UINT_MAX,       //!< No Ref set (invalid Ref)
      ImplicitRef = UINT_MAX - 1  //!< Implicit Ref (i.e. hydrogen, N lone pair, ...).
    };
    /**
     * A list (std::vector) of Ref variables.
     */
    typedef std::vector<Ref> Refs;
    /**
     * Iterator (std::iterator) for a Refs list.
     */
    typedef Refs::iterator RefIter;
    /**
     * Iterator (std::iterator) for a const Refs list.
     */ 
    typedef Refs::const_iterator ConstRefIter;
    //@}

    ///@name Refs utility functions
    //@{
    /**
     * Create a Refs list filled with @p ref1, @p ref2, @p ref3 & @p ref4.
     * @p ref4 is not added to the returned Refs if it is equal to NoRef.
     *
     * @return A Refs list containing the specified Ref values. 
     */
    static Refs MakeRefs(Ref ref1, Ref ref2, Ref ref3, Ref ref4 = NoRef)
    {
      Refs refs(3);
      refs[0] = ref1;
      refs[1] = ref2;
      refs[2] = ref3;
      if (ref4 != NoRef)
        refs.push_back(ref4);
      return refs;
    }
    /**
     * Check if @p refs1 and @p refs2 contain the same Ref values regardless 
     * of their order.
     *
     * @code
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(2, 3, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 2, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 4, 1)) // false
     * @endcode
     *
     * @return True if @p refs1 and @p refs2 contain the same Ref values.
     */
    static bool ContainsSameRefs(const Refs &refs1, const Refs &refs2);
    /**
     * @return True if @p refs contains @p ref.
     */
    static bool ContainsRef(const Refs &refs, unsigned long ref);
    //@}

    ///@name Low-level functions used by implementation.
    //@{
    /**
     * Compute the inversion vector for @p refs and return the sum of it's 
     * elements. The ith element in the inversion vector is the number of 
     * element to the right of element i with a lower value.
     *
     * The number of inversions is the same as the number of interchanges
     * of consecutive elements. 
     *
     * When working with 3 refs from a tetrahedral configuration:
     * @code
     * permutation   inversion vector    sum
     * -------------------------------------
     * 123           0 0 0               0 (even) -> clockwise
     * 132           0 1 0               1 (odd)  -> anti-clockwise
     * 213           1 0 0               1 (odd)  -> anti-clockwise
     * 231           1 1 0               2 (even) -> clockwise
     * 312           2 0 0               2 (even) -> clockwise
     * 321           2 1 0               3 (odd)  -> anti-clockwise
     * @endcode
     */
    static int NumInversions(const Refs &refs);
    /**
     * Permutate element @p i with @p j in @p refs.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element i (0...N-1) will be mutated to j and vice versa.
     * @param j Element j (0...N-1) will be mutated to i and vice versa.
     *
     * @note This method does nothing if i equals j.
     */
    static void Permutate(Refs &refs, int i, int j);
    /**
     * Get @p refs with element @p i and @p j permutated.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element @p i (0...N-1) will be mutated to @p j and vice versa.
     * @param j Element @p j (0...N-1) will be mutated to @p i and vice versa.
     *
     * @return @p refs with elements @p i and @p j permutated.
     *
     * @note This method does nothing if @p i equals @p j.
     */
    static Refs Permutated(const Refs &refs, int i, int j);
    //@}
 
  };

  struct StereogenicUnit
  {
    StereogenicUnit() : type(static_cast<OBStereo::Type>(0)), id(OBStereo::NoRef), para(false)
    {
    }

    StereogenicUnit(OBStereo::Type _type, unsigned long _id, bool _para = false) :
        type(_type), id(_id), para(_para)
    {
    }
    
    OBStereo::Type type; //!< the type for this stereogenic unit
    unsigned long id; //! the atom/bond (depends on type) unique id
    bool para; //! para- (=ressemble) or true-stereocenter
  };


  // fwd decl
  class OBMol;
  /**
   * @brief Base class for all stereochemistry classes.
   *
   * All stereochemistry classes are derived from OBStereoBase. This class
   * inherits from OBGenericData which allows the objects to be stored in
   * the molecule. The attribute (OBGenericData::GetAttribute) is set to
   * "StereoData" and the data type is OBGenericDataType::StereoData. The 
   * pure virtual OBStereoBase::GetType function must be implemented by
   * derived classes to return a type defined in OBStereo::Type.
   *
   * Use the OBStereoFacade for easy access to the derived classes. 
   *
   * OBStereoBase keeps track of the OBMol object. This must always be
   * a valid (not 0 or deleted) pointer and can only be set using the 
   * constructor. Subclasses can use this to get more information on bonding 
   * for example. Finally, OBStereoBase also keeps track of the specified
   * flag. By default, this is always set to true.
   *
   * @sa OBStereo OBStereoFacade
   */
  class OBAPI OBStereoBase : public OBGenericData
  {
    public:
      /**
       * Constructor. By default, the stereochemistry is specified. Use 
       * SetSpecified(false) for unspecified/unknown stereochemistry.
       *
       * @param mol The molecule.
       */
      OBStereoBase(OBMol *mol) : 
        OBGenericData("StereoData", OBGenericDataType::StereoData, perceived),
        m_mol(mol), m_specified(true)
      {
      }
      /**
       * Destructor.
       */
      virtual ~OBStereoBase() { m_mol = 0; }

      ///@name Geniric (for all OBStereo::Type) stereochemistry
      //@{
      /**
       * Get the molecule. This can be used by subclasses when more
       * information is needed (e.g. OBCisTransStereo::GetCisRef, ...).
       */
      OBMol* GetMolecule() const { return m_mol; }
      /**
       * Reimplemented by subclasses to return the type defined in OBStereo::Type.
       */
      virtual OBStereo::Type GetType() const = 0;
      /**
       * Set whether the stereochemistry is specified. Comparing a specified 
       * OBStereoBase derived class (or it's Config struct) with an unspecified 
       * one, always returns true.
       */
      void SetSpecified(bool specified) { m_specified = specified; }
      /**
       * @return True if the stereochemistry is specified.
       */
      bool IsSpecified() const { return m_specified; }
      //@}
    private:
      OBMol *m_mol; //!< The parent molecule.
      bool m_specified; //!< True if the stereochemistry is specified, false if unknown/unspecified.
  };

  // fwd decl
  class OBTetrahedralStereo;
  class OBCisTransStereo;
  /**
   * @brief Facade to simplify retrieval of OBStereoBase derived objects.
   *
   * The OBStereoFacade helps with retrieving OBStereoBase derived objects 
   * (i.e. OBTetrahedralStereo, OBCisTransStereo, ...) from an OBMol. This 
   * is done by iterating over all OBGenericData objects with data type 
   * OBGenericDataType::StereoData and checking the OBStereo::Type using 
   * OBStereoBase::GetType.
   *
   * @sa OBStereo OBStereoBase
   */
  class OBAPI OBStereoFacade
  {
    public:
      /**
       * Constructor with @p mol and @p perceive parameter.
       *
       * @param mol The molecule.
       * @param perceive If true, PerceiveStereo will be called if the 
       * OBMol::HasChiralityPerceived() flag is not set. (default is true)
       */
      OBStereoFacade(OBMol *mol, bool perceive = true) : 
          m_mol(mol), m_init(false), m_perceive(perceive)
      {
      }

      ///@name Tetrahedral stereochemistry
      ///@{
      /**
       * Get the number of tetrahedral stereocenters.
       */
      unsigned int NumTetrahedralStereo();
      /**
       * Check if atom with @p id is a tetrahedral center. 
       * @return True if the atom with @p id has tetrahedral stereochemistry.
       */
      bool HasTetrahedralStereo(unsigned long atomId);
      /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This 
       * function returns 0 if there is no OBTetrahedralStereo object found 
       * with the specified center.
       */
      OBTetrahedralStereo* GetTetrahedralStereo(unsigned long atomId);
      ///@}
 
      ///@name Cis/Trans stereochemistry
      ///@{
      /**
       * Get the number of cis/trans stereocenters.
       */
      unsigned int NumCisTransStereo();
      /**
       * Check if bond with @p id is a stereogenic cis/trans double bond. 
       * @return True if the bond with @p id has cis/trans stereochemistry.
       */
      bool HasCisTransStereo(unsigned long bondId);
      /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This 
       * function returns 0 if there is no OBTetrahedralStereo object found 
       * with the specified center.
       */
      OBCisTransStereo* GetCisTransStereo(unsigned long bondId);
      ///@}

    private:
      /**
       * Ensure the maps are initialized and initialize them only once.
       */
      inline void EnsureInit() { if (!m_init) InitMaps(); }
      /**
       * Initialize @p m_tetrahedralMap and m_cistransMap to contain the 
       * data objects. If @p m_perceive is true and chirality isn't perceived
       * yet, PerceiveStereo will be called.
       */
      void InitMaps();

      OBMol *m_mol;
      bool m_init;
      bool m_perceive;
      std::map<unsigned long, OBTetrahedralStereo*> m_tetrahedralMap;
      std::map<unsigned long, OBCisTransStereo*> m_cistransMap;
  };

  // fwd decl
  class OBBond;
  ///@name High level functions
  ///@{
  /**
   * Convert 0D/2D/3D coordinates to OBStereo objects. The right function will 
   * be selected based on the molecule's dimensionality 
   * (i.e. OBMol::GetDimension()).
   *
   * @sa StereoFrom3D StereoFrom2D StereoFrom0D
   */
  OBAPI void PerceiveStereo(OBMol *mol, bool force = false); 
  /**
   * Convert the 2D depiction of molecule @p mol to OBStereo objects.
   * This function makes use of the lower level functions 
   * TetrahedralFrom2D(), CisTransFrom2D(), SquarePlanarFrom2D(), ...
   *
   * First, symmetry analysis taking stereochemistry into account is 
   * performed iteratively (see OBGraphSym). Next the 2D coordinates, 
   * OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans 
   * are used to construct OBStereoBase derived objects to store the 
   * stereochemistry. These objects will be added to @p mol.
   *
   * Unless perception is forced, this function does nothing if stereochemistry
   * has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
   * doing the actual perception, any data of the OBGenericDataType::StereoData
   * type will be deleted.
   *
     @verbatim
     Reference:
     [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
     Unambiguous Identification of the Stereochemical Characteristics of 
     Compounds During Their Registration in Databases. Molecules 2000, 6,
     915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
     @endverbatim
   *
   * @param mol The molecule containing 2D coordinates.
   * @param force Force to run the perception even if the results are cached.
   *
   * @sa StereoFrom3D StereoFrom0D PerceiveStereo
   */
  OBAPI void StereoFrom2D(OBMol *mol, bool tetfrom0D = false, bool force = false);
  /**
   * Convert the 3D coordinates of molecule @p mol to OBStereo objects. This
   * function makes use of the lower level functions TetrahedralFrom3D(),
   * CisTransFrom3D(), SquarePlanarFrom3D(), ...
   *
   * Unless perception is forced, this function does nothing if stereochemistry
   * has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
   * doing the actual perception, any data of the OBGenericDataType::StereoData
   * type will be deleted.
   *
   * @param mol The molecule containing 3D coordinates.
   * @param force Force to run the perception even if the results are cached.
   *
   * @sa StereoFrom3D StereoFrom0D PerceiveStereo
   */
  OBAPI void StereoFrom3D(OBMol *mol, bool force = false);
  /**
   * Add missing OBStereo objects. Unlike StereoFrom3D() and StereoFrom2D(), this
   * method only adds objects for previously unidentified objects since we 
   * don't want to loose any information. The Config::specified flag for the 
   * newly added structs is always set to false.
   *
   * For example, a smiles is read which has two tetrahedral centers. Only one has
   * stereochemisrty specified using a '@' character. StereoFrom0D() will detect the
   * second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
   *
   * @param mol The molecule.
   * @param updown A map of OBStereo::BondDirection for cis/trans bonds
   *
   * @sa StereoFrom3D StereoFrom2D PerceiveStereo
   */
  OBAPI void StereoFrom0D(OBMol *mol,
      std::map<OBBond*, enum OBStereo::BondDirection> *updown = NULL);
  ///@}

  ///@name Low level functions
  ///@{
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom3D() with the @p addToMol parameter is set 
   * to true. 
   *
   * The algorithm to convert the 3D coordinates to OBTetrahedralStereo object
   * uses the sign of the volume described by the 4 center atom neighbors. Given
   * 4 points \f$a\f$, \f$b\f$, \f$c\f$ and \f$d\f$, the signed volume \f$S_v\f$ 
   * is defined as:
   *
     \f[ S_v = \left| \begin{array}{ccc}
     x_b - x_a & y_b - y_a & z_b - z_a \\
     x_c - x_a & y_c - y_a & z_c - z_a \\
     x_d - x_a & y_d - y_a & z_d - z_a 
     \end{array} \right| \f]
   *
   * The sign of \f$S_v\f$ changes when any of the points cross the plane defined
   * by the other 3 points. To make this less abstract one could say that 
   * a change of sign is equal to inverting the tetrahedral stereochemistry.
   *
   * In case there are only 3 neighbor atoms for the tetrahedral center, the 
   * center atom itself is used as 4th point. This only changes the magnitude 
   * and not the sign of \f$S_v\f$ because the center atom is still on the same
   * side of the plane.
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom3D FindTetrahedralAtoms
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom3D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom2D() with the @p addToMol parameter is set 
   * to true.
   *
   * The algorithm to convert the 2D coordinates and bond properties 
   * (i.e. OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans)
   * uses the sign of a triangle. Given 3 points \f$a\f$, \f$b\f$ and \f$c\f$, the 
   * sign of the trianle \f$S_t\f$ is defined as: 
   *
     \f[ S_t = (x_a - x_c) (y_b - y_c) - (y_a - y_c) (x_b - x_c) \f]
   *
   * This is equation 6 from on the referenced web page. The 3 points used 
   * to calculate the triangle sign always remain in the same plane (i.e. z = 0).
   * The actual meaning of \f$S_t\f$ (i.e. assignment of OBStereo::Winding) depends 
   * on the 4th atom. When the atom is in front of the plane, the sign should be
   * changed to have the same absolute meaning for an atom behind the plane and the
   * same triangle. It is important to note that none of the z coordinates is ever 
   * changed, the molecule always stays 2D (unlike methods which set a pseudo-z 
   * coordinate).
   *
   * @todo document bond property interpretation!
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
     @verbatim
     Reference:
     [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
     Unambiguous Identification of the Stereochemical Characteristics of 
     Compounds During Their Registration in Databases. Molecules 2000, 6,
     915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
     @endverbatim
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom2D FindTetrahedralAtoms
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom0D() with the @p addToMol parameter is set 
   * to true. There is no algorithm used here, all specified flags will be
   * set to false. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom0D FindTetrahedralAtoms
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol = true);
 
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom3D() with the @p addToMol parameter is set 
   * to true.
   *
   * The algorithm to convert the 3D coordinates to OBCisTransStereo objects
   * uses the angles between the single bonds connected to the double bond. 
   * These angles can have 3 values in the ideal case. 60 degrees for cis atoms,
   * 120 degrees for atoms on the same side of the double bond and 180 degrees
   * for trans atoms. A tolerance of 10 degrees is used when comparing these
   * angles. Missing atom coordinates (OBStereo::ImplicitRef) and their bond 
   * vectors will be computed if needed.
   *
   @verbatim
     0      3       0   3
      \    /         \ /      angle 0-*-3 & 1-*-2: 60 degrees (cis)
       C==C    -->    *       angle 0-*-1 & 2-*-3: 120 degrees (same bond atom)
      /    \         / \      angle 0-*-2 & 1-*-3: 180 degrees (trans)
     1      2       1   2
   @endverbatim 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom3D FindCisTransBonds
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom3D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom2D() with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * The algorithm for converting the 2D coordinates uses the same triangle 
   * sign as TetrahedralFrom2D(). Depending on sign of 2 triangles, the right
   * OBStereo::Shape is selected. 
   @verbatim
      0      3       
       \    /        2 triangles: 0-1-b & 2-3-a
        a==b    -->  same sign: U
       /    \        opposite sign: Z
      1      2       
   @endverbatim
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom2D FindCisTransBonds
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom2D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom0D() with the @p addToMol parameter is set 
   * to true. There is no algorithm used here, all specified flags will be
   * set to false.
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom0D FindCisTransBonds
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom0D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits,
      std::map<OBBond*, OBStereo::BondDirection> *updown = NULL,
      bool addToMol = true);
  ///@}


  ///@name Stereogenic unit identification
  ///@{
  OBAPI std::vector<StereogenicUnit> FindStereogenicUnits(OBMol *mol, 
      const std::vector<unsigned int> &symClasses);
  /**
   * Find all tetrahedral centers using the symmetry classes in the 
   * @p symClasses map.
   *
   * There are several criteria for tetrahedral centers:
   * - Cannot be Nitrogen, Phosphorus or Sulfur
   * - Must be a sp3 atom (i.e. OBAtom::GetHyb() == 3)
   * - Must have at least 3 heavy neighbor atoms
   *
   * @c True-stereocenters are identified first. These have atoms 
   * with 4 different symmetry classes attached. When there are only 3 neighbor 
   * atoms, it will suffice that their symmetry classes are different. Since we 
   * already checked that there should be at least 3 heavy atom neighbors, we 
   * can conclude the 4th missing atom is an implicit hydrogen 
   * (OBStereo::ImplicitRef).
   *
   * The term @c para-stereocenters (para = ressemble) is used to denote
   * stereocenters, which altough not true-stereocenters (i.e. have 4 different
   * symmetry classes for neighbours), are still stereogenic.[1] 
   * 
   *
     @verbatim
     Reference:
     [1] M. Razinger, K. Balasubramanian, M. Perdih, M. E. Munk,
     Stereoisomere Generation in Computer-Enhanced Structure Elucidation,
     J. Chem. Inf. Comput. Si. 1993, 33, 812-825
     http://www.mcs.csueastbay.edu/~kbalasub/reprints/282.pdf
     @endverbatim
    *
   * @param mol The molecule.
   * @param symClasses The current symmetry classes (OBGraphSym)
   *
   * @return A vector with a pairs consisting of the unique atom id for each 
   * tetrahedral center. The second bool field in the std::pair is to mark if
   * the stereo center is a true-stereocenter (true) or a para-stereocenter 
   * (false).
   */
  OBAPI std::vector<unsigned long> FindTetrahedralAtoms(OBMol *mol, 
      const std::vector<unsigned int> &symClasses);
  /**
   * Find all cis/trans bonds using the symmetry classes in the 
   * @p symClasses map.
   *
   * There are several criteria for cis/trans bonds:
   * - Must be a double bond (i.e. OBBond::GetBO() == 2)
   * - Cannot be in a ring
   * - The double bond begin and end atom should also have a single bond 
   *   (i.e. OBAtom::HasSingleBond())
   *
   * Next, the double bond begin and end atom are separately checked using 
   * the same logic:
   * - If the atom's valence is 2, we only need to check if the atom 
   *   connected by the single bond is not a hydrogen since we will 
   *   assume the 3th missing atom is an implicit hydrogen.
   * - If there are 3 neighbors, the 2 symmetry classes for the atoms 
   *   connected by a single bond, cannot be the same.
   *
   * The logic above can be simplified: Make sure the begin and end atom 
   * are not connected to 2 atoms with the same symmetry class.
   *
   * @param mol The molecule.
   * @param symClasses The current symmetry classes (OBGraphSym)
   *
   * @return The bond ids for all cis/trans bonds.
   */
  OBAPI std::vector<unsigned long> FindCisTransBonds(OBMol *mol, 
      const std::vector<unsigned int> &symClasses);
  ///@}
 
  /**
   * Create and fill OBCisTransStereo objects using the specified
   * @p ctbonds (bond ids) and map containing directions 
   * (OBStereo::BondDirection). This function is intended to be used
   * by 0D formats. The OBCisTransStereo objects will be stored in the 
   * OBMol.
   *
   * @param mol The molecule.
   * @param ctbonds std::vector containing bond ids for cis/trans bonds
   * @param updown std::map containing up to four bond directions per bond 
   *        id in @p ctbonds.
   */
  OBAPI void CisTransFromUpDown(OBMol *mol,
      const std::vector<unsigned long> &ctbonds,
      std::map<OBBond*, OBStereo::BondDirection> *updown);

  /**
   * @page Stereochemistry
   * @section overview Overview of classes 
   *
   * There are many molecules which contain stereogenic elements. However,
   * certain cases (i.e. tetrahedral, cis/trans) are more common than others
   * (i.e. allene, biphenyl, octrahedral, ...). For the common stereogenic 
   * units, classes are provided. The inheritance of these classes resembles
   * the way they are split into groups.
   *
   * - OBStereoBase
   *   - OBTetraNonPlanarStereo
   *     - OBTetrahedralStereo
   *     - OBExtendedTetrahedralStereo
   *   - OBTetraPlanarStereo
   *     - OBCisTransStereo
   *     - OBExtendedCisTransStereo
   *     - OBSquarePlanarStereo
   *   - OBAxialStereo
   *     - OBTrigonalBipyrimidalStereo
   *     - OBOctahedralStereo
   *
   * @image html tetranonplanar.png
   * @image html tetraplanar.png
   *
   * All specific classes (i.e. OBTetrahedralStereo, ...) have embedded Config 
   * structs which define the actual stereochemistry. All these Config structs 
   * use OBStereo::Ref values to reference or uniquely identify atoms. Make sure
   * to read about OBStereo::Ref and the related functions (in OBStereo). OBStereo
   * is also a placeholder for various enums with predefined values for parameters 
   * etc. These enums are used throughout the different stereo classes but having 
   * these enums in a single location makes it easier to remember. When working 
   * with stereo classes, you normally don't need to use any of the parent classes
   * directly. Only OBStereo and the specific class are needed.
   *
   * @section usage Basic usage
   *
   * The OBStereoFacade hides the complexity of working with stereochemistry. When 
   * using openbabel as a library, this is by far the easiest way to access 
   * stereochemistry information. 
   * The header for the specific OBStereo::Type type is all you need to include.
   * These are:
   * - @em openbabel/stereo/tetrahedral.h
   * - @em openbabel/stereo/cistrans.h
   * - @em openbabel/stereo/squareplanar.h
   *
   * All these headers also include @em openbabel/stereo/stereo.h providing 
   * declarations for OBStereo & OBStereoFacade.
   *
     @code
     #include <iostream>
     #include <openbabel/mol.h>
     #include <openbabel/obconversion.h>

     #include <openbabel/stereo/tetrahedral.h>

     using namespace OpenBabel;

     int main()
     {
       OBMol mol;
       OBConversion conv;
       conv.SetInFormat("smi");
       conv.ReadString(&mol, "C[C@H](Cl)Br");

       OBStereoFacade facade(&mol);

       FOR_ATOMS_OF_MOL(atom, mol) {
         if (facade.HasTetrahedralStereo(atom->GetId()))
           std::cout << facade.GetTetrahedralStereo(atom->GetId()) << std::endl;
       }
     }
     @endcode
   *
   * All specific stereo classes and their embedded Config struct have an
   * operator<< function which allows them to be used with std::ostream objects
   * (e.g. std::cout, std::err, ...). These functions are often useful when 
   * debugging code.
   * 
   * @section details Details on implementation
   *
   * The detection of stereogenic units start with symmetry analysis. However, a 
   * complete symmetry analysis also needs to take stereochemistry into account.
   * In practice, this means stereochemistry will be found iteratively. At each
   * iteration, the current atom symmetry classes are used to identify stereogenic
   * units. The details about how the symmetry classes are used depends on the type
   * (OBStereo::Type) of stereogenic unit. For tetrahedral centers, having 3 heavy
   * atom neighbors with different symmetry classes or 4 neighbors with different
   * symmetry classes means the atom is chiral. See FindTetrahedralAtoms(), 
   * FindCisTransBonds(), FindSquarePlanarAtoms(), ... for details.
   *
   * After identifying the stereogenic units, Config structs with all the 
   * information on the spacial arrangement of the groups still have to be 
   * created. This involves interpreting various ways to represent 
   * stereochemisrty:
   *
   * - 3D coordinates: StereoFrom3D()
   * - 2D coordinates: StereoFrom2D()
   * - 0D coordinates: StereoFrom0D()
   * 
   * Both StereoFrom3D() and StereoFrom2D() delete all existing stereochemistry objects
   * before adding new ones. For molecules with 3D coordinates, it is evident that 
   * all information is specified by the coordinates itself. However, if a file format 
   * uses stereo parity flags, Config structs must be constructed using lower level 
   * functions and StereoFrom3D() should not be called. In these cases information 
   * could be lost by calling StereoFrom3D() after reading the file (the stereo flag might have 
   * indicated the stereochemistry was unspecified or the flag might not match the 
   * coordinates). In the case of 2D molecules, the coordinates together with bond
   * properties (OBBond::Hash, OBBond::Wedge, OBBond::WedgeOrHash and 
   * OBBond::CisOrTrans) define the stereochemistry. Again, lower level functions 
   * can be used when stereo flags need to be used.
   *
   * StereoFrom0D() works slightly different than 3D/2D. Here, deleting the 
   * stereochemistry would always result in lost information. Instead StereoFrom0D()
   * only adds new objects for stereogenic units which were previously not found.
   * For example, a smiles is read which has two tetrahedral centers. Only one has
   * stereochemistry specified using a '@' character. StereoFrom0D() will detect the
   * second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
   * The Config::specified flag for the newly added structs is always set to false.
   * 
   * Assuming the format code has correctly set the molecule dimensions (OBMol::GetDimesions), 
   * PerceiveStereo() will automatically select the correct function to call.
   * When StereoFrom3D(), StereoFrom2D() or StereoFrom0D() are not used, make sure to always
   * set OBMol::HasChiralityPerceived() before returning from the format's ReadMolecule().
   *
   *
   * @section formats Guidelines for formats
   * @subsection input Reading files
   *
   * - Read the section above
   * - The MDL format (mdlformat.cpp) is a good example for 2D/3D formats with or 
   *   without parity flags.
   * - The SMILES format (smilesformat.cpp) is a good example for 0D formats.
   *
   * @subsection output Writing files
   *
   * For many file formats no additional code is needed. For example, if a 3D format
   * doesn't require stereo parity flags, writing the coordinates is enough. For 2D
   * file formats it will often suffice to write the coordinates and bond properties.
   * If parity flags are needed, the OBStereoFacade class can be used to retreive the
   * objects for all types of stereochemistry supported by the file format.
   *
   *
   *
   *
   *
   */
  
  ///@}  addtogroup
}

#endif
