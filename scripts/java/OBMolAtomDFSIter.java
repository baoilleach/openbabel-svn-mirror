/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolAtomDFSIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolAtomDFSIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolAtomDFSIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMolAtomDFSIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolAtomDFSIter() {
    this(openbabelJNI.new_OBMolAtomDFSIter__SWIG_0(), true);
  }

  public OBMolAtomDFSIter(OBMol mol) {
    this(openbabelJNI.new_OBMolAtomDFSIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolAtomDFSIter(OBMolAtomDFSIter ai) {
    this(openbabelJNI.new_OBMolAtomDFSIter__SWIG_2(OBMolAtomDFSIter.getCPtr(ai), ai), true);
  }

  public boolean good() {
    return openbabelJNI.OBMolAtomDFSIter_good(swigCPtr, this);
  }

  public OBMolAtomDFSIter inc() {
    return new OBMolAtomDFSIter(openbabelJNI.OBMolAtomDFSIter_inc__SWIG_0(swigCPtr, this), false);
  }

  public OBMolAtomDFSIter inc(int arg0) {
    return new OBMolAtomDFSIter(openbabelJNI.OBMolAtomDFSIter_inc__SWIG_1(swigCPtr, this, arg0), true);
  }

  public OBAtom deref() {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_deref(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom __ref__() {
    return new OBAtom(openbabelJNI.OBMolAtomDFSIter___ref__(swigCPtr, this), false);
  }

  public void setVisit(boolean value) {
    openbabelJNI.OBMolAtomDFSIter_Visit_set(swigCPtr, this, value);
  }

  public boolean getVisit() {
    return openbabelJNI.OBMolAtomDFSIter_Visit_get(swigCPtr, this);
  }

  public void Clear() {
    openbabelJNI.OBMolAtomDFSIter_Clear(swigCPtr, this);
  }

  public void SetIdx(int idx) {
    openbabelJNI.OBMolAtomDFSIter_SetIdx(swigCPtr, this, idx);
  }

  public void SetHyb(int hyb) {
    openbabelJNI.OBMolAtomDFSIter_SetHyb(swigCPtr, this, hyb);
  }

  public void SetAtomicNum(int atomicnum) {
    openbabelJNI.OBMolAtomDFSIter_SetAtomicNum(swigCPtr, this, atomicnum);
  }

  public void SetIsotope(long iso) {
    openbabelJNI.OBMolAtomDFSIter_SetIsotope(swigCPtr, this, iso);
  }

  public void SetImplicitValence(int val) {
    openbabelJNI.OBMolAtomDFSIter_SetImplicitValence(swigCPtr, this, val);
  }

  public void IncrementImplicitValence() {
    openbabelJNI.OBMolAtomDFSIter_IncrementImplicitValence(swigCPtr, this);
  }

  public void DecrementImplicitValence() {
    openbabelJNI.OBMolAtomDFSIter_DecrementImplicitValence(swigCPtr, this);
  }

  public void SetFormalCharge(int fcharge) {
    openbabelJNI.OBMolAtomDFSIter_SetFormalCharge(swigCPtr, this, fcharge);
  }

  public void SetSpinMultiplicity(short spin) {
    openbabelJNI.OBMolAtomDFSIter_SetSpinMultiplicity(swigCPtr, this, spin);
  }

  public void SetType(String type) {
    openbabelJNI.OBMolAtomDFSIter_SetType__SWIG_0(swigCPtr, this, type);
  }

  public void SetType(SWIGTYPE_p_std__string type) {
    openbabelJNI.OBMolAtomDFSIter_SetType__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__string.getCPtr(type));
  }

  public void SetPartialCharge(double pcharge) {
    openbabelJNI.OBMolAtomDFSIter_SetPartialCharge(swigCPtr, this, pcharge);
  }

  public void SetVector(vector3 v) {
    openbabelJNI.OBMolAtomDFSIter_SetVector__SWIG_0(swigCPtr, this, vector3.getCPtr(v), v);
  }

  public void SetVector(double x, double y, double z) {
    openbabelJNI.OBMolAtomDFSIter_SetVector__SWIG_1(swigCPtr, this, x, y, z);
  }

  public void SetVector() {
    openbabelJNI.OBMolAtomDFSIter_SetVector__SWIG_2(swigCPtr, this);
  }

  public void SetCoordPtr(SWIGTYPE_p_p_double c) {
    openbabelJNI.OBMolAtomDFSIter_SetCoordPtr(swigCPtr, this, SWIGTYPE_p_p_double.getCPtr(c));
  }

  public void SetResidue(OBResidue res) {
    openbabelJNI.OBMolAtomDFSIter_SetResidue(swigCPtr, this, OBResidue.getCPtr(res), res);
  }

  public void SetParent(OBMol ptr) {
    openbabelJNI.OBMolAtomDFSIter_SetParent(swigCPtr, this, OBMol.getCPtr(ptr), ptr);
  }

  public void SetAromatic() {
    openbabelJNI.OBMolAtomDFSIter_SetAromatic(swigCPtr, this);
  }

  public void UnsetAromatic() {
    openbabelJNI.OBMolAtomDFSIter_UnsetAromatic(swigCPtr, this);
  }

  public void SetClockwiseStereo() {
    openbabelJNI.OBMolAtomDFSIter_SetClockwiseStereo(swigCPtr, this);
  }

  public void SetAntiClockwiseStereo() {
    openbabelJNI.OBMolAtomDFSIter_SetAntiClockwiseStereo(swigCPtr, this);
  }

  public void SetPositiveStereo() {
    openbabelJNI.OBMolAtomDFSIter_SetPositiveStereo(swigCPtr, this);
  }

  public void SetNegativeStereo() {
    openbabelJNI.OBMolAtomDFSIter_SetNegativeStereo(swigCPtr, this);
  }

  public void UnsetStereo() {
    openbabelJNI.OBMolAtomDFSIter_UnsetStereo(swigCPtr, this);
  }

  public void SetInRing() {
    openbabelJNI.OBMolAtomDFSIter_SetInRing(swigCPtr, this);
  }

  public void SetChiral() {
    openbabelJNI.OBMolAtomDFSIter_SetChiral(swigCPtr, this);
  }

  public void ClearCoordPtr() {
    openbabelJNI.OBMolAtomDFSIter_ClearCoordPtr(swigCPtr, this);
  }

  public int GetFormalCharge() {
    return openbabelJNI.OBMolAtomDFSIter_GetFormalCharge(swigCPtr, this);
  }

  public long GetAtomicNum() {
    return openbabelJNI.OBMolAtomDFSIter_GetAtomicNum(swigCPtr, this);
  }

  public int GetIsotope() {
    return openbabelJNI.OBMolAtomDFSIter_GetIsotope(swigCPtr, this);
  }

  public int GetSpinMultiplicity() {
    return openbabelJNI.OBMolAtomDFSIter_GetSpinMultiplicity(swigCPtr, this);
  }

  public double GetAtomicMass() {
    return openbabelJNI.OBMolAtomDFSIter_GetAtomicMass(swigCPtr, this);
  }

  public double GetExactMass() {
    return openbabelJNI.OBMolAtomDFSIter_GetExactMass(swigCPtr, this);
  }

  public long GetIdx() {
    return openbabelJNI.OBMolAtomDFSIter_GetIdx(swigCPtr, this);
  }

  public long GetCoordinateIdx() {
    return openbabelJNI.OBMolAtomDFSIter_GetCoordinateIdx(swigCPtr, this);
  }

  public long GetCIdx() {
    return openbabelJNI.OBMolAtomDFSIter_GetCIdx(swigCPtr, this);
  }

  public long GetValence() {
    return openbabelJNI.OBMolAtomDFSIter_GetValence(swigCPtr, this);
  }

  public long GetHyb() {
    return openbabelJNI.OBMolAtomDFSIter_GetHyb(swigCPtr, this);
  }

  public long GetImplicitValence() {
    return openbabelJNI.OBMolAtomDFSIter_GetImplicitValence(swigCPtr, this);
  }

  public long GetHvyValence() {
    return openbabelJNI.OBMolAtomDFSIter_GetHvyValence(swigCPtr, this);
  }

  public long GetHeteroValence() {
    return openbabelJNI.OBMolAtomDFSIter_GetHeteroValence(swigCPtr, this);
  }

  public String GetType() {
    return openbabelJNI.OBMolAtomDFSIter_GetType(swigCPtr, this);
  }

  public double GetX() {
    return openbabelJNI.OBMolAtomDFSIter_GetX(swigCPtr, this);
  }

  public double x() {
    return openbabelJNI.OBMolAtomDFSIter_x(swigCPtr, this);
  }

  public double GetY() {
    return openbabelJNI.OBMolAtomDFSIter_GetY(swigCPtr, this);
  }

  public double y() {
    return openbabelJNI.OBMolAtomDFSIter_y(swigCPtr, this);
  }

  public double GetZ() {
    return openbabelJNI.OBMolAtomDFSIter_GetZ(swigCPtr, this);
  }

  public double z() {
    return openbabelJNI.OBMolAtomDFSIter_z(swigCPtr, this);
  }

  public SWIGTYPE_p_double GetCoordinate() {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetCoordinate(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public vector3 GetVector() {
    return new vector3(openbabelJNI.OBMolAtomDFSIter_GetVector(swigCPtr, this), false);
  }

  public double GetPartialCharge() {
    return openbabelJNI.OBMolAtomDFSIter_GetPartialCharge(swigCPtr, this);
  }

  public OBResidue GetResidue() {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetResidue(swigCPtr, this);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public OBMol GetParent() {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetParent(swigCPtr, this);
    return (cPtr == 0) ? null : new OBMol(cPtr, false);
  }

  public boolean GetNewBondVector(vector3 v, double length) {
    return openbabelJNI.OBMolAtomDFSIter_GetNewBondVector(swigCPtr, this, vector3.getCPtr(v), v, length);
  }

  public OBBond GetBond(OBAtom arg0) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetBond(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom GetNextAtom() {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetNextAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator BeginBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator(openbabelJNI.OBMolAtomDFSIter_BeginBonds(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator EndBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator(openbabelJNI.OBMolAtomDFSIter_EndBonds(swigCPtr, this), true);
  }

  public OBBond BeginBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_BeginBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond NextBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_NextBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom BeginNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator arg0) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_BeginNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom NextNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator arg0) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_NextNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public double GetDistance(int index) {
    return openbabelJNI.OBMolAtomDFSIter_GetDistance__SWIG_0(swigCPtr, this, index);
  }

  public double GetDistance(OBAtom arg0) {
    return openbabelJNI.OBMolAtomDFSIter_GetDistance__SWIG_1(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public double GetAngle(int b, int c) {
    return openbabelJNI.OBMolAtomDFSIter_GetAngle__SWIG_0(swigCPtr, this, b, c);
  }

  public double GetAngle(OBAtom b, OBAtom c) {
    return openbabelJNI.OBMolAtomDFSIter_GetAngle__SWIG_1(swigCPtr, this, OBAtom.getCPtr(b), b, OBAtom.getCPtr(c), c);
  }

  public void NewResidue() {
    openbabelJNI.OBMolAtomDFSIter_NewResidue(swigCPtr, this);
  }

  public void DeleteResidue() {
    openbabelJNI.OBMolAtomDFSIter_DeleteResidue(swigCPtr, this);
  }

  public void AddBond(OBBond bond) {
    openbabelJNI.OBMolAtomDFSIter_AddBond(swigCPtr, this, OBBond.getCPtr(bond), bond);
  }

  public void InsertBond(SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator i, OBBond bond) {
    openbabelJNI.OBMolAtomDFSIter_InsertBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBBond_p_t__iterator.getCPtr(i), OBBond.getCPtr(bond), bond);
  }

  public boolean DeleteBond(OBBond arg0) {
    return openbabelJNI.OBMolAtomDFSIter_DeleteBond(swigCPtr, this, OBBond.getCPtr(arg0), arg0);
  }

  public void ClearBond() {
    openbabelJNI.OBMolAtomDFSIter_ClearBond(swigCPtr, this);
  }

  public long CountFreeOxygens() {
    return openbabelJNI.OBMolAtomDFSIter_CountFreeOxygens(swigCPtr, this);
  }

  public long ImplicitHydrogenCount() {
    return openbabelJNI.OBMolAtomDFSIter_ImplicitHydrogenCount(swigCPtr, this);
  }

  public long ExplicitHydrogenCount(boolean ExcludeIsotopes) {
    return openbabelJNI.OBMolAtomDFSIter_ExplicitHydrogenCount__SWIG_0(swigCPtr, this, ExcludeIsotopes);
  }

  public long ExplicitHydrogenCount() {
    return openbabelJNI.OBMolAtomDFSIter_ExplicitHydrogenCount__SWIG_1(swigCPtr, this);
  }

  public long MemberOfRingCount() {
    return openbabelJNI.OBMolAtomDFSIter_MemberOfRingCount(swigCPtr, this);
  }

  public long MemberOfRingSize() {
    return openbabelJNI.OBMolAtomDFSIter_MemberOfRingSize(swigCPtr, this);
  }

  public long CountRingBonds() {
    return openbabelJNI.OBMolAtomDFSIter_CountRingBonds(swigCPtr, this);
  }

  public double SmallestBondAngle() {
    return openbabelJNI.OBMolAtomDFSIter_SmallestBondAngle(swigCPtr, this);
  }

  public double AverageBondAngle() {
    return openbabelJNI.OBMolAtomDFSIter_AverageBondAngle(swigCPtr, this);
  }

  public long BOSum() {
    return openbabelJNI.OBMolAtomDFSIter_BOSum(swigCPtr, this);
  }

  public long KBOSum() {
    return openbabelJNI.OBMolAtomDFSIter_KBOSum(swigCPtr, this);
  }

  public boolean HtoMethyl() {
    return openbabelJNI.OBMolAtomDFSIter_HtoMethyl(swigCPtr, this);
  }

  public boolean SetHybAndGeom(int arg0) {
    return openbabelJNI.OBMolAtomDFSIter_SetHybAndGeom(swigCPtr, this, arg0);
  }

  public void ForceNoH() {
    openbabelJNI.OBMolAtomDFSIter_ForceNoH(swigCPtr, this);
  }

  public boolean HasNoHForced() {
    return openbabelJNI.OBMolAtomDFSIter_HasNoHForced(swigCPtr, this);
  }

  public boolean HasResidue() {
    return openbabelJNI.OBMolAtomDFSIter_HasResidue(swigCPtr, this);
  }

  public boolean IsHydrogen() {
    return openbabelJNI.OBMolAtomDFSIter_IsHydrogen(swigCPtr, this);
  }

  public boolean IsCarbon() {
    return openbabelJNI.OBMolAtomDFSIter_IsCarbon(swigCPtr, this);
  }

  public boolean IsNitrogen() {
    return openbabelJNI.OBMolAtomDFSIter_IsNitrogen(swigCPtr, this);
  }

  public boolean IsOxygen() {
    return openbabelJNI.OBMolAtomDFSIter_IsOxygen(swigCPtr, this);
  }

  public boolean IsSulfur() {
    return openbabelJNI.OBMolAtomDFSIter_IsSulfur(swigCPtr, this);
  }

  public boolean IsPhosphorus() {
    return openbabelJNI.OBMolAtomDFSIter_IsPhosphorus(swigCPtr, this);
  }

  public boolean IsAromatic() {
    return openbabelJNI.OBMolAtomDFSIter_IsAromatic(swigCPtr, this);
  }

  public boolean IsInRing() {
    return openbabelJNI.OBMolAtomDFSIter_IsInRing(swigCPtr, this);
  }

  public boolean IsInRingSize(int arg0) {
    return openbabelJNI.OBMolAtomDFSIter_IsInRingSize(swigCPtr, this, arg0);
  }

  public boolean IsHeteroatom() {
    return openbabelJNI.OBMolAtomDFSIter_IsHeteroatom(swigCPtr, this);
  }

  public boolean IsNotCorH() {
    return openbabelJNI.OBMolAtomDFSIter_IsNotCorH(swigCPtr, this);
  }

  public boolean IsConnected(OBAtom arg0) {
    return openbabelJNI.OBMolAtomDFSIter_IsConnected(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneThree(OBAtom arg0) {
    return openbabelJNI.OBMolAtomDFSIter_IsOneThree(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneFour(OBAtom arg0) {
    return openbabelJNI.OBMolAtomDFSIter_IsOneFour(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsCarboxylOxygen() {
    return openbabelJNI.OBMolAtomDFSIter_IsCarboxylOxygen(swigCPtr, this);
  }

  public boolean IsPhosphateOxygen() {
    return openbabelJNI.OBMolAtomDFSIter_IsPhosphateOxygen(swigCPtr, this);
  }

  public boolean IsSulfateOxygen() {
    return openbabelJNI.OBMolAtomDFSIter_IsSulfateOxygen(swigCPtr, this);
  }

  public boolean IsNitroOxygen() {
    return openbabelJNI.OBMolAtomDFSIter_IsNitroOxygen(swigCPtr, this);
  }

  public boolean IsAmideNitrogen() {
    return openbabelJNI.OBMolAtomDFSIter_IsAmideNitrogen(swigCPtr, this);
  }

  public boolean IsPolarHydrogen() {
    return openbabelJNI.OBMolAtomDFSIter_IsPolarHydrogen(swigCPtr, this);
  }

  public boolean IsNonPolarHydrogen() {
    return openbabelJNI.OBMolAtomDFSIter_IsNonPolarHydrogen(swigCPtr, this);
  }

  public boolean IsAromaticNOxide() {
    return openbabelJNI.OBMolAtomDFSIter_IsAromaticNOxide(swigCPtr, this);
  }

  public boolean IsChiral() {
    return openbabelJNI.OBMolAtomDFSIter_IsChiral(swigCPtr, this);
  }

  public boolean IsAxial() {
    return openbabelJNI.OBMolAtomDFSIter_IsAxial(swigCPtr, this);
  }

  public boolean IsClockwise() {
    return openbabelJNI.OBMolAtomDFSIter_IsClockwise(swigCPtr, this);
  }

  public boolean IsAntiClockwise() {
    return openbabelJNI.OBMolAtomDFSIter_IsAntiClockwise(swigCPtr, this);
  }

  public boolean IsPositiveStereo() {
    return openbabelJNI.OBMolAtomDFSIter_IsPositiveStereo(swigCPtr, this);
  }

  public boolean IsNegativeStereo() {
    return openbabelJNI.OBMolAtomDFSIter_IsNegativeStereo(swigCPtr, this);
  }

  public boolean HasChiralitySpecified() {
    return openbabelJNI.OBMolAtomDFSIter_HasChiralitySpecified(swigCPtr, this);
  }

  public boolean HasChiralVolume() {
    return openbabelJNI.OBMolAtomDFSIter_HasChiralVolume(swigCPtr, this);
  }

  public boolean IsHbondAcceptor() {
    return openbabelJNI.OBMolAtomDFSIter_IsHbondAcceptor(swigCPtr, this);
  }

  public boolean IsHbondDonor() {
    return openbabelJNI.OBMolAtomDFSIter_IsHbondDonor(swigCPtr, this);
  }

  public boolean IsHbondDonorH() {
    return openbabelJNI.OBMolAtomDFSIter_IsHbondDonorH(swigCPtr, this);
  }

  public boolean HasAlphaBetaUnsat(boolean includePandS) {
    return openbabelJNI.OBMolAtomDFSIter_HasAlphaBetaUnsat__SWIG_0(swigCPtr, this, includePandS);
  }

  public boolean HasAlphaBetaUnsat() {
    return openbabelJNI.OBMolAtomDFSIter_HasAlphaBetaUnsat__SWIG_1(swigCPtr, this);
  }

  public boolean HasBondOfOrder(long arg0) {
    return openbabelJNI.OBMolAtomDFSIter_HasBondOfOrder(swigCPtr, this, arg0);
  }

  public int CountBondsOfOrder(long arg0) {
    return openbabelJNI.OBMolAtomDFSIter_CountBondsOfOrder(swigCPtr, this, arg0);
  }

  public boolean HasNonSingleBond() {
    return openbabelJNI.OBMolAtomDFSIter_HasNonSingleBond(swigCPtr, this);
  }

  public boolean HasSingleBond() {
    return openbabelJNI.OBMolAtomDFSIter_HasSingleBond(swigCPtr, this);
  }

  public boolean HasDoubleBond() {
    return openbabelJNI.OBMolAtomDFSIter_HasDoubleBond(swigCPtr, this);
  }

  public boolean HasAromaticBond() {
    return openbabelJNI.OBMolAtomDFSIter_HasAromaticBond(swigCPtr, this);
  }

  public boolean MatchesSMARTS(String arg0) {
    return openbabelJNI.OBMolAtomDFSIter_MatchesSMARTS(swigCPtr, this, arg0);
  }

  public OBBase DoTransformations(SWIGTYPE_p_std__mapTstd__string_std__string_t arg0) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_DoTransformations(swigCPtr, this, SWIGTYPE_p_std__mapTstd__string_std__string_t.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBBase(cPtr, false);
  }

  public String ClassDescription() {
    return openbabelJNI.OBMolAtomDFSIter_ClassDescription(swigCPtr, this);
  }

  public boolean HasData(long type) {
    return openbabelJNI.OBMolAtomDFSIter_HasData__SWIG_2(swigCPtr, this, type);
  }

  public void DeleteData(long type) {
    openbabelJNI.OBMolAtomDFSIter_DeleteData__SWIG_0(swigCPtr, this, type);
  }

  public void DeleteData(OBGenericData arg0) {
    openbabelJNI.OBMolAtomDFSIter_DeleteData__SWIG_1(swigCPtr, this, OBGenericData.getCPtr(arg0), arg0);
  }

  public void DeleteData(vectorData arg0) {
    openbabelJNI.OBMolAtomDFSIter_DeleteData__SWIG_2(swigCPtr, this, vectorData.getCPtr(arg0), arg0);
  }

  public void SetData(OBGenericData d) {
    openbabelJNI.OBMolAtomDFSIter_SetData(swigCPtr, this, OBGenericData.getCPtr(d), d);
  }

  public long DataSize() {
    return openbabelJNI.OBMolAtomDFSIter_DataSize(swigCPtr, this);
  }

  public OBGenericData GetData(long type) {
    long cPtr = openbabelJNI.OBMolAtomDFSIter_GetData__SWIG_0(swigCPtr, this, type);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public vectorData GetData() {
    return new vectorData(openbabelJNI.OBMolAtomDFSIter_GetData__SWIG_3(swigCPtr, this), false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator BeginData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(openbabelJNI.OBMolAtomDFSIter_BeginData(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator EndData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(openbabelJNI.OBMolAtomDFSIter_EndData(swigCPtr, this), true);
  }

}
