/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class vector3 {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected vector3(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(vector3 obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_vector3(swigCPtr);
    }
    swigCPtr = 0;
  }

  public vector3(double x, double y, double z) {
    this(openbabelJNI.new_vector3__SWIG_0(x, y, z), true);
  }

  public vector3(double x, double y) {
    this(openbabelJNI.new_vector3__SWIG_1(x, y), true);
  }

  public vector3(double x) {
    this(openbabelJNI.new_vector3__SWIG_2(x), true);
  }

  public vector3() {
    this(openbabelJNI.new_vector3__SWIG_3(), true);
  }

  public vector3(vector3 v) {
    this(openbabelJNI.new_vector3__SWIG_4(vector3.getCPtr(v), v), true);
  }

  public void Set(double x, double y, double z) {
    openbabelJNI.vector3_Set__SWIG_0(swigCPtr, this, x, y, z);
  }

  public void Set(SWIGTYPE_p_double c) {
    openbabelJNI.vector3_Set__SWIG_1(swigCPtr, this, SWIGTYPE_p_double.getCPtr(c));
  }

  public void SetX(double x) {
    openbabelJNI.vector3_SetX(swigCPtr, this, x);
  }

  public void SetY(double y) {
    openbabelJNI.vector3_SetY(swigCPtr, this, y);
  }

  public void SetZ(double z) {
    openbabelJNI.vector3_SetZ(swigCPtr, this, z);
  }

  public void Get(SWIGTYPE_p_double c) {
    openbabelJNI.vector3_Get(swigCPtr, this, SWIGTYPE_p_double.getCPtr(c));
  }

  public boolean IsApprox(vector3 arg0, double precision) {
    return openbabelJNI.vector3_IsApprox(swigCPtr, this, vector3.getCPtr(arg0), arg0, precision);
  }

  public SWIGTYPE_p_double AsArray() {
    long cPtr = openbabelJNI.vector3_AsArray(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public void randomUnitVector(OBRandom oeRand) {
    openbabelJNI.vector3_randomUnitVector__SWIG_0(swigCPtr, this, OBRandom.getCPtr(oeRand), oeRand);
  }

  public void randomUnitVector() {
    openbabelJNI.vector3_randomUnitVector__SWIG_1(swigCPtr, this);
  }

  public vector3 normalize() {
    return new vector3(openbabelJNI.vector3_normalize(swigCPtr, this), false);
  }

  public boolean CanBeNormalized() {
    return openbabelJNI.vector3_CanBeNormalized(swigCPtr, this);
  }

  public double length() {
    return openbabelJNI.vector3_length(swigCPtr, this);
  }

  public double length_2() {
    return openbabelJNI.vector3_length_2(swigCPtr, this);
  }

  public double x() {
    return openbabelJNI.vector3_x(swigCPtr, this);
  }

  public double y() {
    return openbabelJNI.vector3_y(swigCPtr, this);
  }

  public double z() {
    return openbabelJNI.vector3_z(swigCPtr, this);
  }

  public double distSq(vector3 vv) {
    return openbabelJNI.vector3_distSq(swigCPtr, this, vector3.getCPtr(vv), vv);
  }

  public boolean createOrthoVector(vector3 v) {
    return openbabelJNI.vector3_createOrthoVector(swigCPtr, this, vector3.getCPtr(v), v);
  }

}
