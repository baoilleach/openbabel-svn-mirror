/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolAngleIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolAngleIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolAngleIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMolAngleIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolAngleIter() {
    this(openbabelJNI.new_OBMolAngleIter__SWIG_0(), true);
  }

  public OBMolAngleIter(OBMol mol) {
    this(openbabelJNI.new_OBMolAngleIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolAngleIter(OBMolAngleIter ai) {
    this(openbabelJNI.new_OBMolAngleIter__SWIG_2(OBMolAngleIter.getCPtr(ai), ai), true);
  }

  public boolean good() {
    return openbabelJNI.OBMolAngleIter_good(swigCPtr, this);
  }

  public OBMolAngleIter inc(int arg0) {
    return new OBMolAngleIter(openbabelJNI.OBMolAngleIter_inc(swigCPtr, this, arg0), true);
  }

  public vectorUnsignedInt __ref__() {
    return new vectorUnsignedInt(openbabelJNI.OBMolAngleIter___ref__(swigCPtr, this), true);
  }

}
