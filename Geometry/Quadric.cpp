/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Quadric.h"

Quadric::Quadric(const Matrix4x4<float> &q) { this->mQuadric = q; }

Quadric::~Quadric() {}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const { 
    Vector4<float> pWorld(x, y, z, 1);
    Matrix4x4<float> M = mWorld2Obj;
    Matrix4x4<float> Q = mQuadric;
    Matrix4x4<float> Qhat = M.Inverse().Transpose() * Q * M.Inverse();
    float val = pWorld * (Qhat * pWorld);
    return val;
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */
Vector3<float> Quadric::GetGradient(float x, float y, float z) const {
    TransformW2O(x, y, z);
    Vector4<float> pObject(x, y, z, 1);
    Vector4<float> val = 2.0f * (mQuadric * pObject);
    Vector3<float> result = {val[0], val[1], val[2]};

  return (result);
}
