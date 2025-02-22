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
#ifndef __operatoradvect_h__
#define __operatoradvect_h__

#include "Levelset/LevelSetOperator.h"
#include "Math/Function3D.h"
#include "Math/Matrix4x4.h"

/*! \brief A level set operator that does external advection
 *
 * This class implements level set advectionr in an external vector field by the
 * PDE
 *
 *  \f$
 *  \dfrac{\partial \phi}{\partial t} + \mathbf{V}(\mathbf{x})\cdot \nabla \phi
 * = 0 \f$
 */
//! \lab4 Implement advection in external vector field
class OperatorAdvect : public LevelSetOperator {
protected:
  Function3D<Vector3<float>> *mVectorField;

public:
  OperatorAdvect(LevelSet *LS, Function3D<Vector3<float>> *vf)
      : LevelSetOperator(LS), mVectorField(vf) {}

  virtual float ComputeTimestep() {
    // Compute and return a stable timestep
    // (Hint: Function3D::GetMaxValue())
      float dX = mLS->GetDx();
      Vector3<float> maxes = mVectorField->GetMaxValue();
      return (dX / std::max(abs(maxes[0]), std::max(abs(maxes[1]), abs(maxes[2]))));
  }

  virtual void Propagate(float time) {
    // Determine timestep for stability
    float dt = ComputeTimestep();

    // Propagate level set with stable timestep dt
    // until requested time is reached
    for (float elapsed = 0; elapsed < time;) {

      if (dt > time - elapsed)
        dt = time - elapsed;
      elapsed += dt;

      IntegrateEuler(dt);
      // IntegrateRungeKutta(dt);
    }
  }

  virtual float Evaluate(size_t i, size_t j, size_t k) {
      // Compute the rate of change (dphi/dt)

      // Remember that the point (i,j,k) is given in grid coordinates, while
      // the velocity field used for advection needs to be sampled in
      // world coordinates (x,y,z). You can use LevelSet::TransformGridToWorld()
      // for this task.
      float x = i;
      float y = j;
      float z = k;

      mLS->TransformGridToWorld(x, y, z);
      Vector3<float> v = mVectorField->GetValue(x, y, z);
      Vector3<float> gradient = {0, 0, 0};
      gradient[0] = (v[0] >= 0) ? mLS->DiffXm(i, j, k) : mLS->DiffXp(i, j, k);
      gradient[1] = (v[1] >= 0) ? mLS->DiffYm(i, j, k) : mLS->DiffYp(i, j, k);
      gradient[2] = (v[2] >= 0) ? mLS->DiffZm(i, j, k) : mLS->DiffZp(i, j, k);

     return (-1 * (v * gradient));
  }
};

#endif
