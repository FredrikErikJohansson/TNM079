/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __operatormeancurvatureflow_h__
#define __operatormeancurvatureflow_h__

#include "Levelset/LevelSetOperator.h"

/*! \brief A level set operator that does mean curvature flow.
 *
 * This class implements level set propagation in the normal direction
 * as defined by the mean curvature flow \f$\kappa\f$ in the following PDE
 *
 *  \f[
 *  \dfrac{\partial \phi}{\partial t} + \alpha \kappa|\nabla \phi| = 0
 *  \f]
 */
//! \lab4 Implement mean curvature flow
class OperatorMeanCurvatureFlow : public LevelSetOperator {
protected:
  //! Scaling parameter, affects time step constraint
  float mAlpha;

public:
  OperatorMeanCurvatureFlow(LevelSet *LS, float alpha = .9f)
      : LevelSetOperator(LS), mAlpha(alpha) {}

  virtual float ComputeTimestep() {
    // Compute and return a stable timestep
    return mLS->GetDx() * mLS->GetDx() / 6 * mAlpha;
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
    float dx = mLS->DiffXpm(i, j, k);
    float dy = mLS->DiffYpm(i, j, k);
    float dz = mLS->DiffZpm(i, j, k);

    float dxx = mLS->Diff2Xpm(i, j, k);
    float dyy = mLS->Diff2Ypm(i, j, k);
    float dzz = mLS->Diff2Zpm(i, j, k);

    float dyz = mLS->Diff2YZpm(i, j, k);
    float dzx = mLS->Diff2ZXpm(i, j, k);
    float dxy = mLS->Diff2XYpm(i, j, k);

    float denominator = (2.0f * std::pow((dx * dx + dy * dy + dz * dz), 3.0f / 2.0f));
    float term1 = ((dx * dx) * (dyy + dzz) - (2.0f * dy * dz * dyz)) / denominator;
    float term2 = ((dy * dy) * (dxx + dzz) - (2.0f * dx * dz * dzx)) / denominator;
    float term3 = ((dz * dz) * (dxx + dyy) - (2.0f * dx * dy * dxy)) / denominator;
    float kappa = term1 + term2 + term3;

    return (mAlpha*kappa*(std::sqrt((dx * dx) + (dy * dy) + (dz * dz))));
  }
};

#endif
