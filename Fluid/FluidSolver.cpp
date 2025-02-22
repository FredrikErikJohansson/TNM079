﻿
#include "FluidSolver.h"

void FluidSolver::AddSolid(Implicit *impl) {
  mSolids.insert(impl);
  std::cerr << "'" << impl->GetName() << "' added as a solid" << std::endl;
}

void FluidSolver::AddFluid(LevelSet *LS) {
  if (mFluids.find(LS) != mFluids.end()) {
    std::cerr << "'" << LS->GetName() << "' already added as a fluid"
              << std::endl;
    return;
  }

  mFluids.insert(LS);
  std::cerr << "'" << LS->GetName() << "' added as a fluid" << std::endl;

  // Update the bounding box (in world coordinates)
  mBox = BoxUnion(mBox, LS->GetBoundingBox());

  // Compute size of the grid (given dx)
  Vector3<float> extent = mBox.pMax - mBox.pMin;
  int dimX = (int)Round(extent[0] / mDx) + 1;
  int dimY = (int)Round(extent[1] / mDx) + 1;
  int dimZ = (int)Round(extent[2] / mDx) + 1;

  std::cerr << "Creating fluid grid of size " << dimX << " x " << dimY << " x "
            << dimZ << std::endl;

  // Create new velocity, voxel grid and solid mask
  mVelocityField = Volume<Vector3<float>>(dimX, dimY, dimZ);
  mVoxels = Volume<float>(dimX, dimY, dimZ);
  mSolidMask = Volume<bool>(dimX, dimY, dimZ);

  // Add the volume
  mInitialVolume += LS->ComputeVolume(LS->GetDx());
}

Vector3<float> FluidSolver::GetValue(float x, float y, float z) const {
  // Transform to grid coordinates
  x = (x - mBox.pMin[0]) / mDx;
  y = (y - mBox.pMin[1]) / mDx;
  z = (z - mBox.pMin[2]) / mDx;

  return mVelocityField.GetValue(x, y, z);
}

Vector3<float> FluidSolver::GetMaxValue() const {
  Vector3<float> maxVal(0.0, 0.0, 0.0);
  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
      for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {

        if (maxVal.Norm() < mVelocityField.GetValue(i, j, k).Norm())
          maxVal = mVelocityField.GetValue(i, j, k);
      }
    }
  }

  return maxVal;
}

int FluidSolver::Solve(float time) {
  // Propagate the solution until requested time is reached
  int iterations = 0;
  for (float elapsed = 0; elapsed < time;) {

    // Determine timestep for stability
    float dt = ComputeTimestep();
    std::cerr << "Propagating solution with dt = " << dt << std::endl;

    if (dt > time - elapsed)
      dt = time - elapsed;
    elapsed += dt;

    // Compute current volume
    mCurrentVolume = 0;
    std::set<LevelSet *>::const_iterator iter = mFluids.begin();
    std::set<LevelSet *>::const_iterator iend = mFluids.end();
    while (iter != iend) {
      mCurrentVolume += (*iter)->ComputeVolume((*iter)->GetDx());
      iter++;
    }
    std::cout << "Current volume: " << mCurrentVolume << std::endl;
    std::cout << "Initial volume: " << mInitialVolume << std::endl;
    std::cout << "Loss of volume: " << mInitialVolume - mCurrentVolume
              << std::endl;

    // Classify all voxels as either solid, fluid or empty
    std::cerr << "Classifying voxels..." << std::endl;
    ClassifyVoxels();

    // Self advection
    std::cerr << "Self advection..." << std::endl;
    SelfAdvection(dt, 6);

    // Add the external forces
    std::cerr << "Adding external forces..." << std::endl;
    ExternalForces(dt);

    // Enforce boundary conditions
    std::cerr << "Boundary conditions..." << std::endl;
    EnforceDirichletBoundaryCondition();

    // Compute the projection for preserving volume
    std::cerr << "Projection..." << std::endl;
    Projection();

    // The projection should not violate the boundary conditions,
    // but it might becacuse of numerical errors (the solution is not
    // exact). Enforce the boundary conditions again to safeguard
    // against this.
    std::cerr << "Boundary conditions..." << std::endl;
    EnforceDirichletBoundaryCondition();

    // Extend the velocites to "air" so the entire level set narrowband
    // is advected by the velocity field
    VelocityExtension();

    std::cerr << "Done with iteration" << std::endl;
    iterations++;
  }

  return iterations;
}

// Compute a stable timestep given the external forces
float FluidSolver::ComputeTimestep() {
  // If there are no external forces the system is unconditionally stable
  if (mExternalForces == NULL)
    return std::numeric_limits<float>::max();

  // Get the max value from the vector field
  Vector3<float> maxVal = mExternalForces->GetMaxValue();

  return 0.7f * mDx / maxVal.Length();
}

float FluidSolver::ComputePotentialEnergy() { return 0; }

float FluidSolver::ComputeKineticEnergy() { return 0; }

// Add the external forces
void FluidSolver::ExternalForces(float dt) {
  if (mExternalForces == NULL)
    return;

  float x, y, z;
  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
        for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {
        // If we're in fluid (see FluidSolver::IsFluid()), sample the external
        // force field (using world coordinates, see
        // FluidSolver::TransformGridToWorld()) and perform the integration to
        // update the velocity field (mVelocityField). The simplest possible
        // integrator is the explicit Euler.
        if (IsFluid(i, j, k)) {
          TransformGridToWorld(i, j, k, x, y, z);
          auto F = mExternalForces->GetValue(x, y, z);
          Vector3<float> V2 = mVelocityField.GetValue(i, j, k) + dt * F;
          mVelocityField.SetValue(i, j, k, V2);
        }
      }
    }
  }
}

// Compute the self advection term
void FluidSolver::SelfAdvection(float dt, int steps) {
  // Copy the current velocity field
  Volume<Vector3<float>> velocities = mVelocityField;

  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
      for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {
        
        // If we're in fluid, sample the current velocity field at (i,j,k).
        // Then, trace a particle at initial position (i,j,k) back in time
        // through the velocity field using 'steps' number of steps. Note that
        // each step is dt/steps time units long. Note also that the velocities
        // are given in world space, but you perform the trace in grid space so
        // you need to scale the velocities accordingly (grid spacing: mDx).
        // When you trace the particle you interpolate the velocities inbetween
        // the grid points, use mVelocityField.GetValue(float i, float j, float
        // k) for trilinear interpolation.
        // TODO: Add code here
        if (IsFluid(i, j, k)) {
        Vector3<float> currentVel = mVelocityField.GetValue(i, j, k);
        Vector3<float> currentPos(i, j, k);
        for (int i = 0; i < steps; i++) {
            currentPos -= (currentVel * (dt / steps))/mDx;
            currentVel = mVelocityField.GetValue(currentPos[0], currentPos[1], currentPos[2]);
        }
        velocities.SetValue(i, j, k, currentVel);
        }
      }
    }
  }

  // Update the current velocity field
  mVelocityField = velocities;
}

// Enforce the Dirichlet boundary conditions
void FluidSolver::EnforceDirichletBoundaryCondition() {
  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
      for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {
        // If we're in fluid, check the neighbors of (i,j,k) to
        // see if it's next to a solid boundary. If so, project
        // the velocity to the boundary plane by setting the
        // velocity to zero along the given dimension.

        //Getvalue may be world-coord
        if (IsFluid(i, j, k)) {
          Vector3<float> field = mVelocityField.GetValue(i, j, k);

          if (IsSolid(i - 1, j, k) && field[0] < 0) {
            field[0] = 0;
          } else if (IsSolid(i + 1, j, k) && field[0] > 0) {
            field[0] = 0;
          }
          
          if (IsSolid(i, j - 1, k) && field[1] < 0) {
            field[1] = 0;
          } else if (IsSolid(i, j + 1, k) && field[1] > 0) {
            field[1] = 0;
          }

          if (IsSolid(i, j, k - 1) && field[2] < 0) {
            field[2] = 0;
          } else if (IsSolid(i, j, k + 1) && field[2] > 0) {
            field[2] = 0;
          }

          mVelocityField.SetValue(i, j, k, field);
        } 
      }
    }
  }
}

// Project the velocity field to preserve the volume
void FluidSolver::Projection() {

  // Compute number of elements in the grid
  int elements = mVoxels.GetDimX() * mVoxels.GetDimY() * mVoxels.GetDimZ();

  // Create sparse matrix and guess that we have 7 non-zero elements
  // per grid point
  CoordMatrix<float, size_t> A(elements, elements);
  A.reserve(elements * 7);
  A.beginPush();

  // Create vectors x, b in the linear system of equations Ax=b
  std::vector<float> x(elements, 0), b(elements, 0);
  
  float dx2 = mDx * mDx;
  float missingVolume = mInitialVolume - mCurrentVolume;

  std::cerr << "Building A matrix and b vector..." << std::endl;
  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
        for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {
        // If we're in fluid...
        if (IsFluid(i, j, k)) {

          
          // Compute the linear indices of (i,j,k) and its neighbors
          // (you need these to index into the A matrix and x,b vectors)
          size_t ind = mVoxels.ComputeLinearIndex(i, j, k);
          size_t ind_ip = mVoxels.ComputeLinearIndex(i + 1, j, k);
          size_t ind_im = mVoxels.ComputeLinearIndex(i - 1, j, k);
          size_t ind_jp = mVoxels.ComputeLinearIndex(i, j + 1, k);
          size_t ind_jm = mVoxels.ComputeLinearIndex(i, j - 1, k);
          size_t ind_kp = mVoxels.ComputeLinearIndex(i, j, k + 1);
          size_t ind_km = mVoxels.ComputeLinearIndex(i, j, k - 1);

          // Compute entry for b vector (divergence of the velocity field:
          // \nabla \dot V_i,j,k)
          // TODO: Add code here
          // u
          auto term1 =
              mVelocityField.GetValue(i + 1, j, k)[0] - mVelocityField.GetValue(i - 1, j, k)[0];
          // v
          auto term2 =
              mVelocityField.GetValue(i, j + 1, k)[1] - mVelocityField.GetValue(i, j - 1, k)[1];
          // w
          auto term3 =
              mVelocityField.GetValue(i , j, k + 1)[2] - mVelocityField.GetValue(i , j, k - 1)[2];

          b[ind] = ((term1 + term2 + term3) / (2.0f * mDx));

          if (missingVolume > 0)
            b[ind] -= missingVolume*100;
          

          // Compute entries for A matrix (discrete Laplacian operator).
          // The A matrix is a sparse matrix but can be used like a regular
          // matrix. That is, you access the elements by A(row, column).
          // However, due to the matrix data structure you cannot read
          // elements at this point (until you do A.endPush(), see below).
          // So, only use A(row, column) = ... to set a value in the matrix,
          // don't use A(row, column) to get a value.
          // Remember to enforce the boundary conditions if we're next to
          // a solid (allow no change of flow in that direction).
          // Remember to treat the boundaries of (i,j,k).
          // TODO: Add code here

          int neighbours[6] = {!IsSolid(i + 1, j, k), !IsSolid(i - 1, j, k), !IsSolid(i, j + 1, k),
                        !IsSolid(i, j - 1, k), !IsSolid(i, j, k + 1), !IsSolid(i, j, k - 1)};
          int current = 0;
          for (int i = 0; i < 6; i++) current += neighbours[i];
          A(ind, ind_ip) = neighbours[0] / dx2;
          A(ind, ind_im) = neighbours[1] / dx2;
          A(ind, ind_jp) = neighbours[2] / dx2;
          A(ind, ind_jm) = neighbours[3] / dx2;
          A(ind, ind_kp) = neighbours[4] / dx2;
          A(ind, ind_km) = neighbours[5] / dx2;
          A(ind, ind) = -current / dx2;
        }
      }
    }
  }

  // Rebuild the sparse matrix structure
  A.endPush();

  // Solve Ax=b using conjugate gradient
  std::cerr << "Conjugate gradient solver... ";
  ConjugateGradient<CoordMatrix<float, size_t>, std::vector<float>, float> CG(
      100, 1e-3f);
  CG.solve(A, x, b);
  std::cerr << "finished with tolerance " << CG.getTolerance() << " in "
            << CG.getNumIter() << " iterations" << std::endl;

  // Subtract the gradient of x to preserve the volume
  for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
    for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
      for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {

        // If we're in fluid...
        if (IsFluid(i, j, k)) {

          // Compute the linear indices of (i,j,k) and its neighbors
          size_t ind_ip = mVoxels.ComputeLinearIndex(i + 1, j, k);
          size_t ind_im = mVoxels.ComputeLinearIndex(i - 1, j, k);
          size_t ind_jp = mVoxels.ComputeLinearIndex(i, j + 1, k);
          size_t ind_jm = mVoxels.ComputeLinearIndex(i, j - 1, k);
          size_t ind_kp = mVoxels.ComputeLinearIndex(i, j, k + 1);
          size_t ind_km = mVoxels.ComputeLinearIndex(i, j, k - 1);

          // Compute the gradient of x at (i,j,k) using central differencing
          // and subtract this gradient from the velocity field.
          // Thereby removing divergence - preserving volume.
          // TODO: Add code here
          Vector3<float> gradX = {x[ind_ip] - x[ind_im], x[ind_jp] - x[ind_jm],
                                  x[ind_kp] - x[ind_km]};
          gradX /= (2.0f * mDx);
          Vector3<float> field = mVelocityField.GetValue(i, j, k) - gradX;
          mVelocityField.SetValue(i, j, k, field);
        }
      }
    }
  }
}

// Extent the velocities to "air"
void FluidSolver::VelocityExtension() {
  auto iterations = 0;
  auto narrowbandWidth = 0.0f;
  auto dt = 0.7f;
  std::set<LevelSet *>::const_iterator iter = mFluids.begin();
  std::set<LevelSet *>::const_iterator iend = mFluids.end();
  while (iter != iend) {
    iterations =
        std::max(iterations, (int)((*iter)->GetNarrowBandWidth() * 0.5f / dt));
    narrowbandWidth = std::max(narrowbandWidth, (*iter)->GetNarrowBandWidth() *
                                                    (*iter)->GetDx() * 0.5f);
    iter++;
  }

  iterations = std::min(iterations, (int)(mVoxels.GetDimX() / dt));
  iterations = std::min(iterations, (int)(mVoxels.GetDimY() / dt));
  iterations = std::min(iterations, (int)(mVoxels.GetDimZ() / dt));

  std::cerr << "Velocity extension (" << iterations << " iterations)..."
            << std::endl;
  for (int iter = 0; iter < iterations; iter++) {

    Volume<Vector3<float>> velocities = mVelocityField;

    for (size_t i = 0; i < mVoxels.GetDimX(); i++) {
      for (size_t j = 0; j < mVoxels.GetDimY(); j++) {
        for (size_t k = 0; k < mVoxels.GetDimZ(); k++) {

          const float val = mVoxels.GetValue(i, j, k);

          // Extend only in air or solid (don't extend beyond narrowband for
          // efficiency)
          if (IsFluid(i, j, k) || val >= narrowbandWidth)
            continue;

          Vector3<float> normal(
              mVoxels.GetValue(i + 1, j, k) - mVoxels.GetValue(i - 1, j, k),
              mVoxels.GetValue(i, j + 1, k) - mVoxels.GetValue(i, j - 1, k),
              mVoxels.GetValue(i, j, k + 1) - mVoxels.GetValue(i, j, k - 1));

          if (normal.Length() > 0)
            normal.Normalize();

          Vector3<float> v = mVelocityField.GetValue(i, j, k);

          Vector3<float> ip = mVelocityField.GetValue(i + 1, j, k);
          Vector3<float> im = mVelocityField.GetValue(i - 1, j, k);

          Vector3<float> jp = mVelocityField.GetValue(i, j + 1, k);
          Vector3<float> jm = mVelocityField.GetValue(i, j - 1, k);

          Vector3<float> kp = mVelocityField.GetValue(i, j, k + 1);
          Vector3<float> km = mVelocityField.GetValue(i, j, k - 1);

          Vector3<float> v0 = v - dt * (std::max(0.0f, normal[0]) * (v - im) +
                                        std::min(0.0f, normal[0]) * (ip - v) +
                                        std::max(0.0f, normal[1]) * (v - jm) +
                                        std::min(0.0f, normal[1]) * (jp - v) +
                                        std::max(0.0f, normal[2]) * (v - km) +
                                        std::min(0.0f, normal[2]) * (kp - v));

          velocities.SetValue(i, j, k, v0);
        }
      }
    }

    mVelocityField = velocities;
  }
}

void FluidSolver::ClassifyVoxels() {
  // Iterate the domain and classify all voxels as either
  // fluid, solid or empty
  for (int i = 0; i < mVoxels.GetDimX(); i++) {
    for (int j = 0; j < mVoxels.GetDimY(); j++) {
      for (int k = 0; k < mVoxels.GetDimZ(); k++) {

        ClassifyVoxel(i, j, k);
      }
    }
  }
}

void FluidSolver::ClassifyVoxel(int i, int j, int k) {
  // Transform grid to world coordinates
  float x, y, z;
  TransformGridToWorld(i, j, k, x, y, z);

  // Check if we're inside a solid and update the solid mask
  bool solid = false;
  auto iterSolid = mSolids.begin();
  auto iendSolid = mSolids.end();
  while (iterSolid != iendSolid && !solid) {
    if ((*iterSolid)->GetValue(x, y, z) <= 0)
      solid = true;
    iterSolid++;
  }
  mSolidMask.SetValue(i, j, k, solid);

  // Compute the closest distance to all fluids
  float distance = std::numeric_limits<float>::max();
  auto iterFluid = mFluids.begin();
  auto iendFluid = mFluids.end();
  while (iterFluid != iendFluid) {
    distance = std::min(distance, (*iterFluid)->GetValue(x, y, z));
    iterFluid++;
  }
  mVoxels.SetValue(i, j, k, distance);
}

bool FluidSolver::IsSolid(size_t i, size_t j, size_t k) const {
  return mSolidMask.GetValue(i, j, k);
}

bool FluidSolver::IsFluid(size_t i, size_t j, size_t k) const {
  return mVoxels.GetValue(i, j, k) <= 0.5f * mDx;
}
