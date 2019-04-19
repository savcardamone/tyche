/**
 * @file n_cody.cu
 * @author Salvatore Cardamone
 * @brief Routines for an N-Body simulation which can be ported to the GPU.
 */
#include "n_body/n_body.hpp"

/**
 * @brief Initialise e
 *
 * @param[in] nParticles Number of particles we're dealing with.
 */
__global__
void NBody::InitialiseDevice(const unsigned int& nParticles)
  : nParticles_(nParticles), particles(nParticles), acc(nParticles) {

}

/**
 * @brief
 *
 */
__global__
void NBody::CalculateForces() {

  extern __shared__ float4[] shared_positions;

  // Get the index of the particle the thread is going to be working on
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto thread_particle = particles[tid];

  //
  for (auto i = 0, iTile = 0; i < nParticles; i += nThreadWarp, ++iTile) {

    // Every thread loads a shared position into the shared array
    int idx = tile * blockDim.x + threadIdx.x;
    shared_positions[threadIdx.x] = particles[idx];
    __syncthreads();

    //
    acc = tile_calculation(myPosition, acc);
    __syncthreads();

  }

  // Save the result in global memory for the integration step.
  float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};
  globalA[gtid] = acc4;

}

/**
 * @brief Update the acceleration of particle i as a result of its interaction
 *        with particle j. 19 FLOPs if we count the inverse square root as a
 *        single FLOP, 20 otherwise.
 * @param i Position and mass of particle i.
 * @param j Position and mass of particle j.
 * @param acc_i Acceleration of particle i, updated on return.
 */
__device__
void NBody::PairwiseInteraction(const float4& i, const float4& j,
                                float3& acc_i) {

  // Vector difference between the two particles
  float3 rij;
  rij.x = j.x - i.x;
  rij.y = j.y - i.y;
  rij.z = j.z - i.z;

  // Compute the distance-dependent component of the force
  // Note that we don't need the product of masses since we're skipping the
  // force calculation and going straight to acceleration
  float r2 = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z + eps2;
  float r6 = r2 * r2 * r2;
  float temp = j.w  / sqrt(r6);

  // Update acceleration: Do FMA so we get two FLOPs for the price of one
  acc_i.x += rij.x * temp;
  acc_i.y += rij.y * temp;
  acc_i.z += rij.z * temp;

}

/**
 * @brief Compute the accelerations of the particles in a tile with all others
 *        in the system.
 * @param i Position and mass of particle i.
 */
__device__
void NBody::CalculateTile(const float4& i) {

  // The tile has a shared set of particles. Each thread involved in computing
  // the tile iterates over this shared set of particles
  extern __shared__ float4[] shared_particles;

  // For every particle in the shared particles, compute interaction with
  // particle assigned to thread
  for (auto j=0; j<blockDim.x; ++j) {
    PariwiseInteraction(i, shared_particles[j], acc[i]);
  }

}
